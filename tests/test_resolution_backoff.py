import unittest
from types import SimpleNamespace
import igraph as ig
import numpy as np
import pandas as pd
from src.resolution_backoff import partition_conflicts, refine_resolution


class BackoffTests(unittest.TestCase):
    def test_cross_group_is_not_necessarily_high_fp(self):
        matrix = pd.DataFrame([[1]*10], columns=list('abcdefghij'))
        groups = {m: int(m == 'j') for m in matrix.columns}
        s, cells = partition_conflicts(matrix, groups)
        self.assertEqual(s['multi_group_fraction'], 1)
        self.assertEqual(s['high_conflict_fraction'], 0)
        self.assertAlmostEqual(cells.iloc[0].fp_proxy, 0.1)

    def test_inclusive_cutoff_missingness_and_population(self):
        matrix = pd.DataFrame([[1,1],[1,np.nan],[0,0],[np.nan,np.nan]], columns=['a','b'])
        s, cells = partition_conflicts(matrix, {'a':0,'b':1})
        self.assertEqual(s['n_informative_cells'], 2)
        self.assertEqual(s['high_conflict_fraction'], 0.5)
        self.assertTrue(cells.iloc[0].high_conflict)
        self.assertTrue(np.isnan(cells.iloc[2].fp_proxy))

    def test_real_leiden_complete_graph_avoids_singletons(self):
        import leidenalg as la
        graph = ig.Graph.Full(12); graph.vs['name'] = ['m'+str(i) for i in range(12)]
        graph.es['weight'] = [1.0]*graph.ecount()
        matrix = pd.DataFrame(np.ones((8,12)), columns=graph.vs['name'])
        initial = la.find_partition(graph, la.RBConfigurationVertexPartition,
                                    weights='weight', resolution_parameter=1.2, seed=42)
        self.assertEqual(len(initial), 12)
        result = refine_resolution(graph, matrix, initial, 1.2, threshold=0.1)
        self.assertTrue(result['threshold_met'])
        self.assertAlmostEqual(result['resolution'], 1.0)
        self.assertEqual(len(result['partition']), 1)
        self.assertEqual(len(result['trace']), 3)

    def test_nonmonotone_path_first_passing_candidate(self):
        graph=ig.Graph.Full(4);graph.vs['name']=list('abcd')
        matrix=pd.DataFrame([[1,1,0,0],[0,0,1,1]],columns=list('abcd'))
        initial=SimpleNamespace(membership=[0,1,2,3])
        calls=[]
        def backend(g,res,seed):
            calls.append(res)
            return SimpleNamespace(membership=[0,0,1,2] if res>1.15 else
                                   [0,1,2,3] if res>1.05 else [0,0,1,1])
        r=refine_resolution(graph,matrix,initial,1.3,0.1,partitioner=backend)
        self.assertEqual(calls,[1.2,1.1,1.0])
        self.assertEqual(r['resolution'],1.0)
        self.assertEqual([s['high_conflict_fraction'] for s in r['trace']],[1.0,0.5,1.0,0.0])
        self.assertEqual([s['n_informative_cells'] for s in r['trace']],[2,2,2,2])

    def test_supported_distinct_groups_are_not_merged(self):
        graph=ig.Graph.Full(4);graph.vs['name']=list('abcd')
        matrix=pd.DataFrame([[1,1,0,0],[0,0,1,1]],columns=list('abcd'))
        p=SimpleNamespace(membership=[0,0,1,1])
        def backend(*args):
            self.fail('No extra optimisation should occur when the initial grouping passes')
        r=refine_resolution(graph,matrix,p,1.2,0.1,partitioner=backend)
        self.assertEqual(r['resolution'],1.2)
        self.assertEqual(len(r['trace']),1)

    def test_unmet_budget_retains_initial_not_unjustified_merging(self):
        graph=ig.Graph.Full(2);graph.vs['name']=['a','b']
        matrix=pd.DataFrame([[1,1]],columns=['a','b'])
        p=SimpleNamespace(membership=[0,1]);calls=[]
        def backend(g,res,seed):calls.append(res);return p
        r=refine_resolution(graph,matrix,p,0.6,0.0,partitioner=backend)
        self.assertFalse(r['threshold_met'])
        self.assertEqual(r['resolution'],0.6)
        self.assertEqual(calls,[0.5])
        self.assertEqual(r['stop_reason'],'resolution_floor_threshold_unmet')

    def test_zero_information_does_not_trigger_backoff(self):
        graph=ig.Graph.Full(2);graph.vs['name']=['a','b']
        p=SimpleNamespace(membership=[0,1])
        r=refine_resolution(graph,pd.DataFrame([[0,np.nan]],columns=['a','b']),p,1.2,0.1)
        self.assertEqual(r['stop_reason'],'no_informative_cells')
        self.assertEqual(len(r['trace']),1)

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(ValueError):
            partition_conflicts(pd.DataFrame([[3]],columns=['a']),{'a':0})
        with self.assertRaises(ValueError):
            partition_conflicts(pd.DataFrame([[1]],columns=['a']),{'missing':0})
        graph=ig.Graph.Full(1);graph.vs['name']=['a']
        with self.assertRaises(ValueError):
            refine_resolution(graph,pd.DataFrame([[1]],columns=['a']),SimpleNamespace(membership=[0]),1.2,float('nan'))


class BackoffCliTests(unittest.TestCase):
    def parse_actual_cli(self, argv):
        import sys
        from pathlib import Path
        from unittest.mock import patch
        source = (Path(__file__).resolve().parents[1] / 'src/run_phylosilid_fullTree_scRNA.py').read_text()
        source = source[source.index('import argparse\n'):source.index('# Set all random seeds')]
        with patch.object(sys, 'argv', ['test'] + argv):
            namespace = {}
            exec(compile(source, 'actual_cli', 'exec'), namespace)
        return namespace['args'].resolution_backoff_threshold

    def test_default_and_explicit_disable(self):
        self.assertEqual(self.parse_actual_cli([]), 0.10)
        self.assertIsNone(self.parse_actual_cli(['--disable_resolution_backoff']))
        self.assertEqual(self.parse_actual_cli(['--resolution_backoff_threshold', '0.2']), 0.2)

    def test_invalid_budgets_and_conflicting_options(self):
        import contextlib
        import io
        for argv in [['--resolution_backoff_threshold', value] for value in ['nan', 'inf', '-0.1', '1.1']] + [['--disable_resolution_backoff', '--resolution_backoff_threshold', '0.1']]:
            with contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    self.parse_actual_cli(argv)
                self.assertEqual(error.exception.code, 2)


if __name__=='__main__':
    unittest.main(verbosity=2)
