"""Regression tests for the selected Step 7/Step 8 export bundles."""
import ast
import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

ROOT = Path(os.environ.get('PHYLOSOLID_ROOT', Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT))
from src.scaffold_builder import TreeNode
from src.phylo_export import export_final_phylo_results, prepare_phylo_matrix


def tree(include_artifact=True):
    root = TreeNode('ROOT')
    branch = TreeNode('A|B')
    root.add_child(branch)
    branch.add_child(TreeNode('G'))
    if include_artifact:
        branch.add_child(TreeNode('C'))
    return root


class StageExportTests(unittest.TestCase):
    def setUp(self):
        self.pre = pd.DataFrame({'A|B':[1,1,0], 'C':[0,1,0], 'G':[0,0,0]}, index=['c1','c2','zero'])
        self.observed = pd.DataFrame({'A':[1,1,0], 'B':[0,1,0], 'C':[0,1,0], 'G':[1,0,0]}, index=self.pre.index)
        self.snapshot = (tree(), self.pre.copy(), [], ['unplaced'])
        self.post = self.pre.loc[['c1'], ['A|B','G']].copy()

    def export(self, temp, fallback=False):
        return export_final_phylo_results(
            temp, self.snapshot[0] if fallback else tree(False),
            prepare_phylo_matrix(self.pre if fallback else self.post, []),
            self.observed, .1, pre_qc_snapshot=self.snapshot,
            output_stage='pre_fallback' if fallback else 'final')

    def test_real_exports_retain_pre_qc_cells_and_mutations_only_in_unpruned(self):
        with tempfile.TemporaryDirectory() as temp:
            result = self.export(temp)
            pre_dir = Path(temp, 'phylo_unpruned')
            post_dir = Path(temp, 'phylo')
            pre_I = pd.read_csv(pre_dir/'I_full_withNA3.txt', sep='\t', index_col=0)
            post_I = pd.read_csv(post_dir/'final_cleaned_I_full_withNA3_for_circosPlot.txt', sep='\t', index_col=0)
            self.assertEqual(list(pre_I.index), ['c1','c2'])
            self.assertEqual(list(pre_I.columns), ['A','B','C'])
            self.assertEqual(list(post_I.index), ['c1'])
            self.assertEqual(list(post_I.columns), ['A','B'])
            for directory, stem, expected, stage in [
                (pre_dir,'tree_node', {'A','B','C'},'pre_qc'),
                (post_dir,'final_cleaned_tree_node',{'A','B'},'final'),
            ]:
                matrix_stem = ('M_full_basedPivots.filtered_sites_inferred' if stage == 'pre_qc'
                               else 'final_cleaned_M_full_basedPivots.filtered_sites_inferred')
                inferred = pd.read_csv(directory/(matrix_stem+'.CFMatrix'), sep='\t', index_col=0)
                observations = pre_I if stage == 'pre_qc' else post_I
                self.assertEqual(list(inferred.index), list(observations.index))
                self.assertEqual(list(inferred.columns), list(observations.columns))
                self.assertFalse(inferred.eq(0).all(axis=0).any())
                self.assertFalse(inferred.eq(0).all(axis=1).any())
                text=(directory/(stem+'.txt')).read_text()
                self.assertNotIn('G', text)
                self.assertEqual('C' in text, 'C' in expected)
                metadata=json.loads((directory/'output_stage.json').read_text())
                self.assertEqual(metadata['output_stage'],stage)
                self.assertTrue(metadata['all_zero_rows_columns_removed'])
                flips=pd.read_csv(directory/'df_total_flipping_count.txt',sep='\t')
                self.assertAlmostEqual(flips.loc[0,'weighted_discordance_index_Omega'], metadata['omega'])
            self.assertAlmostEqual(result['omega_pre_qc'], .1)
            self.assertAlmostEqual(result['omega'], .1)
            pd.testing.assert_frame_equal(self.snapshot[1],self.pre)
            self.assertIn('G', self.snapshot[0].all_names())

    def test_fallback_bundles_equal_and_explicitly_marked(self):
        with tempfile.TemporaryDirectory() as temp:
            result=self.export(temp,True)
            pre=pd.read_csv(Path(temp,'phylo_unpruned/I_full_withNA3.txt'),sep='\t',index_col=0)
            post=pd.read_csv(Path(temp,'phylo/final_cleaned_I_full_withNA3_for_circosPlot.txt'),sep='\t',index_col=0)
            pd.testing.assert_frame_equal(pre,post)
            self.assertEqual(Path(temp,'phylo_unpruned/tree_node.txt').read_text(),Path(temp,'phylo/final_cleaned_tree_node.txt').read_text())
            self.assertEqual(json.loads(Path(temp,'phylo/output_stage.json').read_text())['output_stage'],'pre_fallback')
            self.assertEqual(result['omega'],result['omega_pre_qc'])

    def test_root_and_merged_mutations_share_metric_population(self):
        matrix=pd.DataFrame({'ROOT':[1], 'A|B':[1]},index=['c1'])
        expanded=prepare_phylo_matrix(matrix,['R','R','A'])
        self.assertEqual(list(expanded.columns),['R','A','B'])
        self.assertTrue(expanded.eq(1).all().all())
        root=tree(False)
        root.children[0].children=[]
        root.add_child(TreeNode('R'))
        obs=pd.DataFrame({'R':[0], 'A':[1], 'B':[0]},index=['c1'])
        with tempfile.TemporaryDirectory() as temp:
            result=export_final_phylo_results(temp,root,expanded,obs,.1,pre_qc_snapshot=(root,matrix,['R'],[]))
            self.assertAlmostEqual(result['omega'],.2)
            self.assertEqual(result['omega'],result['omega_pre_qc'])

    def test_missing_pre_observations_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(ValueError,'Missing observations'):
                export_final_phylo_results(temp,tree(False),prepare_phylo_matrix(self.post,[]),self.observed.loc[['c1']],.1,pre_qc_snapshot=self.snapshot)

    def test_empty_pre_qc_cleanup_writes_empty_bundle(self):
        empty_pre = pd.DataFrame({'G':[0]}, index=['zero'])
        empty_tree = TreeNode('ROOT')
        empty_tree.add_child(TreeNode('G'))
        with tempfile.TemporaryDirectory() as temp:
            result = export_final_phylo_results(
                temp, tree(False), prepare_phylo_matrix(self.post, []), self.observed, .1,
                pre_qc_snapshot=(empty_tree, empty_pre, [], []))
            self.assertEqual(result['unpruned']['cells'], 0)
            self.assertEqual(result['unpruned']['muts'], 0)
            self.assertTrue(Path(temp, 'phylo_unpruned/M_full_basedPivots.filtered_sites_inferred.CFMatrix').exists())

    def test_all_modes_wire_selected_snapshot_and_stage(self):
        for mode in ['scRNA','scDNA','spacetracer','binput']:
            module=ast.parse((ROOT/'src'/f'run_phylosilid_fullTree_{mode}.py').read_text())
            calls=[n for n in ast.walk(module) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='export_final_phylo_results']
            self.assertEqual(len(calls),1)
            self.assertTrue({'pre_qc_snapshot','output_stage'}.issubset({kw.arg for kw in calls[0].keywords}),mode)


if __name__=='__main__':
    unittest.main(verbosity=2)
