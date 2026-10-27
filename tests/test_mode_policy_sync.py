"""Exercise actual entry-point policy functions without loading input datasets."""
import argparse
import ast
import contextlib
import io
import itertools
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

ROOT = Path(os.environ.get('PHYLOSOLID_ROOT', Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT))
from src.pruning_policy import pruning_ratio

MODES = ('scRNA', 'scDNA', 'spacetracer', 'binput')


def source(mode):
    return (ROOT / 'src' / f'run_phylosilid_fullTree_{mode}.py').read_text()


def load_function(mode, name, namespace):
    node = next(n for n in ast.parse(source(mode)).body
                if isinstance(n, ast.FunctionDef) and n.name == name)
    exec(compile(ast.Module(body=[node], type_ignores=[]), mode, 'exec'), namespace)
    return namespace[name]


def candidate(cv, pre, final, stage='final'):
    return dict(cv_value=cv, omega_pre_qc=pre, omega_final=final,
                output_stage=stage, scaffold_count=6, retention_rate=.9,
                final_cells=105, final_muts=24)


class ModePolicyTests(unittest.TestCase):
    def parse_cli(self, mode, argv):
        namespace = {'argparse': argparse}
        with patch.object(sys, 'argv', ['test'] + argv):
            if mode == 'binput':
                return load_function(mode, 'parse_args', namespace)()
            text = source(mode)
            text = text[text.index('parser = argparse.ArgumentParser()'):text.index('# Set all random seeds')]
            exec(compile(text, mode, 'exec'), namespace)
            return namespace['args']

    def test_backoff_options_all_modes(self):
        for mode in MODES:
            with self.subTest(mode=mode):
                self.assertEqual(self.parse_cli(mode, []).resolution_backoff_threshold, .1)
                self.assertIsNone(self.parse_cli(mode, ['--disable_resolution_backoff']).resolution_backoff_threshold)
                for value in ['0', '0.2', '1']:
                    self.assertEqual(self.parse_cli(mode, ['--resolution_backoff_threshold', value]).resolution_backoff_threshold, float(value))
                invalid = [['--resolution_backoff_threshold', x] for x in ['nan', 'inf', '-.1', '1.1']]
                invalid.append(['--disable_resolution_backoff', '--resolution_backoff_threshold', '.1'])
                for argv in invalid:
                    with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                        self.parse_cli(mode, argv)

    def test_cv_policy_independent_of_candidate_order(self):
        cases = [
            ([candidate(.3,0,0), candidate(.9,1,1)], .3),
            ([candidate(.3,1,0), candidate(.9,1,1)], .9),
            ([candidate(.3,0,1,'pre_fallback'), candidate(.9,2,1)], .9),
            ([candidate(.3,2,1), candidate(.9,2,1)], .3),
            ([candidate(.3,30,1), candidate(.9,20,1)], .9),
            ([candidate(.3,20,1), candidate(.9,20,1)], .3),
            ([candidate(.3,1,2,'pre_fallback'), candidate(.9,1,0,'pre_fallback')], .3),
        ]
        for mode in MODES[:-1]:
            choose = load_function(mode, 'select_best_cv', dict(pd=pd, compute_pruning_ratio=pruning_ratio, PRUNING_CONFIDENCE_THRESHOLD=10))
            self.assertIsNone(choose(pd.DataFrame())[0])
            for rows, expected in cases:
                for perm in itertools.permutations(rows):
                    with self.subTest(mode=mode, rows=perm):
                        result = choose(pd.DataFrame(perm))
                        self.assertEqual(result[1], expected)
                        if all(row['output_stage'] == 'pre_fallback' for row in rows):
                            self.assertFalse(result[3]['is_confident'].any())

    def test_summary_preserves_fallback_status(self):
        import logging
        import time
        for mode in MODES[:-1]:
            save = load_function(mode, 'save_comprehensive_results', dict(pd=pd, os=os, time=time, logger=logging.getLogger(__name__), PRUNING_CONFIDENCE_THRESHOLD=10))
            row = dict(candidate(.3,0,0,'pre_fallback'), success=True,
                       zero_final_omega=True, is_confident=False, pruning_ratio=0)
            with tempfile.TemporaryDirectory() as temp:
                save([row], temp, 'sample')
                frame = pd.read_csv(Path(temp, 'cv_search_summary.csv'))
                self.assertEqual(frame.loc[0,'output_stage'], 'pre_fallback')
                self.assertTrue(frame.loc[0,'zero_final_omega'])
                self.assertFalse(frame.loc[0,'confident'])

    def test_resolution_defaults_preserved(self):
        for mode in MODES:
            node = next(n for n in ast.parse(source(mode)).body
                        if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='SETTING_PARAMS' for t in n.targets))
            values=ast.literal_eval(node.value)
            self.assertEqual(values['resolution_of_graph'], .1 if mode=='scDNA' else None)
            self.assertEqual((values['min_resolution'],values['max_resolution']),(.5,2.0))


if __name__ == '__main__':
    unittest.main(verbosity=2)
