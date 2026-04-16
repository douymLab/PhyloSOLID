"""
SpaceTracer step implementations - custom tree building step
"""

import logging
import sys
import subprocess
from pathlib import Path
from typing import Dict, Any, Optional

from ..base import PipelineStep
from utils.command import CommandRunner, RScriptRunner


class SpaceTracerTreeBuildingStep(PipelineStep):
    """Tree building step for SpaceTracer"""
    
    def __init__(self, workdir: Path, script_dir: Path, config: Dict[str, Any] = None):
        super().__init__("03_tree_building", workdir, config)
        self.script_dir = Path(script_dir) / "scrna" / "tree_building"
        self.main_script = Path(__file__).parent.parent.parent / 'src' / 'run_phylosilid_fullTree_spacetracer.py'
        self.runner = CommandRunner(workdir=self.workdir)
        self.r_runner = RScriptRunner(workdir=self.workdir)
    
    def _prepare_celltype_file(self, sample_id: str, metadata_file: Path, 
                               barcode_file: Path) -> Optional[Path]:
        output_file = self.workdir / f"celltype_file_for_{sample_id}.txt"
        
        if output_file.exists():
            return output_file
        
        script_path = Path(__file__).parent.parent.parent / 'scripts' / 'scrna' / 'celltype_anno' / 'celltype_annotation_using_scATOMIC.R'
        
        if not script_path.exists():
            self.logger.warning(f"Celltype annotation script not found: {script_path}")
            return None
        
        cmd = [
            'Rscript',
            str(script_path),
            '--metadata', str(metadata_file),
            '--barcode', str(barcode_file),
            '--output', str(output_file),
            '--sample', sample_id
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            if output_file.exists():
                return output_file
            return None
        except subprocess.CalledProcessError:
            return None
    
    def _execute(self, sample_id: str, cellnum: int, 
                 celltype_file: Optional[Path] = None,
                 metadata_file: Optional[Path] = None,
                 barcode_file: Optional[Path] = None,
                 features_file: Optional[Path] = None) -> Dict[str, Any]:
        
        self.logger.info(f"Building tree for {sample_id}")
        
        # Determine final celltype file
        final_celltype_file = None
        
        if celltype_file and celltype_file.exists():
            final_celltype_file = celltype_file
        elif metadata_file and barcode_file:
            final_celltype_file = self._prepare_celltype_file(sample_id, metadata_file, barcode_file)
        
        if final_celltype_file is None:
            self.logger.warning("No celltype information available, proceeding without celltype file")
        
        # Features file is None for spacetracer (skip classifier)
        self.logger.info("No features file provided (classifier will be skipped)")
        
        if not self.main_script.exists():
            raise FileNotFoundError(f"Tree building script not found: {self.main_script}")
        
        # Point to the data directory created by tree_input step (in 02_treeinput, not 03_treebuilding)
        treeinput_dir = self.workdir.parent / '02_treeinput'
        data_dir = treeinput_dir / 'data'
        
        if not data_dir.exists():
            self.logger.error(f"Data directory not found: {data_dir}")
            self.logger.error("Please ensure tree_input step has been run successfully")
            raise FileNotFoundError(f"Data directory not found: {data_dir}")
        
        outputpath = self.workdir.absolute()
        celltype_arg = str(final_celltype_file) if final_celltype_file else "None"
        
        # Build command
        cmd = [
            sys.executable,
            str(self.main_script),
            '-i', str(data_dir),
            '-o', str(outputpath),
            '-s', sample_id,
            '-c', celltype_arg,
            '--features_file', 'None',
            '--is_predict_germ', 'no',
            '--is_detect_passtree_by_dp', 'no',
            '--is_filter_quality', 'yes'
        ]
        
        self.logger.info(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if result.stdout:
                self.logger.info(result.stdout)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Tree building script failed: {e}")
            self.logger.error(f"STDERR: {e.stderr}")
            raise
        
        # Find tree file
        tree_file = self.workdir / 'PhyloSOLID' / 'celltree.newick'
        if tree_file.exists():
            self.outputs['tree_file'] = tree_file
            self.logger.info(f"Tree file generated: {tree_file}")
        
        return {
            'sample_id': sample_id,
            'tree_file': str(tree_file) if tree_file.exists() else None,
            'results_dir': str(self.workdir)
        }
