"""
SpaceTracer pipeline - skips feature extraction, uses custom tree building
"""

import logging
from pathlib import Path
from typing import Dict, Any, Optional, List

from ..base import Pipeline
from ..scrna.steps import SCRNAFeatureExtractionStep, SCRNATreeInputStep
from .steps import SpaceTracerTreeBuildingStep


class SpaceTracerPipeline(Pipeline):
    """
    SpaceTracer Pipeline
    
    Skips feature extraction and uses custom tree building step.
    All other parameters and workflows remain the same as scRNA mode.
    """
    
    def __init__(self, workdir: Path, script_dir: Path, config: Dict[str, Any] = None):
        """Initialize SpaceTracer pipeline"""
        super().__init__(workdir, config)
        self.script_dir = script_dir
        self.logger = logging.getLogger(__name__)
        
        # Add steps - note: feature_extraction is omitted
        self.add_step('tree_input', 
                     SCRNATreeInputStep(self.workdir, script_dir, config))
        self.add_step('tree_building', 
                     SpaceTracerTreeBuildingStep(self.workdir, script_dir, config))
        
        self.logger.info("SpaceTracer pipeline initialized")
        self.logger.info("  - Feature extraction is skipped")
        self.logger.info("  - Using custom tree building step")
    
    def get_tree_file(self) -> Optional[Path]:
        """Get the final tree file path"""
        if 'tree_building' in self.steps:
            return self.steps['tree_building'].get_output('tree_file')
        return None
    
    def run(self, sample_id: str, mutation_list: Path, bam_file: Path,
            barcode_file: Path, celltype_file: Optional[Path] = None,
            metadata_file: Optional[Path] = None,
            read_len: int = 91, cellnum: int = 155, threads: int = 4,
            running_type: str = "benchmark", ase_filepath: Optional[str] = None,
            steps: List[str] = None, parallel: bool = False) -> Dict[str, Any]:
        """
        Run SpaceTracer pipeline (skips feature extraction)
        """
        # Default steps exclude feature_extraction
        if steps is None:
            steps = ['tree_input', 'tree_building']
        
        self.logger.info("=" * 50)
        self.logger.info("Running SpaceTracer pipeline (feature extraction skipped)")
        self.logger.info(f"Sample: {sample_id}, Steps: {steps}")
        self.logger.info("=" * 50)
        
        # Prepare step arguments
        step_kwargs = {
            'tree_input': {
                'sample_id': sample_id,
                'mutation_list': mutation_list,
                'bam_file': bam_file,
                'barcode_file': barcode_file,
                'cellnum': cellnum,
                'threads': min(threads, 2)
            },
            'tree_building': {
                'sample_id': sample_id,
                'cellnum': cellnum,
                'celltype_file': celltype_file,
                'metadata_file': metadata_file,
                'barcode_file': barcode_file,
                'features_file': None  # Skip classifier
            }
        }
        
        # Run steps sequentially
        for step_name in steps:
            if step_name in self.steps:
                self.run_step(step_name, **step_kwargs.get(step_name, {}))
        
        # Add summary
        self.results['summary'] = {
            'sample_id': sample_id,
            'workdir': str(self.workdir),
            'steps_completed': list(self.results.keys()),
            'tree_file': self.get_tree_file()
        }
        
        return self.results