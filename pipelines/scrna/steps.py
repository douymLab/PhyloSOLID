"""
Individual step implementations for scRNA pipeline
"""

import logging
from pathlib import Path
from typing import Dict, Any, Optional
import subprocess
import sys

from ..base import PipelineStep
from utils.command import CommandRunner, RScriptRunner

class SCRNAFeatureExtractionStep(PipelineStep):
    """Feature extraction step - calls external bash script"""
    
    def __init__(self, workdir: Path, script_dir: Path, config: Dict[str, Any] = None):
        """
        Initialize feature extraction step
        
        Args:
            workdir: Working directory
            script_dir: Directory containing external scripts
            config: Configuration dictionary
        """
        super().__init__("01_features", workdir, config)
        self.script_dir = Path(script_dir) / "scrna" / "feature_extraction"
        self.main_script = self.script_dir / "scRNA_sites_grep_all_features_and_filtration_withASE_and_patched.whole_steps.sh"
        self.runner = CommandRunner(workdir=self.workdir)
        self.r_runner = RScriptRunner(workdir=self.workdir)
        
    def _execute(self, sample_id: str, mutation_list: Path, bam_file: Path,
                 barcode_file: Path, read_len: int = 91, 
                 running_type: str = "benchmark",
                 ase_filepath: Optional[str] = None,
                 threads: int = 4) -> Dict[str, Any]:
        """
        Execute feature extraction
        
        Args:
            sample_id: Sample identifier
            mutation_list: Path to mutation list file
            bam_file: Path to BAM file
            barcode_file: Path to barcode file
            read_len: Read length
            running_type: Running type (benchmark/production)
            ase_filepath: ASE file path (optional)
            threads: Number of threads
            
        Returns:
            Dictionary with step results
        """
        self.logger.info(f"Extracting features for {sample_id}")
        
        # Check if script exists
        if not self.main_script.exists():
            raise FileNotFoundError(f"Feature extraction script not found: {self.main_script}")
        
        # Prepare arguments matching the original script's 10 parameters
        args = [
            sample_id,
            str(mutation_list.absolute()),
            str(bam_file.absolute()),
            str(read_len),
            str(self.workdir.absolute()),           # out_dir_name
            str(bam_file.absolute()),               # filtered_bam
            str(threads),
            str(barcode_file.absolute()),
            ase_filepath if ase_filepath else "no",
            running_type
        ]
        
        self.logger.info(f"Running feature extraction script with {len(args)} parameters")
        self.runner.run_script(self.main_script, args)
        
        # Find output files matching real output structure
        patched_file = self.workdir / f"{sample_id}.{running_type}_patched.feature.txt"
        if patched_file.exists():
            self.outputs['features_file'] = patched_file
            self.logger.info(f"Features file generated: {patched_file}")
        else:
            self.logger.warning(f"Features file not found: {patched_file}")
        
        # Record depth_in_spots directory
        depth_dir = self.workdir / "depth_in_spots"
        if depth_dir.exists():
            self.outputs['depth_dir'] = depth_dir
            self.logger.info(f"Depth directory created: {depth_dir}")
        
        return {
            'sample_id': sample_id,
            'features_file': str(patched_file) if patched_file.exists() else None,
        }


class SCRNATreeInputStep(PipelineStep):
    """Tree input generation step"""
    
    def __init__(self, workdir: Path, script_dir: Path, config: Dict[str, Any] = None):
        """
        Initialize tree input step
        
        Args:
            workdir: Working directory
            script_dir: Directory containing external scripts
            config: Configuration dictionary
        """
        super().__init__("02_treeinput", workdir, config)
        self.script_dir = Path(script_dir) / "scrna" / "tree_input"
        self.main_script = self.script_dir / "run_PhyloSOLID.treeinput_data.scRNAmode_using_identifiers.sh"
        self.runner = CommandRunner(workdir=self.workdir)
        
    def _execute(self, sample_id: str, mutation_list: Path, bam_file: Path,
                 barcode_file: Path, cellnum: int, threads: int = 2) -> Dict[str, Any]:
        """
        Execute tree input generation
        
        Args:
            sample_id: Sample identifier
            mutation_list: Path to mutation list
            bam_file: Path to BAM file
            barcode_file: Path to barcode file
            cellnum: Number of cells
            threads: Number of threads
            
        Returns:
            Dictionary with step results
        """
        self.logger.info(f"Generating tree input for {sample_id}")
        
        # Check if script exists
        if not self.main_script.exists():
            raise FileNotFoundError(f"Tree input script not found: {self.main_script}")
        
        # Workpath should be parent directory
        workpath = self.workdir.parent.absolute()
        
        # Prepare arguments
        args = [
            str(workpath),
            sample_id,
            str(barcode_file.absolute()),
            str(bam_file.absolute()),
            str(mutation_list.absolute()),
            str(cellnum),
            str(threads)
        ]
        
        # Run script
        self.logger.info(f"Running tree input script with {len(args)} parameters")
        self.runner.run_script(self.main_script, args)
        
        # Collect output files matching real structure
        treeinput_dir = self.workdir
        outputs = {}
        
        # Check if treeinput directory was created
        if not treeinput_dir.exists():
            self.logger.error(f"Tree input directory not created: {treeinput_dir}")
            return {
                'sample_id': sample_id,
                'status': 'failed',
                'error': 'Tree input directory not created'
            }
        
        # Spot matrix file
        spot_files = list(treeinput_dir.glob(f"treeinput_spot_c_{cellnum}.csv"))
        if spot_files:
            self.outputs['spot_matrix'] = spot_files[0]
            outputs['spot_matrix'] = str(spot_files[0])
            self.logger.info(f"Spot matrix file generated: {spot_files[0]}")
        else:
            self.logger.warning(f"Spot matrix file not found for cellnum {cellnum}")
        
        # SCID barcode file
        scid_file = treeinput_dir / "treeinput_scid_barcode.txt"
        if scid_file.exists():
            self.outputs['scid_barcode'] = scid_file
            outputs['scid_barcode'] = str(scid_file)
            self.logger.info(f"SCID barcode file generated: {scid_file}")
        
        # Features file
        features_file = treeinput_dir / "features_file.txt"
        if features_file.exists():
            self.outputs['features_file'] = features_file
            outputs['features_file'] = str(features_file)
            self.logger.info(f"Features file generated: {features_file}")
        
        # Data directory (for subsequent steps)
        data_dir = self.workdir.parent / 'data'
        if data_dir.exists():
            self.outputs['data_dir'] = data_dir
            outputs['data_dir'] = str(data_dir)
            self.logger.info(f"Data directory created: {data_dir}")
        
        return {
            'sample_id': sample_id,
            'outputs': outputs
        }


class SCRNATreeBuildingStep(PipelineStep):
    """Tree building step"""
    
    def __init__(self, workdir: Path, script_dir: Path, config: Dict[str, Any] = None):
        """
        Initialize tree building step
        
        Args:
            workdir: Working directory
            script_dir: Directory containing external scripts
            config: Configuration dictionary
        """
        super().__init__("03_tree_building", workdir, config)
        self.script_dir = Path(script_dir) / "scrna" / "tree_building"
        self.main_script = self.script_dir / "bash_run_other_methods_for_scRNA_part1.sh"
        self.runner = CommandRunner(workdir=self.workdir)
        self.r_runner = RScriptRunner(workdir=self.workdir)
        
    def _prepare_celltype_file(self, sample_id: str, metadata_file: Path, 
                               barcode_file: Path) -> Optional[Path]:
        """
        Generate cell type file from metadata
        
        Args:
            sample_id: Sample identifier
            metadata_file: Path to metadata file
            barcode_file: Path to barcode file
            
        Returns:
            Path to generated cell type file or None if generation failed
        """
        output_file = self.workdir / f"celltype_file_for_{sample_id}.txt"
        
        # Return existing file if already generated
        if output_file.exists():
            self.logger.info(f"Celltype file already exists: {output_file}")
            return output_file
        
        # Path to cell type annotation script
        script_path = Path(__file__).parent.parent.parent / 'scripts' / 'scrna' / 'celltype_anno' / 'celltype_annotation_using_scATOMIC.R'
        
        if not script_path.exists():
            self.logger.warning(f"Celltype annotation script not found: {script_path}")
            self.logger.warning("Please install scATOMIC or provide celltype file directly")
            return None
        
        # Build command
        cmd = [
            'Rscript',
            str(script_path),
            '--metadata', str(metadata_file),
            '--barcode', str(barcode_file),
            '--output', str(output_file),
            '--sample', sample_id
        ]
        
        try:
            self.logger.info(f"Running celltype annotation: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            if result.stdout:
                self.logger.debug(f"Celltype script stdout: {result.stdout}")
            if result.stderr:
                self.logger.debug(f"Celltype script stderr: {result.stderr}")
            
            if output_file.exists():
                self.logger.info(f"Celltype file generated: {output_file}")
                return output_file
            else:
                self.logger.error("Celltype script ran but output file was not created")
                return None
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to generate celltype file: {e}")
            self.logger.error(f"STDERR: {e.stderr}")
            return None
        except Exception as e:
            self.logger.error(f"Unexpected error generating celltype file: {e}")
            return None
    
    def _execute(self, sample_id: str, cellnum: int, 
                 celltype_file: Optional[Path] = None,
                 metadata_file: Optional[Path] = None,
                 barcode_file: Optional[Path] = None,
                 features_file: Optional[Path] = None) -> Dict[str, Any]:
        """
        Execute tree building
        
        Args:
            sample_id: Sample identifier
            cellnum: Number of cells
            celltype_file: Path to pre-computed cell type file (priority 1)
            metadata_file: Path to metadata file for cell type generation (priority 2)
            barcode_file: Path to barcode file (needed for generation)
            features_file: Path to features file from feature extraction step
            
        Returns:
            Dictionary with step results
        """
        self.logger.info(f"Building tree for {sample_id}")
        
        # Determine the final celltype file
        final_celltype_file = None
        
        # Priority 1: User-provided celltype file
        if celltype_file:
            if celltype_file.exists():
                final_celltype_file = celltype_file
                self.logger.info(f"Using user-provided celltype file: {final_celltype_file}")
            else:
                self.logger.error(f"User-provided celltype file not found: {celltype_file}")
                raise FileNotFoundError(f"Celltype file not found: {celltype_file}")
        
        # Priority 2: Generate from metadata
        elif metadata_file and barcode_file:
            self.logger.info(f"Generating celltype file from metadata: {metadata_file}")
            final_celltype_file = self._prepare_celltype_file(sample_id, metadata_file, barcode_file)
            if final_celltype_file:
                self.logger.info(f"Generated celltype file: {final_celltype_file}")
            else:
                self.logger.error("Failed to generate celltype file from metadata")
                # Continue without celltype file (may affect results)
        
        # Priority 3: No celltype information
        else:
            self.logger.warning("No celltype information available (neither file nor metadata)")
            self.logger.warning("Proceeding without celltype file - this may affect results")
        
        # Check features file
        if features_file and features_file.exists():
            self.logger.info(f"Using features file: {features_file}")
        else:
            self.logger.error(f"Features file not found: {features_file}")
            raise FileNotFoundError(f"Features file not found: {features_file}")
        
        # Check if main script exists
        if not self.main_script.exists():
            raise FileNotFoundError(f"Tree building script not found: {self.main_script}")
        
        # Workpath should be parent directory
        workpath = self.workdir.parent.absolute()
        
        # Prepare arguments for bash script
        # celltype_arg should be the file path or "None" as string
        celltype_arg = str(final_celltype_file) if final_celltype_file else "None"
        
        args = [
            str(workpath),
            sample_id,
            str(cellnum),
            celltype_arg
        ]
        
        # Run the main script
        self.logger.info(f"Running tree building with celltype: {celltype_arg}")
        self.logger.info(f"Command: {self.main_script} {' '.join(args)}")
        
        try:
            self.runner.run_script(self.main_script, args)
        except Exception as e:
            self.logger.error(f"Tree building script failed: {e}")
            raise
        
        # Find output files matching real structure
        results_dir = self.workdir
        
        if not results_dir.exists():
            self.logger.error(f"Results directory not created: {results_dir}")
            return {
                'sample_id': sample_id,
                'status': 'failed',
                'error': 'Results directory not created'
            }
        
        # Check for tree file
        tree_file = results_dir / 'PhyloSOLID' / 'celltree.newick'
        if tree_file.exists():
            self.outputs['tree_file'] = tree_file
            self.logger.info(f"Tree file generated: {tree_file}")
        else:
            self.logger.warning(f"Tree file not found at: {tree_file}")
            # Try alternative location
            alt_tree = results_dir / 'PhyloSOLID' / 'celltree.nwk'
            if alt_tree.exists():
                self.outputs['tree_file'] = alt_tree
                self.logger.info(f"Tree file found at alternative location: {alt_tree}")
        
        # Check for CFMatrix
        cfmatrix = results_dir / 'PhyloSOLID' / 'cell_by_mut.CFMatrix'
        if cfmatrix.exists():
            self.outputs['cfmatrix'] = cfmatrix
            self.logger.info(f"CFMatrix file generated: {cfmatrix}")
        
        # Also check scaffold_builder results
        scaffold_tree = results_dir / 'scaffold_builder' / 'phylo_scaffold_tree' / 'final_cleaned_M_scaffold_basedPivots.filtered_sites_inferred.CFMatrix'
        if scaffold_tree.exists():
            self.outputs['scaffold_cfmatrix'] = scaffold_tree
            self.logger.info(f"Scaffold CFMatrix found: {scaffold_tree}")
        
        return {
            'sample_id': sample_id,
            'tree_file': str(tree_file) if tree_file.exists() else None,
            'results_dir': str(results_dir)
        }