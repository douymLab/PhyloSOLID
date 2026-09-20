#!/usr/bin/env python3
"""
PhyloSOLID command line interface

This module provides the main CLI entry point for running PhyloSOLID pipelines.
"""

import argparse
import logging
import sys
import os
import yaml
import subprocess
from pathlib import Path
from typing import Dict, Optional

# Add project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src import __version__
from pipelines.scrna.pipeline import SCRNAPipeline
from pipelines.scdna.pipeline import scDNAPipeline
from pipelines.spacetracer.pipeline import SpaceTracerPipeline

def setup_logging(verbose: bool = False, log_file: Optional[Path] = None):
    """
    Setup logging configuration
    
    Args:
        verbose: Whether to enable debug logging
        log_file: Path to log file (optional)
    """
    level = logging.DEBUG if verbose else logging.INFO
    
    handlers = [logging.StreamHandler()]
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=handlers
    )

def find_script_dir() -> Optional[Path]:
    """
    Find external scripts directory
    
    Returns:
        Path to scripts directory or None if not found
    """
    # Check environment variable
    if 'PHYLOSOLID_SCRIPTS' in os.environ:
        return Path(os.environ['PHYLOSOLID_SCRIPTS'])
    
    # Check installation directory
    package_dir = Path(__file__).parent.parent
    script_dir = package_dir / 'scripts'
    if script_dir.exists():
        return script_dir
    
    return None

def load_config(config_file: Optional[Path]) -> Dict:
    """
    Load configuration from YAML file
    
    Args:
        config_file: Path to YAML config file
        
    Returns:
        Configuration dictionary
    """
    if config_file and config_file.exists():
        with open(config_file, 'r') as f:
            return yaml.safe_load(f) or {}
    return {}


def _deep_update(base: Dict, extra: Dict) -> Dict:
    for key, value in (extra or {}).items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _deep_update(base[key], value)
        else:
            base[key] = value
    return base


def load_runtime_config(config_file: Optional[Path] = None) -> Dict:
    """Load config/paths.yaml if present, then overlay an explicit --config file."""
    config: Dict = {}
    default_paths = Path(__file__).resolve().parent.parent / "config" / "paths.yaml"
    if default_paths.exists():
        config = load_config(default_paths)
    if config_file:
        config = _deep_update(config, load_config(config_file))
    return config

def load_paths_config(config_file: Optional[Path] = None) -> dict:
    """Load paths configuration from YAML file"""
    default_paths = {
        'conda': {
            'python': None,
            'env_path': None
        },
        'annovar': {
            'script_dir': None,
            'humandb': None,
            'build': 'hg38'
        },
        'reference': {
            'genome_fasta': None,
            'gff3_file': None,
            'mappability_file': None,
            'gnomad_file': None,
            'rna_editing_file': None
        },
        'scdna': {
            'genotyper_dir': None
        }
    }
    
    if config_file and config_file.exists():
        with open(config_file, 'r') as f:
            user_paths = yaml.safe_load(f)
            # Merge with defaults
            if user_paths:
                for section in default_paths:
                    if section in user_paths:
                        default_paths[section].update(user_paths[section])
    
    return default_paths

def validate_annovar_paths(paths_config: dict) -> bool:
    """Validate ANNOVAR paths"""
    annovar_dir = paths_config.get('annovar', {}).get('script_dir')
    humandb = paths_config.get('annovar', {}).get('humandb')
    
    if not annovar_dir or not humandb:
        return False
    
    annovar_script = Path(annovar_dir) / 'annotate_variation.pl'
    humandb_path = Path(humandb)
    
    if not annovar_script.exists():
        logging.error(f"ANNOVAR script not found: {annovar_script}")
        return False
    
    if not humandb_path.exists() or not any(humandb_path.glob('hg38_*.txt')):
        logging.error(f"ANNOVAR humandb not properly set up: {humandb}")
        return False
    
    return True

def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(
        description='PhyloSOLID: Tree building from single-cell sequencing data',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Global options
    parser.add_argument('--version', action='version', version=f'PhyloSOLID {__version__}')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose debug output')
    parser.add_argument('--log-file', type=Path, help='Path to log file')
    parser.add_argument('--script-dir', type=Path, help='Directory containing external scripts')
    parser.add_argument('--workdir', '-w', type=Path, default='./phylosolid_work',
                       help='Working directory (default: ./phylosolid_work)')
    parser.add_argument('--config', '-c', type=Path,
                       help='Configuration YAML (default: config/paths.yaml if present)')
    parser.add_argument('--seed', type=int, default=42, 
                       help='Random seed for reproducibility (default: 42)')    
    
    # Subcommands
    subparsers = parser.add_subparsers(dest='mode', required=True, help='Analysis mode')
    
    # check-annovar subcommand
    check_parser = subparsers.add_parser('check-annovar', 
                                        help='Check ANNOVAR installation and configuration')
    check_parser.add_argument('--config', '-c', type=Path, required=True,
                             help='Configuration file path')
    check_parser.add_argument('--threads', '-t', type=int, default=4,
                             help='Number of threads (default: 4)')
    
    # binary-matrix subcommand
    binary_parser = subparsers.add_parser('binary-matrix', 
                                          help='Build tree directly from binary matrix')
    binary_parser.add_argument('-s', '--sampleid', required=True, help='Sample ID')
    binary_parser.add_argument('-i', '--inputfile', required=True, help='Binary matrix file')
    binary_parser.add_argument('-o', '--outputpath', required=True, help='Output directory')
    
    # spacetracer subcommand
    spacetracer_parser = subparsers.add_parser('spacetracer', 
                                               help='SpaceTracer mode (skip feature extraction)')
    spacetracer_parser.add_argument('--sample', '-s', required=True, help='Sample ID')
    spacetracer_parser.add_argument('--mutation-list', '-m', required=True, type=Path, 
                                   help='Mutation list file')
    spacetracer_parser.add_argument('--bam', '-b', required=True, type=Path, help='BAM file')
    spacetracer_parser.add_argument('--barcode', '-bc', required=True, type=Path, 
                                   help='Barcode file')
    spacetracer_parser.add_argument('--celltype-file', type=Path, 
                                   help='Pre-computed cell type file')
    spacetracer_parser.add_argument('--metadata', type=Path, 
                                   help='Metadata file for cell type extraction')
    spacetracer_parser.add_argument('--read-len', type=int, default=91, help='Read length (default: 91)')
    spacetracer_parser.add_argument('--cellnum', type=int, default=155, help='Number of cells (default: 155)')
    spacetracer_parser.add_argument('--threads', '-t', type=int, default=4, help='Number of threads')
    spacetracer_parser.add_argument('--running-type', default='benchmark',
                                   choices=['benchmark', 'production'], 
                                   help='Running type (default: benchmark)')
    spacetracer_parser.add_argument('--ase-filepath', help='ASE file path (optional)')
    spacetracer_parser.add_argument('--steps', nargs='+',
                                   choices=['tree_input', 'tree_building'],
                                   help='Steps to run (default: both)')
    spacetracer_parser.add_argument('--parallel', action='store_true',
                                   help='Run steps in parallel')
    
    # scRNA subcommand
    scrna_parser = subparsers.add_parser('scrna', help='scRNA-seq mode')
    scrna_parser.add_argument('--sample', '-s', required=True, help='Sample ID')
    scrna_parser.add_argument('--mutation-list', '-m', required=True, type=Path, 
                             help='Mutation list file')
    scrna_parser.add_argument('--bam', '-b', required=True, type=Path, help='BAM file')
    scrna_parser.add_argument('--barcode', '-bc', required=True, type=Path, 
                             help='Barcode file')
    scrna_parser.add_argument('--celltype-file', type=Path, 
                             help='Pre-computed cell type file (if not provided, will be generated from metadata)')
    scrna_parser.add_argument('--metadata', type=Path, 
                             help='Metadata file for cell type extraction (used if celltype-file not provided)')
    scrna_parser.add_argument('--read-len', type=int, default=91, help='Read length (default: 91)')
    scrna_parser.add_argument('--cellnum', type=int, default=155, help='Number of cells (default: 155)')
    scrna_parser.add_argument('--threads', '-t', type=int, default=4, help='Number of threads (default: 4)')
    scrna_parser.add_argument('--running-type', default='benchmark',
                             choices=['benchmark', 'production'], 
                             help='Running type (default: benchmark)')
    scrna_parser.add_argument('--ase-filepath', help='ASE file path (optional)')
    scrna_parser.add_argument('--steps', nargs='+',
                             choices=['feature_extraction', 'tree_input', 'tree_building'],
                             help='Steps to run (default: all)')
    scrna_parser.add_argument('--parallel', action='store_true',
                             help='Run feature extraction and tree input in parallel')
    
    # scDNA subcommand: single-cell BAMs + mutation list; bulk BAM is optional
    scdna_parser = subparsers.add_parser(
        'scdna',
        help='scDNA-seq mode: BAM pileup to tree (bulk BAM optional)',
        description=(
            'End-to-end scDNA pipeline: pileup mutation sites from single-cell BAMs, '
            'unphased genotyping, tree-input matrices, then tree building. '
            'Bulk BAM is optional. Requires scdna.genotyper_dir in config/paths.yaml '
            'or PHYLOSOLID_GENOTYPER_DIR.'
        )
    )
    scdna_parser.add_argument('--sample', '-s', required=True, help='Sample ID')
    scdna_parser.add_argument('--mutation-list', '-m', required=True, type=Path,
                             help='Mutation list (chr_pos_ref_alt, one per line)')
    bam_group = scdna_parser.add_mutually_exclusive_group(required=True)
    bam_group.add_argument('--bam-dir', type=Path,
                           help='Directory of single-cell BAM files named {cell}.bam')
    bam_group.add_argument('--bam-list', type=Path,
                           help='Text file with one single-cell BAM path per line')
    scdna_parser.add_argument('--sample-list', type=Path,
                             help='Optional cell IDs (one per line). If omitted, BAM file stems are used')
    scdna_parser.add_argument('--bulk-bam', type=Path,
                             help='Optional bulk BAM. Not required; omitted by default')
    scdna_parser.add_argument('--keep-bulk', action='store_true',
                             help='Keep a BAM named bulk.bam if it appears in --bam-dir/--bam-list')
    scdna_parser.add_argument('--reference', type=Path,
                             help='Reference FASTA for optional classifier BAM features (or scdna.genome_fasta in config)')
    scdna_parser.add_argument('--mappability', type=Path,
                             help='Mappability bigWig for optional classifier BAM features (or scdna.mappability_file in config)')
    scdna_parser.add_argument('--celltype-file', type=Path, help='Optional cell-type annotation file')
    scdna_parser.add_argument('--threshold', type=float, default=0.9,
                             help='Somatic posterior threshold for tree-input filtering (default: 0.9)')
    scdna_parser.add_argument('--threads', '-t', type=int, default=4, help='Number of threads')
    scdna_parser.add_argument('--steps', nargs='+',
                             choices=['feature_extraction', 'tree_input', 'tree_building'],
                             help='Steps to run (default: all). Use feature_extraction tree_input to stop at tree-input files')
    
    args = parser.parse_args()
    
    # Handle binary-matrix mode
    if args.mode == 'binary-matrix':
        binput_script = Path(__file__).parent.parent / 'src' / 'run_phylosilid_fullTree_binput.py'
        if not binput_script.exists():
            print(f"Error: Binary matrix script not found: {binput_script}")
            sys.exit(1)
        
        cmd = [sys.executable, str(binput_script), 
               '-s', args.sampleid, 
               '-i', args.inputfile, 
               '-o', args.outputpath]
        env = os.environ.copy()
        project_root = Path(__file__).parent.parent
        env['PYTHONPATH'] = str(project_root) + os.pathsep + env.get('PYTHONPATH', '')
        
        try:
            subprocess.run(cmd, check=True, cwd=str(project_root), env=env)
        except subprocess.CalledProcessError as e:
            sys.exit(e.returncode)
        return
    
    # Setup logging
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger('phylosolid')
    
    # Display version
    logger.info(f"PhyloSOLID version {__version__}")
    logger.info(f"Mode: {args.mode}")
    
    # Find scripts directory
    script_dir = args.script_dir or find_script_dir()
    if script_dir:
        logger.info(f"Using scripts from: {script_dir}")
    else:
        logger.warning("Script directory not found - external scripts may not be accessible")
    
    # Load configuration (config/paths.yaml by default, overlay --config if given)
    config = load_runtime_config(args.config)
    
    # Safely update configuration based on available attributes
    if hasattr(args, 'threads'):
        config['threads'] = args.threads
    if hasattr(args, 'verbose'):
        config['verbose'] = args.verbose
    elif args.verbose:
        config['verbose'] = True
    
    # Create working directory only for pipeline modes
    if args.mode not in ['check-annovar']:
        workdir = args.workdir / args.sample
        workdir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Working directory: {workdir}")
    
    # Check input files only for pipeline modes
    if args.mode not in ['check-annovar']:
        file_args = ['mutation_list']
        if args.mode == 'scdna':
            file_args.extend(['bam_dir', 'bam_list', 'sample_list', 'bulk_bam', 'reference', 'mappability', 'celltype_file'])
        else:
            file_args.extend(['bam', 'barcode'])
        for file_arg in file_args:
            if hasattr(args, file_arg):
                file_path = getattr(args, file_arg)
                if file_path and not Path(file_path).exists():
                    logger.error(f"Input file not found: {file_path}")
                    sys.exit(1)
        if args.mode == 'scdna' and not getattr(args, 'bam_dir', None) and not getattr(args, 'bam_list', None):
            logger.error("scDNA mode requires --bam-dir or --bam-list")
            sys.exit(1)
    
    try:
        # Handle check-annovar mode
        if args.mode == 'check-annovar':
            # Load paths configuration
            paths_config = load_paths_config(args.config)
            
            print("\n" + "="*60)
            print("PhyloSOLID ANNOVAR Configuration Check")
            print("="*60)
            
            annovar_config = paths_config.get('annovar', {})
            annovar_dir = annovar_config.get('script_dir')
            humandb = annovar_config.get('humandb')
            build = annovar_config.get('build', 'hg38')
            
            print(f"Genome build: {build}")
            print("-" * 40)
            
            all_good = True
            
            # Check ANNOVAR script directory
            if annovar_dir:
                annovar_path = Path(annovar_dir)
                annotate_script = annovar_path / 'annotate_variation.pl'
                table_script = annovar_path / 'table_annovar.pl'
                convert_script = annovar_path / 'convert2annovar.pl'
                
                print(f"ANNOVAR directory: {annovar_path}")
                
                if annotate_script.exists():
                    print(f"  ✅ annotate_variation.pl found")
                else:
                    print(f"  ❌ annotate_variation.pl not found")
                    all_good = False
                
                if table_script.exists():
                    print(f"  ✅ table_annovar.pl found")
                else:
                    print(f"  ⚠️  table_annovar.pl not found (optional)")
                
                if convert_script.exists():
                    print(f"  ✅ convert2annovar.pl found")
                else:
                    print(f"  ⚠️  convert2annovar.pl not found (optional)")
            else:
                print("❌ ANNOVAR script_dir not configured in paths.yaml")
                all_good = False
            
            print("-" * 40)
            
            # Check humandb directory
            if humandb:
                humandb_path = Path(humandb)
                print(f"humandb directory: {humandb_path}")
                
                if humandb_path.exists():
                    print(f"  ✅ Directory exists")
                    
                    # Check for required database files
                    required_dbs = [
                        f'{build}_refGene.txt',
                        f'{build}_refGeneMrna.fa',
                        f'{build}_refGene.txt.gz',
                        f'{build}_refGeneMrna.fa.gz'
                    ]
                    
                    found_dbs = []
                    for db in required_dbs:
                        db_file = humandb_path / db
                        if db_file.exists():
                            found_dbs.append(db)
                    
                    if found_dbs:
                        print(f"  ✅ Found database files: {', '.join(found_dbs[:3])}")
                        if len(found_dbs) > 3:
                            print(f"     and {len(found_dbs)-3} more")
                    else:
                        print(f"  ❌ No {build} database files found")
                        print(f"     Please run: ./annotate_variation.pl -downdb -buildver {build} -webfrom annovar refGene {humandb}")
                        all_good = False
                else:
                    print(f"  ❌ Directory does not exist")
                    all_good = False
            else:
                print("❌ humandb not configured in paths.yaml")
                all_good = False
            
            print("="*60)
            
            if all_good:
                print("\n✅ ANNOVAR check passed! Your configuration looks good.")
                print("   You can now run PhyloSOLID pipelines.")
                sys.exit(0)
            else:
                print("\n❌ ANNOVAR check failed. Please fix the issues above.")
                print("   See documentation for ANNOVAR installation instructions.")
                sys.exit(1)
        
        # Handle spacetracer mode
        elif args.mode == 'spacetracer':
            pipeline = SpaceTracerPipeline(workdir, script_dir, config)
            
            if args.steps is None:
                steps = ['tree_input', 'tree_building']
            else:
                steps = args.steps
            
            results = pipeline.run(
                sample_id=args.sample,
                mutation_list=args.mutation_list,
                bam_file=args.bam,
                barcode_file=args.barcode,
                celltype_file=getattr(args, 'celltype_file', None),
                metadata_file=getattr(args, 'metadata', None),
                read_len=args.read_len,
                cellnum=args.cellnum,
                threads=args.threads,
                running_type=args.running_type,
                ase_filepath=getattr(args, 'ase_filepath', None),
                steps=steps,
                parallel=args.parallel
            )
            
            logger.info("=" * 50)
            logger.info("SpaceTracer pipeline completed successfully!")
            logger.info(f"Steps completed: {', '.join(results.keys())}")
            
            tree_file = pipeline.get_tree_file()
            if tree_file:
                logger.info(f"Tree file: {tree_file}")
            
            summary_file = workdir / 'pipeline_summary.yaml'
            with open(summary_file, 'w') as f:
                yaml.dump(results, f, default_flow_style=False)
            logger.info(f"Summary saved to: {summary_file}")
        
        # Run scRNA pipeline
        elif args.mode == 'scrna':
            pipeline = SCRNAPipeline(workdir, script_dir, config)
            
            results = pipeline.run(
                sample_id=args.sample,
                mutation_list=args.mutation_list,
                bam_file=args.bam,
                barcode_file=args.barcode,
                celltype_file=getattr(args, 'celltype_file', None),
                metadata_file=getattr(args, 'metadata', None),
                read_len=getattr(args, 'read_len', 91),
                cellnum=getattr(args, 'cellnum', 155),
                threads=args.threads,
                running_type=getattr(args, 'running_type', 'benchmark'),
                ase_filepath=getattr(args, 'ase_filepath', None),
                steps=getattr(args, 'steps', None),
                parallel=getattr(args, 'parallel', False)
            )
            
            # Display results
            logger.info("=" * 50)
            logger.info("Pipeline completed successfully!")
            logger.info(f"Steps completed: {', '.join(results.keys())}")
            
            tree_file = pipeline.get_tree_file()
            if tree_file:
                logger.info(f"Tree file: {tree_file}")
            
            # Save summary
            summary_file = workdir / 'pipeline_summary.yaml'
            with open(summary_file, 'w') as f:
                yaml.dump(results, f, default_flow_style=False)
            logger.info(f"Summary saved to: {summary_file}")
        
        # Run scDNA pipeline
        elif args.mode == 'scdna':
            pipeline = scDNAPipeline(workdir, script_dir, config)
            
            results = pipeline.run(
                sample_id=args.sample,
                mutation_list=args.mutation_list,
                bam_dir=getattr(args, 'bam_dir', None),
                bam_list=getattr(args, 'bam_list', None),
                sample_list=getattr(args, 'sample_list', None),
                bulk_bam=getattr(args, 'bulk_bam', None),
                keep_bulk=getattr(args, 'keep_bulk', False),
                reference=getattr(args, 'reference', None),
                mappability=getattr(args, 'mappability', None),
                celltype_file=getattr(args, 'celltype_file', None),
                threads=args.threads,
                threshold=getattr(args, 'threshold', 0.9),
                steps=getattr(args, 'steps', None)
            )
            
            logger.info("=" * 50)
            logger.info("scDNA pipeline completed")
            logger.info(f"Steps completed: {', '.join(results.keys())}")
            data_dir = results.get('summary', {}).get('data_dir')
            if data_dir:
                logger.info(f"Tree-input matrices: {data_dir}")
            tree_file = pipeline.get_tree_file()
            if tree_file:
                logger.info(f"Tree file: {tree_file}")
            summary_file = workdir / 'pipeline_summary.yaml'
            with open(summary_file, 'w') as f:
                yaml.dump(results, f, default_flow_style=False)
            logger.info(f"Summary saved to: {summary_file}")
        
        else:
            logger.error(f"Mode '{args.mode}' not implemented")
            sys.exit(1)
            
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()