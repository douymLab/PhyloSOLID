"""scDNA pipeline steps: BAM pileup, unphased genotyping, tree-input matrices, tree building."""

from __future__ import annotations

import logging
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, Optional

from ..base import PipelineStep
from utils.command import CommandRunner, RScriptRunner


def _package_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _config_path(config: Dict[str, Any], *keys: str, default: Any = None) -> Any:
    current: Any = config or {}
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


class scDNAFeatureExtractionStep(PipelineStep):
    """Pileup mutation sites from scBAMs, run unphased genotyper, optionally extract classifier features."""

    def __init__(self, workdir: Path, script_dir: Path, config: Dict[str, Any] = None):
        super().__init__("01_features", workdir, config)
        self.script_dir = Path(script_dir) / "scdna"
        self.extract_script = self.script_dir / "feature_extraction" / "extract_unphased_reads_from_bams.py"
        self.format_script = self.script_dir / "tree_input" / "format_genotyper_to_treeinput.py"
        self.runner = CommandRunner(workdir=self.workdir)

    def _genotyper_dir(self) -> Path:
        configured = (
            os.environ.get("PHYLOSOLID_GENOTYPER_DIR")
            or _config_path(self.config, "scdna", "genotyper_dir")
            or _config_path(self.config, "genotyper_dir")
        )
        if not configured:
            raise FileNotFoundError(
                "scDNA unphased genotyper directory is not set. "
                "Add scdna.genotyper_dir to config/paths.yaml or export PHYLOSOLID_GENOTYPER_DIR."
            )
        return Path(configured)

    def _execute(
        self,
        sample_id: str,
        mutation_list: Path,
        bam_dir: Optional[Path] = None,
        bam_list: Optional[Path] = None,
        sample_list: Optional[Path] = None,
        bulk_bam: Optional[Path] = None,
        keep_bulk: bool = False,
        reference: Optional[Path] = None,
        mappability: Optional[Path] = None,
        threads: int = 4,
        **kwargs,
    ) -> Dict[str, Any]:
        self.logger.info(f"Extracting scDNA pileups and genotypes for {sample_id}")
        if not self.extract_script.exists():
            raise FileNotFoundError(f"Pileup script not found: {self.extract_script}")

        reads_bed = self.workdir / f"{sample_id}.unphased_reads.bed"
        sample_name_file = self.workdir / "sample_name.txt"
        scid_file = self.workdir / "treeinput_scid_barcode.txt"
        depth_file = self.workdir / f"{sample_id}.average_depth.txt"
        staged_bam_dir = self.workdir / "bam_links"

        extract_cmd = [
            sys.executable,
            str(self.extract_script),
            "--mutation-list",
            str(Path(mutation_list).absolute()),
            "--output",
            str(reads_bed),
            "--sample-name-file",
            str(sample_name_file),
            "--scid-file",
            str(scid_file),
            "--average-depth-file",
            str(depth_file),
            "--stage-bam-dir",
            str(staged_bam_dir),
        ]
        if bam_dir:
            extract_cmd.extend(["--bam-dir", str(Path(bam_dir).absolute())])
        if bam_list:
            extract_cmd.extend(["--bam-list", str(Path(bam_list).absolute())])
        if sample_list:
            extract_cmd.extend(["--sample-list", str(Path(sample_list).absolute())])
        if bulk_bam:
            extract_cmd.extend(["--bulk-bam", str(Path(bulk_bam).absolute())])
        if keep_bulk:
            extract_cmd.append("--keep-bulk")

        self.logger.info("Running BAM pileup at mutation sites")
        self.runner.run(extract_cmd)

        meta_file = reads_bed.with_suffix(".meta.txt")
        cell_num = 0
        with_bulk = "yes" if bulk_bam or keep_bulk else "no"
        if meta_file.exists():
            meta = {}
            with open(meta_file) as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    key, value = line.rstrip().split("\t", 1)
                    meta[key] = value
            cell_num = int(meta.get("cell_num", 0))
            with_bulk = meta.get("with_bulk", with_bulk)
        if cell_num <= 0:
            raise RuntimeError("Failed to infer single-cell count from pileup metadata")

        genotyper_dir = self._genotyper_dir()
        genotyper_script = genotyper_dir / "run_unphased_prod_whether_bulk.py"
        if not genotyper_script.exists():
            raise FileNotFoundError(
                f"Unphased genotyper not found: {genotyper_script}. "
                "Set scdna.genotyper_dir in the config YAML."
            )

        posterior_file = self.workdir / f"{sample_id}.unphased_posterior.txt"
        env = os.environ.copy()
        env["PYTHONPATH"] = str(genotyper_dir) + os.pathsep + env.get("PYTHONPATH", "")
        genotyper_runner = CommandRunner(workdir=self.workdir, env=env)
        genotyper_cmd = [
            sys.executable,
            str(genotyper_script),
            "-n",
            str(cell_num),
            "-c",
            "no" if with_bulk == "no" else "yes",
            "-b",
            with_bulk,
            "-i",
            str(reads_bed),
            "-o",
            str(posterior_file),
            "-j",
            str(max(1, threads)),
            "-s",
            "1",
        ]
        self.logger.info("Running unphased genotyper (with_bulk=%s, cell_num=%s)", with_bulk, cell_num)
        genotyper_runner.run(genotyper_cmd)
        if not posterior_file.exists() or posterior_file.stat().st_size == 0:
            raise RuntimeError(f"Genotyper produced no output: {posterior_file}")

        tree_input_file = self.workdir / f"{sample_id}.tree_input.txt"
        features_input = self.workdir / f"{sample_id}.read_level_features_input.txt"
        format_cmd = [
            sys.executable,
            str(self.format_script),
            "--genotyper-output",
            str(posterior_file),
            "--reads-bed",
            str(reads_bed),
            "--output",
            str(tree_input_file),
            "--features-input",
            str(features_input),
            "--with-bulk",
            with_bulk,
        ]
        self.runner.run(format_cmd)

        features_output = self.workdir / f"{sample_id}.read_level_features_output.bed"
        classifier_features = None
        ref_value = reference or _config_path(self.config, "scdna", "genome_fasta")
        map_value = mappability or _config_path(self.config, "scdna", "mappability_file")
        reference_path = Path(ref_value) if ref_value else None
        mappability_path = Path(map_value) if map_value else None
        feature_script = genotyper_dir / "read_level_features_extraction_whether_bulk.py"
        can_extract = (
            feature_script.exists()
            and reference_path is not None
            and reference_path.exists()
            and mappability_path is not None
            and mappability_path.exists()
            and shutil.which("bigWigAverageOverBed")
            and features_input.exists()
            and features_input.stat().st_size > 0
        )
        if can_extract:
            self.logger.info("Extracting classifier read-level features without requiring bulk")
            feature_cmd = [
                sys.executable,
                str(feature_script),
                "-i",
                str(features_input),
                "-o",
                str(features_output),
                "-s",
                str(sample_name_file),
                "-wb",
                with_bulk,
                "-b",
                str(staged_bam_dir),
                "-r",
                str(reference_path),
                "-u",
                str(mappability_path),
                "-d",
                str(depth_file),
                "-n",
                str(max(1, min(threads, 4))),
            ]
            try:
                genotyper_runner.run(feature_cmd)
                if features_output.exists() and features_output.stat().st_size > 0:
                    classifier_features = features_output
                    self.outputs["classifier_features"] = classifier_features
            except Exception as exc:
                self.logger.warning("Classifier feature extraction failed and will be skipped: %s", exc)
        else:
            self.logger.warning(
                "Skipping classifier BAM feature extraction (need genotyper script, "
                "scdna.genome_fasta, scdna.mappability_file, and bigWigAverageOverBed). "
                "Tree-input matrices will still be built."
            )

        self.outputs["reads_bed"] = reads_bed
        self.outputs["posterior_file"] = posterior_file
        self.outputs["tree_input_file"] = tree_input_file
        self.outputs["scid_file"] = scid_file
        self.outputs["sample_name_file"] = sample_name_file
        self.outputs["depth_file"] = depth_file
        self.outputs["staged_bam_dir"] = staged_bam_dir

        return {
            "sample_id": sample_id,
            "cell_num": cell_num,
            "with_bulk": with_bulk,
            "reads_bed": str(reads_bed),
            "posterior_file": str(posterior_file),
            "tree_input_file": str(tree_input_file),
            "scid_file": str(scid_file),
            "classifier_features": str(classifier_features) if classifier_features else None,
        }


class scDNATreeInputStep(PipelineStep):
    """Convert genotyper tree_input into PhyloSOLID data matrices."""

    def __init__(self, workdir: Path, script_dir: Path, config: Dict[str, Any] = None):
        super().__init__("02_treeinput", workdir, config)
        self.script_dir = Path(script_dir) / "scdna" / "tree_input"
        self.main_script = self.script_dir / "generate_treeinput_data.R"
        self.r_runner = RScriptRunner(workdir=self.workdir)

    def _execute(
        self,
        sample_id: str,
        tree_input_file: Path,
        scid_file: Optional[Path] = None,
        cellnum: int = 0,
        threshold: float = 0.9,
        **kwargs,
    ) -> Dict[str, Any]:
        self.logger.info(f"Generating scDNA tree-input matrices for {sample_id}")
        if not self.main_script.exists():
            raise FileNotFoundError(f"Tree input R script not found: {self.main_script}")
        tree_input_file = Path(tree_input_file)
        if not tree_input_file.exists():
            raise FileNotFoundError(f"tree_input file not found: {tree_input_file}")

        data_dir = self.workdir / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        staged_input = self.workdir / f"{sample_id}.tree_input.txt"
        if tree_input_file.resolve() != staged_input.resolve():
            shutil.copy2(tree_input_file, staged_input)

        args = [
            "--inputfile",
            str(tree_input_file.absolute()),
            "--cellnum",
            str(cellnum),
            "--outputpath",
            str(data_dir.absolute()),
            "--is_remove_cells",
            "no",
            "--threshold",
            str(threshold),
            "--indid",
            sample_id,
        ]
        if scid_file and Path(scid_file).exists():
            staged_scid = self.workdir / "treeinput_scid_barcode.txt"
            if Path(scid_file).resolve() != staged_scid.resolve():
                shutil.copy2(scid_file, staged_scid)
            args.extend(["--scid_file", str(Path(scid_file).absolute())])

        self.logger.info("Running generate_treeinput_data.R")
        self.r_runner.run(self.main_script, args)

        required = [
            "data.posterior_matrix.txt",
            "data.allele_count.txt",
            "features.preprocess_items.txt",
            "data.likelihood_mut_matrix.txt",
            "data.likelihood_unmut_matrix.txt",
        ]
        missing = [name for name in required if not (data_dir / name).exists()]
        if missing:
            raise FileNotFoundError(
                f"Tree-input matrices missing in {data_dir}: {', '.join(missing)}"
            )

        self.outputs["data_dir"] = data_dir
        for name in required:
            self.outputs[name] = data_dir / name
        self.outputs["tree_input_file"] = staged_input

        return {
            "sample_id": sample_id,
            "cellnum": cellnum,
            "data_dir": str(data_dir),
            "tree_input_file": str(staged_input),
        }


class scDNATreeBuildingStep(PipelineStep):
    """Run the scDNA PhyloSOLID tree builder on prepared data matrices."""

    def __init__(self, workdir: Path, script_dir: Path, config: Dict[str, Any] = None):
        super().__init__("03_tree_building", workdir, config)
        self.script_dir = Path(script_dir) / "scdna" / "tree_building"
        self.main_script = _package_root() / "src" / "run_phylosilid_fullTree_scDNA.py"
        env = os.environ.copy()
        env["PYTHONPATH"] = str(_package_root()) + os.pathsep + env.get("PYTHONPATH", "")
        self.runner = CommandRunner(workdir=_package_root(), env=env)
        self.r_runner = RScriptRunner(workdir=self.workdir)

    @staticmethod
    def _usable_features_file(features_file: Optional[Path]) -> Optional[Path]:
        if not features_file:
            return None
        path = Path(features_file)
        if not path.exists() or path.stat().st_size == 0:
            return None
        with open(path) as handle:
            header = handle.readline()
            first_row = handle.readline()
        if not first_row.strip():
            return None
        columns = header.strip().split("\t")
        required = {"chrom", "end", "ref_allele", "alt_allele"}
        if not required.issubset(set(columns)):
            return None
        return path

    def _execute(
        self,
        sample_id: str,
        data_dir: Path,
        cellnum: int,
        celltype_file: Optional[Path] = None,
        features_file: Optional[Path] = None,
        **kwargs,
    ) -> Dict[str, Any]:
        self.logger.info(f"Building scDNA tree for {sample_id}")
        if not self.main_script.exists():
            raise FileNotFoundError(f"scDNA tree building script not found: {self.main_script}")
        data_dir = Path(data_dir)
        if not data_dir.exists():
            raise FileNotFoundError(f"Prepared data directory not found: {data_dir}")

        features_file = self._usable_features_file(features_file)
        if features_file is None:
            self.logger.info("No usable classifier features file; tree building will keep all mutations as candidates")

        results_dir = self.workdir
        cmd = [
            sys.executable,
            str(self.main_script),
            "--sampleid",
            sample_id,
            "--inputpath",
            str(data_dir.absolute()),
            "--outputpath",
            str(results_dir.absolute()),
            "--celltype_file",
            str(celltype_file) if celltype_file else "None",
            "--features_file",
            str(features_file) if features_file else "None",
        ]
        self.logger.info("Running %s", " ".join(cmd))
        self.runner.run(cmd)

        phylo_dir = results_dir / "05_final_results" / "phylo"
        self.outputs["results_dir"] = results_dir
        if phylo_dir.exists():
            self.outputs["phylo_dir"] = phylo_dir
        cfmatrix_candidates = [
            phylo_dir / "final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix",
            phylo_dir / "final_cleaned_M_scaffold_basedPivots.filtered_sites_inferred.CFMatrix",
            results_dir / "03_scaffold_builder" / "phylo_scaffold_tree" / "final_cleaned_M_scaffold_basedPivots.filtered_sites_inferred.CFMatrix",
        ]
        cfmatrix = next((path for path in cfmatrix_candidates if path.exists()), None)
        convert_script = _package_root() / "scripts" / "scrna" / "tree_building" / "convert_PhyloSOLID_tree.R"
        tree_file = results_dir / "PhyloSOLID" / "celltree.newick"
        if cfmatrix:
            self.outputs["cfmatrix"] = cfmatrix
            tree_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(cfmatrix, tree_file.parent / "cell_by_mut.CFMatrix")
            if convert_script.exists():
                try:
                    self.r_runner.run(convert_script, [str(cfmatrix), str(tree_file)])
                    if tree_file.exists():
                        self.outputs["tree_file"] = tree_file
                except Exception as exc:
                    self.logger.warning("Newick conversion skipped: %s", exc)
            else:
                self.logger.warning("convert_PhyloSOLID_tree.R not found; CFMatrix is available at %s", cfmatrix)

        return {
            "sample_id": sample_id,
            "cellnum": cellnum,
            "results_dir": str(results_dir),
            "tree_file": str(self.outputs.get("tree_file")) if "tree_file" in self.outputs else None,
        }
