"""scDNA pipeline: BAM pileup, unphased genotyping, tree-input matrices, tree building."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..base import Pipeline
from .steps import (
    scDNAFeatureExtractionStep,
    scDNATreeBuildingStep,
    scDNATreeInputStep,
)


class scDNAPipeline(Pipeline):
    """End-to-end scDNA pipeline from single-cell BAMs to a phylogenetic tree."""

    def __init__(self, workdir: Path, script_dir: Path, config: Dict[str, Any] = None):
        super().__init__(workdir, config)
        self.script_dir = script_dir
        self.logger = logging.getLogger(__name__)
        self.add_step("feature_extraction", scDNAFeatureExtractionStep(self.workdir, script_dir, config))
        self.add_step("tree_input", scDNATreeInputStep(self.workdir, script_dir, config))
        self.add_step("tree_building", scDNATreeBuildingStep(self.workdir, script_dir, config))

    def get_step_output(self, step_name: str, output_name: str = None):
        if step_name not in self.steps:
            self.logger.warning(f"Step {step_name} not found")
            return None
        if output_name:
            return self.steps[step_name].get_output(output_name)
        return self.steps[step_name].outputs

    def get_tree_file(self) -> Optional[Path]:
        return self.get_step_output("tree_building", "tree_file")

    def run(
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
        celltype_file: Optional[Path] = None,
        threads: int = 4,
        threshold: float = 0.9,
        steps: List[str] = None,
        **kwargs,
    ) -> Dict[str, Any]:
        self.logger.info(f"Starting scDNA pipeline for sample {sample_id}")
        if not bam_dir and not bam_list:
            raise ValueError("scDNA mode requires --bam-dir or --bam-list of single-cell BAM files")

        steps_to_run = steps or ["feature_extraction", "tree_input", "tree_building"]

        feat_result = {}
        if "feature_extraction" in steps_to_run:
            feat_result = self.run_step(
                "feature_extraction",
                sample_id=sample_id,
                mutation_list=mutation_list,
                bam_dir=bam_dir,
                bam_list=bam_list,
                sample_list=sample_list,
                bulk_bam=bulk_bam,
                keep_bulk=keep_bulk,
                reference=reference,
                mappability=mappability,
                threads=threads,
            )
        else:
            feat_step = self.steps["feature_extraction"]
            meta_file = feat_step.workdir / f"{sample_id}.unphased_reads.meta.txt"
            cell_num = kwargs.get("cellnum", 0)
            if meta_file.exists():
                with open(meta_file) as handle:
                    for line in handle:
                        if line.startswith("cell_num\t"):
                            cell_num = int(line.split("\t", 1)[1].strip())
                            break
            classifier_path = feat_step.workdir / f"{sample_id}.read_level_features_output.bed"
            feat_result = {
                "cell_num": cell_num,
                "tree_input_file": str(feat_step.workdir / f"{sample_id}.tree_input.txt"),
                "scid_file": str(feat_step.workdir / "treeinput_scid_barcode.txt"),
                "classifier_features": str(classifier_path) if classifier_path.exists() and classifier_path.stat().st_size > 0 else None,
            }

        tree_input_file = Path(feat_result["tree_input_file"])
        scid_file = Path(feat_result["scid_file"]) if feat_result.get("scid_file") else None
        cellnum = int(feat_result.get("cell_num") or kwargs.get("cellnum") or 0)

        data_result = {}
        if "tree_input" in steps_to_run:
            data_result = self.run_step(
                "tree_input",
                sample_id=sample_id,
                tree_input_file=tree_input_file,
                scid_file=scid_file,
                cellnum=cellnum,
                threshold=threshold,
            )
        else:
            data_dir = self.steps["tree_input"].workdir / "data"
            data_result = {"data_dir": str(data_dir), "cellnum": cellnum}

        if "tree_building" in steps_to_run:
            features_file = feat_result.get("classifier_features")
            if features_file:
                features_file = Path(features_file)
            self.run_step(
                "tree_building",
                sample_id=sample_id,
                data_dir=Path(data_result["data_dir"]),
                cellnum=cellnum,
                celltype_file=celltype_file,
                features_file=features_file,
            )

        self.results["summary"] = {
            "sample_id": sample_id,
            "workdir": str(self.workdir),
            "steps_completed": [name for name in steps_to_run if name in self.results],
            "tree_file": str(self.get_tree_file()) if self.get_tree_file() else None,
            "data_dir": data_result.get("data_dir"),
        }
        return self.results
