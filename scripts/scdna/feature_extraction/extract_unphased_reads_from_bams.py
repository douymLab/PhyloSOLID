#!/usr/bin/env python3
"""Extract unphased allele/base-quality counts from single-cell BAMs.

Writes the 7-column unphased_reads.bed consumed by scMosaicCaller's
run_unphased_prod_whether_bulk.py. Bulk BAM is optional and is not required.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_SCRIPTS = SCRIPT_DIR.parent
if str(PACKAGE_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(PACKAGE_SCRIPTS))

from preprocess_utils import (  # noqa: E402
    BASES,
    collect_sc_bams,
    default_allele_pool,
    default_popaf_dict,
    load_mutation_list,
    stage_bam_dir,
    write_average_depth_file,
    write_sample_name_file,
    write_scid_file,
)

logger = logging.getLogger("extract_unphased_reads")


def resolve_contig(bam: pysam.AlignmentFile, chrom: str) -> Optional[str]:
    names = set(bam.references)
    candidates = [chrom]
    if chrom.startswith("chr"):
        candidates.append(chrom[3:])
    else:
        candidates.append("chr" + chrom)
    if chrom in {"M", "chrM", "MT", "chrMT"}:
        candidates.extend(["MT", "chrM", "chrMT", "M"])
    for name in candidates:
        if name in names:
            return name
    return None


def pileup_site(bam: pysam.AlignmentFile, contig: str, pos_1based: int) -> Dict[str, object]:
    counts: Counter = Counter()
    start0 = pos_1based - 1
    for column in bam.pileup(
        contig,
        start0,
        pos_1based,
        truncate=True,
        stepper="nofilter",
        min_base_quality=0,
        ignore_overlaps=False,
        ignore_orphans=False,
    ):
        if column.reference_pos != start0:
            continue
        for pileup_read in column.pileups:
            if pileup_read.is_del or pileup_read.is_refskip:
                continue
            aln = pileup_read.alignment
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            qpos = pileup_read.query_position
            if qpos is None or aln.query_sequence is None or aln.query_qualities is None:
                continue
            base = aln.query_sequence[qpos].upper()
            if base not in BASES:
                continue
            key = f"{base}_{int(aln.query_qualities[qpos])}"
            counts[key] += 1
    if not counts:
        return {"no_reads": "*"}
    return dict(counts)


def observed_bases(reads_dict: Dict[str, object]) -> List[str]:
    bases: List[str] = []
    for key in reads_dict:
        if key == "no_reads":
            continue
        base = str(key).split("_", 1)[0]
        if base in BASES and base not in bases:
            bases.append(base)
    return bases


def coverage_from_reads(reads_dict: Dict[str, object]) -> int:
    if "no_reads" in reads_dict:
        return 0
    return int(sum(int(v) for v in reads_dict.values()))


def extract_reads(
    mutations: Sequence[Tuple[str, int, str, str]],
    sc_bams: Sequence[Tuple[str, Path]],
    bulk_bam: Optional[Path] = None,
    min_alt_reads: int = 0,
) -> Tuple[List[str], Dict[str, float]]:
    sc_handles = []
    contig_maps: List[Dict[str, str]] = []
    bulk_handle = None
    try:
        for name, path in sc_bams:
            bam = pysam.AlignmentFile(str(path), "rb")
            sc_handles.append(bam)
            contig_maps.append({})
            logger.info("Opened BAM %s (%s)", name, path)

        bulk_contigs: Dict[str, str] = {}
        if bulk_bam:
            bulk_handle = pysam.AlignmentFile(str(bulk_bam), "rb")
            logger.info("Opened optional bulk BAM %s", bulk_bam)

        depths = {name: [] for name, _ in sc_bams}
        if bulk_bam:
            depths["bulk"] = []

        lines: List[str] = []
        for chrom, pos, ref, alt in mutations:
            cell_reads: List[Dict[str, object]] = []
            seen = set()
            for bam, contig_map, (name, _) in zip(sc_handles, contig_maps, sc_bams):
                contig = contig_map.get(chrom)
                if contig is None:
                    contig = resolve_contig(bam, chrom)
                    if contig is None:
                        logger.warning("Contig %s not found in BAM %s", chrom, name)
                        contig_map[chrom] = ""
                        cell_reads.append({"no_reads": "*"})
                        depths[name].append(0)
                        continue
                    contig_map[chrom] = contig
                if contig_map[chrom] == "":
                    cell_reads.append({"no_reads": "*"})
                    depths[name].append(0)
                    continue
                reads = pileup_site(bam, contig, pos)
                cell_reads.append(reads)
                cov = coverage_from_reads(reads)
                depths[name].append(cov)
                seen.update(observed_bases(reads))

            bulk_reads: Optional[Dict[str, object]] = None
            if bulk_handle is not None:
                contig = bulk_contigs.get(chrom)
                if contig is None:
                    contig = resolve_contig(bulk_handle, chrom) or ""
                    bulk_contigs[chrom] = contig
                if contig:
                    bulk_reads = pileup_site(bulk_handle, contig, pos)
                else:
                    bulk_reads = {"no_reads": "*"}
                depths["bulk"].append(coverage_from_reads(bulk_reads))
                seen.update(observed_bases(bulk_reads))

            if min_alt_reads > 0:
                alt_total = 0
                for reads in ([bulk_reads] if bulk_reads else []) + cell_reads:
                    for key, value in reads.items():
                        if str(key).split("_", 1)[0] == alt:
                            alt_total += int(value)
                if alt_total < min_alt_reads:
                    continue

            allele_pool = default_allele_pool(ref, alt, seen)
            popaf = default_popaf_dict(ref, alt)
            all_reads = ([bulk_reads] if bulk_reads is not None else []) + cell_reads
            start0 = pos - 1
            line = "\t".join(
                [
                    str(chrom),
                    str(start0),
                    str(pos),
                    ref,
                    str(allele_pool),
                    str(all_reads),
                    str(popaf),
                ]
            )
            lines.append(line)

        avg_depth = {}
        for sample, values in depths.items():
            avg_depth[sample] = (sum(values) / len(values)) if values else 0.0
        return lines, avg_depth
    finally:
        for bam in sc_handles:
            bam.close()
        if bulk_handle is not None:
            bulk_handle.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract unphased reads at mutation sites from single-cell BAMs (bulk optional)."
    )
    parser.add_argument("-m", "--mutation-list", required=True, type=Path, help="Mutation list (chr_pos_ref_alt)")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output unphased_reads.bed")
    parser.add_argument("--bam-dir", type=Path, help="Directory of {cell}.bam files")
    parser.add_argument("--bam-list", type=Path, help="Text file with one BAM path per line")
    parser.add_argument("--sample-list", type=Path, help="Optional cell IDs, one per line, matching BAM stems")
    parser.add_argument("--bulk-bam", type=Path, help="Optional bulk BAM; omitted by default")
    parser.add_argument("--keep-bulk", action="store_true", help="Keep a BAM named bulk.bam if present in bam-dir/list")
    parser.add_argument("--stage-bam-dir", type=Path, help="Write named BAM symlinks for feature extraction")
    parser.add_argument("--sample-name-file", type=Path, help="Output sample name list")
    parser.add_argument("--scid-file", type=Path, help="Output scid file for generate_treeinput_data.R")
    parser.add_argument("--average-depth-file", type=Path, help="Output per-cell average depth table")
    parser.add_argument("--min-alt-reads", type=int, default=0, help="Drop sites with fewer alt reads across cells")
    return parser.parse_args()


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    args = parse_args()
    mutations = load_mutation_list(args.mutation_list)
    sc_bams, detected_bulk = collect_sc_bams(
        bam_dir=args.bam_dir,
        bam_list=args.bam_list,
        sample_list=args.sample_list,
        bulk_bam=args.bulk_bam,
        keep_bulk=args.keep_bulk,
    )
    bulk_bam = args.bulk_bam or (detected_bulk if args.keep_bulk else None)
    lines, avg_depth = extract_reads(mutations, sc_bams, bulk_bam=bulk_bam, min_alt_reads=args.min_alt_reads)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        handle.write("\n".join(lines) + ("\n" if lines else ""))
    logger.info("Wrote %d sites to %s", len(lines), args.output)

    sample_names = [name for name, _ in sc_bams]
    if args.sample_name_file:
        write_sample_name_file(args.sample_name_file, sample_names, bulk_name="bulk" if bulk_bam else None)
    if args.scid_file:
        write_scid_file(args.scid_file, sample_names)
    if args.average_depth_file:
        write_average_depth_file(args.average_depth_file, avg_depth)
    if args.stage_bam_dir:
        bulk_entry = ("bulk", Path(bulk_bam)) if bulk_bam else None
        stage_bam_dir(args.stage_bam_dir, sc_bams, bulk=bulk_entry)
        logger.info("Staged BAM symlinks in %s", args.stage_bam_dir)

    meta_path = args.output.with_suffix(".meta.txt")
    with open(meta_path, "w") as handle:
        handle.write(f"cell_num\t{len(sc_bams)}\n")
        handle.write(f"with_bulk\t{'yes' if bulk_bam else 'no'}\n")
        handle.write(f"n_sites\t{len(lines)}\n")
        handle.write("samples\t" + ",".join(sample_names) + "\n")
    logger.info("Wrote metadata %s", meta_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
