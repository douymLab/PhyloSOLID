#!/usr/bin/env python3
"""Convert unphased genotyper output into generate_treeinput_data.R input.

The R script expects 12 tab-separated columns, no header, with quoted Python
literals, and a dummy bulk dict as the first element of the reads list when
bulk was not used during genotyping.
"""

from __future__ import annotations

import argparse
import ast
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_SCRIPTS = SCRIPT_DIR.parent
if str(PACKAGE_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(PACKAGE_SCRIPTS))

from preprocess_utils import quote_r_field  # noqa: E402

logger = logging.getLogger("format_treeinput")

DUMMY_BULK = {"no_reads": "*"}
GENOTYPER_HEADER_PREFIX = "chr"


def split_tsv(line: str) -> List[str]:
    return line.rstrip("\n").split("\t")


def parse_reads_key(fields: Sequence[str]) -> Tuple[str, str, str, str]:
    return fields[0], fields[1], fields[2], fields[3]


def load_reads_bed(path: Path) -> Dict[Tuple[str, str, str, str], str]:
    reads = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = split_tsv(line)
            if len(fields) < 6:
                continue
            reads[parse_reads_key(fields)] = fields[5]
    return reads


def maybe_eval(value: str):
    try:
        return ast.literal_eval(value)
    except Exception:
        return value


def ensure_reads_with_bulk(reads_raw: str, with_bulk: bool) -> str:
    reads = maybe_eval(reads_raw)
    if not isinstance(reads, list):
        raise ValueError(f"reads_unphased is not a list: {reads_raw[:80]}")
    if with_bulk:
        return str(reads)
    if reads and isinstance(reads[0], dict):
        return str([DUMMY_BULK] + reads)
    return str([DUMMY_BULK] + list(reads))


def ensure_whether_mutated(raw: str, with_bulk: bool) -> str:
    values = maybe_eval(raw)
    if not isinstance(values, list):
        return raw
    if with_bulk:
        return str(values)
    # Feature extraction without bulk expects one flag per single cell.
    return str(values)


def convert_line(
    fields: Sequence[str],
    reads_lookup: Dict[Tuple[str, str, str, str], str],
    with_bulk: bool,
) -> Optional[Tuple[str, str]]:
    if len(fields) < 12:
        return None
    key = parse_reads_key(fields)
    reads_raw = reads_lookup.get(key)
    if reads_raw is None:
        logger.warning("No pileup reads for site %s", key)
        return None
    tree_cols = [
        fields[0],
        fields[1],
        fields[2],
        fields[3],
        fields[4],
        fields[5],
        fields[6],
        fields[7],
        fields[8],
        fields[10],  # log likelihood mut
        fields[11],  # log likelihood unmut
        ensure_reads_with_bulk(reads_raw, with_bulk),
    ]
    quoted = []
    for i, col in enumerate(tree_cols):
        if i in {4, 5, 6, 8, 9, 10, 11}:
            quoted.append(quote_r_field(col))
        else:
            quoted.append(str(col))
    tree_line = "\t".join(quoted)

    whether = fields[-1] if len(fields) >= 21 else "[]"
    feature_cols = [
        fields[0],
        fields[1],
        fields[2],
        fields[3],
        fields[4],
        fields[5],
        fields[6],
        fields[7],
        ensure_whether_mutated(whether, with_bulk),
    ]
    feature_quoted = []
    for i, col in enumerate(feature_cols):
        if i in {4, 5, 6, 8}:
            feature_quoted.append(quote_r_field(col))
        else:
            feature_quoted.append(str(col))
    feature_line = "\t".join(feature_quoted)
    return tree_line, feature_line


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Format unphased genotyper output for PhyloSOLID scDNA tree input.")
    parser.add_argument("-g", "--genotyper-output", required=True, type=Path)
    parser.add_argument("-r", "--reads-bed", required=True, type=Path)
    parser.add_argument("-o", "--output", required=True, type=Path, help="12-column tree_input file")
    parser.add_argument("--features-input", type=Path, help="Optional classifier feature-extraction input")
    parser.add_argument("--with-bulk", choices=["yes", "no"], default="no")
    return parser.parse_args()


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    args = parse_args()
    with_bulk = args.with_bulk == "yes"
    reads_lookup = load_reads_bed(args.reads_bed)
    tree_lines: List[str] = []
    feature_lines: List[str] = []
    skipped = 0
    with open(args.genotyper_output) as handle:
        for i, line in enumerate(handle):
            if not line.strip():
                continue
            fields = split_tsv(line)
            if i == 0 and fields and fields[0] == GENOTYPER_HEADER_PREFIX:
                continue
            converted = convert_line(fields, reads_lookup, with_bulk)
            if converted is None:
                skipped += 1
                continue
            tree_line, feature_line = converted
            tree_lines.append(tree_line)
            feature_lines.append(feature_line)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        handle.write("\n".join(tree_lines) + ("\n" if tree_lines else ""))
    logger.info("Wrote %d tree-input sites to %s (skipped %d)", len(tree_lines), args.output, skipped)

    if args.features_input:
        args.features_input.parent.mkdir(parents=True, exist_ok=True)
        with open(args.features_input, "w") as handle:
            handle.write("\n".join(feature_lines) + ("\n" if feature_lines else ""))
        logger.info("Wrote classifier feature input to %s", args.features_input)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
