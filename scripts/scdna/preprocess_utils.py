#!/usr/bin/env python3
"""Shared helpers for scDNA BAM-to-tree-input preprocessing."""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)

BASES = ("A", "T", "G", "C")
DEFAULT_POPAF = 1e-4
BULK_STEMS = {"bulk", "pseudobulk", "pseudo_bulk"}


def parse_mutation_id(raw: str) -> Optional[Tuple[str, int, str, str]]:
    """Parse chr_pos_ref_alt (or whitespace-separated) mutation identifiers."""
    text = raw.strip()
    if not text or text.startswith("#"):
        return None
    if re.match(r"^(chr|chrom|id|mutid)\b", text, flags=re.IGNORECASE):
        return None

    parts = re.split(r"[\s,;]+", text)
    if len(parts) >= 4 and parts[1].isdigit() and len(parts[2]) == 1 and len(parts[3]) == 1:
        chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2].upper(), parts[3].upper()
        chrom = chrom[3:] if chrom.lower().startswith("chr") and not chrom.lower().startswith("chrom") else chrom
        if chrom.lower().startswith("chrom"):
            chrom = chrom[5:]
        return chrom, pos, ref, alt

    token = parts[0]
    fields = token.split("_")
    if len(fields) < 4:
        return None
    ref, alt = fields[-2].upper(), fields[-1].upper()
    if len(ref) != 1 or len(alt) != 1 or ref not in "ACGTN" or alt not in "ACGTN":
        return None
    pos_str = fields[-3]
    if not pos_str.isdigit():
        return None
    chrom = "_".join(fields[:-3])
    if chrom.lower().startswith("chr") and not chrom.lower().startswith("chrom"):
        chrom = chrom[3:]
    return chrom, int(pos_str), ref, alt


def load_mutation_list(path: Path) -> List[Tuple[str, int, str, str]]:
    mutations: List[Tuple[str, int, str, str]] = []
    seen = set()
    with open(path) as handle:
        for line in handle:
            parsed = parse_mutation_id(line)
            if parsed is None:
                continue
            if parsed in seen:
                continue
            seen.add(parsed)
            mutations.append(parsed)
    if not mutations:
        raise ValueError(f"No mutations parsed from {path}")
    logger.info("Loaded %d mutations from %s", len(mutations), path)
    return mutations


def default_popaf_dict(ref: str, alt: str) -> Dict[str, float]:
    popaf = {base: DEFAULT_POPAF for base in BASES}
    remaining = 1.0 - DEFAULT_POPAF * 3
    if ref in popaf:
        popaf[ref] = round(remaining, 8) if remaining > 0 else DEFAULT_POPAF
    if alt in popaf and alt != ref:
        popaf[alt] = DEFAULT_POPAF
    return popaf


def default_allele_pool(ref: str, alt: str, observed: Iterable[str]) -> List[str]:
    pool: List[str] = []
    for base in (ref, alt):
        if base in "ACGT" and base not in pool:
            pool.append(base)
    for base in observed:
        if base in "ACGT" and base not in pool:
            pool.append(base)
    return pool


def is_bulk_name(name: str) -> bool:
    stem = Path(name).stem.lower()
    return stem in BULK_STEMS


def find_bam_index(bam_path: Path) -> Optional[Path]:
    candidates = [
        bam_path.with_suffix(bam_path.suffix + ".bai"),
        bam_path.with_suffix(".bai"),
        Path(str(bam_path) + ".bai"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _read_name_list(path: Path) -> List[str]:
    names: List[str] = []
    with open(path) as handle:
        for i, line in enumerate(handle):
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            cols = re.split(r"[\s,]+", text)
            if i == 0 and cols[0].lower() in {"sample", "sampleid", "scid", "scid_basedtree", "cell", "cellid"}:
                continue
            name = cols[0]
            if name.lower() in {"scid_basedtree", "sampleid"}:
                continue
            names.append(name)
    return names


def _bam_stem(path: Path) -> str:
    name = path.name
    if name.endswith(".bam"):
        return name[:-4]
    return path.stem


def collect_sc_bams(
    bam_dir: Optional[Path] = None,
    bam_list: Optional[Path] = None,
    sample_list: Optional[Path] = None,
    bulk_bam: Optional[Path] = None,
    keep_bulk: bool = False,
) -> Tuple[List[Tuple[str, Path]], Optional[Path]]:
    """
    Resolve single-cell BAM paths.

    Bulk is never required. A BAM named bulk.bam is skipped unless keep_bulk=True
    or --bulk-bam is provided.
    """
    detected_bulk: Optional[Path] = Path(bulk_bam) if bulk_bam else None
    entries: List[Tuple[str, Path]] = []

    if bam_list:
        with open(bam_list) as handle:
            for line in handle:
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                path = Path(text.split()[0])
                stem = _bam_stem(path)
                if is_bulk_name(stem):
                    if detected_bulk is None:
                        detected_bulk = path
                    logger.info("Detected bulk BAM in bam-list: %s", path)
                    continue
                entries.append((stem, path))
    elif bam_dir:
        bam_dir = Path(bam_dir)
        for path in sorted(bam_dir.glob("*.bam")):
            if path.name.endswith(".bai"):
                continue
            stem = _bam_stem(path)
            if is_bulk_name(stem):
                if detected_bulk is None:
                    detected_bulk = path
                logger.info("Detected bulk BAM in bam-dir: %s", path)
                continue
            entries.append((stem, path))
    else:
        raise ValueError("Either --bam-dir or --bam-list is required")

    if sample_list:
        wanted = _read_name_list(Path(sample_list))
        wanted = [name for name in wanted if not is_bulk_name(name)]
        by_name = {name: path for name, path in entries}
        missing = [name for name in wanted if name not in by_name]
        if missing:
            # Allow sample-list names to match BAM filenames in bam_dir.
            if bam_dir:
                for name in missing:
                    candidate = Path(bam_dir) / f"{name}.bam"
                    if candidate.exists():
                        by_name[name] = candidate
            missing = [name for name in wanted if name not in by_name]
        if missing:
            raise FileNotFoundError(
                "BAM files not found for sample-list entries: " + ", ".join(missing[:10])
            )
        entries = [(name, by_name[name]) for name in wanted]

    if not entries:
        raise FileNotFoundError("No single-cell BAM files found")

    for name, path in entries:
        if not path.exists():
            raise FileNotFoundError(f"BAM not found: {path}")
        if find_bam_index(path) is None:
            raise FileNotFoundError(
                f"BAM index not found for {path}. Expected {path}.bai or {path.with_suffix('.bai')}"
            )

    if detected_bulk is not None and not Path(detected_bulk).exists():
        logger.warning("Bulk BAM path does not exist and will be ignored: %s", detected_bulk)
        detected_bulk = None

    logger.info("Using %d single-cell BAM files", len(entries))
    if keep_bulk:
        return entries, detected_bulk
    return entries, None


def write_sample_name_file(path: Path, sample_names: Sequence[str], bulk_name: Optional[str] = None) -> None:
    with open(path, "w") as handle:
        if bulk_name:
            handle.write(f"{bulk_name}\n")
        for name in sample_names:
            handle.write(f"{name}\n")


def write_scid_file(path: Path, sample_names: Sequence[str]) -> None:
    with open(path, "w") as handle:
        handle.write("scid_basedTree\tsampleid\n")
        for name in sample_names:
            handle.write(f"{name}\t{name}\n")


def write_average_depth_file(path: Path, depths: Dict[str, float]) -> None:
    with open(path, "w") as handle:
        handle.write("sample\taverage_depth\n")
        for sample, depth in depths.items():
            handle.write(f"{sample}\t{depth}\n")


def stage_bam_dir(staging_dir: Path, entries: Sequence[Tuple[str, Path]], bulk: Optional[Tuple[str, Path]] = None) -> Path:
    """Create a directory of named BAM symlinks for feature extraction."""
    staging_dir.mkdir(parents=True, exist_ok=True)
    staged: List[Tuple[str, Path]] = []
    if bulk:
        staged.append(bulk)
    staged.extend(entries)
    for name, src in staged:
        dst = staging_dir / f"{name}.bam"
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        os.symlink(src.resolve(), dst)
        index = find_bam_index(src)
        if index is not None:
            for dst_index in (staging_dir / f"{name}.bam.bai", staging_dir / f"{name}.bai"):
                if dst_index.exists() or dst_index.is_symlink():
                    dst_index.unlink()
            os.symlink(index.resolve(), staging_dir / f"{name}.bam.bai")
    return staging_dir


def quote_r_field(value: str) -> str:
    """Quote a field so R read.table(header=FALSE, sep='\\t') keeps Python dicts intact."""
    text = str(value)
    if text.startswith('"') and text.endswith('"'):
        return text
    return '"' + text.replace('"', '\\"') + '"'
