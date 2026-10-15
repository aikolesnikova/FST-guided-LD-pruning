#!/usr/bin/env python3
"""Greedy, FST-prioritised LD pruning using PLINK 1.9 pairwise r²."""
from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    snp: str
    fst: float


def load_candidates(path: Path) -> list[Variant]:
    """Read seven-column PLINK BIM-like input: CHR ID CM BP A1 A2 MaxFst."""
    candidates = []
    seen = set()
    with path.open(encoding="utf-8") as f:
        for n, line in enumerate(f, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 7:
                raise ValueError(f"{path}:{n}: expected 7 columns (BIM + MaxFst), got {len(fields)}")
            chrom, snp, cm, bp, a1, a2, score = fields
            try:
                pos, fst = int(bp), float(score)
            except ValueError as exc:
                raise ValueError(f"{path}:{n}: invalid POS or MaxFst") from exc
            if pos < 1 or not math.isfinite(fst):
                raise ValueError(f"{path}:{n}: POS must be positive and MaxFst finite")
            if snp in seen:
                raise ValueError(f"{path}:{n}: duplicate SNP ID {snp}")
            seen.add(snp)
            candidates.append(Variant(chrom, pos, snp, fst))
    return sorted(candidates, key=lambda v: (-v.fst, v.chrom, v.pos, v.snp))


def parse_ld(path: Path, lead: str) -> dict[str, float]:
    """Parse standard PLINK 1.9 --r2 LD output; return partners and r²."""
    partners = {}
    with path.open(encoding="utf-8") as f:
        header = f.readline().split()
        needed = {"SNP_A", "SNP_B", "R2"}
        if not needed.issubset(header):
            raise ValueError(f"{path}: missing PLINK LD columns: {sorted(needed - set(header))}")
        idx = {key: header.index(key) for key in needed}
        for n, line in enumerate(f, 2):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < len(header):
                raise ValueError(f"{path}:{n}: truncated LD row")
            a, b = fields[idx["SNP_A"]], fields[idx["SNP_B"]]
            if lead not in (a, b):
                raise ValueError(f"{path}:{n}: pair does not include lead {lead}")
            other = b if a == lead else a
            try:
                r2 = float(fields[idx["R2"]])
            except ValueError:
                # PLINK can emit NA for pairs without sufficient data.
                continue
            if not math.isfinite(r2) or not 0 <= r2 <= 1:
                raise ValueError(f"{path}:{n}: invalid r² value {r2}")
            partners[other] = max(r2, partners.get(other, -1.0))
    return partners


def validate_bfile(bfile: Path, candidates: list[Variant]) -> None:
    for suffix in (".bed", ".bim", ".fam"):
        if not Path(f"{bfile}{suffix}").is_file():
            raise FileNotFoundError(f"Missing PLINK input: {bfile}{suffix}")
    bim_ids = set()
    with Path(f"{bfile}.bim").open(encoding="utf-8") as f:
        for n, line in enumerate(f, 1):
            row = line.split()
            if len(row) < 6:
                raise ValueError(f"{bfile}.bim:{n}: incomplete BIM row")
            if row[1] in bim_ids:
                raise ValueError(f"{bfile}.bim:{n}: duplicate variant ID {row[1]}")
            bim_ids.add(row[1])
    missing = set(v.snp for v in candidates) - bim_ids
    if missing:
        raise ValueError(f"{len(missing)} scored SNP IDs are missing from the PLINK BIM (e.g. {sorted(missing)[0]})")


def write_table(path: Path, columns: list[str], rows: list[tuple]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        w.writerows(rows)


def run(args: argparse.Namespace) -> None:
    scores = Path(args.scores)
    bfile = Path(args.bfile)
    out = Path(args.outdir)
    if out.exists() and any(out.iterdir()) and not args.overwrite:
        raise FileExistsError(f"Output directory is not empty: {out}. Use --overwrite to replace outputs.")
    variants = load_candidates(scores)
    if not variants:
        raise ValueError("No candidate variants in the score file")
    validate_bfile(bfile, variants)
    if not shutil.which(args.plink):
        raise FileNotFoundError(f"PLINK executable not found: {args.plink}")
    out.mkdir(parents=True, exist_ok=True)
    remaining = {v.snp for v in variants}
    leads, excluded = [], []
    with tempfile.TemporaryDirectory(prefix="fst_ld_") as temp:
        prefix = str(Path(temp) / "pairwise")
        for variant in variants:
            if variant.snp not in remaining:
                continue
            iteration = len(leads) + 1
            cmd = [args.plink, "--bfile", str(bfile), "--ld-snp", variant.snp,
                   "--ld-window", str(args.ld_window), "--ld-window-kb", str(args.window_kb),
                   "--r2", "--ld-window-r2", "0", "--threads", str(args.threads),
                   "--memory", str(args.memory_mb), "--out", prefix]
            result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            if result.returncode != 0:
                raise RuntimeError(f"PLINK failed for {variant.snp} (exit {result.returncode}):\n{result.stdout[-2000:]}")
            ld_path = Path(prefix + ".ld")
            if not ld_path.is_file():
                raise RuntimeError(f"PLINK did not produce {ld_path} for {variant.snp}")
            partners = parse_ld(ld_path, variant.snp)
            remaining.remove(variant.snp)  # explicit progress, independent of PLINK self-LD rows
            leads.append((iteration, variant.chrom, variant.pos, variant.snp, variant.fst))
            for other in sorted(remaining):
                r2 = partners.get(other)
                if r2 is not None and r2 > args.r2_threshold:
                    remaining.remove(other)
                    excluded.append((other, variant.snp, r2, iteration))
            ld_path.unlink()
            print(f"[lead {iteration}] {variant.snp}: remaining {len(remaining)}", flush=True)
    write_table(out / "lead_variants.tsv", ["Iteration", "CHR", "POS", "SNP", "MaxFst"], leads)
    write_table(out / "excluded_variants.tsv", ["SNP", "Excluded_by", "R2", "Lead_iteration"], excluded)
    write_table(out / "summary.tsv", ["Input_candidates", "Lead_variants", "Excluded_variants", "R2_threshold", "Window_kb"],
                [(len(variants), len(leads), len(excluded), args.r2_threshold, args.window_kb)])
    if len(leads) + len(excluded) != len(variants):
        raise AssertionError("Invariant violated: every candidate must be selected or excluded exactly once")
    print(f"Finished: {len(variants)} candidates -> {len(leads)} leads, {len(excluded)} excluded; {out}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bfile", required=True, help="PLINK 1.9 BED/BIM/FAM prefix")
    p.add_argument("--scores", required=True, help="Seven-column BIM-style candidate file, 7th column = MaxFst")
    p.add_argument("--outdir", required=True)
    p.add_argument("--r2-threshold", type=float, default=0.05, help="Exclude if r² > threshold (default 0.05)")
    p.add_argument("--window-kb", type=int, default=10000, help="LD window in kb (default 10000 = 10 Mb)")
    p.add_argument("--ld-window", type=int, default=99999, help="PLINK count-based LD window")
    p.add_argument("--plink", default="plink", help="PLINK 1.9 executable")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--memory-mb", type=int, default=8000)
    p.add_argument("--overwrite", action="store_true", help="Allow existing output directory")
    a = p.parse_args()
    if not 0 <= a.r2_threshold <= 1 or a.window_kb < 1 or a.ld_window < 1 or a.threads < 1 or a.memory_mb < 1:
        p.error("Check ranges: 0 <= r² <= 1 and positive window/threads/memory")
    try:
        run(a)
    except (FileNotFoundError, FileExistsError, ValueError, RuntimeError) as exc:
        p.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
