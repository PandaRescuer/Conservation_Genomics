#!/usr/bin/env python3

import argparse
import gzip
import sqlite3
import sys
import tempfile
from pathlib import Path


def open_text(path, mode="rt"):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def get_coordinate(line, path, line_no):
    fields = line.split("\t", 2)
    if len(fields) < 3:
        raise ValueError(f"{path}:{line_no}: Invalid VCF record")

    try:
        pos = int(fields[1])
        if pos < 1:
            raise ValueError
    except ValueError as exc:
        raise ValueError(
            f"{path}:{line_no}: Invalid VCF position"
        ) from exc

    return fields[0], pos


def reference_sites(path):
    """Yield coordinates incrementally from the RMA VCF."""
    with open_text(path) as infile:
        for line_no, line in enumerate(infile, 1):
            if not line.strip() or line.startswith("#"):
                continue
            yield get_coordinate(line, path, line_no)


def build_index(db, rma_vcf):
    """Store supported coordinates on disk with a bounded memory cache."""
    db.executescript("""
        PRAGMA journal_mode=OFF;
        PRAGMA synchronous=OFF;
        PRAGMA temp_store=FILE;
        PRAGMA cache_size=-32768;
        PRAGMA mmap_size=0;

        CREATE TABLE sites (
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            PRIMARY KEY (chrom, pos)
        ) WITHOUT ROWID;
    """)

    with db:
        db.executemany(
            "INSERT OR IGNORE INTO sites (chrom, pos) VALUES (?, ?)",
            reference_sites(rma_vcf),
        )

    count = db.execute("SELECT COUNT(*) FROM sites").fetchone()[0]
    print(f"RMA_UNIQUE_SITES\t{count}", file=sys.stderr)


def filter_vcf(db, input_vcf, output_vcf):
    """Retain original records whose coordinates occur in the RMA VCF."""
    total = kept = 0
    cursor = db.cursor()

    with open_text(input_vcf) as infile, open_text(output_vcf, "wt") as outfile:
        for line_no, line in enumerate(infile, 1):
            if line.startswith("#"):
                outfile.write(line)
                continue

            if not line.strip():
                continue

            total += 1
            chrom, pos = get_coordinate(line, input_vcf, line_no)

            match = cursor.execute(
                "SELECT 1 FROM sites WHERE chrom = ? AND pos = ?",
                (chrom, pos),
            ).fetchone()

            if match is not None:
                outfile.write(line)
                kept += 1

    print(f"\nINPUT\t{input_vcf}", file=sys.stderr)
    print(f"OUTPUT\t{output_vcf}", file=sys.stderr)
    print(f"TOTAL\t{total}", file=sys.stderr)
    print(f"KEPT\t{kept}", file=sys.stderr)
    print(f"FILTERED\t{total - kept}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Filter two VCF files by coordinates in an RMA VCF."
    )
    parser.add_argument("rma_vcf")
    parser.add_argument("deleterious_vcf")
    parser.add_argument("tolerated_vcf")
    parser.add_argument("output_deleterious")
    parser.add_argument("output_tolerated")
    parser.add_argument(
        "--tmp-dir",
        help="Disk directory for the temporary SQLite index",
    )
    args = parser.parse_args()

    inputs = {
        Path(path).resolve()
        for path in (
            args.rma_vcf,
            args.deleterious_vcf,
            args.tolerated_vcf,
        )
    }
    outputs = [
        Path(args.output_deleterious).resolve(),
        Path(args.output_tolerated).resolve(),
    ]

    if any(path in inputs for path in outputs):
        parser.error("Output files must differ from all input files")
    if outputs[0] == outputs[1]:
        parser.error("The two output files must differ")

    with tempfile.TemporaryDirectory(
        prefix="filter_rma_", dir=args.tmp_dir
    ) as temp_dir:
        db = sqlite3.connect(str(Path(temp_dir) / "sites.sqlite"))
        try:
            build_index(db, args.rma_vcf)
            filter_vcf(
                db, args.deleterious_vcf, args.output_deleterious
            )
            filter_vcf(
                db, args.tolerated_vcf, args.output_tolerated
            )
        finally:
            db.close()


if __name__ == "__main__":
    main()