
#!/usr/bin/env python3

import sys
import csv
import os
import subprocess
import tempfile
from pathlib import Path
from collections import Counter


BASES = set("ACGT")


def parse_maf_block(block, raw_file, stats):
    """
    Process one MAF alignment block and stream valid
    ancestral allele observations to a temporary file.
    """

    stats["TOTAL_BLOCKS"] += 1

    records = []

    for line in block:
        fields = line.split()

        if fields and fields[0] == "s":
            if len(fields) < 7:
                stats["MALFORMED_S_ROWS"] += 1
                return

            records.append(fields)

    if len(records) != 6:
        stats["SKIPPED_NON_SIX_SPECIES_BLOCKS"] += 1
        return

    ref = records[0]

    ref_src = ref[1]
    ref_start = int(ref[2])
    ref_size = int(ref[3])
    ref_strand = ref[4]

    if not ref_src.startswith("Ailuropoda_melanoleuca."):
        raise ValueError(
            f"Unexpected reference species: {ref_src}"
        )

    if ref_strand != "+":
        raise ValueError(
            f"Negative-strand reference block: {ref_src}"
        )

    chrom = ref_src.split(".", 1)[1]

    sequences = [record[6].upper() for record in records]

    if len(set(map(len, sequences))) != 1:
        raise ValueError(
            f"Inconsistent alignment lengths: {ref_src}"
        )

    if sum(base != "-" for base in sequences[0]) != ref_size:
        raise ValueError(
            f"Reference alignment size mismatch: {ref_src}"
        )

    stats["VALID_SIX_SPECIES_BLOCKS"] += 1

    ref_pos = ref_start

    for column in zip(*sequences):

        stats["TOTAL_ALIGNMENT_COLUMNS"] += 1

        ref_base = column[0]

        if ref_base == "-":
            stats["SKIPPED_REFERENCE_GAP_COLUMNS"] += 1
            continue

        pos0 = ref_pos
        ref_pos += 1

        if any(base not in BASES for base in column):
            stats["SKIPPED_INVALID_BASE_COLUMNS"] += 1
            continue

        outgroup_bases = column[1:]

        if len(set(outgroup_bases)) != 1:
            stats["SKIPPED_OUTGROUP_DISAGREEMENT"] += 1
            continue

        ancestral_base = outgroup_bases[0]

        status = (
            "SAME"
            if ref_base == ancestral_base
            else "DIFFERENT"
        )

        raw_file.write(
            f"{chrom}\t{pos0}\t"
            f"{ref_base}\t{ancestral_base}\t{status}\n"
        )

        stats["VALID_SITE_OBSERVATIONS"] += 1


def stream_maf(input_maf, raw_file, stats):
    """
    Read MAF blocks sequentially without retaining
    previously processed blocks in memory.
    """

    current_block = []

    with open(input_maf, "r") as infile:

        for line in infile:

            fields = line.split(None, 1)

            is_block_start = (
                bool(fields) and fields[0] == "a"
            )

            if is_block_start:

                if current_block:
                    parse_maf_block(
                        current_block,
                        raw_file,
                        stats
                    )

                current_block = [line]

            elif current_block:
                current_block.append(line)

        if current_block:
            parse_maf_block(
                current_block,
                raw_file,
                stats
            )


def sort_observations(raw_path, sorted_path, temp_dir):
    """
    Use GNU sort for memory-limited external sorting.

    Records are sorted by chromosome and numeric position.
    """

    command = [
        "sort",
        "-T", temp_dir,
        "-S", "1G",
        "--parallel=2",
        "-t", "\t",
        "-k1,1",
        "-k2,2n",
        raw_path,
        "-o", sorted_path
    ]

    subprocess.run(
        command,
        check=True,
        env={**os.environ, "LC_ALL": "C"}
    )


def resolve_observations(
    sorted_path,
    different_bed,
    same_bed,
    conflict_file,
    stats
):
    """
    Resolve duplicate reference coordinates in sorted order.

    Identical observations are collapsed.
    Conflicting observations are excluded from both BED files.
    """

    with (
        open(sorted_path, "r") as infile,
        open(different_bed, "w") as fdiff,
        open(same_bed, "w") as fsame,
        open(conflict_file, "w", newline="") as fconf
    ):

        writer = csv.writer(
            fconf,
            delimiter="\t"
        )

        writer.writerow([
            "CHROM",
            "START",
            "END",
            "OBSERVATION_COUNT",
            "OBSERVATIONS"
        ])

        current_key = None
        current_count = 0
        observations = set()

        def flush_group():
            nonlocal current_key
            nonlocal current_count
            nonlocal observations

            if current_key is None:
                return

            chrom, pos0 = current_key

            stats["TOTAL_UNIQUE_COORDINATES"] += 1

            if current_count > 1:
                stats["DUPLICATE_COORDINATES"] += 1
                stats["EXTRA_DUPLICATE_OBSERVATIONS"] += (
                    current_count - 1
                )

            if len(observations) > 1:

                stats["CONFLICTING_COORDINATES"] += 1

                writer.writerow([
                    chrom,
                    pos0,
                    pos0 + 1,
                    current_count,
                    ";".join(
                        f"{ref},{anc},{status}"
                        for ref, anc, status
                        in sorted(observations)
                    )
                ])

            else:

                ref_base, ancestral_base, status = next(
                    iter(observations)
                )

                bed_line = (
                    f"{chrom}\t{pos0}\t{pos0 + 1}\t"
                    f"{ref_base},{ancestral_base}\n"
                )

                if status == "SAME":
                    fsame.write(bed_line)
                    stats["SAME_UNIQUE_SITES"] += 1

                else:
                    fdiff.write(bed_line)
                    stats["DIFFERENT_UNIQUE_SITES"] += 1

        for line in infile:

            chrom, pos, ref, anc, status = line.rstrip(
                "\n"
            ).split("\t")

            pos0 = int(pos)
            key = (chrom, pos0)

            observation = (ref, anc, status)

            if key != current_key:

                flush_group()

                current_key = key
                current_count = 0
                observations = set()

            current_count += 1
            observations.add(observation)

        flush_group()


def process_maf(input_maf, output_bed):

    output_path = Path(output_bed)

    prefix = (
        str(output_path.with_suffix(""))
        if output_path.suffix == ".bed"
        else str(output_path)
    )

    different_bed = output_bed
    same_bed = prefix + ".same.bed"
    conflict_file = prefix + ".conflicts.tsv"
    summary_file = prefix + ".summary.tsv"

    stats = Counter()

    # Temporary files are removed automatically on completion.
    with tempfile.TemporaryDirectory(
        prefix="maf_bed_",
        dir=str(output_path.parent.resolve())
    ) as temp_dir:

        raw_path = os.path.join(
            temp_dir,
            "observations.raw.tsv"
        )

        sorted_path = os.path.join(
            temp_dir,
            "observations.sorted.tsv"
        )

        # Stream observations to disk.
        with open(raw_path, "w") as raw_file:

            stream_maf(
                input_maf,
                raw_file,
                stats
            )

        # Sort observations using bounded memory.
        sort_observations(
            raw_path,
            sorted_path,
            temp_dir
        )

        # Resolve duplicates and write final BED files.
        resolve_observations(
            sorted_path,
            different_bed,
            same_bed,
            conflict_file,
            stats
        )

    with open(summary_file, "w") as outfile:

        outfile.write("METRIC\tCOUNT\n")

        for metric in sorted(stats):
            outfile.write(
                f"{metric}\t{stats[metric]}\n"
            )

    print("\nMAF ancestral allele extraction summary")
    print("=" * 60)

    for metric in [
        "TOTAL_BLOCKS",
        "VALID_SIX_SPECIES_BLOCKS",
        "SKIPPED_NON_SIX_SPECIES_BLOCKS",
        "TOTAL_ALIGNMENT_COLUMNS",
        "SKIPPED_REFERENCE_GAP_COLUMNS",
        "SKIPPED_INVALID_BASE_COLUMNS",
        "SKIPPED_OUTGROUP_DISAGREEMENT",
        "VALID_SITE_OBSERVATIONS",
        "TOTAL_UNIQUE_COORDINATES",
        "DUPLICATE_COORDINATES",
        "EXTRA_DUPLICATE_OBSERVATIONS",
        "CONFLICTING_COORDINATES",
        "SAME_UNIQUE_SITES",
        "DIFFERENT_UNIQUE_SITES"
    ]:
        print(f"{metric}: {stats[metric]:,}")

    print("\nOutput files:")
    print(f"DIFFERENT BED: {different_bed}")
    print(f"SAME BED: {same_bed}")
    print(f"Conflicts: {conflict_file}")
    print(f"Summary: {summary_file}")

    print("\nProcessing completed.")


if __name__ == "__main__":

    if len(sys.argv) != 3:
        sys.exit(
            "Usage: python maf_find_variant_bed.py "
            "input.maf output.bed"
        )

    process_maf(
        sys.argv[1],
        sys.argv[2]
    )
