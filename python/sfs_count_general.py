#!/usr/bin/env python3

import argparse
from collections import defaultdict
from pathlib import Path


POPULATIONS = {
    "SC_WLP": [
        "GP01", "GP02", "GP03", "GP04", "GP05",
        "GP06", "GP07", "GP09", "GP10"
    ],
    "QLI": [
        "QLI05c", "QLI06c", "QLI07c", "QLI04c", "GP55"
    ],
    "SC_CAP": [
        "GP12", "GP13", "GP14", "GP15", "GP16",
        "GP17", "GP18", "GP19", "GP23", "GP24",
        "GP25", "GP27", "GP28", "GP29", "GP30",
        "GP31", "GP32", "GP33", "GP34", "GP35",
        "GP38", "GP40", "GP43", "GP44", "GP46"
    ],
    "MIX_CAP": [
        "GP39", "GP41", "GP42", "GP45", "GP47",
        "GP48", "GP49", "GP50", "GP51", "GP52",
        "GP53"
    ]
}


def read_sample_mapping(pop_file):
    """Map VCF sample names to population IDs."""
    name_dict = {}
    header = None

    with pop_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()

            if not fields:
                continue

            if fields[0] == "ID":
                header = fields
                continue

            if header is None or len(fields) < 3:
                continue

            population_id = fields[0]
            sample_name = fields[2]
            name_dict[sample_name] = population_id

    return name_dict


def initialize_sfs_dict():
    """Initialize alternative allele count bins for each population."""
    sfs_dict = {}

    for population, samples in POPULATIONS.items():
        sfs_dict[population] = defaultdict(int)

        for allele_count in range(1, 2 * len(samples) + 1):
            sfs_dict[population][allele_count] = 0

    return sfs_dict


def classify_genotype(genotype):
    """Return the alternative allele count for a genotype."""
    genotype = genotype.split(":", 1)[0]

    if genotype in {"./.", ".|."}:
        return None

    if genotype in {"0/0", "0|0"}:
        return 0

    if genotype in {"0/1", "0|1", "1/0", "1|0"}:
        return 1

    if genotype in {"1/1", "1|1"}:
        return 2

    return None


def count_sfs(vcf_file, name_dict):
    """Count alternative allele frequencies by population."""
    sfs_dict = initialize_sfs_dict()
    sample_names = None

    with vcf_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                fields = line.rstrip("\r\n").split()
                sample_names = fields[9:]
                continue

            if line.startswith("#"):
                continue

            if sample_names is None:
                raise ValueError(
                    f"VCF header not found before variant records: {vcf_file}"
                )

            fields = line.rstrip("\r\n").split()

            alt_counts = {
                population: 0
                for population in POPULATIONS
            }

            for index, sample_name in enumerate(sample_names):
                genotype = fields[index + 9]
                allele_count = classify_genotype(genotype)

                if allele_count is None:
                    continue

                population_id = name_dict.get(sample_name)

                if population_id is None:
                    continue

                for population, population_samples in POPULATIONS.items():
                    if population_id in population_samples:
                        alt_counts[population] += allele_count
                        break

            for population, allele_count in alt_counts.items():
                if allele_count > 0:
                    sfs_dict[population][allele_count] += 1

    return sfs_dict


def write_sfs(output_file, sfs_dict):
    """Write the site-frequency spectrum table."""
    with output_file.open("w", encoding="utf-8") as handle:
        handle.write("POP\tfreq\tnum\n")

        for population, frequency_dict in sfs_dict.items():
            for frequency, count in sorted(frequency_dict.items()):
                handle.write(
                    f"{population}\t{frequency}\t{count}\n"
                )


def process_vcf(vcf_prefix, pop_file, output_dir):
    """Process one VCF prefix."""
    vcf_file = Path(f"{vcf_prefix}.vcf")
    output_file = output_dir / f"{vcf_file.stem}_SFS_no_MISs.txt"

    name_dict = read_sample_mapping(pop_file)
    sfs_dict = count_sfs(vcf_file, name_dict)
    write_sfs(output_file, sfs_dict)

    print(f"Processed: {vcf_file}")
    print(f"Output: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate population-level alternative allele frequency spectra."
    )

    parser.add_argument("--LOF", required=True, help="LOF VCF prefix")
    parser.add_argument("--DEL", required=True, help="DEL VCF prefix")
    parser.add_argument("--TOL", required=True, help="TOL VCF prefix")
    parser.add_argument("--MIS", required=True, help="MISSENSE VCF prefix")
    parser.add_argument("--SYN", required=True, help="SYN VCF prefix")
    parser.add_argument("--CNE", required=True, help="CNE VCF prefix")
    parser.add_argument("--pop-file", required=True, help="Population metadata file")
    parser.add_argument("--output-dir", required=True, help="Output directory")

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pop_file = Path(args.pop_file)

    vcf_prefixes = [
        args.LOF,
        args.DEL,
        args.SYN,
        args.MIS,
        args.TOL,
        args.CNE,
    ]

    for vcf_prefix in vcf_prefixes:
        process_vcf(
            vcf_prefix=vcf_prefix,
            pop_file=pop_file,
            output_dir=output_dir,
        )


if __name__ == "__main__":
    main()
