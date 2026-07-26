# Classify VCF variants according to SIFT predictions.
# The script reads SIFT annotation results and a VCF file, then separates
# variants into DELETERIOUS and TOLERATED VCF files.


import argparse
from collections import defaultdict


def load_sift_predictions(sift_filename):
    """Load and unify SIFT predictions by chromosome and position."""
    sift_dict = defaultdict(lambda: defaultdict(list))

    with open(sift_filename, "r") as sift_file:
        for line in sift_file:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue

            chromosome = fields[0]
            position = fields[1]
            prediction = fields[-1]

            sift_dict[chromosome][position].append(prediction)

    unified_sift_dict = {}

    for chromosome, positions in sift_dict.items():
        unified_sift_dict[chromosome] = {}

        for position, predictions in positions.items():
            if "DELETERIOUS" in predictions:
                unified_sift_dict[chromosome][position] = "DELETERIOUS"
            elif "DELETERIOUS (*WARNING! Low confidence)" in predictions:
                unified_sift_dict[chromosome][position] = (
                    "DELETERIOUS (*WARNING! Low confidence)"
                )
            elif "TOLERATED" in predictions:
                unified_sift_dict[chromosome][position] = "TOLERATED"

    return unified_sift_dict


def classify_vcf(
    sift_filename,
    input_vcf_filename,
    deleterious_vcf_filename,
    tolerated_vcf_filename
):
    """Split VCF variants into DELETERIOUS and TOLERATED files."""
    sift_predictions = load_sift_predictions(sift_filename)

    with open(input_vcf_filename, "r") as vcf_file, \
            open(deleterious_vcf_filename, "w") as deleterious_file, \
            open(tolerated_vcf_filename, "w") as tolerated_file:

        for line in vcf_file:
            if line.startswith("#"):
                deleterious_file.write(line)
                tolerated_file.write(line)
                continue

            fields = line.split()
            chromosome = fields[0]
            position = fields[1]

            prediction = sift_predictions.get(chromosome, {}).get(position)

            # Write variants according to their SIFT prediction
            if prediction == "DELETERIOUS":
                deleterious_file.write(line)
            elif prediction == "TOLERATED":
                tolerated_file.write(line)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Classify VCF variants into DELETERIOUS and TOLERATED "
            "files using SIFT predictions."
        )
    )
    parser.add_argument("sift_file", help="Input SIFT prediction file")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument(
        "deleterious_vcf",
        help="Output VCF file for DELETERIOUS variants"
    )
    parser.add_argument(
        "tolerated_vcf",
        help="Output VCF file for TOLERATED variants"
    )

    args = parser.parse_args()

    classify_vcf(
        args.sift_file,
        args.input_vcf,
        args.deleterious_vcf,
        args.tolerated_vcf
    )


if __name__ == "__main__":
    main()