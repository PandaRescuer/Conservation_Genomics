#!/usr/bin/env python3

import argparse
from collections import defaultdict


SC_W = [
    "GP01", "GP02", "GP03", "GP04", "GP05",
    "GP06", "GP07", "GP08", "GP09", "GP10"
]

SC_CAP = [
    "GP12", "GP13", "GP14", "GP15", "GP16", "GP17",
    "GP18", "GP19", "GP20", "GP21", "GP22", "GP23",
    "GP24", "GP25", "GP26", "GP27", "GP29", "GP30",
    "GP31", "GP32", "GP33", "GP34", "GP35", "GP36",
    "GP37", "GP38", "GP40", "GP43", "GP44", "GP46"
]

MIX_CAP = [
    "GP41", "GP42", "GP45", "GP47", "GP48",
    "GP49", "GP50", "GP51", "GP52", "GP53"
]

RECIPIENTS = ["GP11"]

RESCUE_POPULATIONS = {
    "SC_WLP": SC_W,
    "SC_CAP": SC_CAP,
    "MIX_CAP": MIX_CAP
}


def read_name_dict(pop_file):
    """Map VCF sample names to panda IDs."""
    name_dict = {}

    with open(pop_file, "r") as handle:
        for line in handle:
            fields = line.split()

            if not fields or fields[0] == "ID":
                continue

            name_dict[fields[2]] = fields[0]

    return name_dict


def count_genotypes(vcf_file, name_dict, recipient, rescue_pop):
    """Count recipient-absent heterozygous and homozygous genotypes."""
    het_counts = defaultdict(int)
    hom_counts = defaultdict(int)

    rescue_indices = []
    recipient_index = None
    header = None

    with open(vcf_file, "r") as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            fields = line.split()

            if line.startswith("#"):
                header = fields

                for index in range(9, len(fields)):
                    sample_id = name_dict[fields[index]]

                    if sample_id in rescue_pop:
                        rescue_indices.append(index)
                        het_counts[sample_id] = 0
                        hom_counts[sample_id] = 0

                    if sample_id == recipient:
                        recipient_index = index

                continue

            recipient_gt = fields[recipient_index][:3]

            if recipient_gt != "0/0":
                continue

            for index in rescue_indices:
                rescue_id = name_dict[header[index]]
                rescue_gt = fields[index][:3]

                if rescue_gt in {"0/1", "0|1", "1|0"}:
                    het_counts[rescue_id] += 1
                elif rescue_gt in {"1/1", "1|1"}:
                    hom_counts[rescue_id] += 1

    return het_counts, hom_counts


def main():
    parser = argparse.ArgumentParser(
        description="Count recipient-absent DEL and LOF genotypes."
    )
    parser.add_argument("del_vcf", help="DEL VCF file")
    parser.add_argument("lof_vcf", help="LOF VCF file")
    parser.add_argument("pop_file", help="Population metadata file")
    parser.add_argument("output_file", help="Output TSV file")
    args = parser.parse_args()

    name_dict = read_name_dict(args.pop_file)

    with open(args.output_file, "w") as output:
        output.write(
            "NAME\thet_DEL\thom_DEL\thet_LOF\thom_LOF"
            "\trecipient\trescue_pop\n"
        )

        for rescue_name, rescue_pop in RESCUE_POPULATIONS.items():
            for recipient in RECIPIENTS:
                del_het, del_hom = count_genotypes(
                    args.del_vcf,
                    name_dict,
                    recipient,
                    rescue_pop
                )

                lof_het, lof_hom = count_genotypes(
                    args.lof_vcf,
                    name_dict,
                    recipient,
                    rescue_pop
                )

                for sample_id in rescue_pop:
                    output.write(
                        "\t".join(
                            [
                                sample_id,
                                str(del_het[sample_id]),
                                str(del_hom[sample_id]),
                                str(lof_het[sample_id]),
                                str(lof_hom[sample_id]),
                                recipient,
                                rescue_name,
                            ]
                        )
                        + "\n"
                    )


if __name__ == "__main__":
    main()
