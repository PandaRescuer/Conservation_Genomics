
import sys
import gzip
import os
from pathlib import Path


def find_RMA(sites_dict, same_sites_dict, input_VCF_filename, output_filename, replacements):
    total = 0
    kept = 0
    corrected = 0
    unique_sites = set()

    with open(input_VCF_filename, 'r') as VCF_infile, open(
        output_filename, 'w'
    ) as outfile:
        for line in VCF_infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
                if not line.strip():
                    continue

                total += 1

                line_info = line.split()
                ref = line_info[3]
                alt = line_info[4]
                annoinfo = line_info[7].split(';')
                chrom = line_info[0]
                pos = line_info[1]
                genotype_list = line_info[9:]
                pos0 = int(pos) - 1

                if pos0 in sites_dict.get(chrom, {}):
                    site_info = sites_dict[chrom][pos0].split(',')
                    maf_ingroup_site = site_info[0]
                    maf_outgroup_site = site_info[1]

                    if alt == maf_outgroup_site:
                        if len(ref) != 1 or len(alt) != 1 or ref.upper() not in "ACGT" or alt.upper() not in "ACGT":
                            raise ValueError(f"Not a biallelic SNP: {chrom}:{pos}")
                        key = (chrom, pos0)
                        change = (ref.upper(), alt.upper())
                        if key in replacements and replacements[key] != change:
                            raise ValueError(f"Conflicting replacements: {chrom}:{pos}")
                        replacements[key] = change
                        new_genotype_list = []

                        for i in range(0, len(genotype_list), 1):
                            genotype_info = genotype_list[i].split(':')
                            genotype = genotype_info[0]
                            depth = genotype_info[1].split(',')
                            ref_depth = depth[0]
                            alt_depth = depth[1]
                            depth = alt_depth + ',' + ref_depth

                            if genotype == './.':
                                pass
                            elif genotype == '0/0' or genotype == '0|0':
                                genotype = '1/1'
                            elif (
                                genotype == '0/1'
                                or genotype == '0|1'
                                or genotype == '1|0'
                            ):
                                genotype = '0/1'
                            elif genotype == '1/1' or genotype == '1|1':
                                genotype = '0/0'

                            new_genotype = (
                                genotype + ':' + depth + ':'
                                + ':'.join(genotype_info[2:])
                            )
                            new_genotype_list.append(new_genotype)

                        new_annoinfo = annoinfo[2:] + ['RMA_corrected']
                        out_list = (
                            line_info[:3]
                            + [alt, ref]
                            + line_info[5:7]
                            + [';'.join(new_annoinfo)]
                            + [line_info[8]]
                            + ['\t'.join(new_genotype_list)]
                        )
                        outfile.write('\t'.join(out_list) + '\n')

                        kept += 1
                        corrected += 1
                        unique_sites.add((chrom, pos0))

                    else:
                        outfile.write(line)
                        kept += 1
                        unique_sites.add((chrom, pos0))

                elif pos0 in same_sites_dict.get(chrom, set()):
                    # Retain sites supported by SAME BED.
                    outfile.write(line)
                    kept += 1
                    unique_sites.add((chrom, pos0))

                # Omit sites absent from both BED files.

    # Generate statistics.
    stats = [
        ("RMA_UNIQUE_SITES", len(unique_sites)),
        ("INPUT", input_VCF_filename),
        ("OUTPUT", output_filename),
        ("TOTAL", total),
        ("KEPT", kept),
        ("FILTERED", total - kept),
        ("RMA_CORRECTED", corrected),
    ]

    # Save statistics to a separate file.
    stats_file = output_filename + ".stats.txt"

    with open(stats_file, "w") as stats_out:
        for key, value in stats:
            stats_out.write(f"{key}\t{value}\n")

    # Also print statistics to terminal.
    print()
    for key, value in stats:
        print(f"{key}\t{value}", flush=True)

    print(f"STATS_FILE\t{stats_file}", flush=True)

    return total, kept, corrected


def write_revised_fasta(input_fasta, output_fasta, replacements):
    # Stream FASTA lines and apply only sites corrected in the input VCFs.
    by_chrom = {}
    for (chrom, pos0), (ref, alt) in replacements.items():
        by_chrom.setdefault(chrom, []).append((pos0, ref, alt))
    for values in by_chrom.values():
        values.sort()

    opener = gzip.open if str(input_fasta).endswith('.gz') else open
    temporary = str(output_fasta) + '.partial'
    chrom = None
    offset = 0
    index = 0
    targets = []
    seen = set()
    changed = 0
    try:
        with opener(input_fasta, 'rt') as src, open(temporary, 'w') as dst:
            for line in src:
                if line.startswith('>'):
                    if index != len(targets):
                        raise ValueError(f"Replacement outside sequence: {chrom}")
                    chrom = line[1:].split()[0]
                    if chrom in seen:
                        raise ValueError(f"Duplicate FASTA sequence: {chrom}")
                    seen.add(chrom)
                    targets = by_chrom.get(chrom, [])
                    offset = 0
                    index = 0
                    dst.write(line.rstrip() + '\n')
                    continue
                seq = ''.join(line.split())
                if not seq:
                    continue
                if chrom is None:
                    raise ValueError('Sequence before FASTA header')
                end = offset + len(seq)
                bases = None
                while index < len(targets) and targets[index][0] < end:
                    pos0, ref, alt = targets[index]
                    local = pos0 - offset
                    if local < 0 or seq[local].upper() != ref:
                        raise ValueError(f"FASTA/VCF REF mismatch: {chrom}:{pos0 + 1}, expected {ref}")
                    if bases is None:
                        bases = list(seq)
                    bases[local] = alt.lower() if seq[local].islower() else alt
                    changed += (ref != alt)
                    index += 1
                dst.write((''.join(bases) if bases is not None else seq) + '\n')
                offset = end
            if not seen:
                raise ValueError('Empty FASTA')
            if index != len(targets):
                raise ValueError(f"Replacement outside sequence: {chrom}")
            missing = set(by_chrom) - seen
            if missing:
                raise ValueError(f"Chromosomes absent from FASTA: {sorted(missing)}")
        os.replace(temporary, output_fasta)
    except Exception:
        if os.path.exists(temporary):
            os.remove(temporary)
        raise
    print(f"FASTA_UNIQUE_CORRECTED_SITES\t{changed}", flush=True)
    with open(str(output_fasta) + '.stats.txt', 'w') as out:
        out.write(f"FASTA_UNIQUE_CORRECTED_SITES\t{changed}\n")


if __name__ == "__main__":
    # Arguments: DIFFERENT BED, SAME BED, followed by input/output VCF pairs.
    if len(sys.argv) < 7 or (len(sys.argv) - 5) % 2 != 0:
        sys.exit(
            "Usage: python maf_find_variant_RMA.py "
            "site.bed site.same.bed original.fa revised.fa "
            "input1.vcf output1.vcf "
            "[input2.vcf output2.vcf ...]"
        )

    input_BED_file = sys.argv[1]
    input_SAME_BED_file = sys.argv[2]
    input_fasta = sys.argv[3]
    output_fasta = sys.argv[4]
    vcf_pairs = list(zip(sys.argv[5::2], sys.argv[6::2]))
    inputs = [input_BED_file, input_SAME_BED_file, input_fasta] + [p[0] for p in vcf_pairs]
    outputs = [output_fasta, output_fasta + '.partial', output_fasta + '.stats.txt', 'RMA_all_VCF_summary.stats.txt']
    outputs += [p[1] for p in vcf_pairs] + [p[1] + '.stats.txt' for p in vcf_pairs]
    input_paths = {Path(p).resolve() for p in inputs}
    output_paths = [Path(p).resolve() for p in outputs]
    if input_paths.intersection(output_paths) or len(set(output_paths)) != len(output_paths):
        sys.exit('Input and output paths must not overlap; output paths must be unique.')
    replacements = {}

    # Read DIFFERENT BED once.
    print(f"Reading DIFFERENT BED: {input_BED_file}", flush=True)
    special_site_dictionaries = {}

    with open(input_BED_file, 'r') as BED_infile:
        for line in BED_infile:
            if not line.strip() or line.startswith('#'):
                continue

            line_info = line.split()
            chrom = line_info[0]
            pos0 = int(line_info[1])

            if chrom not in special_site_dictionaries:
                special_site_dictionaries[chrom] = {}

            special_site_dictionaries[chrom][pos0] = line_info[3]

    # Read SAME BED once.
    print(f"Reading SAME BED: {input_SAME_BED_file}", flush=True)
    same_site_dictionaries = {}

    with open(input_SAME_BED_file, 'r') as BED_infile:
        for line in BED_infile:
            if not line.strip() or line.startswith('#'):
                continue

            line_info = line.split()
            chrom = line_info[0]
            pos0 = int(line_info[1])

            if chrom not in same_site_dictionaries:
                same_site_dictionaries[chrom] = set()

            same_site_dictionaries[chrom].add(pos0)

    print("BED loading completed.", flush=True)

    # Reuse the BED dictionaries for every VCF.
    grand_total = 0
    grand_kept = 0
    grand_corrected = 0

    for index, (input_VCF_file, output_file) in enumerate(vcf_pairs, 1):
        print(
            f"\n[{index}/{len(vcf_pairs)}] Processing: {input_VCF_file}",
            flush=True,
        )

        total, kept, corrected = find_RMA(
            special_site_dictionaries,
            same_site_dictionaries,
            input_VCF_file,
            output_file,
            replacements,
        )

        grand_total += total
        grand_kept += kept
        grand_corrected += corrected

        print(f"Saved: {output_file}", flush=True)

    summary_file = "RMA_all_VCF_summary.stats.txt"

    summary = [
        ("TOTAL", grand_total),
        ("KEPT", grand_kept),
        ("FILTERED", grand_total - grand_kept),
        ("RMA_CORRECTED", grand_corrected),
    ]

    with open(summary_file, "w") as summary_out:
        for key, value in summary:
            summary_out.write(f"{key}\t{value}\n")

    print("\nALL_VCF_SUMMARY")
    for key, value in summary:
        print(f"{key}\t{value}", flush=True)

    print(f"SUMMARY_FILE\t{summary_file}", flush=True)
    # Release BED dictionaries before streaming the original genome.
    del special_site_dictionaries, same_site_dictionaries
    write_revised_fasta(input_fasta, output_fasta, replacements)
    print(f"Saved FASTA: {output_fasta}", flush=True)
