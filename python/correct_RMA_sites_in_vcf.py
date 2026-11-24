import sys


def find_RMA(sites_dict, same_sites_dict, input_VCF_filename, output_filename):
    with open(input_VCF_filename, 'r') as VCF_infile, open(
        output_filename, 'w'
    ) as outfile:
        for line in VCF_infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
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
                    else:
                        outfile.write(line)

                elif pos0 in same_sites_dict.get(chrom, set()):
                    # Retain sites supported by SAME BED.
                    outfile.write(line)

                # Omit sites absent from both BED files.


if __name__ == "__main__":
    # Arguments: DIFFERENT BED, SAME BED, followed by input/output VCF pairs.
    if len(sys.argv) < 5 or (len(sys.argv) - 3) % 2 != 0:
        sys.exit(
            "Usage: python maf_find_variant_RMA.py "
            "site.bed site.same.bed "
            "input1.vcf output1.vcf "
            "[input2.vcf output2.vcf ...]"
        )

    input_BED_file = sys.argv[1]
    input_SAME_BED_file = sys.argv[2]
    vcf_pairs = list(zip(sys.argv[3::2], sys.argv[4::2]))

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
    for index, (input_VCF_file, output_file) in enumerate(vcf_pairs, 1):
        print(
            f"[{index}/{len(vcf_pairs)}] Processing: {input_VCF_file}",
            flush=True,
        )

        find_RMA(
            special_site_dictionaries,
            same_site_dictionaries,
            input_VCF_file,
            output_file,
        )

        print(f"Saved: {output_file}", flush=True)

    print("All VCF files completed.", flush=True)
