# Correct reference-misidentified allele (RMA) sites in a VCF file and
# update the corresponding reference genome based on candidate sites
# provided in a BED file.

import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord

def find_RMA(sites_dict, input_VCF_filename, output_filename):
    """
    Identify and correct RMA sites by swapping reference/alt alleles and genotypes
    when the VCF alt matches the outgroup allele.
    """
    replace_count = 0  # Counter for number of sites corrected
    with open(input_VCF_filename, 'r') as VCF_infile, open(output_filename, 'w') as outfile:
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
                if int(pos)-1 in sites_dict[chrom].keys():
                    site_info = sites_dict[chrom][int(pos)-1].split(',')
                    maf_ingroup_site = site_info[0]
                    maf_outgroup_site = site_info[1]
                    # If VCF alt equals the outgroup allele, reverse alleles and correct genotypes
                    if alt == maf_outgroup_site:
                        replace_count += 1
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
                            elif genotype == '0/1' or genotype == '0|1' or genotype == '1|0':
                                genotype = '0/1'
                            elif genotype == '1/1' or genotype == '1|1':
                                genotype = '0/0'
                            new_genotype = genotype + ':' + depth + ':' + ':'.join(genotype_info[2:])
                            new_genotype_list.append(new_genotype)
                        new_annoinfo = annoinfo[2:] + ['RMA_corrected']
                        out_list = line_info[:3] + [alt, ref] + line_info[5:7] + [';'.join(new_annoinfo)] + [line_info[8]] + [
                            '\t'.join(new_genotype_list)]
                        chr2seq_dict[chrom][int(pos)-1] = alt
                        outfile.write('\t'.join(out_list) + '\n')
                    else:
                        outfile.write(line)
                else:
                    outfile.write(line)
    return replace_count

if __name__ == "__main__":
    input_BED_file = sys.argv[1]    # Input BED file path
    input_VCF_file = sys.argv[2]    # Input VCF file path
    ref_fasta = sys.argv[3]         # Original reference genome
    ref_fasta_new = sys.argv[4]     # New reference genome (output)
    output_file = sys.argv[5]       # Output VCF file path

    # Load original reference genome
    chr2seq_dict = {}
    for seq_record in SeqIO.parse(ref_fasta, "fasta"):
        chrseq = seq_record.seq
        chr2seq_dict[str(seq_record.id)] = MutableSeq(chrseq)

    # Parse BED file to get special sites (RMA candidates)
    special_site_dictionaries = {}
    with open(input_BED_file, 'r') as BED_infile:
        for line in BED_infile:
            line_info = line.split()
            try:
                special_site_dictionaries[line_info[0]][int(line_info[1])] = line_info[3]
            except KeyError:
                special_site_dictionaries[line_info[0]] = {int(line_info[1]): line_info[3]}

    # Run correction and get count of replaced sites
    replaced_count = find_RMA(special_site_dictionaries, input_VCF_file, output_file)

    # Write new reference genome with corrected bases
    my_records = []
    for key in chr2seq_dict.keys():
        records = SeqRecord(Seq(chr2seq_dict[key]), id=key)
        my_records.append(records)
    SeqIO.write(my_records, ref_fasta_new, "fasta")

    # Print summary statistics
    print(f"Processing complete. {replaced_count} sites were replaced. Results saved to {output_file} and new reference genome {ref_fasta_new}")