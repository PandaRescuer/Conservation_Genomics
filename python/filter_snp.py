import sys

f_in = open(sys.argv[1], 'r')
f_out = open(sys.argv[2], 'w')
f_filter = open(sys.argv[3], 'w')
f_depth99 = open(sys.argv[4], 'r')
# read file of depth 99% in each individual in to dict


maxD_dict = {}
for line in f_depth99:
    newline = line.split()
    maxD_dict[newline[0]] = newline[1]

# get list of samples name 
samples = []

for line in f_in:
    if line.startswith('##'):
        pass
    else:
        for i in line.split()[9:]:
            samples.append(i)
        break

f_in.seek(0)

#samples not check missing rate 
not_check_sample_list=['MSH01','MSH02','MSH03','MSH04','MSH05','MSH06','MSH07','MSH08','MSH09','QLA01c','QLA02c','QLA03','QLA04','QLA05','QLA06','QLA07','QLA08','QLI01','QLI02','QLI03','QLI04c','QLI05c','QLI06c','QLI07c']

# keep filters in filter file
# FILTER is PASS or .
for line in f_in:
    # header
    if line.startswith('#'):
        f_out.write(line)
        continue

    newline = line.split()
    if newline[6] not in ('.', 'PASS'):
        f_filter.write('\t'.join(newline[0:6]) + '\t' + 'FAIL_notPASS' + '\t' + '\t'.join(newline[7:]) + '\n')
        continue
    # only keep biallelic SNP
    if len(newline[4]) > 1:
        f_filter.write('\t'.join(newline[0:6]) + '\t' + 'FAIL_notBiaSNP' + '\t' + '\t'.join(newline[7:]) + '\n')
        continue

    missing, excesshet = 0, 0
    GT_list = []

    for i in range(0, len(samples)):
 
        GT = newline[i + 9]
        GT_info_list = GT.split(':')

        sample_AB = GT_info_list[1]
        sample_Depth = GT_info_list[2]
        if samples[i] in not_check_sample_list:
            if GT[:3] == './.':
                GT_list.append(GT)
            else:
                # filter depth > 99% depth in each individual
                if int(sample_Depth) >= int(maxD_dict[samples[i]]):
                        GT_list.append('./.' + GT[3:])
                else:
                    # allele balance
                    if GT[:3] == '0/1' or GT[:3] == '0|1' or GT[:3] == '1|0':
                        sample_AB_list = sample_AB.split(',')
                        sample_AB_REF = int(sample_AB_list[0])
                        sample_AB_ALT = int(sample_AB_list[1])
                        # Filter uninformative genotypes
                        if sample_AB_REF == 0 and sample_AB_ALT == 0:
                            GT_list.append('./.' + GT[3:])
                        else:
                            AB = sample_AB_REF / (sample_AB_REF + sample_AB_ALT)
                            if 0.2 <= AB <= 0.8:
                                GT_list.append(GT)
                                excesshet += 1
                            else:
                                GT_list.append('./.' + GT[3:])
                    else:
                        GT_list.append(GT)
        else:
            if GT[:3] == './.':
                missing += 1
                GT_list.append(GT)
            else:
                # filter depth > 99% depth in each individual
                if int(sample_Depth) >= int(maxD_dict[samples[i]]):
                        GT_list.append('./.' + GT[3:])
                        missing += 1
                else:
                    # allele balance
                    if GT[:3] == '0/1' or GT[:3] == '0|1' or GT[:3] == '1|0':
                        sample_AB_list = sample_AB.split(',')
                        sample_AB_REF = int(sample_AB_list[0])
                        sample_AB_ALT = int(sample_AB_list[1])
                        # Filter uninformative genotypes
                        if sample_AB_REF == 0 and sample_AB_ALT == 0:
                            GT_list.append('./.' + GT[3:])
                            missing += 1
                        else:
                            AB = sample_AB_REF / (sample_AB_REF + sample_AB_ALT)
                            if 0.2 <= AB <= 0.8:
                                GT_list.append(GT)
                                excesshet += 1
                            else:
                                GT_list.append('./.' + GT[3:])
                                missing += 1
                    else:
                        GT_list.append(GT)
    # Filter out sites with more than 20% missingness
    if missing > 0.2 * (len(samples) - len(not_check_sample_list)):
        f_filter.write(
            '\t'.join(newline[0:6]) + '\t' + 'FAIL_missngdata' + '\t' + '\t'.join(newline[7:9])+'\t' + '\t'.join(
                GT_list) + '\n')
        continue
    # Filter out sites without variant after foward steps
    GT_front_list=[]
    for i in GT_list:
        GT_front_list.append(i[:3])
    if len(GT_front_list) == GT_front_list.count('./.') + GT_front_list.count('0/0') + GT_front_list.count('0|0'):
        f_filter.write(
            '\t'.join(newline[0:6]) + '\t' + 'FAIL_novariant' + '\t' + '\t'.join(newline[7:9])+'\t' + '\t'.join(
                GT_list) + '\n')
        continue
    f_out.write('\t'.join(newline[0:9]) + '\t' + '\t'.join(GT_list) + '\n')
f_in.close()
f_out.close()
f_filter.close()
f_depth99.close()

