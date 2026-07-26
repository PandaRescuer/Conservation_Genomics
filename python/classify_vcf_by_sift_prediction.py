file_dir = 'D://analyse'
f_sift = open(file_dir + '//sift//SIFT.txt', 'r')
f_vcf = open(file_dir + '//MAF//missense_snp_RMA.vcf', 'r')
f_DELE = open(file_dir + '//sift//mis_snp_DELETERIOUS.vcf', 'w')
f_TOLE = open(file_dir + '//sift//mis_snp_TOLERATED.vcf', 'w')

siftdict = {}
for line in f_sift:
    newline = line.split('\t')
    if newline[0] in siftdict.keys():
        pass
    else:
        siftdict[newline[0]] = {}
    try:
        siftdict[newline[0]][newline[1]].append(newline[-1][:-1])
    except KeyError:
        siftdict[newline[0]][newline[1]] = [newline[-1][:-1]]

total_withsift = 0

uni_siftdict = {}
for i in siftdict.keys():
    uni_siftdict[i] = {}
    for j in siftdict[i].keys():
        if 'DELETERIOUS' in siftdict[i][j]:
            uni_siftdict[i][j] = 'DELETERIOUS'
        elif 'DELETERIOUS (*WARNING! Low confidence)' in siftdict[i][j]:
            uni_siftdict[i][j] = 'DELETERIOUS (*WARNING! Low confidence)'
        elif 'TOLERATED' in siftdict[i][j]:
            uni_siftdict[i][j] = 'TOLERATED'
        else:
            print(i,j)
            pass

for line in f_vcf:
    if line.startswith('#'):
        f_DELE.write(line)
        f_TOLE.write(line)
    else:
        newline = line.split()
        chrom=newline[0]
        pos=newline[1]
        info=newline[7].split(';')
        if pos in uni_siftdict[chrom].keys():
            # if uni_siftdict[chrom][pos] == 'DELETERIOUS' or uni_siftdict[chrom][pos] == 'DELETERIOUS (*WARNING! Low confidence)':
            if uni_siftdict[chrom][pos] == 'DELETERIOUS':
                f_DELE.write(line)

            elif uni_siftdict[chrom][pos] == 'TOLERATED':
                f_TOLE.write(line)

f_sift.close()
f_vcf.close()
f_TOLE.close()
f_DELE.close()

