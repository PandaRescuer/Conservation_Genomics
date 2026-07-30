LOF_file_dir = 'LOF_snp_RMA'
DEL_file_dir = 'mis_snp_DELETERIOUS'
TOL_file_dir = 'mis_snp_TOLERATED'
Mis_file_dir = 'missense_snp_RMA'
SYN_file_dir = 'synonymous_snp_RMA'
CNE_file_dir = 'CNE_snp_RMA'

filelist = [LOF_file_dir,DEL_file_dir,SYN_file_dir,Mis_file_dir,TOL_file_dir,CNE_file_dir]
for file in filelist:
    f_vcf = open(file + '.vcf', 'r')
    f_SFS = open(file + '_SFS_no_Miss.txt', 'w')
    f_SFS.write('POP\tfreq\tnum\n')
    f_pop = open('pop.txt', 'r')
    pop_dict = {}
    name_dict = {}
    for line in f_pop:
        newline = line.split()
        if newline[0] == 'ID':
            column = newline[1:]
        else:
            pop_dict[newline[0]] = {}
            for i in range(0, len(column)):
                pop_dict[newline[0]][column[i]] = newline[i + 1]
            name_dict[newline[2]]=newline[0]
    samples = []

    # samples
    SC_WLP = ['GP01', 'GP02', 'GP03', 'GP04', 'GP05', 'GP06', 'GP07', 'GP09', 'GP10']
    QLI = ['QLI05c', 'QLI06c', 'QLI07c', 'QLI04c', 'GP55']
    SC_CAP = ['GP12', 'GP13', 'GP14', 'GP15', 'GP16', 
              'GP17', 'GP18', 'GP19', 'GP23', 'GP24',
              'GP25', 'GP27', 'GP28', 'GP29', 'GP30', 
              'GP31', 'GP32', 'GP33', 'GP34', 'GP35', 
              'GP38',  'GP40', 'GP43', 'GP44',  'GP46']
    MIX_CAP = ['GP39','GP41','GP42', 'GP45','GP47', 'GP48', 'GP49', 'GP50', 'GP51', 'GP52', 'GP53']


    Pop_alt_num_dict = {}

    def create_dict(POP_list, keyname):
        Pop_alt_num_dict[keyname] = {}
        for i in range(1, len(POP_list) + 1):
            for j in range(i, i * 2 + 1, 1):
                Pop_alt_num_dict[keyname][j] = 0

    create_dict(SC_WLP, 'SC_WLP')
    create_dict(SC_CAP, 'SC_CAP')
    create_dict(QLI, 'QLI')
    create_dict(MIX_CAP, 'MIX_CAP')

    for line in f_vcf:
        if line[0:2] == '##':
            continue
        if line.startswith('#'):
            sample_list = line.split()[9:]
            for i in range(0, len(sample_list)):
                samples.append(sample_list[i])
            continue
        else:
            # if './.' in line:
                # continue
            newline = line.split()
            GT_list = []
            sample_GT_dict = {'0/0': [], '0/1': [], '1/1': []}
            for i in range(0, len(samples)):
                # if samples[i] in not_include_pop:
                #     continue
                GT = newline[i + 9]
                if GT[:3] == './.':
                    pass
                elif GT[:3] == '0/0' or GT[:3] == '0|0':
                    sample_GT_dict['0/0'].append(name_dict[samples[i]])
                elif GT[:3] == '1/1' or GT[:3] == '1|1':
                    sample_GT_dict['1/1'].append(name_dict[samples[i]])
                elif GT[:3] == '0/1' or GT[:3] == '0|1' or GT[:3] == '1|0':
                    sample_GT_dict['0/1'].append(name_dict[samples[i]])
                else:
                    print(line)

            # count the alt num in each populations
            alt_num_dict = {'SC_WLP': 0, 'QLI': 0, 'SC_CAP': 0, 'MIX_CAP': 0}
            for i in sample_GT_dict['0/1']:
                if i in SC_WLP:
                    alt_num_dict['SC_WLP'] += 1
                elif i in QLI:
                    alt_num_dict['QLI'] += 1
                elif i in SC_CAP:
                    alt_num_dict['SC_CAP'] += 1
                elif i in MIX_CAP:
                    alt_num_dict['MIX_CAP'] += 1
            for i in sample_GT_dict['1/1']:
                if i in SC_WLP:
                    alt_num_dict['SC_WLP'] += 2
                elif i in QLI:
                    alt_num_dict['QLI'] += 2
                elif i in SC_CAP:
                    alt_num_dict['SC_CAP'] += 2
                elif i in MIX_CAP:
                    alt_num_dict['MIX_CAP'] += 2
            alt_num_list = list(alt_num_dict.values())
            # SFS
            # print(alt_num_dict)
            for pop in alt_num_dict.keys():
                if alt_num_dict[pop]>0:
                    Pop_alt_num_dict[pop][alt_num_dict[pop]] += 1

    print(Pop_alt_num_dict)
    for pop in Pop_alt_num_dict.keys():
        for freq in Pop_alt_num_dict[pop].keys():
            f_SFS.write(pop + '\t' + str(freq) + '\t' + str(Pop_alt_num_dict[pop][freq]) + '\n')
    f_vcf.close()
    f_pop.close()
    f_SFS.close()
