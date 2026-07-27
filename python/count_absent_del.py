#f_out = open('DEL.txt', 'w')
#f_out2 = open('DEL_2.txt', 'w')
#file_in = 'mis_snp_DELETERIOUS.vcf'

f_out = open('LOF.txt', 'w')
f_out2 = open('LOF_2.txt', 'w')
file_in = 'LOF_snp_RMA.vcf'

f_out.write('NAME\tnum\ttype\trecipient\trescue_pop\n')
f_out2.write('NAME\thet\thom\trecipient\trescue_pop\n')

SC_W = ['GP01', 'GP02', 'GP03', 'GP04', 'GP05', 'GP06', 'GP07', 'GP08', 'GP09', 'GP10']
SC_CAP = ['GP12', 'GP13', 'GP14', 'GP15', 'GP16', 'GP17', 'GP18', 'GP19', 'GP20', 'GP21', 'GP22', 'GP23', 'GP24',
              'GP25', 'GP26', 'GP27', 'GP29', 'GP30', 'GP31', 'GP32', 'GP33', 'GP34', 'GP35', 'GP36', 'GP37',
              'GP38', 'GP40', 'GP43', 'GP44', 'GP46']
MIX_CAP = ['GP41', 'GP42', 'GP45', 'GP47', 'GP48', 'GP49', 'GP50', 'GP51', 'GP52', 'GP53']

f_pop = open('pop.txt', 'r')

name_dict = {}
for line in f_pop:
    newline = line.split()
    if newline[0] == 'ID':
        column = newline[1:]
    else:
        name_dict[newline[2]] = newline[0]

recp_pop_all = ['GP11']

rescue_pop_ALL = SC_W

# for SC_W
for i in range(len(recp_pop_all)):
    recp_pop = recp_pop_all[i]
    rescue_pop = rescue_pop_ALL
    f_in = open(file_in, 'r')
    rescue_pop_index_list = []
    rescue_pop_Dict_het = {}
    rescue_pop_Dict_hom = {}
    for line in f_in:
        if line.startswith('##'):
            pass
        else:
            newline = line.split()
            if line.startswith('#'):
                info_line = newline
                for i in range(9, len(newline)):
                    sample_id = name_dict[newline[i]]
                    if sample_id in rescue_pop:
                        rescue_pop_Dict_het[sample_id] = 0
                        rescue_pop_Dict_hom[sample_id] = 0
                        rescue_pop_index_list.append(i)
                    if sample_id == recp_pop:
                        recp_pop_index = i
            else:
                recp_GT = newline[recp_pop_index][:3]
                if recp_GT == '0/0':
                    for i in rescue_pop_index_list:
                        rescue_GT = newline[i][:3]
                        if rescue_GT == '0/1' or rescue_GT == '0|1' or rescue_GT == '1|0':
                            rescue_pop_Dict_het[name_dict[info_line[i]]] += 1
                        elif rescue_GT == '1/1' or rescue_GT == '1|1':
                            rescue_pop_Dict_hom[name_dict[info_line[i]]] += 1
    sum_list = []

    for key in rescue_pop_Dict_het:
        sum_list.append(rescue_pop_Dict_het[key])
        f_out.write('\t'.join([str(key), str(rescue_pop_Dict_het[key]), 'het', recp_pop,'SC_WLP']) + '\n')
    print(recp_pop, sum(sum_list) / len(sum_list))
    for key in rescue_pop_Dict_hom:
        f_out.write('\t'.join([str(key), str(rescue_pop_Dict_hom[key]), 'hom', recp_pop,'SC_WLP']) + '\n')
        f_out2.write(
            '\t'.join([str(key), str(rescue_pop_Dict_het[key]), str(rescue_pop_Dict_hom[key]), recp_pop,'SC_WLP']) + '\n')

    f_in.close()

rescue_pop_ALL = SC_CAP

# for SC_CAP
for i in range(len(recp_pop_all)):
    recp_pop = recp_pop_all[i]
    rescue_pop = rescue_pop_ALL
    f_in = open(file_in, 'r')
    rescue_pop_index_list = []
    rescue_pop_Dict_het = {}
    rescue_pop_Dict_hom = {}
    for line in f_in:
        if line.startswith('##'):
            pass
        else:
            newline = line.split()
            if line.startswith('#'):
                info_line = newline
                for i in range(9, len(newline)):
                    sample_id = name_dict[newline[i]]
                    if sample_id in rescue_pop:
                        rescue_pop_Dict_het[sample_id] = 0
                        rescue_pop_Dict_hom[sample_id] = 0
                        rescue_pop_index_list.append(i)
                    if sample_id == recp_pop:
                        recp_pop_index = i
            else:
                recp_GT = newline[recp_pop_index][:3]
                if recp_GT == '0/0':
                    for i in rescue_pop_index_list:
                        rescue_GT = newline[i][:3]
                        if rescue_GT == '0/1' or rescue_GT == '0|1' or rescue_GT == '1|0':
                            rescue_pop_Dict_het[name_dict[info_line[i]]] += 1
                        elif rescue_GT == '1/1' or rescue_GT == '1|1':
                            rescue_pop_Dict_hom[name_dict[info_line[i]]] += 1
    sum_list = []

    for key in rescue_pop_Dict_het:
        sum_list.append(rescue_pop_Dict_het[key])
        f_out.write('\t'.join([str(key), str(rescue_pop_Dict_het[key]), 'het', recp_pop,'SC_CAP']) + '\n')
    print(recp_pop, sum(sum_list) / len(sum_list))
    for key in rescue_pop_Dict_hom:
        f_out.write('\t'.join([str(key), str(rescue_pop_Dict_hom[key]), 'hom', recp_pop,'SC_CAP']) + '\n')
        f_out2.write(
            '\t'.join([str(key), str(rescue_pop_Dict_het[key]), str(rescue_pop_Dict_hom[key]), recp_pop,'SC_CAP']) + '\n')

    f_in.close()

#for MIX_CAP
rescue_pop_ALL = MIX_CAP

for i in range(len(recp_pop_all)):
    recp_pop = recp_pop_all[i]
    rescue_pop = rescue_pop_ALL
    f_in = open(file_in, 'r')
    rescue_pop_index_list = []
    rescue_pop_Dict_het = {}
    rescue_pop_Dict_hom = {}
    for line in f_in:
        if line.startswith('##'):
            pass
        else:
            newline = line.split()
            if line.startswith('#'):
                info_line = newline
                for i in range(9, len(newline)):
                    sample_id = name_dict[newline[i]]
                    if sample_id in rescue_pop:
                        rescue_pop_Dict_het[sample_id] = 0
                        rescue_pop_Dict_hom[sample_id] = 0
                        rescue_pop_index_list.append(i)
                    if sample_id == recp_pop:
                        recp_pop_index = i
            else:
                recp_GT = newline[recp_pop_index][:3]
                if recp_GT == '0/0':
                    for i in rescue_pop_index_list:
                        rescue_GT = newline[i][:3]
                        if rescue_GT == '0/1' or rescue_GT == '0|1' or rescue_GT == '1|0':
                            rescue_pop_Dict_het[name_dict[info_line[i]]] += 1
                        elif rescue_GT == '1/1' or rescue_GT == '1|1':
                            rescue_pop_Dict_hom[name_dict[info_line[i]]] += 1
    sum_list = []

    for key in rescue_pop_Dict_het:
        sum_list.append(rescue_pop_Dict_het[key])
        f_out.write('\t'.join([str(key), str(rescue_pop_Dict_het[key]), 'het', recp_pop,'MIX_CAP']) + '\n')
    print(recp_pop, sum(sum_list) / len(sum_list))
    for key in rescue_pop_Dict_hom:
        f_out.write('\t'.join([str(key), str(rescue_pop_Dict_hom[key]), 'hom', recp_pop,'MIX_CAP']) + '\n')
        f_out2.write(
            '\t'.join([str(key), str(rescue_pop_Dict_het[key]), str(rescue_pop_Dict_hom[key]), recp_pop,'MIX_CAP']) + '\n')

    f_in.close()
f_out.close()
f_out2.close()
