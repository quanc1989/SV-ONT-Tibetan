import configargparse
import numpy as np


def config_opts():
    parser = configargparse.ArgumentParser(
        description='func_select_SEQ_fromVCFs',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add('-f', '--file', required=True,
               help='sample list file')
    parser.add('-s', '--save', required=True,
               help='save path prefix')
    return parser


if __name__ == '__main__':

    parser = config_opts()
    opt = parser.parse_args()

    with open(opt.file, 'r') as file_filter:
        list_samples = file_filter.read().strip().split('\n')

    list_ins = {}
    list_rnames = []

    with open(opt.save + '.seq', 'w') as file_write, open(opt.save + '.rnames', 'w') as file_write_rnames:
        for sample in list_samples:

            print(sample)

            with open(sample + '.files', 'r') as file_filter:
                list_vcfs = file_filter.read().strip().split('\n')

            for file_vcf in list_vcfs:

                print(file_vcf)

                if 'svim' in file_vcf:
                    str_seq = 'SEQS'
                    str_rnames = 'READS'
                elif 'sniffles' in file_vcf:
                    str_seq = 'SEQ'
                    str_rnames = 'RNAMES'
                else:
                    continue

                # index_line = 1
                with open(file_vcf, 'r') as file_read:
                    while True:
                        line = file_read.readline().strip()
                        if line:
                            if line.startswith('#'):
                                continue
                            else:

                                # print(str(index_line))
                                # index_line += 1
                                strs = line.split('\t')
                                id = strs[2]
                                chr = strs[0]
                                pos = strs[1]
                                info = strs[7].split(';')

                                svtype = [x.split('=')[1] for x in info if 'SVTYPE=' in x][0]
                                # end = [x.split('=')[1] for x in info if 'END=' in x][0]
                                seq = [x.split('=')[1] for x in info if str_seq + '=' in x]
                                rnames = [x.split('=')[1] for x in info if str_rnames + '=' in x]

                                if svtype.__contains__('INS'):
                                    if len(seq) > 0 and seq[0] != '':
                                        if chr+'_'+pos not in list_ins:
                                            list_ins[chr+'_'+pos] = []

                                        if seq[0].__contains__(','):
                                            list_ins[chr + '_' + pos].extend(seq[0].split(','))
                                        else:
                                            list_ins[chr+'_'+pos].append(seq[0])
                                        # if chr+'_'+pos == 'chr2_92314292':
                                        #     print(seq[0])
                                        #     print(line)
                                        #     print(list_ins[chr+'_'+pos])
                                    else:

                                        if not strs[4].__contains__('INS'):

                                            if chr + '_' + pos not in list_ins:
                                                list_ins[chr + '_' + pos] = []
                                            if strs[4].__contains__(','):
                                                list_ins[chr + '_' + pos].extend(strs[4].split(','))
                                            else:
                                                list_ins[chr + '_' + pos].append(strs[4])


                                if len(rnames) > 0:
                                    for rname in rnames[0].split(','):
                                        if rname not in list_rnames and rname != 'input':
                                            list_rnames.append(rname)

                        else:
                            break

        for id_ins, arry_seq in list_ins.items():

            arry_seq_len = [len(x) for x in arry_seq]
            index_seq = arry_seq_len.index(min(arry_seq_len))

            # if id_ins == 'chr2_92314292':
            #     print(arry_seq)
            #     print(arry_seq_len)
            #     print(index_seq)
            #     print(arry_seq[index_seq])
            #     print(id_ins + '\t' + arry_seq[index_seq] + '\t' + str(len(arry_seq)) +'\t' + str(np.mean(np.array(arry_seq_len))) + '\t' + str(np.std(np.array(arry_seq_len))) + '\n')

            file_write.writelines(id_ins + '\t' + arry_seq[index_seq] + '\t' + str(len(arry_seq)) +'\t' + str(np.mean(np.array(arry_seq_len))) + '\t' + str(np.std(np.array(arry_seq_len))) + '\n')

        for rname in list_rnames:
            file_write_rnames.writelines(rnames+'\n')









