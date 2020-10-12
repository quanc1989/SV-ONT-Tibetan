import configargparse
import numpy as np


def config_opts():
    parser = configargparse.ArgumentParser(
        description='func_select_SEQ_fromVCFs',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add('-f', '--file', required=True,
               help='vcf path file')
    parser.add('-s', '--save', required=True,
               help='save path')
    return parser


if __name__ == '__main__':

    parser = config_opts()
    opt = parser.parse_args()

    with open(opt.file, 'r') as file_filter:
        list_vcfs = file_filter.read().strip().split('\n')

    list_ins = {}
    with open(opt.save, 'w') as file_write:
        for file_vcf in list_vcfs:
            with open(file_vcf, 'r') as file_read:
                while True:
                    line = file_read.readline().strip()
                    if line:
                        if line.startswith('#'):
                            continue
                        else:
                            strs = line.split('\t')
                            id = strs[2]
                            chr = strs[0]
                            pos = strs[1]
                            info = strs[7].split(';')
                            # print(info)
                            svtype = [x.split('=')[1] for x in info if 'SVTYPE=' in x][0]
                            end = [x.split('=')[1] for x in info if 'END=' in x][0]
                            seq = [x.split('=')[1] for x in info if 'SEQ=' in x]

                            if svtype.__contains__('INS'):
                                if len(seq) > 0 and seq[0] != '':
                                    if chr+'_'+pos not in list_ins:
                                        list_ins[chr+'_'+pos] = []
                                    # if pos == 'chr2_92314292':
                                    #     print(list_ins[chr+'_'+pos])

                                    if seq[0].__contains__(','):
                                        list_ins[chr + '_' + pos].extend(seq[0].split(','))
                                    else:
                                        list_ins[chr + '_' + pos].append(seq[0])

                    else:
                        break

        for id_ins, arry_seq in list_ins.items():
            arry_seq_len = [len(x) for x in arry_seq]
            index_seq = arry_seq_len.index(min(arry_seq_len))

            file_write.writelines(id_ins + '\t' + arry_seq[index_seq] + '\t' + str(len(arry_seq)) +'\t' + str(np.mean(np.array(arry_seq_len))) + '\t' + str(np.std(np.array(arry_seq_len))) + '\n')











