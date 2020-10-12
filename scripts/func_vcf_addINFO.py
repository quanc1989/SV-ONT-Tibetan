import configargparse


def config_opts():
    parser = configargparse.ArgumentParser(
        description='func_vcf_filter.py',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add('-v', '--vcf', required=True,
               help='vcf source path')
    parser.add('-f', '--filter', required=True,
               help='id filter path')
    parser.add('-c', '--column', type=int, default=3,
               help='column for info')
    parser.add('-t', '--type', default=int,
               help='info type int or float')
    parser.add('-i', '--info', required=True,
               help='add info name')
    parser.add('-s', '--save', required=True,
               help='prefix for save path')
    return parser


if __name__ == '__main__':

    parser = config_opts()
    opt = parser.parse_args()

    path_write_vcf = opt.vcf

    dict_pval = {}
    dict_fst = {}

    arry_filter = opt.filter.split(',')
    arry_info = opt.info.split(',')

    index_filter = 0
    while index_filter < len(arry_filter):
        with open(arry_filter[index_filter], 'r') as file_filter:
            list_ids = file_filter.read().strip().split('\n')
            for pinfo in list_ids[1:]:
                plist = pinfo.split('\t')
                if plist[0]+'_'+plist[1] not in dict_fst:
                    dict_fst[plist[0] + '_' + plist[1]] = {}
                if arry_info[index_filter] not in dict_fst[plist[0] + '_' + plist[1]]:
                    dict_fst[plist[0] + '_' + plist[1]][arry_info[index_filter]] = []
                dict_fst[plist[0] + '_' + plist[1]][arry_info[index_filter]].append(plist[opt.column - 1])

        index_filter += 1

    index_vcf = 1
    with open(opt.vcf, 'r') as file_read, open(opt.save, 'w') as file_write:
        while True:
            line = file_read.readline().strip()

            if line:
                if line.startswith('#'):
                    if line.startswith('#CHROM'):

                        index_filter = 0
                        while index_filter < len(arry_filter):
                            if opt.type == 'int':
                                file_write.writelines(
                                    '##INFO=<ID='+arry_info[index_filter]+',Number=1,Type=Integer,Description="'+arry_info[index_filter]+'">' + '\n')
                            elif opt.type == 'float':
                                file_write.writelines(
                                    '##INFO=<ID='+ arry_info[index_filter]+',Number=1,Type=Float,Description="'+arry_info[index_filter]+'">' + '\n')
                            index_filter += 1

                    file_write.writelines(line + '\n')
                else:

                    strs = line.split('\t')
                    id = strs[2]
                    chr = strs[0]
                    pos = strs[1]

                    info = strs[7].split(';')

                    # svtype = [x.split('=')[1] for x in info if 'SVTYPE=' in x][0]
                    # end = [x.split('=')[1] for x in info if 'END=' in x][0]

                    tag_tmp = chr + '_' + pos

                    index_filter = 0
                    while index_filter < len(arry_filter):
                        # if tag_tmp not in dict_fst or len(dict_fst[tag_tmp]) == 0:
                        if tag_tmp not in dict_fst or arry_info[index_filter] not in dict_fst[tag_tmp]:
                            strs[7] += ';'+arry_info[index_filter]+'=0'
                        else:
                            # print(tag_tmp, dict_fst[tag_tmp])
                            if 'nan' in dict_fst[tag_tmp][arry_info[index_filter]][0]:
                                strs[7] += ';' + arry_info[index_filter] + '=0'
                            else:
                                strs[7] += ';' + arry_info[index_filter] + '=' + dict_fst[tag_tmp][arry_info[index_filter]][0]
                            if len(dict_fst[tag_tmp][arry_info[index_filter]]) > 1:
                                dict_fst[tag_tmp][arry_info[index_filter]].pop(0)
                        index_filter += 1

                    file_write.writelines('\t'.join(strs) + '\n')

            else:
                break












