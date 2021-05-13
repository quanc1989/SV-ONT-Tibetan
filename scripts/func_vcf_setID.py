import configargparse


def config_opts():
    parser = configargparse.ArgumentParser(
        description='func_vcf_filter.py',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add('-v', '--vcf', required=True,
               help='vcf source path')
    parser.add('-c', '--chimp', required=False, action='store_true', default=False,
               help='if add chimp genotype from info')
    parser.add('-s', '--save', required=True,
               help='save path')
    return parser


if __name__ == '__main__':

    parser = config_opts()
    opt = parser.parse_args()

    dict_ids = {}
    list_ids = []


    with open(opt.vcf, 'r') as file_read, open(opt.save, 'w') as file_write:
        while True:
            line = file_read.readline().strip()

            if line:
                if line.startswith('##'):
                    file_write.writelines(line + '\n')
                elif line.startswith('#'):
                    if opt.chimp:
                        file_write.writelines(line + '\tCHIMP' +'\n')
                    else:
                        file_write.writelines(line + '\n')
                else:
                    strs = line.split('\t')

                    chrom = strs[0]
                    id = strs[2]
                    pos = strs[1]
                    ref = strs[3]
                    alt = strs[4]

                    if id == '.':
                        snpid = chrom + ':' + pos

                        if snpid not in dict_ids:
                            dict_ids[snpid] = 0

                        dict_ids[snpid] += 1
                        snpid = snpid + ':' + str(dict_ids[snpid])
                        strs[2] = snpid

                    if strs[2] in list_ids:
                        continue
                    else:
                        list_ids.append(strs[2])


                    if opt.chimp:
                        info = strs[7].split(';')
                        chimp = [x.split('=')[1] for x in info if 'CHIMP=' in x][0]
                        chimp = chimp.upper()
                        if chimp == ref.upper():
                            gp_chimp = '0|0'
                        elif chimp == alt.upper():
                            gp_chimp = '1|1'
                        else:
                            gp_chimp = '.|.'

                        file_write.writelines('\t'.join(strs) + '\t' + gp_chimp + '\n')
                    else:
                        file_write.writelines('\t'.join(strs) + '\n')

            else:
                break












