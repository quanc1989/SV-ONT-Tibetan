import configargparse


def config_opts():
    parser = configargparse.ArgumentParser(
        description='func_vcf_filter.py',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add('-v', '--vcf', required=True,
               help='vcf source path')
    parser.add('-s', '--save', required=True,
               help='save path')
    return parser


if __name__ == '__main__':

    parser = config_opts()
    opt = parser.parse_args()

    list_chrom = []
    for x in range(23):
        list_chrom.append(str(x+1))
        list_chrom.append('chr' + str(x+1))
    list_chrom.append('x')
    list_chrom.append('chrx')

    with open(opt.vcf, 'r') as file_read, open(opt.save, 'w') as file_write:
        while True:
            line = file_read.readline().strip()

            if line:
                if line.startswith('#'):
                    file_write.writelines(line+'\n')
                elif line:
                    vcf_array = line.split('\t')
                    chrom = vcf_array[0]

                    # if not chrom.lower().__contains__('chr'):
                    #     chrom =
                    svtype = [x.split('=')[1] for x in vcf_array[7].split(';') if 'SVTYPE=' in x][0]
                    chr2 = [x.split('=')[1] for x in vcf_array[7].split(';') if 'CHR2=' in x]
                    alt = vcf_array[4]

                    if chrom.lower() in list_chrom:
                        if len(chr2) > 0:
                            chr2 = chr2[0]
                            if chr2.lower() not in list_chrom:
                                continue
                        elif alt.__contains__('[') or alt.__contains__(']'):
                            tmp = 'chr' + alt.lower().split(':')[0]

                            if tmp.__contains__('['):
                                chr2 = tmp.split('[')[1]
                            else:
                                chr2 = tmp.split(']')[1]

                            if chr2 not in list_chrom:
                                continue
                        else:
                            chr2='NA'

                        vcf_array[0] = chrom.replace('CHR', 'chr')
                        vcf_array[4] = alt.replace('CHR','chr')

                        if chr2 != 'NA':
                            vcf_array[7] = vcf_array[7].replace('CHR2='+chr2,'CHR2='+chr2.replace('CHR','chr'))

                        file_write.writelines('\t'.join(vcf_array)+'\n')

            else:
                break










