import configargparse


def config_opts():
    parser = configargparse.ArgumentParser(
        description='func_trans_vcf_format_INDELtoSYMBOL.py',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add('-v', '--vcf', required=True,
               help='vcf path')
    parser.add('-s', '--save', required=True,
               help='save path')
    return parser


if __name__ == '__main__':

    parser = config_opts()
    opt = parser.parse_args()

    with open(opt.vcf, 'r') as file_read, open(opt.save, 'w') as file_write:
        lines = file_read.read().strip().split('\n')
        for line in lines:

            if line.startswith('#'):
                file_write.writelines(line+'\n')
            else:
                repline = line

                infolist = repline.split('\t')
                info = infolist[7]

                svtype = [x.split('=')[1] for x in info.split(';') if 'SVTYPE=' in x][0]

                if svtype.__contains__('INS'):
                    if info.__contains__('SEQS'):
                        seq = [x.split('=')[1] for x in info.split(';') if 'SEQS=' in x][0]
                        info = info.replace('SEQS=' + seq + ';', '')
                    else:
                        seq = [x.split('=')[1] for x in info.split(';') if 'SEQ=' in x][0]
                        info = info.replace('SEQ=' + seq + ';', '')

                    infolist[7] = info
                    infolist[3] = 'N'
                    infolist[4] = seq

                file_write.writelines('\t'.join(infolist) + '\n')
