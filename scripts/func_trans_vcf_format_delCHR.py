import configargparse


def config_opts():
    parser = configargparse.ArgumentParser(
        description='func_selectColFromBed.py',
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
                if 'chr' in line:
                    file_write.writelines(line.replace('chr', '') + '\n')
                else:
                    file_write.writelines(line+'\n')
            else:
                repline = line

                if line.startswith('chr'):
                    repline = line.replace('chr', '')

                file_write.writelines(repline + '\n')

