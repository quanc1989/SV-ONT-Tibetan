import os
import configargparse
import os

import configargparse


def config_opts():
    parser = configargparse.ArgumentParser(
        description='run_inference_summary.py',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add('-p', '--path', required=True,
               help='path for result of inference')
    parser.add('-s', '--save', required=True,
               help='save path')
    return parser


if __name__ == '__main__':

    parser = config_opts()
    opt = parser.parse_args()

    list_result = []
    list_params = []
    list_last_line = []

    for fname in os.listdir(opt.path):
        if fname.endswith('.optimized.txt'):
            with open(os.path.join(opt.path, fname), 'r') as file_read:
                lines = file_read.read().strip().split('\n')
                max_ll = 0
                tmp_line = ''
                tmp_params = ''
                last_line = ''
                for line in lines:
                    if line.startswith('Model'):
                        if tmp_params == '':
                            tmp_params = line.split('\t')[6].split('(')[1].split(')')[0]
                            list_params.append(tmp_params)
                    else:
                        ll = float(line.split('\t')[2])
                        last_line = line
                        if ll > max_ll or max_ll == 0:
                            max_ll = ll
                            tmp_line = line

                list_result.append(tmp_line)
                list_last_line.append(last_line)

    with open(opt.save, 'w') as file_write:
        index_models = 0
        while index_models < len(list_result):
            file_write.writelines(list_result[index_models] + '\t' + list_params[index_models] + '\n')
            index_models += 1

    with open(opt.save + '.last', 'w') as file_write:
        index_models = 0
        while index_models < len(list_last_line):
            file_write.writelines(list_last_line[index_models] + '\n')
            index_models += 1
