import json
import math
import os
import random

import configargparse

import dadi


def config_opts():
    parser = configargparse.ArgumentParser(
        description='run_convert_vcf_to_dadi.py',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add('-d', '--data', required=True,
               help='datadict file')
    parser.add('-p', '--prefix', required=True,
               help='save prefix')
    parser.add('-l', '--poplist', required=True,
               help='population string')
    parser.add('--projections', required=True,
               help='projections down, same size to poplist')
    parser.add('-w', '--window', type=int, required=False,
              help='window width')
    parser.add('-u', '--unfold', required=False, action='store_true', default= False,
              help='whether it is unfold')
    parser.add('-r', '--random', required=False, action='store_true', default= False,
              help='generate random bootstrap sets')
    return parser

if __name__ == '__main__':

    parser = config_opts()
    opt = parser.parse_args()

    datadict = {}
    total_variants_chrom = {}

    with open(opt.data, 'r') as file_read:
        lines = file_read.read().strip().split('\n')
        # limit_variants = 10000
        # index_variant = 0

        for line in lines:
            # print(line)
            key_line = line.split('{')[0]
            val_line = line[line.index('{'):]

            # print(key_line,'--',val_line)
            info = key_line.split('-')
            chrom = info[0]
            pos = int(info[1])

            if chrom not in datadict:
                datadict[chrom] = {}
                total_variants_chrom[chrom] = 0

            if opt.window:
                order_pos = int(pos / opt.window)
                if order_pos not in datadict[chrom]:
                    datadict[chrom][order_pos] = {}

                datadict[chrom][order_pos][key_line] = json.loads(val_line.replace("\'", "\"").replace("(", "[").replace(")", "]"))
            else:
                datadict[chrom][key_line] = json.loads(val_line.replace("\'", "\"").replace("(", "[").replace(")", "]"))
                total_variants_chrom[chrom] += 1

    if not os.path.exists(opt.prefix):
        os.makedirs(opt.prefix)

    # list_pops = ['TIB', 'HANN', 'HANS']
    list_pops = [x.strip() for x in opt.poplist.split(',')]
    list_pojections = [int(x.strip()) for x in opt.projections.split(',')]

    if not opt.window: datadict_bootstrap={}

    with open(os.path.join(opt.prefix, 'generate_sfs_segments.log'), 'w') as file_write:
        for chrom in datadict:
            if opt.window:
                print('Simulation #%s start' % chrom)
                for order_pos in datadict[chrom]:
                    fs = dadi.Spectrum.from_data_dict(data_dict=datadict[chrom][order_pos],
                                                      pop_ids=list_pops,
                                                      projections=list_pojections,
                                                      polarized=opt.unfold)
                    if not os.path.exists(os.path.join(opt.prefix, chrom)):
                        os.makedirs(os.path.join(opt.prefix, chrom))
                    fs.to_file(os.path.join(opt.prefix, chrom, chrom + '-' + str(order_pos) + '.fs'))
                    file_write.writelines('%s\t%d\t%d\n' %(chrom, order_pos, len(datadict[chrom][order_pos])))
            else:
                window = total_variants_chrom[chrom]/100
                list_keys = list(datadict[chrom].keys())

                if opt.random:
                    random.shuffle(list_keys)

                index_variant = 0
                for key_variant in list_keys:

                    # pos = int(key_line.split('-')[1])
                    index_variant += 1
                    order_pos = math.floor(index_variant / window)

                    if order_pos not in datadict_bootstrap:
                        datadict_bootstrap[order_pos] = {}

                    # print(key_variant,'--',datadict[chrom][key_variant],'--',str(order_pos))
                    datadict_bootstrap[order_pos][key_variant] = datadict[chrom][key_variant]

        if not opt.window:

            index_bootstrap = 0
            list_keys = list(datadict_bootstrap.keys())

            while index_bootstrap < 100:

                index_bootstrap += 1
                print('Bootstrap #%d start' % index_bootstrap)

                list_keys_bootstrap = random.choices(list_keys, k=100)
                set_keys_bootstrap = set(list_keys_bootstrap)

                dict_to_spectrum = {}

                for order_pos in set_keys_bootstrap:
                    group_order_pos = list_keys_bootstrap.count(order_pos)
                    for index_group in range(group_order_pos):
                        for key_variant in datadict_bootstrap[order_pos]:
                            dict_to_spectrum[key_variant + '-' + str(index_group)] = datadict_bootstrap[order_pos][key_variant]

                with open(os.path.join(opt.prefix, "datadict."+str(index_bootstrap)+".txt"), 'w') as outfile:
                    for x, y in dict_to_spectrum.items():
                        outfile.write(x + str(y) + "\n")

                fs = dadi.Spectrum.from_data_dict(data_dict=dict_to_spectrum,
                                                  pop_ids=list_pops,
                                                  projections=list_pojections,
                                                  polarized=opt.unfold)
                fs.to_file(os.path.join(opt.prefix, str(index_bootstrap) + '.fs'))
                file_write.writelines('%d\t%d\n' %(index_bootstrap, len(dict_to_spectrum)))
                print(str(len(dict_to_spectrum)))





