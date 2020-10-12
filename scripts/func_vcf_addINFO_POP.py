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
    parser.add('--file_samples', required=True,
               help='file_samples')
    return parser

def init_info_dict(list_groups):
    res = {}
    for group in list_groups:
        res[group] = 0.0
    return res

def clear(dict_info):
    for group in dict_info:
        dict_info[group] = 0.0

if __name__ == '__main__':

    parser = config_opts()
    opt = parser.parse_args()

    list_svid = []
    num_samples = 0
    dict_samples = {}
    list_groups = []

    with open(opt.file_samples, 'r') as file_filter:
        list_samples = file_filter.read().strip().split('\n')
        for line in list_samples:
            info = line.split('\t')
            if info[0] in dict_samples:
                print('ERROR : duplicated sample ' + info[0])
            else:
                dict_samples[info[0]] = info[1]
                if info[1] not in list_groups:
                    list_groups.append(info[1])

    with open(opt.vcf, 'r') as file_read, open(opt.save, 'w') as file_write:
        lines = file_read.read().strip().split('\n')
        for line in lines:

            if line.startswith('#'):
                if line.startswith('#CHROM'):

                    for group in list_groups:
                        file_write.writelines(
                            '##INFO=<ID=SUPP_'+group+'_NGS,Number=1,Type=Integer,Description="NGS supp '+group+'">' + '\n')
                        file_write.writelines(
                            '##INFO=<ID=SUPP_RATE_'+group+'_NGS,Number=1,Type=Float,Description="NGS supp rate '+group+'">' + '\n')
                        file_write.writelines(
                            '##INFO=<ID=AF_'+group+'_NGS,Number=1,Type=Float,Description="NGS '+group+' sv frequency">' + '\n')
                        file_write.writelines(
                            '##INFO=<ID=DP_'+group+'_NGS,Number=1,Type=Integer,Description="NGS DP '+group+'">' + '\n')
                        file_write.writelines(
                            '##INFO=<ID=MR_'+group+'_NGS,Number=1,Type=Float,Description="NGS Missing Rate for '+group+'">' + '\n')

                    infolist = line.split('\t')[:9]
                    samplelist = line.split('\t')[9:]
                    num_samples = len(samplelist)

                file_write.writelines(line+'\n')
            else:
                repline = line

                infolist = repline.split('\t')
                info = infolist[7]

                info_supp = init_info_dict(list_groups)
                info_num_samples = init_info_dict(list_groups)
                info_af = init_info_dict(list_groups)
                info_dp = init_info_dict(list_groups)
                info_missing = init_info_dict(list_groups)
                info_supp_rate = init_info_dict(list_groups)

                missing = 0
                supp_rewrite = 0

                for index_sample in range(num_samples):

                    label_sample = samplelist[index_sample]

                    info_sample = infolist[index_sample+9].split(':')

                    dp_sample = int(info_sample[1])

                    if label_sample in dict_samples:

                        group_sample = dict_samples[label_sample]
                        info_dp[group_sample] += dp_sample
                        info_num_samples[group_sample] += 1

                        if info_sample[0] == './.' or info_sample[0] == '.':
                            info_sample[0] = './.'

                            info_missing[group_sample] += 1

                        if info_sample[0] == '0/1' or info_sample[0] == '1/1':

                            supp_rewrite += 1

                            info_supp[group_sample] += 1
                            info_af[group_sample] += sum([int(x) for x in info_sample[0].split('/')])

                for subgroup in info_missing:
                    missing += info_missing[subgroup]

                if num_samples == missing:

                    supp_rate = 0
                    clear(info_supp_rate)

                    af_all = 0
                    clear(info_af)

                    dp_all = 0
                    clear(info_dp)

                else:

                    for subgroup in list_groups:
                        if info_missing[subgroup] == info_num_samples[subgroup]:
                            info_supp_rate[subgroup] = 0
                            info_af[subgroup] = 0
                            info_dp[subgroup] = 0
                        else:
                            info_supp_rate[subgroup] = info_supp[subgroup] * 1.0 / (info_num_samples[subgroup] - info_missing[subgroup])
                            info_af[subgroup] = info_af[subgroup] * 1.0 / ((info_num_samples[subgroup] - info_missing[subgroup]) * 2)
                            info_dp[subgroup] = info_dp[subgroup] * 1.0 / (info_num_samples[subgroup] - info_missing[subgroup])

                missing = missing * 1.0 / num_samples
                for subgroup in list_groups:
                    info_missing[subgroup] = info_missing[subgroup] * 1.0/ info_num_samples[subgroup]

                for group in list_groups:
                    infolist[7] = infolist[7] + (';SUPP_'+group+'_NGS=%d;SUPP_RATE_'+group+'_NGS=%f;AF_'+group+'_NGS=%f;DP_'+group+'_NGS=%d;MR_'+group+'_NGS=%f') % (info_supp[group], info_supp_rate[group], info_af[group], info_dp[group], info_missing[group])

                # infolist[7] = infolist[7][:-1]

                # infolist[8] = 'GT:DR:DV'
                file_write.writelines('\t'.join(infolist) + '\n')
