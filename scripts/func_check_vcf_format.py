import configargparse


def config_opts():
    parser = configargparse.ArgumentParser(
        description='func_vcf_filter.py',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add('-v', '--vcf', required=True,
               help='vcf source path')
    parser.add('-p', '--prefix', required=True,
               help='invalid vcf save prefix')
    parser.add('-s', '--seq', required=True,
               help='seq insertion file array')
    return parser


if __name__ == '__main__':

    parser = config_opts()
    opt = parser.parse_args()

    list_svtype = ['DEL', 'INS', 'DUP', 'INV', 'TRA']
    list_svtype_extra = []

    dict_seq = {}
    dict_pos = {}
    for path_seq in opt.seq.split(','):
        with open(path_seq) as file_read:
            array_seq = file_read.read().strip().split('\n')
            for line in array_seq:
                info = line.split('\t')
                chrom = info[0].split('_')[0]
                pos = int(info[0].split('_')[1])

                if chrom+'_'+str(pos) not in dict_seq:
                    dict_seq[chrom+'_'+str(pos)] = info[1]

                if chrom not in dict_pos:
                    dict_pos[chrom] = []

                if pos not in dict_pos[chrom]:
                    dict_pos[chrom].append(pos)

    for chrom in dict_pos:
        dict_pos[chrom].sort()

    dict_res = {}
    for svtype in list_svtype:
        dict_res[svtype] = {'num': 0, 'invalid': 0, 'corrected':0}

    dict_cache_ins = []
    dict_cache_dup = []
    dict_cache_inv = []
    dict_cache_dupins = []
    dict_cache_invdup = []

    with open(opt.vcf, 'r') as file_read, open(opt.prefix+'.correct.log', 'w') as file_write, open(opt.prefix+'.corrected.vcf', 'w') as file_write_correct:
        while True:
            line = file_read.readline().strip()
            if line:
                if line.startswith('#'):
                    file_write_correct.writelines(line + '\n')
                    continue
                elif line:
                    vcf_array = line.split('\t')

                    chrom = vcf_array[0]
                    pos = vcf_array[1]
                    svtype_symbol = vcf_array[4]

                    svtype = [x.split('=')[1] for x in vcf_array[7].split(';') if 'SVTYPE=' in x][0]
                    chr2 = [x.split('=')[1] for x in vcf_array[7].split(';') if 'CHR2=' in x][0]
                    array_seq = [x.split('=')[1] for x in vcf_array[7].split(';') if 'SEQ=' in x]
                    pos_end = int([x.split('=')[1] for x in vcf_array[7].split(';') if 'END=' in x][0])
                    svlen = float([x.split('=')[1] for x in vcf_array[7].split(';') if 'SVLEN=' in x][0])

                    if svtype not in list_svtype:
                        if svtype not in list_svtype_extra:
                            file_write.writelines('_'.join([chrom, pos, svtype_symbol])+ '\t not in svtype\n')
                            list_svtype_extra.append(svtype)
                    else:
                        dict_res[svtype]['num'] += 1

                        flag_multi = 0
                        tmp_arry = svtype_symbol.split(',')
                        # if len(tmp_arry) > 1:
                        #     flag_multi += 1
                            # print('-'.join([chrom, pos, svtype_symbol]))

                        for symbol in tmp_arry:
                            if symbol.__contains__('<'):
                                tmp_type = svtype_symbol.split('<')[1].split('>')[0]
                                if tmp_type not in list_svtype and tmp_type not in list_svtype_extra:
                                    list_svtype_extra.append(tmp_type)

                        # dict_res[svtype]['multi'] += flag_multi

                        if svtype == 'INS':
                            if len(array_seq) == 0 or \
                                    array_seq[0] == '':
                                dict_res[svtype]['invalid'] += 1
                                # file_write.writelines(line + '\n')

                                if int(pos) in dict_pos[chrom]:
                                    vcf_array[7] = vcf_array[7].replace('SEQ=', '').replace('SEQ', '') + ';SEQ=' + dict_seq[chrom + '_' + pos]
                                    dict_res['INS']['corrected'] += 1
                                else:
                                    # print(dict_pos[chrom])
                                    # dict_pos[chrom] - int(pos)

                                    arry_pos = [str(x) for x in dict_pos[chrom] if abs(x-int(pos)) < 1000]
                                    # if len(arry_pos) > 3:
                                    #     arry_pos = arry_pos[:3]
                                    # arry_seq = [dict_seq[chrom + '_' + x] for x in arry_pos]

                                    if len(arry_pos) > 0:
                                        vcf_array[7] = vcf_array[7].replace('SEQ=', '').replace('SEQ', '') + ';SEQ=' + \
                                                       dict_seq[chrom + '_' + arry_pos[0]]
                                        dict_res['INS']['corrected'] += 1

                                    else:
                                        file_write.writelines('_'.join([chrom, pos, svtype_symbol])+ '\t no sequence\n')


                            vcf_array[7] = vcf_array[7].replace('END=' + str(pos_end),
                                                                'END=' + str(int(pos) + 1))
                            file_write_correct.writelines('\t'.join(vcf_array) + '\n')

                            dict_cache_ins.append(chrom+'_'+pos)

                        elif svtype == 'DEL':
                            if pos_end - int(pos) < 50:
                                dict_res[svtype]['invalid'] += 1
                                # file_write.writelines(line + '\n')
                            # else:
                            file_write_correct.writelines(line + '\n')

                        elif svtype == 'DUP':
                            if pos_end - int(pos) < 50:
                                # file_write.writelines(line + '\n')
                                vcf_array[7] = vcf_array[7].replace('END=' + str(pos_end), 'END=' + str(int(pos) + int(svlen)))

                            if svtype_symbol.__contains__('INS'):
                                # if chrom + '_' + pos in dict_seq:
                                #     vcf_array[7] = vcf_array[7].replace('SEQ=', '').replace('SEQ', '') + ';SEQ=' + dict_seq[
                                #         chrom + '_' + pos]
                                #     vcf_array[7] = vcf_array[7].replace('END=' + str(pos_end),
                                #                                         'END=' + str(int(pos) + 1))
                                #     vcf_array[4] = '<INS>'
                                #     file_write_correct.writelines('\t'.join(vcf_array) + '\n')
                                # else:
                                # print('DUP/INS: ' + chrom + '_' + pos + '\n')
                                dict_res[svtype]['invalid'] += 1
                                dict_cache_dupins.append('\t'.join(vcf_array))
                            else:
                                dict_cache_dup.append(chrom+'_'+pos)
                                vcf_array[4] = '<DUP>'
                                file_write_correct.writelines('\t'.join(vcf_array) + '\n')

                        elif svtype == 'INV':
                            if pos_end - int(pos) < 50:
                                # file_write.writelines(line + '\n')
                                vcf_array[7] = vcf_array[7].replace('END='+str(pos_end), 'END='+str(int(pos) + 1))

                            if svtype_symbol.__contains__('DUP'):
                                # print('DUPINV: ' + chrom + '_' + pos + '_' + svtype_symbol +'\n')
                                dict_cache_invdup.append('\t'.join(vcf_array))
                                dict_res[svtype]['invalid'] += 1
                            else:
                                file_write_correct.writelines('\t'.join(vcf_array) + '\n')
                                dict_cache_inv.append(chrom+'_'+pos)

                        elif svtype == 'TRA':
                            file_write_correct.writelines(line + '\n')
                        #     if pos_end - int(pos) <= 1:
                        #         dict_res[svtype]['invalid'] += 1
                        #         file_write.writelines(line + '\n')
                        else:
                            print('error : ' + svtype)

            else:
                break

        for line in dict_cache_dupins:
            vcf_array = line.split('\t')

            chrom = vcf_array[0]
            pos = vcf_array[1]
            svtype_symbol = vcf_array[4]

            svtype = [x.split('=')[1] for x in vcf_array[7].split(';') if 'SVTYPE=' in x][0]
            chr2 = [x.split('=')[1] for x in vcf_array[7].split(';') if 'CHR2=' in x][0]
            array_seq = [x.split('=')[1] for x in vcf_array[7].split(';') if 'SEQ=' in x]
            pos_end = int([x.split('=')[1] for x in vcf_array[7].split(';') if 'END=' in x][0])
            svlen = float([x.split('=')[1] for x in vcf_array[7].split(';') if 'SVLEN=' in x][0])

            if chrom + '_' + pos in dict_cache_ins:
                if chrom + '_' + pos not in dict_cache_dup:
                    vcf_array[4] = '<DUP>'
                    dict_cache_dup.append(chrom + '_' + pos)
                    file_write_correct.writelines('\t'.join(vcf_array) + '\n')
                    dict_res['DUP']['corrected'] += 1
                else:
                    file_write.writelines('_'.join([chrom, pos, svtype_symbol])+ '\t both in ins and dup\n')

            elif chrom + '_' + pos in dict_cache_dup:
                if int(pos) in dict_pos[chrom]:
                    vcf_array[7] = vcf_array[7].replace('SEQ=', '').replace('SEQ', '') + ';SEQ=' + dict_seq[
                        chrom + '_' + pos]
                    dict_res['DUP']['corrected'] += 1
                    dict_cache_ins.append(chrom + '_' + pos)
                else:
                    # print(dict_pos[chrom])
                    # dict_pos[chrom] - int(pos)

                    arry_pos = [str(x) for x in dict_pos[chrom] if abs(x - int(pos)) < 1000]

                    if len(arry_pos) > 0:
                        vcf_array[7] = vcf_array[7].replace('SEQ=', '').replace('SEQ', '') + ';SEQ=' + \
                                       dict_seq[chrom + '_' + arry_pos[0]]
                        dict_res['DUP']['corrected'] += 1
                        dict_cache_ins.append(chrom + '_' + pos)
                    else:
                        file_write.writelines('_'.join([chrom, pos, svtype_symbol])+ '\t no sequence\n')

                vcf_array[4] = '<INS>'
                vcf_array[7] = vcf_array[7].replace('END=' + str(pos_end),
                                                    'END=' + str(int(pos) + 1))
                file_write_correct.writelines('\t'.join(vcf_array) + '\n')
            else:
                file_write.writelines('_'.join([chrom, pos, svtype_symbol]) + '\t split to DUP and INS\n')
                dict_res['DUP']['corrected'] += 1

                vcf_array[4] = '<DUP>'
                dict_cache_dup.append(chrom + '_' + pos)
                file_write_correct.writelines('\t'.join(vcf_array) + '\n')

                if int(pos) in dict_pos[chrom]:
                    vcf_array[7] = vcf_array[7].replace('SEQ=', '').replace('SEQ', '') + ';SEQ=' + dict_seq[
                        chrom + '_' + pos]
                    dict_cache_ins.append(chrom + '_' + pos)
                else:
                    # print(dict_pos[chrom])
                    # dict_pos[chrom] - int(pos)

                    arry_pos = [str(x) for x in dict_pos[chrom] if abs(x - int(pos)) < 1000]

                    if len(arry_pos) > 0:
                        vcf_array[7] = vcf_array[7].replace('SEQ=', '').replace('SEQ', '') + ';SEQ=' + \
                                       dict_seq[chrom + '_' + arry_pos[0]]
                        dict_cache_ins.append(chrom + '_' + pos)
                    else:
                        file_write.writelines('_'.join([chrom, pos, svtype_symbol])+ '\t no sequence\n')

                vcf_array[4] = '<INS>'
                vcf_array[7] = vcf_array[7].replace('END=' + str(pos_end),
                                                    'END=' + str(int(pos) + 1))
                file_write_correct.writelines('\t'.join(vcf_array) + '\n')

        for line in dict_cache_invdup:
            vcf_array = line.split('\t')

            chrom = vcf_array[0]
            pos = vcf_array[1]
            svtype_symbol = vcf_array[4]

            svtype = [x.split('=')[1] for x in vcf_array[7].split(';') if 'SVTYPE=' in x][0]
            chr2 = [x.split('=')[1] for x in vcf_array[7].split(';') if 'CHR2=' in x][0]
            array_seq = [x.split('=')[1] for x in vcf_array[7].split(';') if 'SEQ=' in x]
            pos_end = int([x.split('=')[1] for x in vcf_array[7].split(';') if 'END=' in x][0])
            svlen = float([x.split('=')[1] for x in vcf_array[7].split(';') if 'SVLEN=' in x][0])

            if chrom + '_' + pos in dict_cache_inv:
                if chrom + '_' + pos not in dict_cache_dup:
                    vcf_array[4] = '<DUP>'
                    dict_cache_dup.append(chrom + '_' + pos)
                    file_write_correct.writelines('\t'.join(vcf_array) + '\n')
                    dict_res['INV']['corrected'] += 1
                else:
                    file_write.writelines('_'.join([chrom, pos, svtype_symbol])+ '\t both in inv and dup\n')

            elif chrom + '_' + pos in dict_cache_dup:
                vcf_array[4] = '<INV>'
                dict_cache_inv.append(chrom + '_' + pos)
                file_write_correct.writelines('\t'.join(vcf_array) + '\n')
                dict_res['INV']['corrected'] += 1
            else:
                file_write.writelines('_'.join([chrom, pos, svtype_symbol])+ '\t split into INV and DUP\n')
                dict_res['INV']['corrected'] += 1

                vcf_array[4] = '<DUP>'
                dict_cache_dup.append(chrom + '_' + pos)
                file_write_correct.writelines('\t'.join(vcf_array) + '\n')

                vcf_array[4] = '<INV>'
                dict_cache_inv.append(chrom + '_' + pos)
                file_write_correct.writelines('\t'.join(vcf_array) + '\n')


        for svtype in list_svtype:
            print(svtype + '\n')
            for key, val in dict_res[svtype].items():
                print(key + ' : ' + str(val) + '\n')

        print('EXtra types: \n')
        for svtype in list_svtype_extra:
            print(svtype + '\t')








