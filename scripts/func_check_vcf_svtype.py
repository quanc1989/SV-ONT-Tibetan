import configargparse


def config_opts():
    parser = configargparse.ArgumentParser(
        description='func_vcf_filter.py',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add('-v', '--vcf', required=True,
               help='vcf source path')
    parser.add('-s', '--save', required=True,
               help='vcf save path')
    return parser

def throwERR(label, err_type):
    print('error -> ', label, err_type)

if __name__ == '__main__':

    parser = config_opts()
    opt = parser.parse_args()

    list_svtype = ['DEL', 'INS', 'DUP', 'INV', 'BND', 'TRA']

    with open(opt.vcf, 'r') as file_read, open(opt.save, 'w') as file_write:
        while True:
            line = file_read.readline().strip()
            if line:
                if line.startswith('#'):
                    file_write.writelines(line + '\n')
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
                    svlen = [x.split('=')[1] for x in vcf_array[7].split(';') if 'SVLEN=' in x][0]

                    if float(svlen) < 0:
                        vcf_array[7] = vcf_array[7].replace('SVLEN=' + svlen, 'SVLEN=' + str(int(-float(svlen))))
                        svlen = str(int(-float(svlen)))

                    tmp_svid = chrom + '_' + pos + '_' + str(pos_end) + '_' + svtype
                    if svtype not in list_svtype:
                        throwERR(tmp_svid, 'not in svtypelist')
                        break
                    else:

                        if svtype not in list_svtype or not svtype_symbol.__contains__(svtype):
                            vcf_array[7] = vcf_array[7].replace('SVTYPE=' + svtype, 'SVTYPE='+svtype_symbol.split('<')[1].split('>')[0])
                            svtype = svtype_symbol.split('<')[1].split('>')[0]

                        if svtype == 'INS':
                            if len(array_seq) == 0 or \
                                    array_seq[0] == '':
                                throwERR(tmp_svid, 'invalid INS seq NO SEQ')
                                # break

                            if pos_end - int(pos) != 1:
                                throwERR(tmp_svid, 'END-POS != 1')
                                break

                        elif svtype == 'DEL':
                            if abs(float(svlen)) < 0 or str(int(int(pos) + float(svlen))) != str(pos_end):
                                throwERR(tmp_svid, 'SVLEN change ' + svlen + ' to ' + str(abs(pos_end - int(pos))))
                                vcf_array[7] = vcf_array[7].replace('SVLEN=' + svlen, 'SVLEN=' + str(abs(pos_end - int(pos))))

                        elif svtype == 'DUP':

                            if pos_end - int(pos) < 50:
                                throwERR(tmp_svid,
                                         'change end ' + str(pos_end) + ' -> ' + str(int(int(pos) + float(svlen))))
                                if abs(float(svlen)) > 50:
                                    vcf_array[7] = vcf_array[7].replace('END=' + str(pos_end),
                                                                    'END=' + str(int(int(pos) + float(svlen))))
                                    vcf_array[2] = chrom + '_' + pos + '_' +str(int(int(pos) + float(svlen))) + '_' + 'DUP'
                                else:
                                    break
                            vcf_array[7] = vcf_array[7].replace('SVLEN=' + svlen,
                                                                'SVLEN=' + str(abs(pos_end - int(pos))))

                        elif svtype == 'INV':
                            if pos_end - int(pos) < 50:
                                throwERR(tmp_svid, 'change end ' + str(pos_end) + ' -> ' + str(int(int(pos) + float(svlen))))
                                if abs(float(svlen)) > 50:
                                    vcf_array[7] = vcf_array[7].replace('END=' + str(pos_end),
                                                                        'END=' + str(int(int(pos) + float(svlen))))
                                    vcf_array[2] = chrom + '_' + pos + '_' + str(int(int(pos) + float(svlen))) + '_' + 'INV'
                                else:
                                    break
                                # break
                            vcf_array[7] = vcf_array[7].replace('SVLEN=' + svlen,
                                                                'SVLEN=' + str(abs(pos_end - int(pos))))

                        elif svtype == 'BND' or svtype == 'TRA':
                            vcf_array[4] = '<BND>'

                            chr2 = [x.split('=')[1] for x in vcf_array[7].split(';') if 'CHR2=' in x][0]
                            vcf_array[2] = chrom + '_' + pos + '_' + chr2 + '_' + str(pos_end) + '_' + 'BND'
                            vcf_array[7] = vcf_array[7].replace('SVTYPE=TRA',
                                                                'SVTYPE=BND')
                    file_write.writelines('\t'.join(vcf_array) + '\n')

            else:
                break








