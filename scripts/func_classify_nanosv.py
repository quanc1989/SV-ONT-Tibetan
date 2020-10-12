#!/usr/bin/env python
import os
import re
import sys


def classify(line, ALT_INDEX):
    #get alt, chrom1, chrom2, position (pos), id, old SVTYPE (should be BND if virgin svaba vcf) from line
    s = line.split("\t")
    alt = s[ALT_INDEX]
    chrom1 = s[0].lower()
    pos = int(s[1])
    id=s[2]

    oldType = [ x.split('=')[1] for x in line.split(';') if 'SVTYPE' in x][0]

    # print('\n' + oldType + '\n')

    # get new type
    if oldType == 'BND':

        mateChrom = 'chr' + alt.lower().split(':')[0].split('chr')[1]

        # print(chrom1, mateChrom)

        if chrom1 == mateChrom:
            INV_PATTERN_1 = re.compile(r'\D\].+\]')
            INV_PATTERN_2 = re.compile(r'\[.+\[\D')
            if INV_PATTERN_1.match(alt) or INV_PATTERN_2.match(alt):
                return "INV"
            # DEL
            DEL_PATTERN_THIS = re.compile(r'\D\[.+\[')
            if DEL_PATTERN_THIS.match(alt):
                return "DEL"
            # DUP
            INS_PATTERN_THIS = re.compile(r'\].+\]\D')
            if INS_PATTERN_THIS.match(alt):
                # return "DUP/INS"
                return "DUP"
        else:
            return 'BND'
    else:
        return oldType
    # return "INS"

if __name__ == "__main__":
    file = sys.argv[1]
    if not os.path.exists(file):
        raise IOError(file)
    alt_index = -1
    #generate mate:mate dictionary
    #load file into ram
    vcf_file=[]
    with open (file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            vcf_file.append(line)
    with open(file, "r") as f:
        for line in f:
            # print comments
            if line.startswith("##"):
                sys.stdout.write(line)
                continue
            # header contains indexes
            if line.startswith('#'):
                split = line.split("\t")
                for index, val in enumerate(split):
                    if val == "ALT":
                        alt_index = index
                        break
                sys.stdout.write(line)
                continue
            if alt_index == -1:
                print("ERROR: NO ALT INDEX FOUND")
                exit(1)
            newType = classify(line, alt_index)
            if newType != "NONE":
                newLine = re.sub(r'SVTYPE=BND',"SVTYPE="+newType,line)
                sys.stdout.write(newLine)