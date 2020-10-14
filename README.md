# SV-ONT-Tibetan
This repository includes scripts to analyze structural variations of 25 Chinese samples using nanopore sequencing. 

## Pipeline for Multi-sample SV-calling and annotation.

### Requirements
- [NGMLR](https://github.com/philres/ngmlr)
- [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
- [SVIM](https://github.com/eldariont/svim)
- [NanoSV](https://github.com/mroosmalen/nanosv)
- [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)
- [SVTK](https://github.com/talkowski-lab/svtk)
- [ANNOTSV](https://github.com/lgmgeo/AnnotSV)
- python3
- numpy

### Summary
In the bash file ```pipeline.sv-calling.sh```, we use sample data to demonstrate the complete process of detecting and annotating structural variations based on nanopore sequencing technology.  

1. Firstly, long reads were mapped to GRCh37 human reference from NCBI without alternate sequences. Mapping was performed with NGMLR with ONT default parameters. 

2. Then SV calling was performed on each sample using Sniffles, NanoSV, and SVIM. These tools have been reported to be compatible with NGMLR and show better accuracy and sensitivity than others. Five minimum supporting reads with at least 50 bp length was required. The insertion sequence and read ID was required for each method, and the rest are all default parameters. SURVIVOR was used to merge the SVs, which are supported by at least two methods with a maximum allowed pairwise distance of 1,000 bp between breakpoints. Meanwhile, SVs obtained from different tools are not necessary to agree on the SV-type or the strand, so that we could capture as many potential breakpoints as possible. Finally, we merged the SVs obtained from all the samples as long as one sample supports it. 

3. So far, we have obtained potential regions for all samples, and we need to get a fully genotyped multi-samples dataset. We re-ran Sniffles across all the samples with all these potential regions (--Ivcf) and finally combined SVs with SURVIVOR. This time, we asked SURVIVOR only to report SVs supported by at least one sample, and they have to agree on the SV-type. Furthermore, we used a hard threshold with five minimum supporting reads, and all non-missing genotypes less than this threshold were modified to reference. 

    ![模型示意图](pipeline-sv-calling.png)