### 1. mapping ###

path_data='0_raw_data/ONT'
path_save='1_mapping'
path_ref='0_raw_data/hg19.fa'

arrray_sample=('CQ132' 'CQ144')
suffix_fastq='.pass.fastq'
num_thread=60

if [ ! -d $path_save ];then
	mkdir $path_save
fi

for sampleID in ${arrray_sample[@]};do
	ngmlr -t $num_thread -r $path_ref \
	-q $path_data/$sampleID$suffix_fastq -o $path_save/$sampleID'.sam' -x ont
	samtools view -@ $num_thread -S -b $path_save/$sampleID'.sam' \
	| samtools sort -@ $num_thread -o $path_save/$sampleID'.bam'
	samtools index -@ $num_thread $path_save/$sampleID'.bam'
done

### 2. sv calling ###

path_data='1_mapping'
path_save='2_sv_calling'

command_sniffles='sniffles'

command_nanosv='NanoSV.py'
path_scripts='scripts'
path_example='example'

command_svim='svim'

command_SURVIVOR='SURVIVOR'

minimum_support_read=5
minimum_length=50
num_thread=60

if [ ! -d $path_save ];then
	mkdir $path_save
	mkdir $path_save/sniffles
	mkdir $path_save/nanosv
	mkdir $path_save/svim

	mkdir $path_save/result_denovo
	mkdir $path_save/result_genotyped
fi

for sampleID in ${arrray_sample[@]};do
	
	#sniffles
	$command_sniffles -s $minimum_support_read -t $num_thread \
	-m $path_data/$sampleID'.bam'\
	-v $path_save/sniffles/$sampleID'.sniffles.vcf' \
	--genotype --report_BND -n -1 -l $minimum_length
	bcftools sort -o $path_save/sniffles/$sampleID'.sniffles.sorted.vcf' \
	-O v $path_save/sniffles/$sampleID'.sniffles.vcf';
	python $path_scripts/func_vcf_filter_byCHR.py -v $path_save/sniffles/$sampleID'.sniffles.sorted.vcf' \
	-s $path_save/sniffles/$sampleID'.sniffles.sorted.primary.vcf'

	#nanosv
	python $command_nanosv -t $num_thread -o $path_save/nanosv/$sampleID'.nanosv.vcf' \
	-c $path_example/config_nanosv.ini -s samtools $path_data/$sampleID'.bam'
	vcftools --vcf $path_save/nanosv/$sampleID'.nanosv.vcf' --remove-filtered LowQual \
	--recode --recode-INFO-all --stdout > $path_save/nanosv/$sampleID'.nanosv.q20.vcf'
  bcftools sort -o $path_save/nanosv/$sampleID'.nanosv.q20.sorted.vcf' \
  -O v $path_save/nanosv/$sampleID'.nanosv.q20.vcf'
  python $path_scripts/func_classify_nanosv.py $path_save/nanosv/$sampleID'.nanosv.q20.sorted.vcf' > $path_save/nanosv/$sampleID'.nanosv.q20.sorted.classified.vcf'
  python $path_scripts/func_vcf_filter_byCHR.py \
  -v $path_save/nanosv/$sampleID'.nanosv.q20.sorted.classified.vcf' \
  -s $path_save/nanosv/$sampleID'.nanosv.q20.sorted.classified.primary.vcf'

	#svim
	$command_svim alignment --min_sv_size $minimum_length \
	--sample $sampleID --insertion_sequence \
	--read_names $path_save/svim/$sampleID $path_data/$sampleID'.bam' $path_ref
	awk '{ if($1 ~ /^#/) { print $0 } else { if($6>=5) { print $0 } } }' $path_save/svim/$sampleID/final_results.vcf > $path_save/svim/$sampleID/$sampleID'.q5.vcf'
	bcftools sort -o $path_save/svim/$sampleID/$sampleID'.q5.sorted.vcf' \
	-O v $path_save/svim/$sampleID/$sampleID'.q5.vcf'
  python $path_scripts/func_vcf_filter_byCHR.py -v $path_save/svim/$sampleID/$sampleID'.q5.sorted.vcf' \
  -s $path_save/svim/$sampleID/$sampleID'.q5.sorted.primary.vcf'

  #merge
  echo -e $path_save/sniffles/$sampleID'.sniffles.sorted.primary.vcf\n'$path_save/nanosv/$sampleID'.nanosv.q20.sorted.classified.primary.vcf\n'$path_save/svim/$sampleID/$sampleID'.q5.sorted.primary.vcf' > $path_save/$sampleID'.files'
  $command_SURVIVOR merge $path_save/$sampleID'.files' 1000 2 -1 0 0 50 $path_save/$sampleID'.vcf'
done

ls $path_save/*.vcf > $path_save/merge.denovo.files
$command_SURVIVOR merge $path_save/merge.denovo.files 1000 1 -1 0 0 50 $path_save/result_denovo/merge.denovo.vcf
bcftools sort -o $path_save/result_denovo/merge.denovo.sorted.vcf -O v $path_save/result_denovo/merge.denovo.vcf

for sampleID in ${arrray_sample[@]};do
	$command_sniffles -s $minimum_support_read -t $num_thread \
	-m $path_data/$sampleID'.bam' -v $path_save/$sampleID'.genotyped.vcf' \
	--genotype --report_BND -n -1 -l 50 --Ivcf $path_save/result_denovo/merge.denovo.sorted.vcf
	bcftools sort -o $path_save/$sampleID'.genotyped.sorted.vcf' \
	-O v $path_save/$sampleID'.genotyped.vcf';
	rm $path_save/$sampleID'.genotyped.vcf'

	echo $sampleID >> $path_save/vcf_samples.txt
	echo $path_save/$sampleID'.genotyped.vcf'
done

ls $path_save/*.genotyped.sorted.vcf > $path_save/merge.genotyped.files
$command_SURVIVOR $path_save/merge.genotyped.files 1000 -1 1 0 0 50 $path_save/result_genotyped/merge.genotyped.vcf
bcftools sort -o $path_save/result_genotyped/merge.genotyped.sorted.vcf -O v $path_save/result_genotyped/merge.genotyped.vcf


### 3. stats ###
path_data='2_sv_calling'
path_save='2_sv_calling'

python $path_scripts/func_select_SEQ_and_Rnames_fromVCFs.py -f $path_data/vcf_samples.txt -s $path_save/merge.files
python $path_scripts/func_select_SEQ_fromVCFs.py -f $path_data/merge.genotyped.files -s $path_save/merge.genotyped.files.seq

python $path_scripts/func_check_vcf_format.py \
-v $path_save/result_genotyped/merge.genotyped.sorted.vcf \
-p $path_save/result_genotyped/merge.genotyped.sorted \
-s $path_save/merge.genotyped.files.seq,$path_save/merge.files.seq

# trans to elite format
python $path_scripts/func_trans_vcf_format_INDELtoSYMBOL.py \
-v $path_save/result_genotyped/merge.genotyped.sorted.corrected.vcf \
-s $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.vcf
bcftools sort -o $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted.vcf \
-O v $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.vcf

python $path_scripts/func_trans_vcf_format_INDELtoSYMBOL.py \
-v $path_save/result_genotyped/merge.genotyped.sorted.corrected.vcf \
-s $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.elite.vcf -e
grep -v -E 'SVTYPE=TRA' $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.elite.vcf > $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.elite.excludeTRA.vcf

vcftools --vcf $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted.vcf \
--missing-site --out $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted
vcftools --vcf $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted.vcf \
--weir-fst-pop $path_example/vcf_samples_TIB_ordered.txt \
--weir-fst-pop $path_example/vcf_samples_HAN.txt \
--out $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted

python $path_scripts/func_vcf_addINFO.py \
-v $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted.vcf \
-f $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted.lmiss \
-c 6 -t float -i MR -d 'Missing Rate' -s $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted.addMR.vcf
python $path_scripts/func_vcf_addINFO.py -v $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted.addMR.vcf \
-f $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted.weir.fst \
-c 3 -t float -i FST -d 'Fst' -s $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted.addMR.addFST.vcf
bgzip -c $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted.addMR.addFST.vcf > $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted.addMR.addFST.vcf.gz
bcftools index $path_save/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted.addMR.addFST.vcf.gz

### 4. genotyping ###

path_data='0_raw_data/NGS'
path_save='3_genotyping_NGS'
path_ref='0_raw_data/hg19.fa'
path_scripts='scripts'
path_example='example'
path_sv_calling='2_sv_calling'

command_paragraph='python multigrmpy.py'

arrray_sample=('WGC025244D' 'WGC025246D')
suffix_bam='_dedup_realign_recal.bam'

mean_coverage=40
max_read=800
length_pair_read=150

num_thread=60


if [ ! -d $path_save ];then
	mkdir $path_save
fi

for sampleID in ${arrray_sample[@]};do
	
	mkdir $path_save/$sampleID
	echo -e 'id,path,depth,read length\n'$sampleID','$path_data/$sampleID$suffix_bam','$mean_coverage','$length_pair_read > $path_save/sampleID/sampleID.txt
	$command_paragraph -i $path_sv_calling/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.elite.excludeTRA.vcf \
	-m $path_save/sampleID/sampleID.txt -o $path_save/sampleID -r $path_ref -t $num_thread -M $max_read

	bcftools sort -o $path_save/sampleID/genotypes.sorted.vcf.gz \
	-O z $path_save/sampleID/genotypes.vcf.gz
	tabix -p vcf $path_save/sampleID/genotypes.sorted.vcf.gz

  bcftools filter -i 'FILTER="PASS"' \
  -S . -O z \
  -o $path_save/sampleID/genotypes.sorted.filter.vcf.gz $path_save/sampleID/genotypes.sorted.vcf.gz
  tabix -p vcf $path_save/sampleID/genotypes.sorted.filter.vcf.gz

	echo $sampleID >> $path_save/samples.NGS.txt
done

### 5. NGS states ###
ls $path_save/*/genotypes.sorted.filter.vcf.gz > $path_save/vcf.NGS.files
bcftools merge -O v -o $path_save/merge.paragraph.genotypes.vcf \
-m id -l $path_save/vcf.NGS.files

vcftools --vcf $path_save/merge.paragraph.genotypes.vcf \
--weir-fst-pop $path_example/samples_TIB_NGS.txt \
--weir-fst-pop $path_example/samples_HAN_NGS.txt \
--out $path_save/merge.paragraph.genotypes.TIBvsHAN
vcftools --vcf $path_save/merge.paragraph.genotypes.vcf \
--weir-fst-pop $path_example/samples_TIB_NGS.txt \
--weir-fst-pop $path_example/samples_HAN_NGS.txt \
--out $path_save/merge.paragraph.genotypes.TIBvsHAN.20k_5k \
--fst-window-size 20000 --fst-window-step 5000

python $path_scripts/func_trans_vcf_format_paragraph_general.py \
-v $path_save/merge.paragraph.genotypes.vcf \
-s $path_save/merge.paragraph.genotypes.reformat.vcf \
--file_samples $path_save/samples_group_NGS.bed
python $path_scripts/func_vcf_addINFO_POP.py \
-v $path_save/merge.paragraph.genotypes.reformat.vcf \
-s $path_save/merge.paragraph.genotypes.reformat.addPopInfo.vcf \
--file_samples $path_save/samples_subgroup_NGS.bed

python $path_scripts/func_vcf_addINFO.py \
-v $path_save/merge.paragraph.genotypes.reformat.addPopInfo.vcf \
-f $path_save/merge.paragraph.genotypes.TIBvsHAN.weir.fst \
-s $path_save/merge.paragraph.genotypes.reformat.addPopInfo.addFST.vcf

bgzip -c $path_save/merge.paragraph.genotypes.reformat.addPopInfo.addFST.vcf > $path_save/merge.paragraph.genotypes.reformat.addPopInfo.addFST.vcf.gz
bcftools index $path_save/merge.paragraph.genotypes.reformat.addPopInfo.addFST.vcf.gz


### 6. merge ONT and NGS calls ###
path_data='0_raw_data'
path_scripts='scripts'
path_save='.'

command_svtk='svtk'
command_ANNOTSV='AnnotSV.tcl'

bcftools merge -m id -O v \
-o $path_save/merge.genotypes.vcf \
2_sv_calling/result_genotyped/merge.genotyped.sorted.corrected.INDELtoSYMBOL.sorted.addMR.addFST.vcf.gz \
6_genotyping_NGS_paragraph_update/merge.paragraph.genotypes.reformat.addPopInfo.addFST.vcf.gz

# correct
python $path_scripts/func_check_vcf_svtype.py \
-v $path_save/merge.genotypes.vcf \
-s $path_save/merge.genotypes.corrected.vcf > $path_save/merge.genotypes.corrected.log

# del chrom
python $path_scripts/func_trans_vcf_format_delCHR.py \
-v $path_save/merge.genotypes.corrected.vcf \
-s $path_save/merge.genotypes.corrected.delCHR.vcf

### 7. annotation ###
# svtk annotate
$command_svtk annotate \
--gencode $path_data/gencode.v32lift37.canonical_annotation.gtf.gz \
$path_save/merge.genotypes.corrected.delCHR.vcf \
$path_save/merge.genotypes.corrected.delCHR.svtk.vcf
bgzip -c $path_save/merge.genotypes.corrected.delCHR.svtk.vcf > $path_save/merge.genotypes.corrected.delCHR.svtk.vcf.gz
bcftools index $path_save/merge.genotypes.corrected.delCHR.svtk.vcf.gz

#ANNOTSV
$command_ANNOTSV -SVinputFile $path_save/merge.genotypes.corrected.delCHR.svtk.vcf \
-outputDir $path_save/ \
-outputFile $path_save/merge.genotypes.corrected.delCHR.svtk.annotsv.tsv
