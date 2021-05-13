
# 1. preprocess simulation data
root_script=scripts
path_simulation=YRI_TIB_HAN/split_symmig_all.simulation
for chr in `seq 1 22`;do
  root_script=scripts
  for f in $(ls $path_simulation/chr$chr/variants.chr$chr.*.vcf.gz);do
    label_simulate=$(echo $f | cut -d / -f 4 | cut -d . -f 3);
    echo $label_simulate
    prefix=$path_simulation/chr$chr/variants.chr$chr.$label_simulate.filter
    if [ -d $prefix/chr$chr/$label_simulate ];then
      for f_sim in $(ls $prefix/chr$chr/$label_simulate/*.filter.vcf.gz);do
          label_ext=$(echo $f_sim | cut -d / -f 5 | cut -d . -f 4);
          prefix=$prefix/chr$chr/$label_simulate/variants.chr$chr.$label_simulate.$label_ext.filter
          echo $label_simulate.$label_ext
          if [ ! -f $prefix.target.setID.sorted.vcf.gz ];then
            echo $prefix
            bcftools index $prefix.vcf.gz
            bcftools view -R example/sv_fst_ngs_gte_0.1.surround.bed --threads 20 -O v -o $prefix.target.vcf $prefix.vcf.gz
            python $root_script/func_vcf_setID.py -v $prefix.target.vcf -s $prefix.target.setID.vcf
            bcftools sort -O z -o $prefix.target.setID.sorted.vcf.gz $prefix.target.setID.vcf
            rm $prefix.target.vcf
            rm $prefix.target.setID.vcf
          fi;
      done
    fi
  done
done;

# 2. ABBA-BABA test for each simulation
root_script_genetics=genomics_general

for chr in `seq 1 22`;do
  for LABEL in `seq 1 1000`;do
    # Denisovan
    ANC=DEN
    prefix_raw=$path_simulation/chr$chr
    path_bed=example/chr'$chr'.recode.target.setID.sorted.bed
    prefix=$path_simulation/chr$chr

    if [ ! -d $prefix ];then
      mkdir $prefix
    fi

    if [ -d $prefix_raw/$LABEL ];then
      if [ -f $prefix_raw/$LABEL/variants.chr$chr.$LABEL.filter.target.setID.sorted.vcf.gz ];then
        if [ ! -f $prefix_raw/$LABEL/variants.chr$chr.$LABEL.filter.target.setID.sorted.geno ];then
           python $root_script_genetics/VCF_processing/parseVCF.py \
           -i $prefix_raw/$LABEL/variants.chr$chr.$LABEL.filter.target.setID.sorted.vcf.gz \
           -o $prefix_raw/$LABEL/variants.chr$chr.$LABEL.filter.target.setID.sorted.geno
           sed -i 's/0|0/C|C/g' $prefix_raw/$LABEL/variants.chr$chr.$LABEL.filter.target.setID.sorted.geno
           sed -i 's/0|1/C|T/g' $prefix_raw/$LABEL/variants.chr$chr.$LABEL.filter.target.setID.sorted.geno
           sed -i 's/1|0/T|C/g' $prefix_raw/$LABEL/variants.chr$chr.$LABEL.filter.target.setID.sorted.geno
           sed -i 's/1|1/T|T/g' $prefix_raw/$LABEL/variants.chr$chr.$LABEL.filter.target.setID.sorted.geno
        fi;
        if [ ! -d $prefix/$LABEL ];then
          mkdir $prefix/$LABEL
        fi
        echo $prefix_raw/$LABEL/variants.chr$chr.$LABEL.filter.target.setID.sorted.geno
        python $root_script_genetics/ABBABABAwindows.py \
        -g $prefix_raw/$LABEL/variants.chr$chr.$LABEL.filter.target.setID.sorted.geno \
        -m 100 --windType predefined -P1 TIB -P2 YRI -P3 $ANC -O CHIMP \
        --popsFile example/samples.simulation.clust -f phased \
        -o $prefix/$LABEL/variants.chr$chr.$LABEL.filter.target.setID.sorted.stats.TIB.$ANC.csv \
        --windCoords $path_bed
      fi;
    fi;
  done
done


