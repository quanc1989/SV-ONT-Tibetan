# pipeline for demographic inference and simulation

path_software_easySFS=easySFS.py
path_vcf_include_YRI=example/merge.genotypes.corrected.delCHR.svtk.chr1.vcf
path_samples_cluster=example/samples.YRIandTIBandHAN.clust
root_script=scripts

# 1. construct sfs and proj to a smaller set
$path_software_easySFS \
-p $path_samples_cluster \
-i $path_vcf_include_YRI \
-a -f \
-o YRI_TIB_HAN/sfs.pruned \
--proj 50,40,40 \
--unfolded

# 2. generate bootstrap dataset and calculate GIM
arry_trio=('YRI_TIB_HAN')
arry_model=("sfs.pruned")

for trio in ${arry_trio[@]};do
    echo $trio
    fname_fs=$(echo $trio | sed 's/_/-/g');
    poplist=$(echo $trio | sed 's/_/,/g');
    proj='50,40,40'
    echo $poplist
    echo $proj
    for model in ${arry_model[@]};do
        echo $model
        python $root_script/run_generate_sfs_segments.py \
        -d $trio/$model/datadict.txt \
        -p $trio/$model/bootstrap \
        -l $poplist -u --random \
        --projections $proj
    done;
done;

# 3. demographic inference
model_list=('split_symmig_all')
for model in ${model_list[@]};do
  python $root_script/run_inference.py \
  -s YRI_TIB_HAN/YRI-TIB-HAN.sfs \
  -p YRI_TIB_HAN/YRI-TIB-HAN -m $model \
  --unfolded
done

# 4. simulation
model='sfs.pruned'
path_hapmap=example/hapmap
for label_simulate in `seq 1 1000`; do
  for ((i=22;i>=1;i--));do
          chrom='chr'$i
    echo $chrom
    python $root_script/run_simulate_pipeline.py $chrom $label_simulate YRI_TIB_HAN/YRI-TIB-HAN.sfs YRI_TIB_HAN/$model/segments 'split_symmig_all' YRI_TIB_HAN/$model/ normal $path_hapmap
  done;
done;