#!/bin/bash 
phenotype=$1

mkdir output/GTEx_HP/${phenotype}

for chr in {1..22}; do
  echo ${chr}
  Rscript FUSION/fusion_twas-master/FUSION.assoc_test.R \
  --sumstats data/${phenotype}_TWASInput.txt \
  --weights data/GTExv8.EUR.Brain_Hippocampus/GTExv8.EUR.Brain_Hippocampus.nofilter.pos \
  --weights_dir data/GTExv8.EUR.Brain_Hippocampus/ \
  --ref_ld_chr FUSION/fusion_twas-master/LDREF/1000G.EUR. \
  --chr ${chr} \
  --out output/GTEx_HP/${phenotype}/${phenotype}.${chr}.dat
done
