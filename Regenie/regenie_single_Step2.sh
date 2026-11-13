#!/bin/bash 

modality=$1
chr=$2
mkdir output_${modality}_All_Step2/chr${chr}
regenie \
  --step 2 \
  --chr ${chr} \
  --keep Data/${modality}_ID.txt \
  --bed Data/Q1_c${chr} \
  --phenoFile Data/Average_${modality}_All.txt \
  --covarFile Data/Average_${modality}_All_Cov.txt \
  --pred output_${modality}_All_Step1/All_Step1_pred.list \
  --qt \
  --apply-rint \
  --minMAC 20 \
  --threads 8 \
  --bsize 1000 \
  --write-samples \
  --print-pheno \
  --out output_${modality}_All_Step2/chr${chr}/All_Step2_c${chr}

