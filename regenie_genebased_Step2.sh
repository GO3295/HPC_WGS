#!/bin/bash

modality=$1
chr=$2
anno_type=$3

mkdir output_${modality}_All_Step2/${anno_type}/chr${chr}

regenie \
  --step 2 \
  --chr ${chr} \
  --keep Data/${modality}_ID.txt \
  --bed Data/Q1_c${chr} \
  --phenoFile Data/Average_${modality}_All.txt \
  --covarFile Data/Average_${modality}_All_Cov.txt \
  --pred output_${modality}_All_Step1/All_Step1_pred.list \
  --anno-file Anno/chr${chr}/Main/${anno_type}_chr${chr}.txt \
  --set-list  Anno/chr${chr}/Main/chr${chr}_${anno_type}.setlist \
  --mask-def Anno/Mask/Mask_${anno_type}.txt \
  --aaf-bins 0.01 \
  --minMAC 1 \
  --vc-tests skato,acato-full \
  --vc-maxAAF 0.01 \
  --write-mask \
  --write-mask-snplist \
  --qt \
  --apply-rint \
  --bsize 1000 \
  --threads 8 \
  --write-samples \
  --print-pheno \
  --rgc-gene-p \
  --out output_${modality}_All_Step2/${anno_type}/chr${chr}/c${chr}_${anno_type}


# score
score=$4

mkdir output_${modality}_All_Step2/${anno_type}_${score}/chr${chr}
regenie \
  --step 2 \
  --chr ${chr} \
  --keep Data/${modality}_ID.txt \
  --bed Data/Q1_c${chr} \
  --phenoFile Data/Average_${modality}_All.txt \
  --covarFile Data/Average_${modality}_All_Cov.txt \
  --pred output_${modality}_All_Step1/All_Step1_pred.list \
  --anno-file Anno/chr${chr}/Main/${anno_type}_chr${chr}.txt \
  --set-list  Anno/chr${chr}/Main/chr${chr}_${anno_type}.setlist \
  --mask-def Anno/Mask/Mask_${anno_type}.txt \
  --extract Anno/chr${chr}/Subset/chr${chr}_${score}.txt \
  --aaf-bins 0.01 \
  --minMAC 1 \
  --vc-tests skato,acato-full \
  --vc-maxAAF 0.01 \
  --write-mask \
  --write-mask-snplist \
  --qt \
  --apply-rint \
  --bsize 1000 \
  --threads 8 \
  --write-samples \
  --print-pheno \
  --rgc-gene-p \
  --out output_${modality}_All_Step2/${anno_type}_${score}/chr${chr}/c${chr}_${anno_type}_${score}
  
