#!/bin/bash 

modality=$1

regenie \
  --step 1 \
  --bed Data/ukb_cal_allChrs_hg38 \
  --extract Data/qc_pass.snplist \
  --keep Data/${modality}_ID.txt \
  --phenoFile Data/Average_${modality}_All.txt \
  --covarFile Data/Average_${modality}_All_Cov.txt \
  --qt \
  --apply-rint \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix tmp/tmp_${modality} \
  --threads 8 \
  --out output_${modality}_All_Step1/All_Step1
