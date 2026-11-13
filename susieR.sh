#!/bin/bash 

sed -n '2,100p' landmark.txt | while IFS=$'\t' read -r CHR POS LEAD PHENO landmark MODALITY; do
    
    echo "Process: CHR=$CHR, LANDMARK=$landmark, PHENO=$PHENO, MODALITY=$MODALITY"
    
    sumstats_file=tmp_data/loci$landmark.txt
    sample_file=Average_All_ID.txt
    
    # variant id
    cut -f3 -d$'\t' "$sumstats_file" | tail -n+2 > tmp/test_id.txt

    # get LD
    plink2 \
     --bfile Data/Q1_c${CHR} \
     --threads 8 \
     --memory 100000 \
     --r-unphased square yes-really \
     --extract tmp/test_id.txt \
     --keep ${sample_file} \
     --out tmp/tmp

    # # susieR fine mapping
    Rscript susieR_sumstats.R $sumstats_file $landmark

    rm -f tmp/*
    
    echo "$landmark done"
    
done