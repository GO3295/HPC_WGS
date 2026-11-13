library(data.table)
library(coloc)

hippo_trait <- c("HP_tail","Whole_HP_body","Whole_HP_head","Whole_HP","HP_fissure","Subiculum_body","Subiculum_head",
                 "Presubiculum_body","Presubiculum_head","CA1_body","CA1_head","CA3_body","CA3_head","CA4_body","CA4_head",
                 "Parasubiculum","MLHP_body","MLHP_head","GCMLDG_body","GCMLDG_head","HATA","Fimbria","FA","ICVF","ISOVF",
                 "L1","L2","L3","MD","MO","OD","MT2star","MMS")
T1_trait <- c("HP_tail","Whole_HP_body","Whole_HP_head","Whole_HP","HP_fissure","Subiculum_body","Subiculum_head",
              "Presubiculum_body","Presubiculum_head","CA1_body","CA1_head","CA3_body","CA3_head","CA4_body","CA4_head",
              "Parasubiculum","MLHP_body","MLHP_head","GCMLDG_body","GCMLDG_head","HATA","Fimbria")
DTI_trait <- c("FA","ICVF","ISOVF","L1","L2","L3","MD","MO","OD")
SWI_trait <- c("MT2star","MMS")

AllSigVarPair <- fread("SingleVar_Anno.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
AllSigVarPair <- as.data.frame(AllSigVarPair)
Gene_evidence <- readRDS("Variant_GeneEvidence_v2.rds")

# coloc
qtl_sumstats <- fread("SigGeneSumStats/Brain_pQTL_SigGeneSumStats_Processed.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
qtl_sumstats <- qtl_sumstats[!is.na(qtl_sumstats$BETA),]

# Initialize result list
coloc_res_list <- list()
res_idx <- 1

# Iterate over genes
for (g in unique(qtl_sumstats$GENE)) {
  
  message("Processing gene: ", g)
  
  gene_qtl <- qtl_sumstats[qtl_sumstats$GENE == g, ]
  
  gene_pos <- as.numeric(gene_qtl$pos)
  min_pos <- min(gene_pos, na.rm = TRUE)
  max_pos <- max(gene_pos, na.rm = TRUE)
  
  gene_pattern <- paste0(",", g, ",")
  gene_sig <- AllSigVarPair[
    grepl(gene_pattern, AllSigVarPair$Brain_pQTL, fixed = TRUE),
    c("id", "chr", "trait")
  ]
  
  if (nrow(gene_sig) == 0) {
    next
  }
  
  gene_sig <- unique(gene_sig[, c("chr", "trait")])
  
  # Iterate over chromosome-trait combinations
  for (i in seq_len(nrow(gene_sig))) {
    
    chr <- gene_sig$chr[i]
    trait <- gene_sig$trait[i]
    
    message("  Processing trait: ", trait, " (chr", chr, ")")
    
    file_path <- file.path(
      "tmp_res",
      trait,
      paste0("chr", chr, "_", trait, ".txt")
    )
    
    # Skip missing files
    if (!file.exists(file_path)) {
      warning("File not found: ", file_path)
      next
    }
    
    # Read trait summary statistics
    trait_sumstats <- fread(file_path, sep = "\t", header = TRUE)
    
    # Restrict variants to the gene's genomic region
    trait_pos <- as.numeric(trait_sumstats$GENPOS)
    
    trait_sumstats <- trait_sumstats[
      !is.na(trait_pos) &
        trait_pos >= min_pos &
        trait_pos <= max_pos
    ]
    
    # Match variants shared by both datasets
    common_snps <- intersect(
      gene_qtl$processed_id,
      trait_sumstats$processed_id
    )
    
    # Skip if no shared variants
    if (length(common_snps) == 0) {
      next
    }
    
    # Subset both datasets and align SNP order
    gene_qtl_sub <- gene_qtl[
      match(common_snps, gene_qtl$processed_id),
    ]
    
    trait_sumstats_sub <- trait_sumstats[
      match(common_snps, trait_sumstats$processed_id),
    ]
    
    valid_idx <- complete.cases(
      gene_qtl_sub[, c("BETA", "varbeta", "pos")],
      trait_sumstats_sub[, c("BETA", "varbeta", "GENPOS")]
    ) &
      is.finite(gene_qtl_sub$BETA) &
      is.finite(gene_qtl_sub$varbeta) &
      is.finite(trait_sumstats_sub$BETA) &
      is.finite(trait_sumstats_sub$varbeta) &
      gene_qtl_sub$varbeta > 0 &
      trait_sumstats_sub$varbeta > 0
    
    gene_qtl_sub <- gene_qtl_sub[valid_idx, ]
    trait_sumstats_sub <- trait_sumstats_sub[valid_idx, ]
    
    if (nrow(gene_qtl_sub) < 2) {
      warning(
        "Insufficient valid shared SNPs for gene ",
        g, " and trait ", trait
      )
      next
    }
    
    # Prepare datasets for colocalization analysis
    coloc_eGene <- list(
      beta = gene_qtl_sub$BETA,
      varbeta = gene_qtl_sub$varbeta,
      snp = gene_qtl_sub$processed_id,
      position = as.numeric(gene_qtl_sub$pos),
      type = "quant",
      sdY = 1
    )
    
    coloc_trait <- list(
      beta = trait_sumstats_sub$BETA,
      varbeta = trait_sumstats_sub$varbeta,
      snp = trait_sumstats_sub$processed_id,
      position = as.numeric(trait_sumstats_sub$GENPOS),
      type = "quant",
      sdY = 1
    )
    
    # Run colocalization analysis
    res <- coloc.abf(
      dataset1 = coloc_eGene,
      dataset2 = coloc_trait
    )
    
    # Store results
    res_df <- as.data.frame(t(res$summary))
    res_df$Gene <- g
    res_df$trait <- trait
    res_df$chr <- chr
    res_df$source <- "Brain_pQTL"
    
    coloc_res_list[[res_idx]] <- res_df
    res_idx <- res_idx + 1
  }
}

# Combine all results
if (length(coloc_res_list) > 0) {
  coloc_res <- rbindlist(coloc_res_list, fill = TRUE)
} else {
  coloc_res <- data.table()
}

fwrite(coloc_res,"coloc_output/coloc_Brain_pQTL.txt",sep="\t",quote=FALSE)
print("done")