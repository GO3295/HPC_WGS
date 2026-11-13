library(data.table)
library(TwoSampleMR)

# -----------------------------------------------------------------------------
# preprocessing
# -----------------------------------------------------------------------------
PATH_SUMSTATS   <- "processed_sumstats"
PATH_TMP       <- "tmp/HP_eQTL"
OUT_DIR        <- "output" 
PLINK_SCRIPT   <- "Clumping.sh"
PLINK_R2       <- "0.1" #0.01/0.001                  

# Trait groups
hippo_trait <- c("HP_tail","Whole_HP_body","Whole_HP_head","Whole_HP","HP_fissure",
                 "Subiculum_body","Subiculum_head","Presubiculum_body","Presubiculum_head",
                 "CA1_body","CA1_head","CA3_body","CA3_head","CA4_body","CA4_head",
                 "Parasubiculum","MLHP_body","MLHP_head","GCMLDG_body","GCMLDG_head",
                 "HATA","Fimbria","FA","ICVF","ISOVF","L1","L2","L3","MD","MO","OD",
                 "MT2star","MMS")
T1_trait  <- c("HP_tail","Whole_HP_body","Whole_HP_head","Whole_HP","HP_fissure",
               "Subiculum_body","Subiculum_head","Presubiculum_body","Presubiculum_head",
               "CA1_body","CA1_head","CA3_body","CA3_head","CA4_body","CA4_head",
               "Parasubiculum","MLHP_body","MLHP_head","GCMLDG_body","GCMLDG_head",
               "HATA","Fimbria")
DTI_trait <- c("FA","ICVF","ISOVF","L1","L2","L3","MD","MO","OD")
SWI_trait <- c("MT2star","MMS")

# -----------------------------------------------------------------------------
# 2. Load data
# -----------------------------------------------------------------------------
AllSigVarPair      <- fread("SingleVar.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

qtl_sumstats <- fread("eGeneSigSumStats_HP_eQTL.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
qtl_sumstats <- qtl_sumstats[!is.na(slope)]

qtl_sumstats[, c("chr","pos","ref","alt") := tstrsplit(variant_id, "_", fixed = TRUE)]
qtl_sumstats[, chr := sub("^chr", "", chr)]
qtl_sumstats[, pos := as.integer(pos)]

qtl_sumstats[, chrpos := paste0(chr, ":", pos)]
qtl_sumstats[, id := paste0(chr, ":", pos, ":", ref, ":", alt)]

# -----------------------------------------------------------------------------
# 3. loop over eGenes
# -----------------------------------------------------------------------------
mr_res <- data.frame() 

for (g in unique(qtl_sumstats$gene_id)) {
  message("Processing gene: ", g)
  g_qtl_sumstats <- qtl_sumstats[gene_id == g]
  
  min_pos <- min(g_qtl_sumstats$pos)
  max_pos <- max(g_qtl_sumstats$pos)
  
  # Identify all traits significantly associated with this eGene
  tmp_g   <- paste0(",", g, ",")
  tmp_sig <- AllSigVarPair[grepl(tmp_g, AllSigVarPair$GTEx_HP_eQTL), c("id","chr","trait")]
  tmp_chr <- unique(tmp_sig$chr)
  tmp_trait <- unique(tmp_sig$trait)
  
  stopifnot(length(tmp_chr) == 1L)
  
  tmp_trait1 <- tmp_trait[1]
  tmp_trait1_name <- sprintf("%s/chr%s/All_c%s_%s.txt",
                             PATH_SUMSTATS, tmp_chr, tmp_chr, tmp_trait1)
  tmp_trait1_sumstats <- fread(tmp_trait1_name, sep = "\t", header = TRUE,
                               stringsAsFactors = FALSE)
  tmp_trait1_sumstats <- tmp_trait1_sumstats[GENPOS >= min_pos & GENPOS <= max_pos]
  
  # Keep only variants present in both exposure and outcome
  g_qtl_sumstats <- g_qtl_sumstats[
    id %in% unique(c(tmp_trait1_sumstats$simple_flipped_id,
                     tmp_trait1_sumstats$simple_id))
  ]
  tmp_trait1_sumstats <- tmp_trait1_sumstats[
    simple_id %in% g_qtl_sumstats$id | simple_flipped_id %in% g_qtl_sumstats$id
  ]
  
  flip_idx <- g_qtl_sumstats$id %in% tmp_trait1_sumstats$simple_flipped_id
  if (any(flip_idx)) {
    g_qtl_sumstats[flip_idx, slope := -slope]
    g_qtl_sumstats[flip_idx, id := paste0(chr, ":", pos, ":", alt, ":", ref)]
  }
  g_qtl_sumstats[, c("ref","alt") := tstrsplit(id, ":", fixed = TRUE)[3:4] |> as.list() |> setNames(c("ref","alt"))]
  
  g_qtl_sumstats[, ref := tstrsplit(id, ":", fixed = TRUE)[[3]]]
  g_qtl_sumstats[, alt := tstrsplit(id, ":", fixed = TRUE)[[4]]]
  
  expo <- as.data.frame(g_qtl_sumstats)
  expo$Phenotype <- g; expo$gene <- g
  expo <- expo[, .(SNP = id, beta = slope, se = slope_se,
                   effect_allele = alt, other_allele = ref,
                   Phenotype, chr, position = pos,
                   pval = pval_nominal, gene)]
  
  # ---  outcome ---
  out1 <- as.data.frame(tmp_trait1_sumstats)
  out1$Phenotype <- tmp_trait1; out1$gene <- g
  out1 <- out1[, c("simple_id","BETA","SE","alt","ref","Phenotype",
                   "CHROM","GENPOS","N","pval","gene","ID")]
  colnames(out1) <- c("SNP","beta","se","effect_allele","other_allele",
                      "Phenotype","chr","position","samplesize","pval","gene","original_id")
  
  stopifnot(all(expo$SNP %in% rownames(setNames(seq_len(nrow(out1)), out1$SNP))))
  rownames(expo) <- expo$SNP
  rownames(out1) <- out1$SNP
  expo <- expo[order(expo$position), ]
  out1 <- out1[expo$SNP, ]
  expo$original_id <- out1$original_id
  
  # --- clumping ---
  tmp_clump <- data.table(SNP = expo$original_id, P = expo$pval)
  fwrite(tmp_clump[, .(variant = SNP)], file.path(PATH_TMP, "tmpVar_forclump.txt"),
         sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  fwrite(tmp_clump, file.path(PATH_TMP, "tmp_forclump.txt"),
         sep = "\t", quote = FALSE, row.names = FALSE)
  
  status <- system2("sh", args = c(PLINK_SCRIPT, tmp_chr, PLINK_R2, "HP_eQTL"))
  if (status != 0) stop("clumping failed for", g, " chr ", tmp_chr)
  
  clumped_file <- file.path(PATH_TMP, paste0("clump", PLINK_R2, ".clumps"))
  clumped_data <- fread(clumped_file, stringsAsFactors = FALSE)
  
  qtl_clumped <- expo[expo$original_id %in% clumped_data$ID, ]
  
  for (x in tmp_trait) {
    message("  Trait: ", x)
    
    x_name <- sprintf("%s/chr%s/All_c%s_%s.txt",
                      PATH_SUMSTATS,tmp_chr, tmp_chr, x)
    out_x <- fread(x_name, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    out_x <- out_x[ID %in% expo$original_id]
    out_x <- as.data.frame(out_x)
    out_x$Phenotype <- x; out_x$gene <- g
    out_x <- out_x[, c("simple_id","BETA","SE","alt","ref","Phenotype",
                       "CHROM","GENPOS","N","pval","gene","ID")]
    colnames(out_x) <- c("SNP","beta","se","effect_allele","other_allele",
                         "Phenotype","chr","position","samplesize","pval","gene","original_id")
    rownames(out_x) <- out_x$SNP
    out_x_clumped <- out_x[qtl_clumped$SNP, ]
    
    # --- Run MR ---
    # method_list = c("mr_wald_ratio","mr_egger_regression","mr_weighted_median","mr_ivw","mr_simple_mode",
    # "mr_weighted_mode","mr_ivw_radial","mr_raps")
    expo_fmt  <- format_data(qtl_clumped,    type = "exposure")
    out_fmt   <- format_data(out_x_clumped,  type = "outcome")
    dat       <- harmonise_data(expo_fmt, out_fmt, action = 1)
    res       <- mr(dat)
    
    mr_res <- rbind(mr_res, res)
  }
}

fwrite(mr_res, file.path(OUT_DIR, "MR_HP_eQTL.txt"),
       sep = "\t", quote = FALSE)