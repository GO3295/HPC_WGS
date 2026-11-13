library(susieR)
library(data.table)

sumstats_file <- commandArgs(trailingOnly=TRUE)[1]
landmark <- commandArgs(trailingOnly=TRUE)[2]

# read sumstats data
sumstats <- fread(sumstats_file,header=TRUE,stringsAsFactors=FALSE)
sumstats <- as.data.frame(sumstats)

# read LD data
LD_data <- fread(paste0("tmp/tmp.unphased.vcor1"), sep="\t", header=FALSE, stringsAsFactors=FALSE)
LD_data <- as.matrix(LD_data)
LD_id <- fread(paste0("tmp/tmp.unphased.vcor1.vars"), sep="\t", header=FALSE, stringsAsFactors=FALSE)
if(sum(LD_id$V1==sumstats$ID)==nrow(LD_id)){
  rownames(LD_data) <- sumstats$ID
  colnames(LD_data) <- sumstats$ID
}else{
  stop("Error: Order/Number incorrect.")
}
LD_data[is.na(LD_data)] <- 0

fitted_rss <- susie_rss(bhat=sumstats$BETA,shat=sumstats$SE,R=LD_data,n=median(sumstats$N),L=10)
saveRDS(fitted_rss,paste0("output/rds/susieR_",landmark,".rds"))

