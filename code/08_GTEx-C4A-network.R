rm(list = ls())
options(stringsAsFactors = FALSE)
library(tidyverse)
path <- "/Users/minsookim/Desktop/C4A-network"


# 01. Construct C4A-seeded network in GTEx frontal cortex samples
# ----------------------------------

# Load data and metadata
load(paste0(path, "/data/GTEx_regressed_lmm_RSEM_final.RData"))

length(unique(datMeta$sub_num)) # 153 unique samples
length(unique(datMeta$region)) # 10 unique brain regions
sum(colnames(datExpr.reg) == datMeta$SAMPID)

gencode <- read.csv(paste0(path, "/data/annotation.gene.gencodeV19.csv"))
gencode <- gencode[match(rownames(datExpr.reg), gencode$gene_id), ]

# Load C4 imputation data
imputed <- read.table(paste0(path, "/data/GTEx-C4-imputed.txt"), header = TRUE)

datMeta$C4A_CN <- rep(10, nrow(datMeta))

for (i in 1:nrow(datMeta)) {
  if (length(grep(datMeta$SUBJID[i], imputed$sample)) > 0) {
    datMeta$C4A_CN[i] <- imputed$C4A_CN[grep(datMeta$SUBJID[i], imputed$sample)]
  }
}

datMeta <- datMeta[datMeta$C4A_CN != 10,] # 540 samples, 20765 features
datExpr.reg <- datExpr.reg[, match(datMeta$SAMPID, colnames(datExpr.reg))]
datExpr.norm <- datExpr.norm[, match(colnames(datExpr.reg), colnames(datExpr.norm))]

table(datMeta$region)
xtabs(~ region + C4A_CN, datMeta)

# Calculate C4A co-expression
prsCor = function(i, gene, datExpr){
  c <- cor.test(datExpr[i,], datExpr[gene,], use = 'pairwise.complete.obs')
  dfPrs <- data.frame(Gene = rownames(datExpr)[i], R = c$estimate, P = c$p.value)
  return(dfPrs)
}

dfBA9 <- datExpr.reg[, datMeta$region == "Frontal_Cortex_BA9" & datMeta$C4A_CN == 2] # 36 samples
dfcorBA9 = data.frame()

for (i in 1:nrow(datExpr.reg)) {
  if (i%%100 == 0) {print(i)}
  dfcorBA9 = rbind(dfcorBA9, prsCor(i, "ENSG00000244731", dfBA9))
}

hist(dfcorBA9$P)
dfcorBA9$FDR = p.adjust(dfcorBA9$P, method = "fdr")
dfcorBA9$Gene = gencode$gene_name[match(dfcorBA9$Gene, gencode$gene_id)]
sum(dfcorBA9$FDR < 0.05) # 6032 genes

