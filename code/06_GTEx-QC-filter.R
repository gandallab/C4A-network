rm(list = ls())
options(stringsAsFactors = FALSE)

library(tidyverse)
library(data.table)



# 01. Initial filter for GTEx v7 samples and genes that were used for eQTL analyses ----

# Read raw, un-normalized gene read counts from GTEx
gtex.read <- as.data.frame(fread("./data/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz",
                                 header = TRUE)) # 56202 gene features, 11688 samples

# Read GTEx sample- and subject-level information
gtex.sample.attrib <- as.data.frame(fread("./data/GTEx/phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt", 
                                          header = TRUE))
gtex.pheno <- as.data.frame(fread("./data/GTEx/phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt", 
                                  header = TRUE))
colnames(gtex.pheno) <- as.character(gtex.pheno[1, ])
gtex.pheno <- gtex.pheno[-1, ]

# Identify brain samples
brain.regions <- c("Brain - Hippocampus", "Brain - Anterior cingulate cortex (BA24)", "Brain - Cerebellar Hemisphere", "Brain - Hypothalamus", 
                   "Brain - Caudate (basal ganglia)", "Brain - Nucleus accumbens (basal ganglia)", "Brain - Putamen (basal ganglia)",
                   "Brain - Amygdala", "Brain - Cerebellum", "Brain - Frontal Cortex (BA9)", "Brain - Substantia nigra",
                   "Brain - Cortex", "Brain - Spinal cord (cervical c-1)")
brain.reg.short <- c("Hippocampus", "Anterior_cingulate_cortex_BA24", "Cerebellar_Hemisphere", "Hypothalamus", "Caudate_basal_ganglia", 
                     "Nucleus_accumbens_basal_ganglia", "Putamen_basal_ganglia", "Amygdala", "Cerebellum", "Frontal_Cortex_BA9", 
                     "Substantia_nigra", "Cortex", "Spinal_cord_cervical_c-1")
gtex.brain.attrib <- gtex.sample.attrib[which(gtex.sample.attrib$SMTSD %in% brain.regions), ] # 2076 brain samples out of 15598

# Determine which genes are in all brain regions
count.df <- NULL
for (reg.idx in 1:length(brain.regions)) {
  # Read fully normalized expression values that were used in the FastQTL analyses to subset the unnormalized count values
  reg <- brain.reg.short[reg.idx]
  print(reg)
  tissue.file <- paste0("./data/GTEx/GTEx_Analysis_v7_eQTL_expression_matrices/Brain_", reg, ".v7.normalized_expression.bed")
  tissue.dat <- read.table(tissue.file, comment.char = "", header = TRUE)
  tissue.expr <- tissue.dat[, grep("GTEX", colnames(tissue.dat))]
  rownames(tissue.expr) <- tissue.dat$gene_id
  count.df[[reg]] <- tissue.expr
}

# Only keep genes that are present in all brain regions
all.genes <- NULL
for (i in names(count.df)) {
  all.genes <- c(all.genes, rownames(count.df[[i]]))
}
gene.table <- table(all.genes) # create table of gene counts
genes.keep <- names(gene.table)[which(gene.table == length(brain.regions))]
length(genes.keep) # 20980 genes present in all brain regions

# Extract count data for samples passing GTEx QC thresholds and genes in all brain regions
covariate.df <- NULL
count.df <- matrix()

for (reg.idx in 1:length(brain.regions)) { 
  reg <- brain.reg.short[reg.idx]
  print(reg)
  tissue.file <- paste0("./data/GTEx/GTEx_Analysis_v7_eQTL_expression_matrices/Brain_", reg, ".v7.normalized_expression.bed")
  tissue.dat <- read.table(tissue.file, comment.char = "", header = TRUE)
  tissue.expr <- tissue.dat[, grep("GTEX", colnames(tissue.dat))]
  rownames(tissue.expr) <- tissue.dat$gene_id
  
  # Subset count data to only include samples from the current region
  region <- brain.regions[reg.idx]
  region.sampids <- gtex.sample.attrib$SAMPID[which(gtex.sample.attrib$SMTSD %in% region)]
  region.samp.att <- gtex.sample.attrib[which(gtex.sample.attrib$SMTSD %in% region), ]
  
  # Pull the reads for this region (from the "gtex.read" variable)
  region.counts <- gtex.read[, which(colnames(gtex.read) %in% region.sampids)]
  rownames(region.counts) <- gtex.read$Name

  # Identify the genes that passed GTEx QC pipeline
  tpm.use.genes <- gtex.read$Name[which(gtex.read$Name %in% tissue.dat$gene_id)]
  # Identify samples that passed GTEx QC pipeline
  gtex.use.samples <- colnames(tissue.dat)[grep("GTEX", colnames(tissue.dat))] # gtex processed
  gtex.read.refs <- gsub("-", ".", colnames(region.counts)) # count data
  use.samp.idxs <- unlist(lapply(gtex.use.samples, grep, x = gtex.read.refs)) # idxs in the count data for usable samples
  
  # Subset count data so that they reflects GTEx QC methods
  # Genes that pass GTEx QC and are present in all brain samples
  use.genes <- intersect(tpm.use.genes, genes.keep)
  USE.reg.counts <- region.counts[which(rownames(region.counts) %in% use.genes), use.samp.idxs]
  
  # Append the samples as columns in the data frame
  count.df <- cbind(count.df, as.data.frame(USE.reg.counts)) # First column is NA
  
  # Get region-wise covariates
  covar.file <- paste0("./data/GTEx/GTEx_Analysis_v7_eQTL_covariates/Brain_", reg, ".v7.covariates.txt")
  tissue.covars.tmp <- read.table(covar.file, comment.char = "", header = TRUE)  
  tissue.covar <- tissue.covars.tmp[, grep("GTEX", colnames(tissue.covars.tmp))]
  rownames(tissue.covar) <- tissue.covars.tmp$ID
  if (dim(tissue.covar)[1] > 20) {
    tissue.covar = tissue.covar[-(19:33), ]
  } 
  region <- rep(reg, dim(tissue.covar)[2])
  
  # Phenotype information for the current tissue
  region.pheno <- gtex.pheno[gsub("-", ".", gtex.pheno$SUBJID) %in% colnames(tissue.covar), ]
  
  # Sequencing metric
  seq <- region.samp.att[which(region.samp.att$SAMPID %in% colnames(USE.reg.counts)), ]

  # Check the order of each data frame: count / covariate / phenotype 
  # (make sure subject order is maintained before concatenating into a larger dataframe)
  print(length(which(gsub("-", ".", region.pheno$SUBJID) == colnames(tissue.covar))))
  count.sub.names <- unlist(lapply(colnames(USE.reg.counts), function(x) paste(strsplit(x, '-')[[1]][1:2], collapse = ".")))
  length(which(count.sub.names == colnames(tissue.covar)))
  seq <- seq[match(count.sub.names, colnames(tissue.covar)), ]
  
  # Transpose covariates
  t.tissue.covar <- as.data.frame(t(tissue.covar))
  
  # Append regional covariates to the overall covariate dataframe
  sub_num <- rownames(t.tissue.covar)
  covar.out <- cbind(t.tissue.covar, sub_num, region, region.pheno, seq) # 89 covariates
  covariate.df <- rbind(covariate.df, covar.out)
}

count.df$count.df <- NULL # get rid of column header with NAs

datExpr.counts <- count.df # 20980 features
rownames(datExpr.counts) <- gsub("\\..*", "", rownames(datExpr.counts))
datMeta <- covariate.df # 1497 brain samples



# 02. Convert transcript from RSEM to gene-level counts and filter for above features and samples ----

gtex.rsem <- as.data.frame(fread("./data/GTEx/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_expected_count.txt.gz", 
                                 header = TRUE)) 
# 196520 features, 11688 samples
gtex.rsem[1:5,1:5]

tx2gene <- data.frame(gtex.rsem[,1:2])
rownames(gtex.rsem) <- gtex.rsem$transcript_id 
gtex.rsem <- gtex.rsem[, c(-1,-2)]
gtex.rsem <- as.matrix(gtex.rsem)
gtex.rsem <- rowsum(gtex.rsem, tx2gene$gene_id) # 57820 features, 11688 samples
rownames(gtex.rsem) <- gsub("\\..*", "", rownames(gtex.rsem))
gtex.rsem <- gtex.rsem[match(rownames(datExpr.counts), rownames(gtex.rsem)), ]
gtex.rsem <- gtex.rsem[, match(colnames(datExpr.counts), colnames(gtex.rsem))]

save(gtex.rsem, datMeta, file = "./data/GTEx/GTEX_v7_datExpr_counts_brain_RSEM.RData")


