rm(list = ls())
options(stringsAsFactors = FALSE)
library(tidyverse)
path <- "/Users/minsookim/Desktop/C4A-network"


# 01. Further exclude samples that are processed differently or suffer from brain-related disorders
# ----------------------------------
load(paste0(path, "/data/GTEX_v7_datExpr_counts_brain_RSEM.RData")) # Output from 04_GTEx-QC-filter.R script
datExpr.counts <- gtex.rsem
rm(gtex.rsem)

datMeta <- datMeta %>% 
  filter(SMMTRLTP != "Tissue:PAXgene Preserved", # Remove cerebellum, cortex samples w/ different sample preparation
         region != "Spinal_cord_cervical_c-1", # Remove spinal cord
         !is.na(datMeta$SMRDLGTH),
         INCEXC == "True",
         MHALS != "1",
         MHALZDMT != "1",
         MHDMNTIA != "1",
         MHENCEPHA != "1",
         MHFLU != "1",
         MHJAKOB != "1",
         MHMS != "1",
         MHPRKNSN != "1",
         MHREYES != "1",
         MHSEPSIS != "1",
         MHLUPUS != "1",
         MHCVD != "1",
         MHHIVCT != "1",
         MHALZHMR != "1") # 938 samples remaining

datExpr.counts <- datExpr.counts[, match(datMeta$SAMPID, colnames(datExpr.counts))] 

datMeta$region <- as.factor(datMeta$region)
datMeta$DTHCODD <- as.numeric(datMeta$DTHCODD)
datMeta$DTHCODD[datMeta$DTHCODD >= 0 & datMeta$DTHCODD < 2] <- "0to2h"
datMeta$DTHCODD[datMeta$DTHCODD >= 2 & datMeta$DTHCODD < 10] <- "2to10h"
datMeta$DTHCODD[datMeta$DTHCODD >= 10 & datMeta$DTHCODD < 72] <- "10hto3d"
datMeta$DTHCODD[datMeta$DTHCODD >= 72 & datMeta$DTHCODD < 504] <- "3dto3w"
datMeta$DTHCODD[datMeta$DTHCODD >= 504] <- "3wplus"
datMeta$DTHCODD[is.na(datMeta$DTHCODD)] <- "unknown"


# 02. Visualize raw count data
# ----------------------------------
boxplot(datExpr.counts, range = 0, main = "Raw Counts", xlab = "GTEx samples")
boxplot(log2(1 + datExpr.counts), range = 0, main = "log2 (counts + 1)", xlab = "GTEx samples")

i = 1; plot(density(log2(.001 + datExpr.counts[,i])), col = as.factor(datMeta$region)[i], main = "Gene read counts", 
          xlab = "log2(raw_counts + .001)", xlim = c(-15,30), ylim = c(0,0.45))

for (i in 2:ncol(datExpr.counts)) {
    lines(density(log2(.001 + datExpr.counts[,i])), col = as.factor(datMeta$region)[i])
}

# Calculate and visualize top 10 principle components
mds = cmdscale(dist(t(log2(datExpr.counts + 1))), k = 10)
colnames(mds) = paste0("PC", 1:ncol(mds))

pairs(mds, col = factor(datMeta$region), main="Region", pch = 19, cex = 0.5)
par(xpd = TRUE, oma = c(1,1,1,1)); legend("topright", levels(factor(datMeta$region)), fill = 1:10, cex = .4)

pairs(mds, col = factor(datMeta$sex), main = "Sex", pch = 19, cex = 0.5)
par(xpd = TRUE, oma = c(1,1,1,1)); legend("topright", levels(factor(datMeta$sex)), fill = 1:2, cex = .4)

pairs(mds, col = factor(datMeta$AGE), main = "Age", pch = 19, cex = 0.5)
par(xpd = TRUE, oma = c(1,1,1,1)); legend("topright", levels(factor(datMeta$AGE)), fill = 1:6, cex = .4)

pairs(mds, col = factor(datMeta$DTHHRDY), main = "DeathProcess", pch = 19, cex = 0.5)
par(xpd = TRUE, oma = c(1,1,1,1)); legend("topright", levels(factor(datMeta$DTHHRDY)), fill = 1:5, cex = .4)

library(WGCNA)
tree = hclust(as.dist(1 - bicor(log2(datExpr.counts + 1)), "average"))

plotDendroAndColors(tree, colors = cbind(labels2colors(datMeta$region), labels2colors(datMeta$sex), 
                                         labels2colors(datMeta$AGE), labels2colors(datMeta$DTHHRDY)), 
                    groupLabels = c("Region", "Sex", "Age", "DeathProcess"))


# 03. Visualize gene-filtered & normalized data
# ----------------------------------
library(edgeR)

# Filter genes using edgeR package
genes_to_keep <- apply(cpm(datExpr.counts) > 0.1, 1, sum) > 0.25 * ncol(datExpr.counts)
table(genes_to_keep)
# FALSE  TRUE 
#   215 20765 

# Account for library composition, normalize by library size, and log2 transform
datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr.counts[genes_to_keep, ]), method = "TMM"), log = TRUE)

boxplot(datExpr.norm, range = 0, main = "Filtered, normalized counts")

i = 1; plot(density(datExpr.norm[,i]), col = as.factor(datMeta$region)[i], main = "Gene read counts", 
          xlab = "log2(raw_counts + .001)", xlim = c(-15,30), ylim = c(0,0.45))

for (i in 2:ncol(datExpr.norm)){
    lines(density(datExpr.norm[, i]), col = as.factor(datMeta$region)[i])
}

# Calculate and visualize top 10 principle components
mds <- cmdscale(dist(t(datExpr.norm)), k = 10)
colnames(mds) <- paste0("PC", 1:ncol(mds))

pairs(mds, col = factor(datMeta$region), main = "Region", pch = 19, cex = 0.5)
par(xpd = TRUE, oma = c(1,1,1,1)); legend("topright", levels(factor(datMeta$region)), fill = 1:13, cex = .4)

pairs(mds, col = factor(datMeta$sex), main = "Sex", pch = 19, cex = 0.5)
par(xpd = TRUE, oma = c(1,1,1,1)); legend("topright", levels(factor(datMeta$sex)), fill = 1:2, cex = .4)

pairs(mds, col = factor(datMeta$AGE), main = "Age", pch = 19, cex = 0.5)
par(xpd = TRUE, oma = c(1,1,1,1)); legend("topright", levels(factor(datMeta$AGE)), fill = 1:6, cex = .5)

pairs(mds, col = factor(datMeta$DTHHRDY), main = "DeathProcess", pch = 19, cex = 0.5)
par(xpd = TRUE, oma = c(1,1,1,1)); legend("topright", levels(factor(datMeta$DTHHRDY)), fill = 1:5, cex = .5)

tree <- hclust(as.dist(1 - bicor(datExpr.norm)), "average")
plotDendroAndColors(tree, colors = cbind(labels2colors(datMeta$region), labels2colors(datMeta$sex), 
                                         labels2colors(datMeta$AGE), labels2colors(datMeta$DTHHRDY)), 
                    groupLabels = c("Region", "Sex", "Age", "DeathProcess"))

View(cbind(as.character(datMeta$region), labels2colors(datMeta$region)))


# 04. Remove outliers by connectivity Z-score
# ----------------------------------
# Calculate connectivity z-scores from the normalized data
excludesampleID <- data.frame()

par(mfrow = c(1,1))
for (i in 1:length(unique(datMeta$region))) {
  normadj <- bicor(datExpr.norm[, datMeta$region == unique(datMeta$region)[i]])
  netsummary <- fundamentalNetworkConcepts(normadj)
  C <- netsummary$Connectivity
  Z.C <- (C - mean(C)) / sqrt(var(C))
  
  # Plot the connectivity Z-score for the normalized data and set a cut-off at 3 sd
  plot(1:length(Z.C), Z.C, main = "Outlier Plot of log2(normalized_counts)", 
       xlab = "Samples", ylab = "Connectivity Z Score", col = factor(datMeta$region))
  abline(h = -3, col="red")
  
  # Determine which samples fail the "outlier threshold" via Z-score
  excludesampleID <- rbind(excludesampleID, as.data.frame(names(which(Z.C < -3))))
}

datMeta <- datMeta[-match(t(excludesampleID), datMeta$SAMPID), ] # 920 samples remaining

# Recalculate normalized expression data
datExpr.norm <- cpm(calcNormFactors(DGEList(datExpr.counts[genes_to_keep, match(datMeta$SAMPID,colnames(datExpr.norm))]), 
                                   method = "TMM"), log = TRUE)

table(datMeta$region)


# 05. Construct principal components from sequencing metrics
# ----------------------------------
seqMet <- c('SME2MPRT', 'SMCHMPRS', 'SMNTRART', 'SMMAPRT', 'SMEXNCRT',
            'SMGNSDTC', 'SME1MMRT', 'SMSFLGTH', 'SMMPPD', 
            'SMNTERRT', 'SMRRNANM', 'SMRDTTL', 'SMVQCFL', 'SMTRSCPT', 
            'SMMPPDPR', 'SMNTRNRT', 'SMMPUNRT', 'SMEXPEFF', 'SMMPPDUN', 
            'SME2MMRT', 'SME2ANTI', 'SMALTALG', 'SME2SNSE', 'SMMFLGTH', 
            'SME1ANTI', 'SMSPLTRD', 'SMBSMMRT', 'SME1SNSE', 'SME1PCTS', 
            'SMRRNART', 'SME1MPRT', 'SME2PCTS') # 'SMRDLGTH' skipped due to zero variance
datSeq <- datMeta[, seqMet] # 32 columns

seqPCs <- cmdscale(dist((scale(datSeq))), k=15, eig = TRUE)
colnames(seqPCs$points) <- paste0('seqPC', 1:15)

plot(seqPCs$points[,1], seqPCs$points[,2], main = 'PCA of scaled sequencing metric')
plot(seqPCs$points[,3], seqPCs$points[,4], main = 'PCA of scaled sequencing metric')
plot(seqPCs$points[,4], seqPCs$points[,5], main = 'PCA of scaled sequencing metric')

plot(seqPCs$eig / sum(seqPCs$eig), main = 'Variance Explained') 
# First 13 PCs explain > 99% variance, first 5 PCs 83% variance

plot(datMeta$SMRIN, datMeta$SMMPUNRT)
summary(lm(datMeta$SMMPUNRT ~ datMeta$SMRIN))

library(corrplot)
corrplot(cor(cbind(seqPCs$points, datSeq), use = "pairwise.complete.obs"), tl.cex = .5)


# 06. Regress out unwanted biological and technical covariates
# ----------------------------------
library(lme4)

# Generate covariate matrix
cov <- cbind(seqPCs$points[,1], seqPCs$points[,2], seqPCs$points[,3], seqPCs$points[,4], seqPCs$points[,5],
            datMeta[, c("SMRIN","TRISCHD","DTHCODD","DTHHRDY","DTHRFG","SEX","AGE","region","SUBJID")]) # "SMCENTER","DTHVNT" not considered
cov <- cbind(cov, seqPCs$points[,6:13])

cov[, 1] <- as.numeric(cov[, 1])
cov[, 2] <- as.numeric(cov[, 2])
cov[, 3] <- as.numeric(cov[, 3])
cov[, 4] <- as.numeric(cov[, 4])
cov[, 5] <- as.numeric(cov[, 5])
cov[, 6] <- as.numeric(cov[, 6])
cov[, 7] <- as.numeric(cov[, 7])
cov[, 8] <- as.factor(cov[, 8])
cov[, 9] <- as.factor(cov[, 9])
cov[, 10] <- as.factor(cov[, 10])
cov[, 11] <- as.factor(cov[, 11])
cov[, 12] <- as.numeric(cov[, 12])
cov[, 13] <- as.factor(cov[, 13])
cov[, 14] <- as.factor(cov[, 14])
cov[, 15] <- as.numeric(cov[, 15])
cov[, 16] <- as.numeric(cov[, 16])
cov[, 17] <- as.numeric(cov[, 17])
cov[, 18] <- as.numeric(cov[, 18])
cov[, 19] <- as.numeric(cov[, 19])
cov[, 20] <- as.numeric(cov[, 20])
cov[, 21] <- as.numeric(cov[, 21])
cov[, 22] <- as.numeric(cov[, 22])

colnames(cov) <- c("seqPC1", "seqPC2", "seqPC3", "seqPC4", "seqPC5", "SMRIN", "TRISCHD", "DTHCODD", "DTHHRDY",
                  "DTHRFG", "SEX", "AGE", "region", "SUBJID", "seqPC6", "seqPC7", "seqPC8", "seqPC9", "seqPC10",
                  "seqPC11", "seqPC12", "seqPC13")

# Regress out covariates using linear mixed model w/ a random subject intercept term
datExpr.reg <- datExpr.norm
form_dat <- "y ~ region + seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + SMRIN + TRISCHD + DTHCODD + DTHHRDY + DTHRFG + SEX + AGE + (1|SUBJID)"
mod_remove <- model.matrix(~seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + SMRIN + TRISCHD + DTHCODD + DTHHRDY + DTHRFG + SEX + AGE, cov)[, -1]

d <- c(1:nrow(datExpr.norm))
max <- 2000
d1 <- split(d, ceiling(d/max))

for (j in c(1:length(names(d1)))) {
    lmmod <- apply(as.matrix(datExpr.norm[d1[[j]], ]), 1, function(y) lmer(as.formula(form_dat), data = cov))
    
    for (i in c(1:length(d1[[j]]))) {
        if (i%%100 == 0) print(paste(i, "; Set ", j, sep=""))
        summary <- summary(lmmod[[i]])
        datExpr.reg[d1[[j]][i], ] <- datExpr.norm[d1[[j]][i], ] - mod_remove %*% coef(summary)[11:nrow(coef(summary)), 1]
    }
}

l_converge <- c(261, 307, 312, 1271, 1445, 3613, 20624, 19841, 7940, 7808, 
               7256, 12742, 13552, 12179, 12364, 10707, 10734, 8419, 8865, 9592,
               9359, 9399, 11158, 11398, 15536, 15296, 15324, 16612, 16027, 16158,
               18754, 18866) # 32 genes

# Increase iterations for genes whose estimates did not converge
for (j in l_converge) {
    print(j)
    lmmod <- lmer(datExpr.norm[j, ] ~ region + seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + 
                    seqPC9 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + SMRIN + TRISCHD + DTHCODD + DTHHRDY + DTHRFG + 
                    SEX + AGE + (1|SUBJID), data = cov)
    ss <- getME(lmmod, c("theta", "fixef"))
    lmmod <- update(lmmod, start = ss, control = lmerControl(optCtrl = list(maxfun = 1e4)))
    summary <- summary(lmmod)
    datExpr.reg[j, ] <- datExpr.norm[j, ] - mod_remove %*% coef(summary)[11:nrow(coef(summary)), 1]
}

# Save normalized and regressed expression data for downstream analyses
save(datMeta, datExpr.reg, datExpr.norm, file = paste0(path, "/data/GTEx_regressed_lmm_RSEM_final.RData"))


# 07. Spatial distribution of C4A gene expression
# ----------------------------------
datMeta <- datMeta %>% filter(!region %in% c("Substantia_nigra", "Amygdala"))
datMeta <- datMeta %>% 
  mutate(region = ifelse(region == "Anterior_cingulate_cortex_BA24", "ACC",
                        ifelse(region == "Caudate_basal_ganglia", "CDT",
                               ifelse(region == "Cerebellar_Hemisphere", "CBH",
                                      ifelse(region == "Frontal_Cortex_BA9", "BA9",
                                             ifelse(region == "Hippocampus", "HIP",
                                                    ifelse(region == "Hypothalamus", "HYP",
                                                           ifelse(region == "Nucleus_accumbens_basal_ganglia", "NAc", "PUT"))))))))

datMeta$region <- factor(datMeta$region, levels = c("BA9", "ACC", "HIP", "CDT", "PUT", "HYP", "CBH", "NAc"))

datExpr.reg <- datExpr.reg %>% as.data.frame() %>% dplyr::select(datMeta$SAMPID) %>% as.matrix()

df <- data.frame(C4A = datExpr.reg["ENSG00000244731", ],
                region = datMeta$region)

pdf(paste0(path, "/results/GTEx-expr.pdf"), width = 6.5, height = 3)
# Figure 5B

ggplot(df, aes(x = region, y = C4A)) + 
  geom_violin(aes(col = region, fill = region), trim = FALSE) + 
  geom_boxplot(outlier.size = 0, width = 0.1) + labs(x = "", y = "Residualized C4A expression") + 
  theme_bw() + theme(legend.position = "none") + ylim(-0.1, 10.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 10), 
        axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 8))

dev.off()

