rm(list = ls())
options(stringsAsFactors = FALSE)

library(tidyverse)
library(lme4)
library(pbkrtest)



# 01. Effect of smoking and other relevant non-genetic factors on C4A gene expression ----

# Load data and metadata
load("./data/GTEx/GTEx_regressed_lmm_RSEM_final.RData")

# Load metadata from GTEx v8
datMeta_v8 <- read_csv("./data/GTEx/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.csv")
# phs000424.v8.pht002742.v8.GTEx_Subject_Phenotypes.data_dict.csv is the corresponding key file

# Merge metadata v7 and v8
datMeta <- bind_cols(datMeta, 
                     datMeta_v8[match(datMeta$SUBJID, datMeta_v8$SUBJID), 
                        is.na(match(colnames(datMeta_v8), colnames(datMeta)))])
rm(datMeta_v8)

# Construct principal components from sequencing metrics
seqMet <- c('SME2MPRT', 'SMCHMPRS', 'SMNTRART', 'SMMAPRT', 'SMEXNCRT',
            'SMGNSDTC', 'SME1MMRT', 'SMSFLGTH', 'SMMPPD', 
            'SMNTERRT', 'SMRRNANM', 'SMRDTTL', 'SMVQCFL', 'SMTRSCPT', 
            'SMMPPDPR', 'SMNTRNRT', 'SMMPUNRT', 'SMEXPEFF', 'SMMPPDUN', 
            'SME2MMRT', 'SME2ANTI', 'SMALTALG', 'SME2SNSE', 'SMMFLGTH', 
            'SME1ANTI', 'SMSPLTRD', 'SMBSMMRT', 'SME1SNSE', 'SME1PCTS', 
            'SMRRNART', 'SME1MPRT', 'SME2PCTS') # 'SMRDLGTH' skipped due to zero variance
datSeq <- datMeta[, seqMet] # 32 columns

seqPCs <- cmdscale(dist((scale(datSeq))), k = 15, eig = TRUE)
colnames(seqPCs$points) <- paste0('seqPC', 1:15)

plot(seqPCs$eig / sum(seqPCs$eig), main = 'Variance Explained') 
# First 13 PCs explain > 99% variance, first 5 PCs 83% variance

# Generate covariate matrix
cov <- cbind(seqPCs$points[, 1], seqPCs$points[, 2], seqPCs$points[, 3], seqPCs$points[, 4], seqPCs$points[, 5],
             datMeta[, c("SMRIN", "TRISCHD", "DTHCODD", "DTHHRDY", "DTHRFG", "SEX", "AGE", "region", "SUBJID")],
             seqPCs$points[, 6:13]) # "SMCENTER","DTHVNT" not considered

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

cov <- cov %>% as_tibble() %>% bind_cols(MHSMKSTS = as.factor(datMeta$MHSMKSTS),
                                         MHDRNKSTS = as.factor(datMeta$MHDRNKSTS),
                                         WGHT = as.numeric(datMeta$WGHT),
                                         HGHT = as.numeric(datMeta$HGHT),
                                         BMI = as.numeric(datMeta$BMI),
                                         C4A = datExpr.norm["ENSG00000244731", ],
                                         C4B = datExpr.norm["ENSG00000224389", ],
                                         C1 = as.numeric(datMeta$C1),
                                         C2 = as.numeric(datMeta$C2),
                                         C3 = as.numeric(datMeta$C3),
                                         MHBRNPH = as.numeric(datMeta$MHBRNPH))

# Linear mixed model using lme4 package 

# LRT 
cov1 <- cov[!is.na(cov$MHDRNKSTS), ]

mod1 <- lmer(C4A ~ region + seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + 
               seqPC10 + seqPC11 + seqPC12 + seqPC13 + SMRIN + TRISCHD + DTHCODD + DTHHRDY + DTHRFG + SEX + 
               AGE + C1 + C2 + C3 + MHBRNPH + (1|SUBJID), data = cov1, REML = FALSE)

mod2 <- lmer(C4A ~ MHSMKSTS + region + seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + 
               seqPC10 + seqPC11 + seqPC12 + seqPC13 + SMRIN + TRISCHD + DTHCODD + DTHHRDY + DTHRFG + SEX + 
               AGE + C1 + C2 + C3 + MHBRNPH + (1 |SUBJID), data = cov1, REML = FALSE)

anova(mod1, mod2) # p = 0.00639

as.numeric(2 * (logLik(mod2, REML = FALSE) - logLik(mod1, REML = FALSE)))
pchisq(7.4368, 1, lower.tail = FALSE) # p = 0.00639

# Kenward-Roger adjusted F-tests for REML
mod1 <- lmer(C4A ~ region + seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + 
               seqPC10 + seqPC11 + seqPC12 + seqPC13 + SMRIN + TRISCHD + DTHCODD + DTHHRDY + DTHRFG + SEX + 
               AGE + C1 + C2 + C3 + MHBRNPH + (1|SUBJID), data = cov1, REML = TRUE)

mod2 <- lmer(C4A ~ MHSMKSTS + region + seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + 
               seqPC10 + seqPC11 + seqPC12 + seqPC13 + SMRIN + TRISCHD + DTHCODD + DTHHRDY + DTHRFG + SEX + 
               AGE + C1 + C2 + C3 + MHBRNPH + (1 |SUBJID), data = cov1, REML = TRUE)

KRmodcomp(mod2, mod1) # p = 0.01136

# Parametric bootstrap
set.seed(123)
B <- 1000
lrstat <- numeric(B)
for (i in 1:B) {
  ryield <- unlist(simulate(mod1))
  nmodr <- suppressMessages(lmer(ryield ~ region + seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + 
                                   seqPC10 + seqPC11 + seqPC12 + seqPC13 + SMRIN + TRISCHD + DTHCODD + DTHHRDY + DTHRFG + SEX + 
                                   AGE + C1 + C2 + C3 + MHBRNPH + (1|SUBJID), data = cov1, REML = FALSE))
  amodr <- suppressMessages(lmer(ryield ~ MHSMKSTS + region + seqPC1 + seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + seqPC9 + 
                                   seqPC10 + seqPC11 + seqPC12 + seqPC13 + SMRIN + TRISCHD + DTHCODD + DTHHRDY + DTHRFG + SEX + 
                                   AGE + C1 + C2 + C3 + MHBRNPH + (1|SUBJID), data = cov1, REML = FALSE))
  lrstat[i] <- 2 * (logLik(amodr, REML = FALSE) - logLik(nmodr, REML = FALSE))
}

(pval <- mean(lrstat > 7.4368)) # p = 0.011

# Three ways of assessing significance in lmm all give the same consistent results. 
# Note that smoking status itself is not an ideal covariate that tracks the degree of smoking well (i.e. pack-years preferable). 
