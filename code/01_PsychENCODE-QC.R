rm(list = ls())
options(stringsAsFactors = FALSE)
library(tidyverse)
path <- "/Users/minsookim/Desktop/C4A-network"


# 01. Load PsychENCODE data for primary analysis
# ----------------------------------
library(SummarizedExperiment)
load(paste0(path, "/data/se.CombinedDataForPrimaryAnalysis.RData"))

datMeta <- colData(se.Primary)

datMeta$RIN.squared <- se.Primary$RIN^2
datMeta$seqPC2.squared <- se.Primary$seqPC2^2
datMeta$seqPC3.squared <- se.Primary$seqPC3^2
datMeta$age.squared <- se.Primary$age^2

datExpr <- assays(se.Primary)$counts
datMeta$sex <- factor(datMeta$sex, levels = c("M","F"))


# 02. Construct model matrix for covariates
# ----------------------------------
mod1 <- model.matrix(~ Group + age + age.squared + study + sex + PMI + RIN + RIN.squared + 
                      individualIDSource + tissue + seqPC3.squared + seqPC1 + seqPC10 + seqPC11 + 
                      seqPC12 + seqPC13 + seqPC14 + seqPC16 + seqPC18 + seqPC19 + seqPC2 + seqPC20 + 
                      seqPC21 + seqPC22 + seqPC23 + seqPC24 + seqPC25 + seqPC27 + seqPC28 + seqPC29 + 
                      seqPC3 + seqPC5 + seqPC6 + seqPC7 + seqPC8, data = datMeta)
colnames(mod1) <- make.names(colnames(mod1))

# Remove redundant terms
mod1 <- mod1[, !colnames(mod1) %in% c("individualIDSourceNICHD", "individualIDSourcePitt", 
                                     "individualIDSourceSMRI..New.", "individualIDSourceYale.Pathology")] 

mod0 <- mod1[, !grepl("Group", colnames(mod1))] # Make a baseline model with all terms but group

library(edgeR)
dge.voom <- voom(calcNormFactors(DGEList(datExpr), method = "TMM"), mod1, plot = TRUE)

datExpr <- dge.voom$E

datMeta$individualIDSource <- as.character(datMeta$individualIDSource)
datMeta$individualIDSource[datMeta$individualIDSource %in% 
                             c("NICHD", "Pitt", "SMRI \"New\"", "Yale Pathology")] <- ""
datMeta$individualIDSource <- factor(datMeta$individualIDSource)


# 03. Regress out covariates except the diagnosis effect
# ----------------------------------
Y <- as.matrix(dge.voom$E)
X <- mod1
beta <- (solve(t(X) %*% X) %*% t(X)) %*% t(Y)
datExpr.AllRegressed <- Y - t(X[, c(5:ncol(X))] %*% beta[c(5:ncol(X)), ])

datMeta <- as.data.frame(datMeta)
datMeta$X <- rownames(datMeta)

save(datMeta, datExpr.AllRegressed, file = paste0(path, "/data/PsychENCODE_regressed.RData"))


# 04. Effect of age and sex on C4A gene expression
# ----------------------------------
summary(lm(datExpr["ENSG00000244731", ] ~ Group + age + age.squared + study + sex + PMI + RIN + RIN.squared + 
             individualIDSource + tissue + seqPC3.squared + seqPC1 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + 
             seqPC14 + seqPC16 + seqPC18 + seqPC19 + seqPC2 + seqPC20 + seqPC21 + seqPC22 + seqPC23 + seqPC24 + 
             seqPC25 + seqPC27 + seqPC28 + seqPC29 + seqPC3 + seqPC5 + seqPC6 + seqPC7 + seqPC8, data = datMeta))

datExpr.AllRegressed.NotAge <- Y - t(X[, c(7:ncol(X))] %*% beta[c(7:ncol(X)), ])
df <- data.frame(C4A = datExpr.AllRegressed.NotAge["ENSG00000244731", datMeta$tissue == "frontal cortex"], 
                 datMeta[datMeta$tissue == "frontal cortex", ])

df <- df[df$diagnosis == "Control" | df$diagnosis == "Schizophrenia", ]
df$diagnosis <- factor(df$diagnosis, levels = c("Schizophrenia","Control")) # 593 SCZ, 1137 CTL samples

pdf(paste0(path, "/results/age-expr.pdf"), width = 9.0, height = 3.7)
# Figure 5D

ggplot(df, aes(x = age, y = C4A, col = diagnosis)) + 
  geom_point(size = 1, alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()

