rm(list = ls())
options(stringsAsFactors = FALSE)
library(tidyverse)
path <- "/Users/minsookim/Desktop/C4A-network"


# 01. Effect of C4 copy number variation on gene expression
# ----------------------------------
imputed <- read.table(paste0(path, "/data/PsychENCODE-C4-imputed.txt"), header = TRUE)
imputed$C4_CN <- imputed$C4A_CN + imputed$C4B_CN
imputed <- imputed %>% select(-c(diagnosis, sex, study)) # datMeta already contains this information

# Load PsychENCODE expression + metadata files (see README file for ways to download these) 
load(paste0(path,"/data/PsychENCODE.RData"))
sum(colnames(datExpr) == datMeta$X) # 2160 total samples

datMeta <- datMeta %>% filter(tissue == "frontal cortex") # 2026 frontal cortical samples
datExpr <- as.data.frame(datExpr)
datExpr <- datExpr %>% select(datMeta$X)

datMeta <- datMeta %>% filter(individualID %in% imputed$sample) # 916 samples with C4 imputation
datExpr <- datExpr %>% select(datMeta$X)
datMeta <- cbind(datMeta, imputed[match(datMeta$individualID, imputed$sample),])
datExpr <- as.matrix(datExpr)

# Data visualization
library(cowplot)

df <- datMeta %>% 
  select(sample, diagnosis, structure_allele_best:C4_CN) %>% 
  bind_cols(C4A = datExpr["ENSG00000244731",],
            C4B = datExpr["ENSG00000224389",])

df <- df %>% 
  mutate(diagnosis = ifelse(diagnosis == "Autism Spectrum Disorder", "ASD",
                            ifelse(diagnosis == "Bipolar Disorder", "BD",
                                   ifelse(diagnosis == "Control", "CTL", "SCZ"))))

p1 <- ggplot(df, aes(x = as.factor(C4A_CN), y = C4A)) + 
  geom_boxplot(outlier.size = 0, width = 0.6) + 
  geom_jitter(width = 0.2, alpha = 0.7, size = 0.6, aes(col = diagnosis)) + 
  labs(x = expression(paste(italic("C4A"), " copy number")), 
       y = expression(paste("Residualized ", italic("C4A"), " expression"))) + theme_classic() + 
  theme(axis.title.x = element_blank(), panel.grid = element_blank(), legend.position = "none")

p2 <- ggplot(df, aes(x = as.factor(C4HERV_CN), y = C4A)) + 
  geom_boxplot(outlier.size = 0, width = 0.6) + 
  geom_jitter(width = 0.2, alpha = 0.7, size = 0.6, aes(col = diagnosis)) + 
  labs(x = expression(paste(italic("C4"), "-HERV copy number")), 
       y = expression(paste("Residualized ", italic("C4A"), " expression"))) + theme_classic() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")

p3 <- ggplot(df, aes(x = as.factor(C4B_CN), y = C4A)) + 
  geom_boxplot(outlier.size = 0, width = 0.6) + 
  geom_jitter(width = 0.2, alpha = 0.7, size = 0.6, aes(col = diagnosis)) + 
  labs(x = expression(paste(italic("C4B"), " copy number")), 
       y = expression(paste("Residualized ", italic("C4A"), " expression"))) + theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")

p4 <- ggplot(df, aes(x = as.factor(C4_CN), y = C4A)) + 
  geom_boxplot(outlier.size = 0, width = 0.6) + 
  geom_jitter(width = 0.2, alpha = 0.7, size = 0.6, aes(col = diagnosis)) + 
  labs(x = expression(paste("Total ", italic("C4")," copy number")), 
       y = expression(paste("Residualized ", italic("C4A"), " expression"))) + theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")

p5 <- ggplot(df, aes(x = as.factor(C4A_CN), y = C4B)) + 
  geom_boxplot(outlier.size = 0, width = 0.6) + 
  geom_jitter(width = 0.2, alpha = 0.7, size = 0.6, aes(col = diagnosis)) + 
  labs(x = expression(paste(italic("C4A"), " copy number")), 
       y = expression(paste("Residualized ", italic("C4B"), " expression"))) + theme_classic() +
  theme(legend.position = "none")

p6 <- ggplot(df, aes(x = as.factor(C4HERV_CN), y = C4B)) + 
  geom_boxplot(outlier.size = 0, width = 0.6) + 
  geom_jitter(width = 0.2, alpha = 0.7, size = 0.6, aes(col = diagnosis)) + 
  labs(x = expression(paste(italic("C4"), "-HERV copy number")), 
       y = expression(paste("Residualized ", italic("C4B"), " expression"))) + theme_classic() + 
  theme(axis.title.y = element_blank(), legend.position = "none")

p7 <- ggplot(df, aes(x = as.factor(C4B_CN), y = C4B)) + 
  geom_boxplot(outlier.size = 0, width = 0.6) + 
  geom_jitter(width = 0.2, alpha = 0.7, size = 0.6, aes(col = diagnosis)) + 
  labs(x = expression(paste(italic("C4B"), " copy number")), 
       y = expression(paste("Residualized ", italic("C4B"), " expression"))) + theme_classic() + 
  theme(axis.title.y = element_blank(), legend.position = "none")

p8 <- ggplot(df, aes(x = as.factor(C4_CN), y = C4B)) + 
  geom_boxplot(outlier.size = 0, width = 0.6) + 
  geom_jitter(width = 0.2, alpha = 0.7, size = 0.6, aes(col = diagnosis)) + 
  labs(x = expression(paste("Total ", italic("C4"), " copy number")), 
       y = expression(paste("Residualized ", italic("C4B"), " expression"))) + theme_classic() +
  theme(axis.title.y = element_blank())

library(cowplot)
p <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8 + theme(legend.position="none"), 
          nrow=2, align="v", axis="lr",
          rel_heights = c(8,7.5))

legend <- get_legend(p8)

pdf(paste0(path, "/results/C4-expression.pdf"), width = 10, height = 6)
# Supplementary Figure 3

plot_grid(p, legend, nrow=1, rel_widths=c(1, .1))

dev.off()

