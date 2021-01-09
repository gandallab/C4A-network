rm(list = ls())
options(stringsAsFactors = FALSE)

library(tidyverse)
library(cowplot)
library(gProfileR)



# 01. Effect of C4 copy number variation on gene expression ----

imputed <- read.table("./data/PsychENCODE/PsychENCODE-C4-imputed.txt", header = TRUE)
imputed$C4_CN <- imputed$C4A_CN + imputed$C4B_CN
imputed <- imputed %>% dplyr::select(-c(diagnosis, sex, study)) # datMeta already contains this information

# Load PsychENCODE expression + metadata
load("./data/PsychENCODE/PsychENCODE_regressed.RData")
datExpr <- datExpr.AllRegressed
rm(datExpr.AllRegressed)

sum(colnames(datExpr) == datMeta$X) # 2160 total samples

datMeta <- datMeta %>% filter(tissue == "frontal cortex") # 2026 frontal cortical samples
datExpr <- as.data.frame(datExpr)
datExpr <- datExpr %>% dplyr::select(datMeta$X)

datMeta <- datMeta %>% filter(individualID %in% imputed$sample) # 916 samples with C4 imputation
datExpr <- datExpr %>% dplyr::select(datMeta$X)
datMeta <- cbind(datMeta, imputed[match(datMeta$individualID, imputed$sample),])
datExpr <- as.matrix(datExpr)

# Data visualization
df <- datMeta %>% 
  dplyr::select(sample, diagnosis, structure_allele_best:C4_CN) %>% 
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

p <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8 + theme(legend.position = "none"), 
          nrow = 2, align = "v", axis = "lr",
          rel_heights = c(8,7.5))

legend <- get_legend(p8)

pdf("./results/C4-expression.pdf", width = 10, height = 6)
# Supplementary Figure 3

plot_grid(p, legend, nrow = 1, rel_widths = c(1, .1))

dev.off()



# 02. Moderation of C4 gene co-expression by C4 copy number variation ----

datMeta <- datMeta %>% mutate(pd_avg = (dose1 + dose2)/2) %>% filter(pd_avg >= 0.7) 
# 552 samples w/ high-confidence C4 imputation results
datExpr <- datExpr %>% as.data.frame() %>% dplyr::select(datMeta$X) %>% as.matrix()

table(datMeta$C4A_CN)
xtabs(~ diagnosis + C4A_CN, datMeta)
xtabs(~ study + C4A_CN, datMeta)

df <- data.frame(C4A = datExpr["ENSG00000244731",],
                 C4B = datExpr["ENSG00000224389",],
                 SLC39A10 = datExpr["ENSG00000196950",], datMeta)

lm.fit <- lm(SLC39A10 ~ C4A * (pd_c4a + pd_c4b), df)
summary(lm.fit)

library(emmeans)
mylist <- list(C4A = seq(-2, 4, by = 0.01), pd_c4a = c(0, 1, 2, 3, 4))
predict <- emmip(lm.fit, pd_c4a ~ C4A, at = mylist, CIs = TRUE, plotit = FALSE)

df$C4A_CN <- as.factor(df$C4A_CN)
predict$pd_c4a <- as.factor(predict$pd_c4a)

ggplot(df, aes(x = C4A, y = SLC39A10)) + 
  geom_point(size = 1, alpha = 0.4, aes(col = C4A_CN)) + 
  geom_line(data = predict, aes(x = C4A, y = yvar, col = pd_c4a)) + 
  geom_ribbon(data = predict, aes(x = C4A, y = yvar, ymax = UCL, ymin = LCL, fill = pd_c4a), alpha = 0.4) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = expression(paste(italic("C4A"), " expression")), 
       y = expression(paste(italic("SLC39A10"), " expression"))) +
  guides(fill = FALSE) + ylim(5, 8.5) 



# 03. Constructing C4A gene co-expression network using control samples ----

table(datMeta$C4A_CN)
# 0   1   2   3   4 
# 9 110 324  99  10 

datExpr_mid <- datExpr[, datMeta$C4A_CN == 2 & datMeta$Group == "CTL"] # 145 samples

prsCor = function(i, gene, datExpr) {
  c = cor.test(datExpr[i, ], datExpr[gene, ], use = "pairwise.complete.obs")
  dfPrs = data.frame(Seed = gene, Gene = rownames(datExpr)[i], R = c$estimate, P = c$p.value)
  return(dfPrs)
}

dfC4A_mid <- data.frame()

for (i in 1:nrow(datExpr)) {
  if (i%%100 == 0) {print(i)}
  dfC4A_mid <- rbind(dfC4A_mid, prsCor(i, "ENSG00000244731", datExpr_mid))
}

hist(dfC4A_mid$P, breaks = 50)

dfC4A_mid$FDR <- p.adjust(dfC4A_mid$P, "fdr")

gencode <- read.csv("./data/annotation.gene.gencodeV19.csv")
gencode <- gencode[match(rownames(datExpr), gencode$gene_id),]

dfC4A_mid <- cbind(gene = gencode$gene_name, dfC4A_mid)
xtabs(~ (FDR < 0.05) + (R > 0), dfC4A_mid)
#           R > 0
#FDR < 0.05 FALSE  TRUE
#     FALSE 12335 10418
#     TRUE   1152  1869



# 04. Enrichment patterns ----

source("./code/04_Fisher-exact-test.R")

complement = read.table("./results/complement-gene-sets.tsv", sep = "\t", header = TRUE)
complement = complement$gene[complement$set == "Complement (57 genes)"]

# Over-representation of complement system
ORA(dfC4A_mid$gene[dfC4A_mid$FDR < .05 & dfC4A_mid$R > 0], complement,
    dfC4A_mid$gene, dfC4A_mid$gene)
#OR               Fisher p              
#17.2141244653363 4.03102279425738e-17

ORA(dfC4A_mid$gene[dfC4A_mid$FDR < .05 & dfC4A_mid$R < 0], complement,
    dfC4A_mid$gene, dfC4A_mid$gene)
#OR Fisher p              
#0  0.26

# SynGO enrichment
synGO = readxl::read_xlsx("./data/GSEA-annotations/syngo_annotations.xlsx")

ORA(dfC4A_mid$gene[dfC4A_mid$FDR < .05 & dfC4A_mid$R > 0], unique(synGO$`human ortholog gene symbol`),
    dfC4A_mid$gene, dfC4A_mid$gene)
#OR               Fisher p
#1.14730972759516 0.224729415660388 

ORA(dfC4A_mid$gene[dfC4A_mid$FDR < .05 & dfC4A_mid$R < 0], unique(synGO$`human ortholog gene symbol`),
    dfC4A_mid$gene, dfC4A_mid$gene)
#OR              Fisher p
#3.7725740859901 1.23567031237039e-34 

# GO enrichment
df <- dfC4A_mid[dfC4A_mid$FDR < .05 & dfC4A_mid$R < 0, ]

go.down <- gprofiler(query = df$gene[order(df$R)], 
                     correction_method = "fdr", ordered_query = TRUE, custom_bg = dfC4A_mid$gene,
                     min_set_size = 10, max_set_size = 1000, src_filter = "GO")

go.down <- data.frame(Set = "C4A-negative", go.down)

df <- dfC4A_mid[dfC4A_mid$FDR < .05 & dfC4A_mid$R > 0, ]

go.up <- gprofiler(query = df$gene[order(df$R, decreasing = TRUE)], 
                   correction_method = "fdr", ordered_query = TRUE, custom_bg = dfC4A_mid$gene,
                   min_set_size = 10, max_set_size = 1000, src_filter = "GO")
                   
go.up <- data.frame(Set = "C4A-positive", go.up)

#write.csv(file = "./results/C4A-CN2-CTL-GO.csv", rbind(go.down, go.up))


