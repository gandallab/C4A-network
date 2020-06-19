source("code/04_Fisher-exact-test.R")
options(stringsAsFactors = FALSE)
library(tidyverse)
path <- "~/Github//C4A-network"


# 01. Effect of C4 copy number variation on gene expression
# ----------------------------------
imputed <- read.table(paste0(path, "/data/PsychENCODE-C4-imputed.txt"), header = TRUE)
imputed$C4_CN <- imputed$C4A_CN + imputed$C4B_CN
imputed <- imputed %>% dplyr::select(-c(diagnosis, sex, study)) # datMeta already contains this information

# Load PsychENCODE expression + metadata
load(paste0(path,"/data/PsychENCODE_regressed.RData"))
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

p <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8 + theme(legend.position="none"), 
          nrow=2, align="v", axis="lr",
          rel_heights = c(8,7.5))

legend <- get_legend(p8)

pdf(paste0(path, "/results/C4-expression.pdf"), width = 10, height = 6)
# Supplementary Figure 3

plot_grid(p, legend, nrow=1, rel_widths=c(1, .1))

dev.off()


# 02. Moderation of C4 gene co-expression by C4 copy number variation
# ----------------------------------
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


# 03. Constructing C4A gene co-expression using all samples
# ----------------------------------
table(datMeta$C4A_CN)
# 0   1   2   3   4 
# 9 110 324  99  10 

datExpr_low <- datExpr[, datMeta$C4A_CN < 2]
datExpr_mid.ctl <- datExpr[, datMeta$C4A_CN == 2 & datMeta$Group == "CTL"]
datExpr_mid <- datExpr[, datMeta$C4A_CN == 2]
datExpr_hig <- datExpr[, datMeta$C4A_CN > 2]

prsCor = function(i, gene, datExpr) {
  c = cor.test(datExpr[i, ], datExpr[gene, ], use = "pairwise.complete.obs")
  dfPrs = data.frame(Seed = gene, Gene = rownames(datExpr)[i], R = c$estimate, P = c$p.value)
  return(dfPrs)
}

dfC4A_low <- data.frame()
dfC4A_mid <- data.frame()
dfC4A_mid.ctl = data.frame()
dfC4A_hig <- data.frame()

for (i in 1:nrow(datExpr)) {
  if (i%%100 == 0) {print(i)}
  
  dfC4A_low <- rbind(dfC4A_low, prsCor(i, "ENSG00000244731", datExpr_low))
  dfC4A_mid <- rbind(dfC4A_mid, prsCor(i, "ENSG00000244731", datExpr_mid))
  dfC4A_mid.ctl <- rbind(dfC4A_mid.ctl, prsCor(i, "ENSG00000244731", datExpr_mid.ctl))
  dfC4A_hig <- rbind(dfC4A_hig, prsCor(i, "ENSG00000244731", datExpr_hig))
}

hist(dfC4A_low$P)
hist(dfC4A_mid$P)
hist(dfC4A_mid.ctl$P)
hist(dfC4A_hig$P)

dfC4A_low$FDR = p.adjust(dfC4A_low$P, "fdr")
dfC4A_mid$FDR = p.adjust(dfC4A_mid$P, "fdr")
dfC4A_mid.ctl$FDR = p.adjust(dfC4A_mid.ctl$P, "fdr")
dfC4A_hig$FDR = p.adjust(dfC4A_hig$P, "fdr")


gencode = read.csv("data/annotation.gene.gencodeV19.csv")
gencode = gencode[match(rownames(datExpr), gencode$gene_id),]

# Control only C4A seeded networks
this_df = cbind(gene=gencode$gene_name, dfC4A_mid.ctl)
complement_system = readxl::read_xlsx("results/manuscript/TableS1.xlsx",sheet = 1)

xtabs(~ (FDR < 0.05) + (R > 0), this_df)


# Complement Enrichment
ORA(this_df$gene[this_df$FDR < .05 & this_df$R > 0], complement_system$`Approved symbol`,
    this_df$gene, this_df$gene)
#OR               Fisher p              
# "17.2141244653363" "4.03102279425738e-17"

ORA(this_df$gene[this_df$FDR < .05 & this_df$R < 0], complement_system$`Approved symbol`,
    this_df$gene, this_df$gene)
#OR               Fisher p              
# 0   0.26

# SynGO enrichment
synGO = readxl::read_xlsx("data/GSEA-annotations/syngo_annotations.xlsx")

ORA(this_df$gene[this_df$FDR < .05 & this_df$R > 0], unique(synGO$`human ortholog gene symbol`),
    this_df$gene, this_df$gene)
ORA(this_df$gene[this_df$FDR < .05 & this_df$R < 0], unique(synGO$`human ortholog gene symbol`),
    this_df$gene, this_df$gene)
#          OR               Fisher p                 -95%CI                 +95%CI 
# "3.7725740859901" "1.23567031237039e-34"     "3.10861053939515"     "4.55368635278887" 

write.csv(file="results/C4A-coexpressed-CN2-ctl-145.csv", this_df)

library(gProfileR)
go.down = gprofiler(query = this_df$gene[this_df$FDR < .05 & this_df$R < 0], correction_method = 'fdr', max_set_size = 1000,
                           organism = 'hsapiens', custom_bg = this_df$gene, src_filter = c("GO", "KEGG", "REAC"))
go.down = data.frame(Set="C4A-negative", go.down)
go.up = gprofiler(query = this_df$gene[this_df$FDR < .05 & this_df$R > 0],  correction_method = 'fdr', max_set_size = 1000, 
                           organism = 'hsapiens', custom_bg = this_df$gene, src_filter = c("GO", "KEGG", "REAC"))
go.up = data.frame(Set="C4A-positive", go.up)

write.csv(file="results/C4A-coexpressed-CN2-ctl-145-GO.csv", rbind(go.down, go.up))
save(file="data/data_for_network_permutation.RData", datExpr_mid.ctl, gencode, complement_system, synGO)



# 04. Test Boostrap
# ----------------------------------
library(parallel)
runOneBootstrap = function(gene) {
  this_net= do.call("rbind", mclapply(1:nrow(datExpr), prsCor, gene, datExpr_mid.ctl,mc.cores = 6))
  this_net$Gene = gencode$gene_name[match(this_net$Gene, gencode$gene_id)]
  this_net$FDR = p.adjust(this_net$P,'fdr')
  up_fdr = this_net$Gene[this_net$FDR < .05 & this_net$R > 0]
  down_fdr = this_net$Gene[this_net$FDR < .05 & this_net$R < 0]
  fisher_up_complement = ORA(up_fdr, complement_system$`Approved symbol`,
                             gencode$gene_name, gencode$gene_name)
  fisher_down_synGO = ORA(down_fdr, unique(synGO$`human ortholog gene symbol`),
                             gencode$gene_name, gencode$gene_name)
  dfBootstrap = data.frame(seed = gencode$gene_name[match(gene, gencode$gene_id)],
                           up_fdr= length(up_fdr), down_fdr = length(down_fdr),
                           up_complement_OR = as.numeric(fisher_up_complement[[1]]),
                           up_complement_P = as.numeric(fisher_up_complement[[2]]),
                           down_synGO_OR = as.numeric(fisher_down_synGO[[1]]),
                           down_synGO_P = as.numeric(fisher_down_synGO[[2]]))
  return(dfBootstrap)
}


df_boot = do.call("rbind", lapply(as.list(rownames(datExpr_mid.ctl)[1:100]), runOneBootstrap))
df_boot100to1000 = do.call("rbind", lapply(as.list(rownames(datExpr_mid.ctl)[101:1000]), runOneBootstrap))

