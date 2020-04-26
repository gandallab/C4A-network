rm(list = ls())
options(stringsAsFactors = FALSE)
library(tidyverse)
path <- "/Users/minsookim/Desktop/C4A-network"


# 01. Examine seeded network size wrt C4A copy number
# ----------------------------------
load(paste0(path, "/data/C4A-network.RData"))
network <- network %>% as_tibble() %>% print(width = Inf)

df <- tibble()
for (i in seq(0.3, 0.99, by = 0.01)) {
  df <- network %>% 
    summarise(PCC = i,
              "CN < 2" = sum(abs(`01.R`) > i),
              "CN = 1" = sum(abs(`1.R`) > i),
              "CN = 2" = sum(abs(`2.R`) > i),
              "CN = 3" = sum(abs(`3.R`) > i),
              "CN > 2" = sum(abs(`34.R`) > i)) %>% 
    bind_rows(df, .)
}

df <- df %>% 
  gather(`CN < 2`, `CN = 1`, `CN = 2`, `CN = 3`, `CN > 2`, key = "CN", value = "value")

df$CN <- factor(df$CN, levels = c("CN > 2", "CN = 3", "CN = 2", "CN = 1", "CN < 2"))

jpeg(paste0(path, "/results/network-size.jpeg"), units = "in", width = 5.5, height = 5, res = 300)
# Supplementary Figure 5

df %>% 
  ggplot(aes(x = PCC, y = (value))) + 
  geom_point(aes(color = CN), size = 0.6, shape = 19) + # geom_line(aes(0.5), col = "black") + 
  labs(x = "PCC cutoff", y = "# of genes passing threshold") +
  theme_classic() + theme(legend.position = c(0.85,0.5)) + 
  guides(color = guide_legend(expression(paste(italic("C4A"), " CN"))))

dev.off()

cor.test(network$`01.R`, network$`2.R`) # 0.74
cor.test(network$`01.R`, network$`34.R`) # 0.65
cor.test(network$`2.R`, network$`34.R`) # 0.83


# 02. PsychENCODE WGCNA module enrichment
# ----------------------------------
source(paste0(path, "/code/04_Fisher-exact-test.R"))

PEmodules <- read.table(paste0(path, "/data/PsychENCODE-WGCNA-modules.txt"), header = TRUE)

# Test top 500 positively co-expressed genes
geneset <- network %>% arrange(desc(`34.R`)) %>% select(Gene) %>% dplyr::slice(1:500) %>% t()

df <- tibble()
for (i in unique(PEmodules$Module)) {
  f <- ORA(testpath = PEmodules$Gene[PEmodules$Module == i], refpath = geneset, 
          testbackground = network$Gene, refbackground = network$Gene)
  df <- bind_rows(df, tibble(Module = i, OR = as.numeric(f[[1]]), P = as.numeric(f[[2]])))
}

df$P.adj <- p.adjust(df$P, method = "bonferroni")
df <- df %>% 
  filter(Module != "geneM0")

if (nrow(df) > 5) {
  df <- df %>% 
    arrange(desc(OR)) %>% dplyr::slice(1:5)
}

ast_q <- rep("", nrow(df))
ast_q[df$P.adj < 0.05 & df$OR > 1] <- "*"
df$ast_q <- ast_q

jpeg(paste0(path, "/results/module-pos.jpeg"), units = "in", width = 3.5, height = 3, res = 300)
# Figure 2C

df %>% 
  ggplot(aes(x = reorder(Module, -log10(P)), y = OR)) + 
  geom_bar(stat = "identity", fill = "#F8766D") + coord_flip() + xlab("") + theme_bw() + 
  geom_text(aes(label = ast_q), size = 6, vjust = 0.75, hjust = -0.1,
            position = position_dodge(width = 1))

dev.off()


# 03. gProfiler pathway enrichment
# ----------------------------------
library(gProfileR)

# Test top 500 negatively co-expressed genes
geneset <- network %>% arrange(`34.R`) %>% dplyr::slice(1:500)

go <- gprofiler(query = geneset$Gene, correction_method = "fdr", ordered_query = TRUE, 
               custom_bg = network$Gene, min_set_size = 10, max_set_size = 1000, src_filter="GO", 
               hier_filtering = "strong")

if (nrow(go) > 5) {
  go <- go %>% 
    arrange(p.value) %>% dplyr::slice(1:5)
}

jpeg(paste0(path, "/results/pathway-neg.jpeg"), units = "in", width = 4, height = 3, res = 300)
# Supplementary Figure 8

go %>% 
  ggplot(aes(x = reorder(term.name, -log10(p.value)), y = -log10(p.value))) + 
  geom_bar(stat = "identity", fill = "royalblue") + coord_flip() + xlab("") + theme_bw() + 
  geom_hline(yintercept = -log10(0.05), lty = 2, color = "red")

dev.off()


# 04. Expression-weighted cell-type enrichment (EWCE)
# ----------------------------------
library(EWCE)

# Pre-calculated expression specificity metrics downloaded from www.hjerling-leffler-lab.org/data/scz_singlecell
load(paste0(path, "/data/ctd_allKI.rda"))

# Tidy column names
colnames(ctd[[1]]$mean_exp) = colnames(ctd[[1]]$specificity) = 
  c("Astrocyte / Ependymal", "Dopaminergic Neuron", "Dopaminergic Neuroblast", "Embryonic Dopaminergic Neuron",
    "Embryonic GABAergic Neuron", "Embryonic Midbrain Nucleus Neuron", "Endothelial-Mural", 
    "Hypothalamic Dopaminergic Neuron", "Hypothalamic GABAergic Neuron", "Hypothalamic Glutamatergic Neuron",
    "Cortical Interneuron", "Medium Spiny Neuron", "Microglia", "Neuronal Progenitor", "Neuroblast",
    "Oligodendrocyte Precursor", "Oligodendrocyte", "Oxytocin/Vasopressin Expressing Neuron", "Pyramidal (CA1)",
    "Pyramidal (SS)", "Radial glia-like", "Serotonergic Neuron", "Striatal Interneuron", 
    "Vascular Leptomeningeal Cell")

# Test top 500 negatively co-expressed genes
geneset <- network %>% arrange(`34.R`) %>% select(Gene) %>% dplyr::slice(1:500) %>% t()

data("mouse_to_human_homologs")
m2h <- unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])
mouse.hits <- unique(m2h[m2h$HGNC.symbol %in% geneset,"MGI.symbol"]) # 435 genes
mouse.bg <- unique(m2h$MGI.symbol) # 15604 genes in background set

reps <- 10000
level <- 1

full_results <- bootstrap.enrichment.test(sct_data = ctd, hits = mouse.hits, bg = mouse.bg, 
                                         reps = reps, annotLevel = level)

jpeg(paste0(path, "/results/cell-type-neg.jpeg"), units = "in", width = 20, height = 6, res = 300)
# Supplementary Figure 10

ewce.plot(full_results$results, mtc_method = "BH")

dev.off()

