rm(list = ls())
options(stringsAsFactors = FALSE)

library(tidyverse)
library(gProfileR)
library(EWCE)
library(WGCNA)
library(igraph)
library(ggrepel)



# 01. Examine seeded network size wrt C4A copy number ----

load("./data/PsychENCODE/C4A-network.RData")
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

jpeg("./results/network-size.jpeg", units = "in", width = 5.5, height = 5, res = 300)
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



# 02. PsychENCODE WGCNA module enrichment ----

source("./code/04_Fisher-exact-test.R")

PEmodules <- read.table("./data/PsychENCODE/PsychENCODE-WGCNA-modules.txt", header = TRUE)

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

jpeg("./results/module-pos.jpeg", units = "in", width = 3.5, height = 3, res = 300)
# Figure 2C

df %>% 
  ggplot(aes(x = reorder(Module, -log10(P)), y = OR)) + 
  geom_bar(stat = "identity", fill = "#F8766D") + coord_flip() + xlab("") + theme_bw() + 
  geom_text(aes(label = ast_q), size = 6, vjust = 0.75, hjust = -0.1,
            position = position_dodge(width = 1))

dev.off()



# 03. gProfiler pathway enrichment ----

# Test top 500 negatively co-expressed genes
geneset <- network %>% arrange(`34.R`) %>% dplyr::slice(1:500)

go <- gprofiler(query = geneset$Gene, correction_method = "fdr", ordered_query = TRUE, 
               custom_bg = network$Gene, min_set_size = 10, max_set_size = 1000, src_filter="GO", 
               hier_filtering = "strong")

if (nrow(go) > 5) {
  go <- go %>% 
    arrange(p.value) %>% dplyr::slice(1:5)
}

jpeg("./results/pathway-neg.jpeg", units = "in", width = 4, height = 3, res = 300)
# Supplementary Figure 8

go %>% 
  ggplot(aes(x = reorder(term.name, -log10(p.value)), y = -log10(p.value))) + 
  geom_bar(stat = "identity", fill = "royalblue") + coord_flip() + xlab("") + theme_bw() + 
  geom_hline(yintercept = -log10(0.05), lty = 2, color = "red")

dev.off()



# 04. Expression-weighted cell-type enrichment (EWCE) ----

# Pre-calculated expression specificity metrics downloaded from www.hjerling-leffler-lab.org/data/scz_singlecell
load("./data/ctd_allKI.rda")

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

jpeg("./results/cell-type-neg.jpeg", units = "in", width = 20, height = 6, res = 300)
# Supplementary Figure 10

ewce.plot(full_results$results, mtc_method = "BH")

dev.off()



# 05. Visualize seeded networks wrt C4A copy number ----

# Load C4 imputation result
imputed <- read.table("./data/PsychENCODE/PsychENCODE-C4-imputed.txt", header = TRUE)
imputed$C4_CN <- imputed$C4A_CN + imputed$C4B_CN
imputed <- imputed %>% dplyr::select(-c(diagnosis, sex, study))  # datMeta already contains this information

# Load PsychENCODE expression + metadata
load("./data/PsychENCODE/PsychENCODE_regressed.RData")
datExpr <- datExpr.AllRegressed
rm(datExpr.AllRegressed)
sum(colnames(datExpr) == datMeta$X)  # 2160 total samples

datMeta <- datMeta %>% filter(tissue == "frontal cortex")  # 2026 frontal cortical samples
datExpr <- as.data.frame(datExpr)
datExpr <- datExpr %>% dplyr::select(datMeta$X)

datMeta <- datMeta %>% filter(individualID %in% imputed$sample)  # 916 samples with C4 imputation
datExpr <- datExpr %>% dplyr::select(datMeta$X)
datMeta <- cbind(datMeta, imputed[match(datMeta$individualID, imputed$sample),])
datExpr <- as.matrix(datExpr)

datMeta <- datMeta %>% mutate(pd_avg = (dose1 + dose2)/2) %>% filter(pd_avg >= 0.7) 
# 552 samples w/ high-confidence C4 imputation results
datExpr <- datExpr %>% as.data.frame() %>% dplyr::select(datMeta$X) %>% as.matrix()

datExpr_low <- datExpr[, datMeta$C4A_CN < 2]  # 119 samples
datExpr_mid <- datExpr[, datMeta$C4A_CN == 2]  # 324 samples
datExpr_hig <- datExpr[, datMeta$C4A_CN > 2]  # 109 samples

idx1 <- abs(network$`01.R`) > 0.5 & network$`01.FDR` < 0.05  # 4 genes for C4A 
idx2 <- abs(network$`2.R`) > 0.5 & network$`2.FDR` < 0.05  # 141 genes for C4A
idx3 <- abs(network$`34.R`) > 0.5 & network$`34.FDR` < 0.05  # 330 genes for C4A

idx <- idx1 | idx2 | idx3
sum(idx)  # 362 genes to plot

special <- c("MPPED2", "ZBTB18", "GRIA3", "MVP", "CTSS", "RELA",
             "EIF4EBP1", "NAGA", "SIPA1")

# Set up coordinates that the seeded networks will share
ADJ1 <- adjacency(t(datExpr_low[idx, ]), type = "unsigned", power = 1)
colnames(ADJ1) = rownames(ADJ1) = network$Gene[idx]
ADJ1[ADJ1 < 0.5] <- 0
graph.gene1 <- graph.adjacency(as.matrix(ADJ1), mode = "undirected", 
                               weighted = TRUE, diag = FALSE)
coords <- data.frame(layout.fruchterman.reingold(graph.gene1), 
                     gene = colnames(ADJ1))

# Plot CN < 2 seeded network
ADJ1 <- adjacency(t(datExpr_low[idx1, ]), type = "unsigned", power = 1)
ADJ1.neg <- cor(t(datExpr_low[idx1, ]))
colnames(ADJ1) = rownames(ADJ1) = network$Gene[idx1]
colnames(ADJ1.neg) = rownames(ADJ1.neg) = network$Gene[idx1]
ADJ1[ADJ1 < 0.5] <- 0
ADJ1.neg[abs(ADJ1.neg) < 0.5] <- 0
colnames(ADJ1)[order(colSums(ADJ1 > 0.5))]
colnames(ADJ1)[order(colSums(ADJ1))]

graph.gene1 <- graph.adjacency(as.matrix(ADJ1), mode = "undirected", 
                               weighted = TRUE, diag = FALSE)
edgelist1 <- as_edgelist(graph.gene1)
edgelist.temp <- as_edgelist(graph.gene1)
edgelist1[, 1] <- match(edgelist1[, 1], coords$gene)
edgelist1[, 2] <- match(edgelist1[, 2], coords$gene)
edges1 <- data.frame(coords[edgelist1[, 1], 1:2], 
                     coords[edgelist1[, 2], 1:2], 
                     Network = "Gene", 
                     PCC = NA)
for (i in 1:nrow(edges1)) {
  edges1$PCC[i] <- ADJ1[edgelist.temp[i, 1], edgelist.temp[i, 2]]
  edges1$PCC.neg[i] <- ADJ1.neg[edgelist.temp[i, 1], edgelist.temp[i, 2]]
}
edges1$PCC.neg[edges1$PCC.neg > 0] <- 0
edges1$PCC[edges1$PCC.neg < 0] <- 0

jpeg("./results/network01.jpeg", units = "in", width = 5.5, height = 4, res = 300)
ggplot() + 
  geom_segment(aes(x = X1, y = X2, xend = X1.1, yend = X2.1), color = 'red', 
               data = edges1, size = (edges1$PCC)^5, alpha = 0.8) +
  geom_segment(aes(x = X1, y = X2, xend = X1.1, yend = X2.1), color = 'blue', 
               data = edges1, size = abs(edges1$PCC.neg)^5, alpha = 0.8) +
  geom_point(aes(X1, X2), data = coords, alpha = 0, color = 'black', 
             fill = '#EEEEEE', shape = 21, size = 1.5) +
  geom_point(aes(X1, X2), data = coords[match(colnames(ADJ1), coords$gene), ], 
             alpha = 0.9, color = 'black', fill = '#EEEEEE', shape = 21, size = 1.2) +
  geom_point(aes(X1, X2), data = coords[coords$gene == "C4A", ], alpha = 1, 
             color = 'black', size = 3.5) + 
  theme_classic() + labs(x = "", y = "") + 
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(1,1,1,1), "cm"), 
        rect = element_rect(fill = "transparent", colour = NA)) +
  geom_text(aes(x = X1, y = X2 - 2.0, label = gene), 
            data = coords[coords$gene == "C4A", ], size = 5.0) + 
  geom_text_repel(aes(x = X1, y = X2, label = gene), 
                  data = coords[match(c("C4B", "GFAP", "APLNR"), coords$gene), ], 
                  size = 3.0, color = "#444444")
dev.off()

# Plot CN = 2 seeded network
ADJ2 <- adjacency(t(datExpr_mid[idx2, ]), type = "unsigned", power = 1)
ADJ2.neg <- cor(t(datExpr_mid[idx2, ]))
colnames(ADJ2) = rownames(ADJ2) = network$Gene[idx2]
colnames(ADJ2.neg) = rownames(ADJ2.neg) = network$Gene[idx2]
ADJ2[ADJ2 < 0.5] <- 0
ADJ2.neg[abs(ADJ2.neg) < 0.5] <- 0
colnames(ADJ2)[order(colSums(ADJ2 > 0.5))]
colnames(ADJ2)[order(colSums(ADJ2))]

graph.gene2 <- graph.adjacency(as.matrix(ADJ2), mode = "undirected", 
                               weighted = TRUE, diag = FALSE)
edgelist2 <- as_edgelist(graph.gene2)
edgelist.temp <- as_edgelist(graph.gene2)
edgelist2[, 1] <- match(edgelist2[, 1], coords$gene)
edgelist2[, 2] <- match(edgelist2[, 2], coords$gene)
edges2 <- data.frame(coords[edgelist2[, 1], 1:2], 
                     coords[edgelist2[, 2], 1:2], 
                     Network = "Gene", 
                     PCC = NA)
for (i in 1:nrow(edges2)) {
  edges2$PCC[i] <- ADJ2[edgelist.temp[i, 1], edgelist.temp[i, 2]]
  edges2$PCC.neg[i] <- ADJ2.neg[edgelist.temp[i, 1], edgelist.temp[i, 2]]
}
edges2$PCC.neg[edges2$PCC.neg > 0] <- 0
edges2$PCC[edges2$PCC.neg < 0] <- 0

astro1 <- c("ALDH1L1", "PLSCR4", "AQP4")
astro2 <- c("AHCYL1", "ITGB4")
emphasize <- c("LIMK2", "MVP")
micro <- c("C1QC", "C1QB", "ITGB2")

jpeg("./results/network2.jpeg", units = "in", width = 5.5, height = 4, res = 300)
ggplot() + 
  geom_segment(aes(x = X1, y = X2, xend = X1.1, yend = X2.1), color = 'red', 
               data = edges2, size = (edges2$PCC)^5, alpha = 0.8) +
  geom_segment(aes(x = X1, y = X2, xend = X1.1, yend = X2.1), color = 'blue', 
               data = edges2, size = abs(edges2$PCC.neg)^5, alpha = 0.8) +
  geom_point(aes(X1, X2), data = coords, alpha = 0, color = 'black', 
             fill = '#EEEEEE', shape = 21, size = 1.5) +
  geom_point(aes(X1, X2), data = coords[match(colnames(ADJ2), coords$gene), ], 
             alpha = 0.9, color = 'black', fill = '#EEEEEE', shape = 21, size = 1.2) +
  geom_point(aes(X1, X2), data = coords[coords$gene == "C4A", ], alpha = 1, 
             color = 'black', size = 3.5) + 
  theme_classic() + labs(x = "", y = "") + 
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(1,1,1,1), "cm"), 
        rect = element_rect(fill = "transparent", colour = NA)) +
  geom_text(aes(x = X1, y = X2 - 2.0, label = gene), 
            data = coords[coords$gene == "C4A",], size = 5.0) + 
  geom_text_repel(aes(x = X1, y = X2, label = gene), 
                  data = coords[match(c("C4B", "GFAP", "APLNR"), coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = -2, nudge_x = -0.25, 
                  segment.size = 0.2, segment.color = "#999999") +
  geom_text_repel(aes(x = X1, y = X2, label = gene), 
                  data = coords[match(astro1, coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 5, nudge_x = -3.5, 
                  segment.size = 0.2, segment.color = "#999999") +
  geom_text_repel(aes(x = X1, y = X2, label = gene), 
                  data = coords[match(astro2, coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 0, nudge_x = -5, 
                  segment.size = 0.2, segment.color = "#999999") +
  geom_text_repel(aes(x = X1, y = X2, label = gene), 
                  data = coords[match(micro, coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 1.75, nudge_x = 4.0, 
                  segment.size = 0.2, segment.color = "#999999") +
  geom_text_repel(aes(x = X1, y = X2, label = gene, fontface = "bold"), 
                  data = coords[match("MVP", coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 5.0, nudge_x = -3.0, 
                  segment.size = 0.2, segment.color = "#999999") +
  geom_text_repel(aes(x = X1, y = X2, label = gene), 
                  data = coords[match("LIMK2", coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 7.0, nudge_x = -1, 
                  segment.size = 0.2, segment.color = "#999999")
dev.off()

# Plot CN > 2 seeded network
ADJ3 <- adjacency(t(datExpr_hig[idx3, ]), type = "unsigned", power = 1)
ADJ3.neg <- cor(t(datExpr_hig[idx3, ]))
colnames(ADJ3) = rownames(ADJ3) = network$Gene[idx3]
colnames(ADJ3.neg) = rownames(ADJ3.neg) = network$Gene[idx3]
ADJ3[ADJ3 < 0.5] <- 0
ADJ3.neg[abs(ADJ3.neg) < 0.5] <- 0
colnames(ADJ3)[order(colSums(ADJ3 > 0.5))]
colnames(ADJ3)[order(colSums(ADJ3))]

graph.gene3 <- graph.adjacency(as.matrix(ADJ3), mode = "undirected", 
                               weighted = TRUE, diag = FALSE)
edgelist3 <- as_edgelist(graph.gene3)
edgelist.temp <- as_edgelist(graph.gene3)
edgelist3[, 1] <- match(edgelist3[, 1], coords$gene)
edgelist3[, 2] <- match(edgelist3[, 2], coords$gene)
edges3 <- data.frame(coords[edgelist3[, 1], 1:2], 
                     coords[edgelist3[, 2], 1:2], 
                     Network = "Gene", 
                     PCC = NA)
for (i in 1:nrow(edges3)) {
  edges3$PCC[i] <- ADJ3[edgelist.temp[i, 1], edgelist.temp[i, 2]]
  edges3$PCC.neg[i] <- ADJ3.neg[edgelist.temp[i, 1], edgelist.temp[i, 2]]
}
edges3$PCC.neg[edges3$PCC.neg > 0] <- 0
edges3$PCC[edges3$PCC.neg < 0] <- 0

neg <- c("SLC39A10", "KCNK1", "RGS4", "DCAF6", "PLCB4")
neg.point <- network$Gene[network$`34.R` < -0.5 & network$`34.FDR` < 0.05]
astro1 <- c("ALDH1L1", "PLSCR4", "AQP4")
astro2 <- c("AHCYL1", "ITGB4", "CD99")
emphasize <- c("MVP", "RELA")
micro <- c("C1QC", "C1QB", "ITGB2", "FCGR3A")

jpeg("./results/network34.jpeg", units = "in", width = 5.5, height = 4, res = 300)
ggplot() + 
  geom_segment(aes(x = X1, y = X2, xend = X1.1, yend = X2.1), color = 'red', 
               data = edges3, size = (edges3$PCC)^5, alpha = 0.8) +
  geom_segment(aes(x = X1, y = X2, xend = X1.1, yend = X2.1), color = 'blue', 
               data = edges3, size = abs(edges3$PCC.neg)^5, alpha = 0.8) +
  geom_point(aes(X1, X2), data = coords, alpha = 0, color = 'black', 
             fill = '#EEEEEE', shape = 21, size = 1.5) +
  geom_point(aes(X1, X2), data = coords[match(colnames(ADJ3), coords$gene), ], 
             alpha = 0.9, color = 'black', fill = '#EEEEEE', shape = 21, size = 1.2) +
  geom_point(aes(X1, X2), data = coords[coords$gene == "C4A", ], alpha = 1, 
             color = 'black', size = 3.5) + 
  geom_point(aes(X1, X2), data=coords[match(neg.point, coords$gene), ], 
             alpha = 0.7, color = 'black', fill = 'royalblue1', shape = 21, size = 1.2) + 
  theme_classic() + labs(x = "", y = "") + 
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(1,1,1,1), "cm"), 
        rect = element_rect(fill = "transparent", colour = NA)) +
  geom_text(aes(x = X1, y = X2 - 2.0, label = gene), 
            data = coords[coords$gene == "C4A",], size = 5.0) + 
  geom_text_repel(aes(x = X1, y = X2, label = gene), 
                  data = coords[match(c("C4B", "GFAP"), coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = -2.5, nudge_x = 0, 
                  segment.size = 0.2, segment.color = "#999999") + 
  geom_text_repel(aes(x = X1, y = X2, label = gene), 
                  data = coords[match(neg, coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 3.5, nudge_x = -7.0, 
                  segment.size = 0.2, segment.color = "#999999") + 
  geom_text_repel(aes(x = X1, y = X2, label = gene, fontface = "bold"), 
                  data = coords[match("GRIA3", coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = -1.5, nudge_x = -7.0, 
                  segment.size = 0.2, segment.color = "#999999") + 
  geom_text_repel(aes(x = X1, y = X2, label = gene, fontface = "bold"), 
                  data = coords[match("ZBTB18", coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 3.0, nudge_x = -0.5, 
                  segment.size = 0.2, segment.color = "#999999") + 
  geom_text_repel(aes(x = X1, y = X2, label = gene, fontface = "bold"), 
                  data = coords[match("MPPED2", coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 1.5, nudge_x = -2.0, 
                  segment.size = 0.2, segment.color = "#999999") + 
  geom_text_repel(aes(x = X1, y = X2, label = gene), 
                  data = coords[match(micro, coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 4.0, nudge_x = 5.5, 
                  segment.size = 0.2, segment.color = "#999999") + 
  geom_text_repel(aes(x = X1, y = X2, label = gene, fontface = "bold"), 
                  data = coords[match(emphasize, coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 7.5, nudge_x = -0.75, 
                  segment.size = 0.2, segment.color = "#999999") + 
  geom_text_repel(aes(x = X1, y = X2, label = gene), 
                  data = coords[match("LIMK2", coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 7.5, nudge_x = 0.25, 
                  segment.size = 0.2, segment.color = "#999999") + 
  geom_text_repel(aes(x = X1, y = X2, label = gene, fontface = "bold"), 
                  data = coords[match("CTSS", coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 11.5, nudge_x = 1.5, 
                  segment.size = 0.2, segment.color = "#999999") + 
  geom_text_repel(aes(x = X1, y = X2, label = gene, fontface = "bold"), 
                  data = coords[match("EIF4EBP1", coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 11.5, nudge_x = -6.5, 
                  segment.size = 0.2, segment.color = "#999999") + 
  geom_text_repel(aes(x = X1, y = X2, label = gene, fontface = "bold"), 
                  data = coords[match("NAGA", coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = -1.5, nudge_x = 1.0, 
                  segment.size = 0.2, segment.color = "#999999") + 
  geom_text_repel(aes(x = X1, y = X2, label = gene, fontface = "bold"), 
                  data = coords[match("SIPA1", coords$gene), ], 
                  size = 3.0, color = "#444444", nudge_y = 11.5, nudge_x = -4.0, 
                  segment.size = 0.2, segment.color = "#999999")
dev.off()