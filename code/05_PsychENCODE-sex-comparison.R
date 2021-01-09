rm(list = ls())
options(stringsAsFactors = FALSE)

library(tidyverse)
library(clusterProfiler)
library(ggrepel)
library(gridExtra)
library(grid)



# 01. Load data and annotations ----

# Load sex-specific C4A coexpression data from PsychENCODE
load("./data/PsychENCODE/C4A-network-PsychENCODE-by-sex.RData")

# Load MSigDB gene sets -- GO and Hallmark
geneSets <- data.frame()
geneSets <- rbind(geneSets, read.gmt("./data/GSEA-annotations/c5.all.v7.1.symbols.gmt"))
geneSets <- rbind(geneSets, read.gmt("./data/GSEA-annotations/h.all.v7.1.symbols.gmt"))

# Load SynGO gene sets
syngo <- as.data.frame(readxl::read_xlsx("./data/GSEA-annotations/syngo_annotations.xlsx")[, c(7,4)])
colnames(syngo) <- colnames(geneSets)
geneSets <- rbind(geneSets, syngo)

# Load Lake et al. cell type clusters
lake <- read.csv("./data/GSEA-annotations/Lake_scClusters_fromSupplement.csv")[, c(2,1)]
colnames(lake) <- colnames(geneSets)
lake$term <- paste0("Lake_", substr(lake$term, 0, 3))
geneSets <- rbind(geneSets, lake)

# Combine Male and Female dataframes
df = do.call("cbind", list(M = df.m, F = df.f))



# 02. Run GSEA on M, F, and M vs F separately ----

# Run GSEA on Female results -- rank genes by Pearson R
geneList.F <- df$F.R
names(geneList.F) <- df$F.Gene
geneList.F <- geneList.F[order(geneList.F, decreasing = T)]
go.F <- GSEA(geneList.F, TERM2GENE = geneSets, pvalueCutoff = 1, nPerm = 10000)
#View(as.data.frame(go.F))
gseaplot(go.F, "GO_COMPLEMENT_ACTIVATION")  # Sanity check

# Run GSEA on Male results -- rank genes by Pearson R
geneList.M = df$M.R
names(geneList.M) <- df$M.Gene
geneList.M <- geneList.M[order(geneList.M, decreasing = T)]
go.M <- GSEA(geneList.M, TERM2GENE = geneSets, pvalueCutoff = 1, nPerm = 10000,)
#View(as.data.frame(go.M))
gseaplot(go.M,"GO_COMPLEMENT_ACTIVATION")  # Sanity check

# Run GSEA on Male vs Female results -- rank genes by Pearson R difference
geneList.MvsF <- df$M.R - df$F.R
names(geneList.MvsF) <- df$M.Gene
geneList.MvsF <- geneList.MvsF[order(geneList.MvsF, decreasing = T)]
go.MvsF <- GSEA(geneList.MvsF, TERM2GENE = geneSets, pvalueCutoff = 1, nPerm = 10000)
#View(as.data.frame(go.MvsF))
gseaplot(go.MvsF,"GO_GABA_ERGIC_SYNAPSE")



# 03. Compare result from M and F GSEA ----

df.combined <- merge(go.M@result, go.F@result, by.x = "ID", by.y = "ID", all = T)
df.combined$NES_difference <- df.combined$NES.x - df.combined$NES.y # NES = normalized enrichment statistic
df.combined$NES_mean <- .5*(df.combined$NES.x + df.combined$NES.y)

# Keep terms that are (nominally) significant in at least M or F
df.combined <- df.combined[!is.na(df.combined$NES_mean) & 
                            (df.combined$pvalue.x < .05 | df.combined$pvalue.y < .05), ]
#View(df.combined[,c('ID', 'NES_difference', 'NES_mean', 'NES.x', 'NES.y', 'pvalue.x', 'pvalue.y')])

df.combined$Y_significant <- F
df.combined$Y_significant[df.combined$qvalues.y < .1] <- T
df.combined$logPdiff <- with(df.combined, -log10(pvalue.x)*sign(NES.x) - -log10(pvalue.y)*sign(NES.y))



# 04. Plot top 25 concordant and discordant enrichmenbt terms ----

to_plot1 <- with(df.combined, c(order(NES.x,decreasing = T)[1:25], 
                                order(NES.x,decreasing = F)[1:25], 
                                order(NES.y,decreasing = T)[1:25],
                                order(NES.y,decreasing = T)[1:25]))
to_plot2 <- with(df.combined, c(order(NES_difference,decreasing = T)[1:25], 
                                order(NES_difference,decreasing = F)[1:25]))

df.combined$Ontology <- "HALLMARK & SynGO"
df.combined$Ontology[grep("BP_", df.combined$ID)] <- "GO:BP"
df.combined$Ontology[grep("CC_", df.combined$ID)] <- "GO:CC+MF"
df.combined$Ontology[grep("MF_", df.combined$ID)] <- "GO:CC+MF"
df.combined$Ontology[grep("Lake_", df.combined$ID)] <- "Lake_scSeq"
df.combined$Ontology[grep("HALLMARK_", df.combined$ID)] <- "HALLMARK & SynGO"

df.combined$ID <- df.combined$ID %>% 
  str_replace("CC_", "") %>% 
  str_replace("BP_", "")  %>% 
  str_replace("MF_", "") %>% 
  str_replace("HALLMARK_", "") 

df.combined$qvalues.x <- p.adjust(df.combined$p.adjust.x, method = "bonferroni")

g1 <- ggplot(df.combined[sign(df.combined$NES.x) == sign(df.combined$NES.y) & 
                           (df.combined$p.adjust.x < .05 & df.combined$p.adjust.y<.05), ], 
             aes(x = -log10(pvalue.x) * sign(NES.x), y = -log10(pvalue.y) * sign(NES.y), 
                 label = tolower(ID), color = Ontology)) +
  geom_point(aes(size = setSize.x)) + 
  geom_abline(slope = 1, lty = 2, size = .2) + 
  geom_hline(yintercept = 0, lty = 2, size = .2) + 
  geom_vline(xintercept = 0, lty = 2, size = .2) + 
  theme_bw() + ylim(-5, 5) + xlim(-5, 5) + 
  geom_text_repel(size = 3, segment.size = .1, force = 2) + 
  facet_wrap(~ Ontology) + 
  ggtitle("Sex concordant pathways")

g1.big <- ggplot(df.combined[sign(df.combined$NES.x) == sign(df.combined$NES.y) & 
                               (df.combined$p.adjust.x < .05 & df.combined$p.adjust.y<.05), ], 
                 aes(x = -log10(pvalue.x) * sign(NES.x), y = -log10(pvalue.y) * sign(NES.y), 
                     label = tolower(ID), color = Ontology)) +   
  geom_point(aes(size = setSize.x)) + 
  geom_abline(slope = 1, lty = 2,size = .2) + 
  geom_hline(yintercept = 0, lty = 2, size = .2) + 
  geom_vline(xintercept = 0, lty = 2, size = .2) + 
  theme_bw() + ylim(-5, 5) + xlim(-5, 5) + 
  geom_text_repel(size = 3,segment.size = .1,force = 2) + 
  ggtitle("Sex concordant pathways")

g2 <- ggplot(df.combined[sign(df.combined$NES.x) == -sign(df.combined$NES.y) & 
                           abs(df.combined$NES_difference) > 2 & 
                           (df.combined$p.adjust.y < .05 | df.combined$p.adjust.x<.05), ], 
             aes(x = -log10(pvalue.x) * sign(NES.x), y = -log10(pvalue.y) * sign(NES.y),  
                 label = tolower(ID), color = Ontology)) +    
  geom_point(aes(size = setSize.x)) + 
  geom_abline(slope = 1, lty = 2, size = .2) + 
  geom_hline(yintercept = 0, lty = 2, size = .2) + 
  geom_vline(xintercept = 0, lty=2, size = .2) + 
  theme_bw() + ylim(-5, 5) + xlim(-5, 5) + 
  geom_text_repel(size = 3, segment.size = .1, force = 10) + 
  facet_wrap(~ Ontology) + ggtitle("Sex specific pathways")

g2.big <- ggplot(df.combined[sign(df.combined$NES.x) == -sign(df.combined$NES.y) & 
                               abs(df.combined$NES_difference) > 2 & 
                               (df.combined$p.adjust.y < .05 | df.combined$p.adjust.x<.05), ], 
                 aes(x = -log10(pvalue.x) * sign(NES.x), y = -log10(pvalue.y) * sign(NES.y),  
                     label = tolower(ID), color = Ontology)) + 
  geom_point(aes(size = setSize.x))  +
  geom_abline(slope = 1, lty = 2, size = .2) + 
  geom_hline(yintercept = 0, lty = 2, size = .2) + 
  geom_vline(xintercept = 0, lty = 2, size = .2) + 
  theme_bw() + ylim(-5, 5) + xlim(-5, 5) + 
  geom_text_repel(size = 3,segment.size = .1, force = 10) +  
  ggtitle("Sex specific pathways")

# ggsave(g1,file = './results/C4A_sexCondordant.pdf', device = cairo_pdf, width = 10, height=10, units = 'in')
# ggsave(g2,file = './results/C4A_sexSpecific.pdf', device = cairo_pdf, width = 14, height=14, units = 'in')
# ggsave(g1.big,file = './results/C4A_sexCondordant_big.pdf', device = cairo_pdf, width = 14, height=14, units = 'in')
# ggsave(g2.big,file = './results/C4A_sexSpecific_big.pdf', device = cairo_pdf, width = 14 , height=14, units = 'in')



# 05. Plot Individual GO Terms for M and F ----

term <- "GO_AXONEME_ASSEMBLY"
gseaplot(go.M, geneSetID = term, title = "Male", ylim = c(0, 1))
gseaplot(go.F, geneSetID = term, title = "Female", ylim = c(0, 1))

g1 <- gseaplot(go.F, geneSetID = term, title = "Female", ylim = c(0, 1), by = "runningScore", 
               color = "pink", color.line = "pink", color.vline = "pink") + ylim(-1, 1)
g2 <- gseaplot(go.M, geneSetID = term, title = "Male", ylim = c(0, 1), by = "runningScore",
               color = "blue", color.line = "blue", color.vline = "blue") + ylim(-1, 1)
grid.arrange(grobs = list(g1, g2), ncol = 1, top = textGrob(term, gp = gpar(fontsize = 20, font = 3)))

term <- "HALLMARK_MTORC1_SIGNALING"
g1 <- gseaplot(go.F, geneSetID = term, title = "Female", ylim = c(0, 1), by = "runningScore",
               color = "#FFC0CB", color.line = "#FFC0CB", color.vline = "#FFC0CB") + ylim(-0.5, .5)
g2 <- gseaplot(go.M, geneSetID = term, title = "Male", ylim = c(0, 1), by = "runningScore", 
               color = "#8EB2E6", color.line = "#8EB2E6", color.vline = "#8EB2E6") + ylim(-0.5, .5)
grid.arrange(grobs = list(g1,g2), ncol = 1, top = textGrob(term, gp = gpar(fontsize = 20, font = 3)))



# 06. Barplot of top concordant and discoard GO term enrichments ----

to_plot <- c("GO_COMPLEMENT_ACTIVATION", "GO_HUMORAL_IMMUNE_RESPONSE", "Lake_Ast", "Lake_Mic", 
             "GO_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY", "GO_TRANSFORMING_GROWTH_FACTOR_BETA_BINDING",
             "Lake_Ex5", "GO_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT", "GO_RESPIRASOME", "Lake_Ex4", 
             "GO_PROTEASOME_ACCESSORY_COMPLEX", "GO_NADH_DEHYDROGENASE_COMPLEX_ASSEMBLY",
             "Lake_Oli", "INFLAMMATORY_RESPONSE", "TNFA_SIGNALING_VIA_NFKB", "CHOLESTEROL_HOMEOSTASIS", 
             "MTORC1_SIGNALING", "GO_OUTER_DYNEIN_ARM_ASSEMBLY", "Lake_Ex1", "Lake_Ex2", "Lake_Ex6",  
             "Lake_In8", "GO_INTRACILIARY_TRANSPORT_INVOLVED_IN_CILIUM_ASSEMBLY", 
             "GO_INTRINSIC_COMPONENT_OF_POSTSYNAPTIC_DENSITY_MEMBRANE")

idx <- match(to_plot, df.combined$ID)

this_df <- rbind(data.frame(ID = df.combined$ID[idx], NES = df.combined$NES.x[idx], NES.x = df.combined$NES.x[idx], 
                            P = df.combined$pvalue.x[idx], Sex = "M", NES.diff = df.combined$NES_difference[idx], 
                            NES.mean = df.combined$NES_mean[idx]),
                 data.frame(ID = df.combined$ID[idx], NES = df.combined$NES.y[idx], NES.x = df.combined$NES.x[idx], 
                            P = df.combined$pvalue.y[idx], Sex = "F", NES.diff = df.combined$NES_difference[idx], 
                            NES.mean = df.combined$NES_mean[idx]))
this_df$Group <- "Condordant"
this_df$Group[abs(this_df$NES.diff) > 3.2] <- "Discordant"

# Figure 6C
g1 <- ggplot(this_df, aes(x = reorder(tolower(ID), NES.x), y = NES, fill = Sex)) + 
  geom_bar(stat = "identity", position = position_dodge()) + coord_flip() + labs(x = "") +
  facet_wrap(Group ~ . , scales = "free_y") + theme_bw() + 
  scale_fill_manual(values = list("M" = "lightblue", "F" = "pink"))

ggsave(g1, file = "./results/PsychENCODE-MvsF.pdf", width = 10, height = 4)

# Look at cell type
idx = grep("Lake_Ex", df.combined$ID)
this_df = rbind(data.frame(ID = df.combined$ID[idx], NES = df.combined$NES.x[idx], NES.x = df.combined$NES.x[idx], 
                           P = df.combined$pvalue.x[idx], Sex = "M", NES.diff = df.combined$NES_difference[idx], 
                           NES.mean = df.combined$NES_mean[idx]),
                data.frame(ID = df.combined$ID[idx], NES = df.combined$NES.y[idx], NES.x = df.combined$NES.x[idx], 
                           P = df.combined$pvalue.y[idx], Sex = "F", NES.diff = df.combined$NES_difference[idx], 
                           NES.mean = df.combined$NES_mean[idx]))

ggplot(this_df, aes(x = ID, y = NES, fill = Sex)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  coord_flip() + theme_bw() 


