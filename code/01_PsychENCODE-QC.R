#01_PsychENCODE-QC.R

options(stringsAsFactors = FALSE)
library(tidyverse)
path <- "~/Github//C4A-network"


# 01. Load PsychENCODE data for primary analysis
# ----------------------------------
library(SummarizedExperiment)
load(paste0(path, '/data/se.CombinedDataForPrimaryAnalysis.RData'))

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



# 05. DGE in Complement pathway
# ----------------------------
imputed <- read.table(paste0(path, "/data/PsychENCODE-C4-imputed.txt"), header = TRUE)
imputed$C4_CN <- imputed$C4A_CN + imputed$C4B_CN
datMeta$X = rownames(datMeta)

datMeta_withCN <- as.tibble(datMeta) %>% filter(individualID %in% imputed$sample) # 916 samples with C4 imputation
datMeta_withCN <- cbind(datMeta_withCN, imputed[match(datMeta_withCN$individualID, imputed$sample),])
datExpr_withCN <- as.data.frame(datExpr)[,datMeta_withCN$X]


complement_genes = readxl::read_xlsx('results/manuscript/TableS1.xlsx',sheet=1)
table(complement_genes$`Ensembl gene ID` %in% rownames(datExpr))
these_complement_genes = complement_genes$`Approved symbol`[complement_genes$`Ensembl gene ID` %in% rownames(datExpr)]

tt = data.frame()
C4A_expr = datExpr[which(rowData(se.Primary)$gene_name=="C4A"),]
C4A_expr_subset = as.numeric(datExpr_withCN[which(rowData(se.Primary)$gene_name=="C4A"),])

library(nlme)

for(this_gene in these_complement_genes[-8]) {
  print(this_gene)
    expr = datExpr[which(rowData(se.Primary)$gene_name==this_gene),]

  s = summary(lme(expr~Group + age + age.squared + study + sex + PMI + RIN + RIN.squared + individualIDSource + tissue + seqPC3.squared +
                           seqPC1 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + seqPC14 + seqPC16 + seqPC18 + seqPC19 + seqPC2 + seqPC20 + seqPC21 + 
                                             seqPC22 + seqPC23 + seqPC24 + seqPC25 + seqPC27 + seqPC28 + seqPC29 + seqPC3 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + SV1 + SV2 + SV3 + SV4, data=datMeta, random=~1 | individualID))$tTable
  tt = rbind(tt, cbind(data.frame(Gene = this_gene, Group="SCZ_full", beta= s["GroupSCZ",1], se =s["GroupSCZ",2], tstat=s["GroupSCZ",4], p=s["GroupSCZ",5])))
  tt = rbind(tt, cbind(data.frame(Gene = this_gene, Group="BD_full", beta= s["GroupBD",1], se =s["GroupBD",2], tstat=s["GroupBD",4], p=s["GroupBD",5])))
  
  s2 = summary(lme(expr~Group + C4A_expr+age + age.squared + study + sex + PMI + RIN + RIN.squared + individualIDSource + tissue + seqPC3.squared +
                    seqPC1 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + seqPC14 + seqPC16 + seqPC18 + seqPC19 + seqPC2 + seqPC20 + seqPC21 + 
                    seqPC22 + seqPC23 + seqPC24 + seqPC25 + seqPC27 + seqPC28 + seqPC29 + seqPC3 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + SV1 + SV2 + SV3 + SV4, data=datMeta, random=~1 | individualID))$tTable
  
  tt = rbind(tt, cbind(data.frame(Gene = this_gene, Group="SCZ_full|C4expr", beta= s2["GroupSCZ",1], se =s2["GroupSCZ",2], tstat=s2["GroupSCZ",4], p=s2["GroupSCZ",5])))
  
  
  expr_subset = as.numeric(datExpr_withCN[which(rowData(se.Primary)$gene_name==this_gene),])
  s3 = summary(lme(expr_subset~Group + age + age.squared + study + sex + PMI + RIN + RIN.squared  + tissue + seqPC3.squared +
                    seqPC1 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + seqPC14 + seqPC16 + seqPC18 + seqPC19 + seqPC2 + seqPC20 + seqPC21 + 
                    seqPC22 + seqPC23 + seqPC24 + seqPC25 + seqPC27 + seqPC28 + seqPC29 + seqPC3 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + SV1 + SV2 + SV3 + SV4, data=datMeta_withCN, random=~1 | individualID))$tTable
  
  tt = rbind(tt, cbind(data.frame(Gene = this_gene, Group="SCZ_EUR", beta= s3["GroupSCZ",1], se =s3["GroupSCZ",2], tstat=s3["GroupSCZ",4], p=s3["GroupSCZ",5])))
  
  s4 = summary(lme(expr_subset~Group + C4A_expr_subset +  age + age.squared + study + sex + PMI + RIN + RIN.squared  + tissue + seqPC3.squared +
                     seqPC1 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + seqPC14 + seqPC16 + seqPC18 + seqPC19 + seqPC2 + seqPC20 + seqPC21 + 
                     seqPC22 + seqPC23 + seqPC24 + seqPC25 + seqPC27 + seqPC28 + seqPC29 + seqPC3 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + SV1 + SV2 + SV3 + SV4, data=datMeta_withCN, random=~1 | individualID))$tTable
  
  tt = rbind(tt, cbind(data.frame(Gene = this_gene, Group="SCZ_EUR|C4expr", beta= s4["GroupSCZ",1], se =s4["GroupSCZ",2], tstat=s4["GroupSCZ",4], p=s4["GroupSCZ",5])))
  
  s5 = summary(lme(expr_subset~Group + C4A_CN +  age + age.squared + study + sex + PMI + RIN + RIN.squared  + tissue + seqPC3.squared +
                     seqPC1 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + seqPC14 + seqPC16 + seqPC18 + seqPC19 + seqPC2 + seqPC20 + seqPC21 + 
                     seqPC22 + seqPC23 + seqPC24 + seqPC25 + seqPC27 + seqPC28 + seqPC29 + seqPC3 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + SV1 + SV2 + SV3 + SV4, data=datMeta_withCN, random=~1 | individualID))$tTable
  tt = rbind(tt, data.frame(Gene = this_gene, Group="SCZ_EUR|C4A_CN", beta= s5["GroupSCZ",1], se =s5["GroupSCZ",2], tstat=s5["GroupSCZ",4], p=s5["GroupSCZ",5]))
  
  s6 = summary(lme(expr_subset~Group + C4A_expr_subset + C4A_CN +  age + age.squared + study + sex + PMI + RIN + RIN.squared  + tissue + seqPC3.squared +
                     seqPC1 + seqPC10 + seqPC11 + seqPC12 + seqPC13 + seqPC14 + seqPC16 + seqPC18 + seqPC19 + seqPC2 + seqPC20 + seqPC21 + 
                     seqPC22 + seqPC23 + seqPC24 + seqPC25 + seqPC27 + seqPC28 + seqPC29 + seqPC3 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + SV1 + SV2 + SV3 + SV4, data=datMeta_withCN, random=~1 | individualID))$tTable
  tt = rbind(tt, data.frame(Gene = this_gene, Group="SCZ_EUR|C4A_CN+C4A_expr", beta= s6["GroupSCZ",1], se =s6["GroupSCZ",2], tstat=s6["GroupSCZ",4], p=s6["GroupSCZ",5]))
}



tt$FDR = p.adjust(tt$p,'fdr')
tt$beta_star = as.character(signif(tt$beta,1))
tt$beta_star[tt$FDR < .1] = paste0(tt$beta_star[tt$FDR < .1], " *")
tt$Gene = factor(tt$Gene, levels=sort(unique(tt$Gene),decreasing = T))
tt$Group = factor(tt$Group, levels=c("BD_full", "SCZ_full", "SCZ_EUR", "SCZ_full|C4expr", "SCZ_EUR|C4expr", "SCZ_EUR|C4A_CN", "SCZ_EUR|C4A_CN+C4A_expr"))

ggplot(tt,aes(x=Group, y=Gene, label=beta_star, fill=-log10(FDR))) + 
  geom_tile() + geom_text(size=3) + scale_fill_continuous(low='white', high='red')



