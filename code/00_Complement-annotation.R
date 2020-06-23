# Compile Complement Gene Sets

library(data.table)
library(STRINGdb)
library(biomaRt)
library(Matrix)

# Load gene annotations
gencode = read.csv("data/annotation.gene.gencodeV19.csv")

# Download complement system genes from HGNC
complement_system = as.data.frame(fread('https://www.genenames.org/cgi-bin/genegroup/download?id=492&type=branch'))
complement_system[57, "Approved symbol"] = 'C5AR2'
complement_system$`Ensembl gene ID`[57] = gencode$gene_id[gencode$gene_name=="C5AR2"]

## String DB PPIs
C4A.direct_stringDB = c('ADAM10','AFP','AHSG','ALB','ANO8','APLP2','APOA1','APOA2','APOA5','APOB','APOE','APOL1','APP','BPIFB2','CALU','CHGB','CSF1','CST3','CYR61','DMP1','EVA1A','F5','FAM20A','FAM20C','FBN1','FGA','FGG','FN1','FSTL3','FUCA2','GAS6','HSP90B1','IGFBP1','IGFBP4','IGFBP5','IGFBP7','IGFPB3','IGLL5','IL6','ITIH2','KNG1','LAMB1','LAMB2','LAMC1','LGALS1','LTBP1','MATN3','MFGE8','MFI2','MIA3','MSLN','MXRA8','NOTUM','NUCB1','P4HB','PCSK9','PDIA6','PNPLA2','PRKCSH','PROC','PRSS23','SDC2','SERPINA1','SERPINA10','SERPINC1','SERPIND1','SPARCL1','SPP1','SPP2','TF','TIMP1','TNC','VCAN','VGF','VWA1')

# Get PPIs
# https://inbio-discover.intomics.com/map.html#downloads
ppi = fread('data/InBio_Map_core_2016_09_12/core.psimitab',head=F)
ppi$Uniprot1=gsub("uniprotkb:", "", ppi$V1)
ppi$Uniprot2=gsub("uniprotkb:", "", ppi$V2)
ppi$Confidence = as.numeric(unlist(lapply(strsplit(ppi$V15, "[|]"), '[',1)))
min(ppi$Confidence)

bm = useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = 'grch37.ensembl.org'); #a = listAttributes(bm)
annot = getBM(attributes = c("uniprotsptrembl", "uniprotswissprot", "ensembl_gene_id", "external_gene_name"), mart=bm)

idx = match(ppi$Uniprot1, annot$uniprotswissprot)
ppi$ENSG1 = annot$ensembl_gene_id[idx]
ppi$HGNC1 = annot$external_gene_name[idx]
idx = match(ppi$Uniprot2, annot$uniprotswissprot)
ppi$ENSG2 = annot$ensembl_gene_id[idx]
ppi$HGNC2 = annot$external_gene_name[idx]
genes.symbol = sort(unique(c(ppi$HGNC1, ppi$HGNC2)))
genes.ensg = sort(unique(c(ppi$ENSG1, ppi$ENSG2)))
idx.symbol = data.frame(idx1=match(ppi$HGNC1, genes.symbol), idx2 = match(ppi$HGNC2, genes.symbol), val = ppi$Confidence)
idx.symbol = idx.symbol[!apply(is.na(idx.symbol), 1, any),]
non_symmetric = which(idx.symbol$idx2 > idx.symbol$idx1)
temp = idx.symbol$idx1[non_symmetric]
idx.symbol$idx1[non_symmetric] = idx.symbol$idx2[non_symmetric]
idx.symbol$idx2[non_symmetric] = temp

ppiDirect.symbol = sparseMatrix(i=idx.symbol$idx1, j=idx.symbol$idx2, x=idx.symbol$val, symmetric = T)
colnames(ppiDirect.symbol) = rownames(ppiDirect.symbol) = genes.symbol


# Get C4A and complement direct PPIs
C4A.direct = names(which(ppiDirect.symbol[,"C4A"] != 0 ))
C4A.direct_highConfidence = names(which(ppiDirect.symbol[,"C4A"] >= 0.7 ))

idx = match(complement_system$`Approved symbol`, rownames(ppiDirect.symbol))
Complement.direct = names(which(apply(ppiDirect.symbol[na.omit(idx),]>0,2,any)))
Complement.direct_highConfidence = names(which(apply(ppiDirect.symbol[na.omit(idx),]>=0.7,2,any)))


csg1 = unique(complement_system$`Approved symbol`)
csg2 = unique(c(csg1, C4A.direct_stringDB))
csg3 = unique(c(csg1, C4A.direct_highConfidence))
csg4 = unique(c(csg1, C4A.direct))
csg5 = unique(c(csg1, Complement.direct_highConfidence))
csg6 = unique(c(csg1, Complement.direct))


geneSets = rbind(data.frame(gene = csg1, set='complement_57genes'),
                 data.frame(gene = csg2, set='complement_and_C4A_PPI_stringDB_132genes'),
                 data.frame(gene = csg3, set='complement_and_C4A_PPI_InWeb_highconf_67genes'),
                 data.frame(gene = csg4, set='complement_and_C4A_PPI_InWeb_86genes'),
                 data.frame(gene = csg5, set='complement_and_complement_PPI_InWeb_highconf_545genes'), 
                 data.frame(gene = csg6, set='complement_and_complement_PPI_InWeb_1066genes'))

geneSets.ensg = geneSets
geneSets.ensg$gene = annot$ensembl_gene_id[match(geneSets.ensg$gene, annot$external_gene_name)]
geneSets.ensg = geneSets.ensg[!is.na(geneSets.ensg$gene),]

write.table(file="results/complement_gene_sets.tsv", geneSets, row.names = F, quote=F, sep='\t')
write.table(file="results/complement_gene_sets_ensg.tsv", geneSets.ensg, row.names = F, quote=F, sep='\t')


library(writexl)
library(reshape)
tableS1 = list(Complement_system = complement_system,
               C4A_PPI_stringDB = data.frame(Gene = C4A.direct_stringDB),
               C4A_PPI_InWeb = data.frame(Gene1 =  "C4A", Gene2 = C4A.direct, Confidence=ppiDirect.symbol[C4A.direct,"C4A"]))

Complement_PPI_Inweb = data.frame(melt(as.matrix(ppiDirect.symbol[complement_system$`Approved symbol`[complement_system$`Approved symbol` %in% rownames(ppiDirect.symbol)], Complement.direct])), stringsAsFactors = F)
colnames(Complement_PPI_Inweb) = c("Complement_gene", "PPI_target", "Confidence") 
tableS1[[4]] = Complement_PPI_Inweb
names(tableS1)[4] = "Complement_PPI_Inweb"
write_xlsx(tableS1, path="results/manuscript/TableS1.xlsx")
