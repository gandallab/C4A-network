# 02b_network-permutation
library(parallel)
load("data/data_for_network_permutation.RData")

## Over-representation analysis functions
## Odds-ratio estimator
OR <- function(q,k,m,t) {
  ## 2 x 2 table:
  ##         inTest   !inTest
  ## inRef     q        k
  ## !inRef    m        t
  
  fisher.out <- fisher.test(matrix(c(q, k-q, m-q, t-m-k+q), 2, 2),conf.int=TRUE)
  OR <- fisher.out$estimate
  pval <- fisher.out$p.value
  upCI <- fisher.out$conf.int[1]
  downCI <- fisher.out$conf.int[2]
  
  output <- c(OR,pval,upCI,downCI)
  names(output) <- c("OR","Fisher p","-95%CI","+95%CI")
  return(output)
}

## count overlaps and run the analysis
ORA <- function(testpath,refpath,testbackground,refbackground) {
  testpath = testpath[testpath %in% testbackground]
  refpath = refpath[refpath %in% refbackground]
  q <- length(intersect(testpath,refpath)) ## overlapped pathway size
  k <- length(intersect(refpath,testbackground))  ## input gene set
  m <- length(intersect(testpath,refbackground)) ## input module
  t <- length(intersect(testbackground,refbackground)) ## Total assessed background (intersect reference and test backgrounds)
  
  empvals <- OR(q,k,m,t)
  
  tmpnames <- names(empvals)
  empvals <- as.character(c(empvals,q,k,m,t,100*signif(q/k,3)))
  names(empvals) <- c(tmpnames,"Overlap","Reference List","Input List","Background","% List Overlap")
  return(empvals)
}


prsCor = function(i, gene, datExpr) {
  c = cor.test(datExpr[i, ], datExpr[gene, ], use = "pairwise.complete.obs")
  dfPrs = data.frame(Seed = gene, Gene = rownames(datExpr)[i], R = c$estimate, P = c$p.value)
  return(dfPrs)
}

runOneBootstrap = function(gene) {
  this_net= do.call("rbind", mclapply(1:nrow(datExpr_mid.ctl), prsCor, gene, datExpr_mid.ctl,mc.cores = 16,mc.allow.recursive = F))
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

for(i in seq(4501,25774,by=250)) {
  print(i)
  out_file = paste0("boot_", i, ".csv")
  if(!file.exists(out_file)) {
    df_boot = do.call("rbind", lapply(as.list(rownames(datExpr_mid.ctl)[i:(i+99)]), runOneBootstrap));   
    write.csv(df_boot, file=out_file)
  }
}
  

