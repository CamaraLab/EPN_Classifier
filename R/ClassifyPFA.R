#'
#' Classify into PF A subgroups (PFA_1 and PFA_2)
#'
#' @param bulk a matrix of bulk transcriptomics data (genes by samples)
#' @param permutations the number permutations for a permutation test
#'
#' @export
#'
ClassifyPFA <- function(bulk, permutations = 100000){

  if (is.null(bulk)){
    message("bulk must be a matrix (genes by samples)")
    return(NULL)
  }
  if (is.null(row.names(bulk))){
    message("bulk must have row names (gene names)")
    return(NULL)
  }
  if (is.null(colnames(bulk))){
    colnames(bulk) <- 1:ncol(bulk)
  }

  if (class(permutations) != "numeric" | permutations <= 0 | permutations != round(permutations)){
    message("permutations must be a numeric whole value greater than zero")
    return(NULL)
  }

  #List of differentially expressed genes for each PF A subgroup
  gene_set <- NULL
  gene_set <- vector(mode = "list")
  gene_set[["PFA_1"]] <- c("HOTAIRM1", "HOXB3", "HOXB2", "HOXB4", "HOXA2", "HOXA3", "HOXA1", "HOXA-AS2", "MIR10A", "HOXB-AS1", "HOXB-AS2", "AC106786.1", "ZNF503", "LPAR3", "PRELP", "SEMA3C", "PRDM6", "ZNF503-AS2", "SKAP2", "SLC35F2", "HOXC4", "HOXD4", "PPP4R4", "AJAP1", "HOXA4", "FGFR2", "FGFRL1", "TICAM1", "LRRTM4", "GPR133", "CXXC5", "CCDC85A", "SEMA3B", "CNKSR3", "DPT", "ZNF703", "HCN1", "LRP2", "LINC00626", "HSPB1", "STAC2", "VLDLR", "AC007743.1", "COTL1", "RBMS3", "WIF1", "INPP4B", "RHOBTB3", "TNFRSF11A", "EGF", "LAMA2", "PAX8", "HOXD-AS1", "SERPINA12", "PLSCR4", "STXBP5L", "PIK3R1", "ARHGEF40")
  gene_set[["PFA_2"]] <- c("EN2", "MPPED2", "CNPY1", "PTPN3", "ST6GAL2", "AC008060.5", "B3GALTL", "AC005235.1", "LRRC8D", "SYT10", "WBSCR17", "PCBP3", "EFCC1", "LPO", "SMCO2", "NEDD4L", "EYA4", "LOXL4", "RNF43", "SERPINI1", "ITIH5", "GPR39", "PROB1", "TOX3", "ARNTL2", "BZRAP1","RFX3", "NTS", "ADAMTS3", "TCTEX1D1", "NR1H4", "KRT81", "TEX15", "ASB18", "SH3GL3", "FAM65C", "LINC00907", "LMX1B", "TRPC6", "LECT1", "AL592528.1", "SEMA3D", "TP73", "ELMOD1", "LGR4", "TMEM254", "CHL1", "GREB1L", "SPATA42", "TDRD1", "ITGA8", "KCNG2", "SVOPL", "NRSN1", "ERC2", "CATSPERD", "SPATA6", "CSPP1")

  #Gene set in bulk data
  gene_set <- lapply(gene_set, function(x){
    x[x %in% row.names(bulk)]
  })

  #Same gene set length
  gene_set <- lapply(gene_set,function(x){
    x[1:min(sapply(gene_set, length))]
  })

  cat(paste0(round(length(gene_set[[1]])/58,2)*100,"% of the PF A marker genes are expressed in your transcriptomic data"))

  if (length(gene_set$PFA_1) <= 5){
    message("\n5 or fewer PFA subgroup marker genes are expressed in your data. This classification might be inaccurate")
  }
  if (length(gene_set$PFA_1) == 0){
    message("\nNone of the PFA subgroup marker genes are expressed in your data. Cannot classify data")
    return(NULL)
  }

  #Estimate the sampling distribution of enrichment scores
  null_dist <- Make_null(permutations, bulk[,sample(1:ncol(bulk), 1), drop=F], length(gene_set$PFA_1))

  #Preallocate list for each bulk sample
  classified_samples <- NULL
  classified_samples <- vector(mode="list", length=ncol(bulk))
  names(classified_samples) <- colnames(bulk)


  for (i in 1:ncol(bulk)){
    #prep bulk sample
    bulk_sample <- bulk[,colnames(bulk)[i], drop = FALSE]
    bulk_sample <- bulk_sample[order(bulk_sample, decreasing = TRUE),, drop = FALSE]

    #index of gene set in bulk sample
    indx_gene_set <- lapply(gene_set, function(x){
      c(0, which(row.names((bulk_sample)) %in% x))
    })

    #Find enrichment score for each gene set
    final_es <- internal_EnrichmentScore(nrow(bulk), length(gene_set[[1]]), indx_gene_set)
    classified_samples[[colnames(bulk)[i]]]$es <- final_es

    #Pvalue for each enrichment score
    sample_pvalues <- FindPvalue(final_es, null_dist, gene_set)
    classified_samples[[colnames(bulk)[i]]]$pvalue <- sample_pvalues
  }
  return(classified_samples)
}


