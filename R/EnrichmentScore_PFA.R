#'
#' Calculate the enrichment score for a single bulk_sample sample and the PFA gene sets
#'
#' @param bulk_sample a matrix of a single bulk_sample RNA-seq sample (genes by 1)
#'
#' @export
#'
EnrichmentScore_PFA <- function(bulk_sample){

  #List of differentially expressed genes for each PF A subgroup
  gene_set <- NULL
  gene_set <- vector(mode = "list")
  gene_set[["PFA_1"]] <- c("HOTAIRM1", "HOXB3", "HOXB2", "HOXB4", "HOXA2", "HOXA3", "HOXA1", "HOXA-AS2", "MIR10A", "HOXB-AS1", "HOXB-AS2", "AC106786.1", "ZNF503", "LPAR3", "PRELP", "SEMA3C", "PRDM6", "ZNF503-AS2", "SKAP2", "SLC35F2", "HOXC4", "HOXD4", "PPP4R4", "AJAP1", "HOXA4", "FGFR2", "FGFRL1", "TICAM1", "LRRTM4", "GPR133", "CXXC5", "CCDC85A", "SEMA3B", "CNKSR3", "DPT", "ZNF703", "HCN1", "LRP2", "LINC00626", "HSPB1", "STAC2", "VLDLR", "AC007743.1", "COTL1", "RBMS3", "WIF1", "INPP4B", "RHOBTB3", "TNFRSF11A", "EGF", "LAMA2", "PAX8", "HOXD-AS1", "SERPINA12", "PLSCR4", "STXBP5L", "PIK3R1", "ARHGEF40")
  gene_set[["PFA_2"]] <- c("EN2", "MPPED2", "CNPY1", "PTPN3", "ST6GAL2", "AC008060.5", "B3GALTL", "AC005235.1", "LRRC8D", "SYT10", "WBSCR17", "PCBP3", "EFCC1", "LPO", "SMCO2", "NEDD4L", "EYA4", "LOXL4", "RNF43", "SERPINI1", "ITIH5", "GPR39", "PROB1", "TOX3", "ARNTL2", "BZRAP1","RFX3", "NTS", "ADAMTS3", "TCTEX1D1", "NR1H4", "KRT81", "TEX15", "ASB18", "SH3GL3", "FAM65C", "LINC00907", "LMX1B", "TRPC6", "LECT1", "AL592528.1", "SEMA3D", "TP73", "ELMOD1", "LGR4", "TMEM254", "CHL1", "GREB1L", "SPATA42", "TDRD1", "ITGA8", "KCNG2", "SVOPL", "NRSN1", "ERC2", "CATSPERD", "SPATA6", "CSPP1")

  #Gene set in bulk_sample data
  gene_set <- lapply(gene_set, function(x){
    x[x %in% row.names(bulk_sample)]
  })

  #Same gene set length
  gene_set <- lapply(gene_set,function(x){
    x[1:min(sapply(gene_set, length))]
  })

  cat(paste0(round(length(gene_set[[1]])/58,2)*100,"% of the PFA subtype marker genes are expressed in your bulk_sample data"))

  if (length(gene_set[[1]]) <= 5){
    message("\n5 or fewer PFA subtype marker genes are expressed in your data. This classification might be inaccurate")
  }
  if (length(gene_set[[1]]) == 0){
    message("\nNone of the PFA subtype marker genes are expressed in your data. Cannot classify data")
    return(NULL)
  }


  #Preallocate list of each subgroup
  final_es <- NULL
  final_es <- vector(mode="list", length=length(gene_set))
  names(final_es) <- names(gene_set)

  #Order bulk_sample RNA genes by expression
  bulk_sample <- bulk_sample[order(bulk_sample, decreasing = TRUE),,drop = FALSE]

  for (i in 1:length(gene_set)){
    ES <- 0
    best_ES <- 0
    #Number of hits
    Nh <- length(gene_set[[names(gene_set)[i]]])
    up_ES <- 1/Nh
    #Number of misses
    Nm <- dim(bulk_sample)[1] - length(gene_set[[names(gene_set)[i]]])
    down_ES <- 1 /Nm
    x <- NULL
    y <- NULL
    for (j in 1:nrow(bulk_sample)){
      if (row.names(bulk_sample)[j] %in% gene_set[[names(gene_set)[i]]]){
        ES <- ES + up_ES
      } else{
        ES <- ES - down_ES
      }
      x <- c(x, j)
      y <- c(y, ES)
      if (abs(ES) > best_ES){
        best_sign <- ES/abs(ES)
        best_ES <- abs(ES)
      }
    }
    final_es[[names(gene_set)[i]]]$es <- best_ES*best_sign
    final_es[[names(gene_set)[i]]]$x <- x
    final_es[[names(gene_set)[i]]]$y <- y
  }

  return(final_es)
}
