#'
#' Classify into eight (of the nine) molecular subgroups
#'
#' @param bulk a matrix of bulk transcriptomics data (genes by samples)
#' @param permutations the number permutations for a permutation test
#'
#' @export
#'
ClassifyEPN <- function(bulk, permutations = 10000){

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

  # Subset of genes expressed in microarray data
  bulk <- bulk[row.names(bulk) %in% microarray_genes,]

  if (class(permutations) != "numeric" | permutations <= 0 | permutations != round(permutations)){
    message("permutations must be a numeric whole value greater than zero")
    return(NULL)
  }

  #List of differentially expressed genes for each PF A subgroup
  gene_set <- NULL
  gene_set <- vector(mode = "list")
  # gene_set[["ST_EPN_RELA"]] <- c("L1CAM", "GPR153", "ADAP1", "PCP4L1", "FBXO31", "ST8SIA2", "ELAVL3", "CYP27C1", "INPP5A", "P2RX5", "LHX2", "ANGPTL6", "FHDC1", "SLC12A7", "PLCH2", "GNG3", "KCNQ2", "GABRA3", "ADAMTS2", "ADAMTSL4", "WNT7B", "ZNF185", "EPHB2", "ELL3", "CELSR3", "SLC35E4", "COL9A3", "RAP1GAP", "CACNA1I", "HES4", "SRPK3", "KCNK3", "ELFN1", "VIPR1", "KCNK10", "SDC1", "TLX1", "RCOR2", "GAREM2", "HDAC4", "ATP6V1C2", "NOTUM", "LINC00982", "COL18A1", "CACNA2D2", "GUCY1B2", "SP5", "LOC729870", "UNC119", "PRDM16")
  # gene_set[["PF_EPN_B"]] <- c("SLC28A3", "RIPK4", "FNDC7", "LRRC7", "C15orf56", "NRG4", "SNTN", "NXNL2", "PTGDR", "ATP4B", "GAL3ST3", "AGR3", "KIAA1522", "ENTPD2", "CCDC65", "RASSF6", "CYP11A1", "CCNO", "C4orf22", "KCNA5", "C9orf135", "SHANK2", "LHX8", "PITX1", "RASSF10", "OSR2", "TEKT3", "TWIST1", "F2R", "KIAA1683", "SHOX2", "FBXO15", "SLFN13", "CCDC78", "STOML3", "CYP4B1", "SCGB1D1", "C1orf87", "HOXD1", "CRIP1", "CCDC39", "LOC153684", "USP29", "FAM19A2", "NPR3", "SVOPL", "DNAH11", "CDHR3", "TSGA10", "NELL2")
  # gene_set[["SP_EPN"]] <- c("CFTR", "VTN", "ARHGEF28", "GPD1", "MYH2", "HOXC6", "VEPH1", "HOXA6", "WNT16", "HOXB8", "HOXD-AS2", "WT1-AS", "HOXD8", "KCNMB1", "HOXC8", "JPH2", "RASGRF2", "ERMN", "ZNF423", "NKX6-1", "EDDM3A", "MAP3K7CL", "PON3", "FOXA1", "CADPS2", "LRRK2", "CD47", "VSTM1", "SLN", "TRPC5OS", "RAB3C", "HOXB-AS1", "CHRNA3", "DIRAS2", "HOXA5", "SCN9A", "HOXC9", "PLD5", "C20orf85", "HOXA-AS3", "RS1", "MARCKSL1", "AK5", "NR4A3", "LOC100506314", "SCN3B", "ZIC2", "FOXA2", "SLC8A1-AS1","TFF3")
  # gene_set[["PF_EPN_A"]] <- c("CXorf67", "PLAG1", "EPOP", "ZNF556", "NSG1", "TKTL1", "MYO16", "ALDH1L1", "PPARG", "LDOC1", "IGF2BP3", "HCN1", "RYR3", "MAB21L2", "SERPINE2",  "SLC22A3", "IGSF1", "LAMA2", "STK26", "RELN", "PAX3", "GRIA4", "ANGPTL1", "DSC3", "EPHA7", "KCND2",  "VIPR2", "ADGRD1", "SSTR1", "CACNA2D1", "FAM69C", "MECOM", "VGLL3", "BMP5", "FXYD1", "CHODL", "SLC6A13", "ANKRD36BP2", "RORB", "FAP", "ARHGAP36", "ALX1", "CPNE6", "NKAIN4", "FAM110C", "EYA2", "CXCL5", "TNC", "SHISA6", "C16orf89")
  # gene_set[["PF_SE"]] <- c("C2orf82", "LINC01539", "KAZALD1", "SLC38A1", "PRKAG2", "ACTC1", "ST3GAL6", "LOC100507477", "BHMT2", "DBH-AS1", "GRIN3A", "MSX2", "OSGIN2", "SCGB1D2", "METTL7B", "OTOGL", "RUNX3", "KIT", "THBS4", "C20orf197", "AKR1C1", "WSCD2", "PNPLA7", "SLC7A5", "TSHR", "MYOT", "C6orf141", "FSTL4", "LOC100128239", "HOXA-AS2", "KLHL34", "WNK4", "SFRP2", "SEMA3B", "PAH", "NRXN3", "PCDHB12", "UST", "CYP2J2", "LCAT", "ITIH1", "PRODH", "TMEM176B", "OLFML2B", "LOC90246", "ADARB2", "ITGBL1", "FERMT1", "JAKMIP1", "LINC01094")
  # gene_set[["ST_EPN_YAP1"]] <-c("DMRT3", "LGALSL", "SHISA9", "ARL4D", "TUNAR", "CLDN1", "COL26A1", "ANKRD1", "CPZ", "CXCL14", "KRT7", "LIPG", "GPR39", "ARHGAP22", "PBX3", "GRM7", "FNDC1", "EP300-AS1", "LOC101929122", "RBM20", "FABP7", "DOCK3", "HCG22", "FEZF2", "P3H2", "C1QTNF4", "DMRT1", "DIRAS3", "C21orf62", "LRRC3B", "NRXN2", "LINC01551", "DRD4", "SULF1", "SYNDIG1", "SHROOM3", "WBSCR17", "NLGN4X", "HEPH", "OGN", "LRRC2", "SCML1", "KIRREL3", "RAMP3", "ZNF204P", "FBLN5", "TNMD", "LINGO2", "RHOU", "ZIC1")
  # gene_set[["SP_MPE"]] <- c("HOXB13", "PRAC1", "NEFL", "HOXA13", "KMO", "SLC39A2", "ARL15", "HOTAIR", "HNF1B", "C14orf105", "LOC101928731", "LOC101060400", "HOXC10", "CYTL1", "HSPB3", "SCGN", "HOXC13", "ACSM3", "LEFTY1", "TM4SF4", "DPP4", "ZNF385B", "NEFM", "DKK1", "MUC3A", "PRAC2", "HMGCS2", "SLC37A4", "SHC4", "CHL1", "DRAIC", "HOXB9", "MTCP1", "LOXL4", "ARMCX2", "CILP", "SLC24A2", "CHST7", "GNA14", "CFD", "MUC12", "LINC00643", "HOXD10", "SOCS2-AS1", "CACHD1", "IGFBP5", "ITGB3", "SLC35F4", "DLGAP1", "RBPMS-AS1")
  # gene_set[["ST_SE"]] <- c( "SLC7A3", "GALNT13", "SPX", "CCK", "NXPH3", "HSD17B6", "UNC5D", "HTR2C", "NTSR2", "GJB6", "BHLHE22", "SORCS3", "SETD7", "MYH6", "SST", "NKX2-1", "GJB2", "ACVR1C", "CNGA3", "ARPP21", "STON1", "LINC00672", "PTER", "BRINP1", "KCNK2", "STAC", "RYR1", "ATRNL1", "MYZAP", "GABRA2", "HPCAL4", "CA4", "AK4", "MT1M", "PI15", "MT1G", "ETNPPL", "TSHZ2", "PRR16", "GABRB1", "MYBPC1", "LDLRAD3", "LRAT", "GLDN", "FGF7", "SIM2", "IP6K3", "NXPH1", "SIX3-AS1", "IL33" )
  gene_set[["ST_EPN_RELA"]] <- c("L1CAM", "GPR153", "ADAP1", "PCP4L1", "FBXO31", "ST8SIA2", "ELAVL3", "CYP27C1", "INPP5A", "P2RX5", "LHX2", "ANGPTL6", "FHDC1", "SLC12A7", "PLCH2", "GNG3", "KCNQ2", "GABRA3", "ADAMTS2", "ADAMTSL4", "WNT7B", "ZNF185", "EPHB2", "ELL3", "CELSR3", "SLC35E4", "COL9A3", "RAP1GAP", "CACNA1I", "HES4", "SRPK3", "KCNK3", "ELFN1", "VIPR1", "KCNK10", "SDC1", "TLX1", "RCOR2", "GAREM2", "HDAC4", "ATP6V1C2", "NOTUM", "LINC00982", "COL18A1", "CACNA2D2", "GUCY1B2", "SP5", "LOC729870", "UNC119", "PRDM16", "STMN4", "IFI30")
  gene_set[["PF_EPN_B"]] <- c("SLC28A3", "RIPK4", "FNDC7", "LRRC7", "C15orf56", "NRG4", "SNTN", "NXNL2", "PTGDR", "ATP4B", "GAL3ST3", "AGR3", "KIAA1522", "ENTPD2", "CCDC65", "RASSF6", "CYP11A1", "CCNO", "C4orf22", "KCNA5", "C9orf135", "SHANK2", "LHX8", "PITX1", "RASSF10", "OSR2", "TEKT3", "TWIST1", "F2R", "KIAA1683", "SHOX2", "FBXO15", "SLFN13", "CCDC78", "STOML3", "CYP4B1", "SCGB1D1", "C1orf87", "HOXD1", "CRIP1", "CCDC39", "LOC153684", "USP29", "FAM19A2", "NPR3", "SVOPL", "DNAH11", "CDHR3", "TSGA10", "NELL2", "SLC27A2", "AGBL2")
  gene_set[["SP_EPN"]] <- c("CFTR", "VTN", "ARHGEF28", "GPD1", "MYH2", "HOXC6", "VEPH1", "HOXA6", "WNT16", "HOXB8", "HOXD-AS2", "WT1-AS", "HOXD8", "KCNMB1", "HOXC8", "JPH2", "RASGRF2", "ERMN", "ZNF423", "NKX6-1", "EDDM3A", "MAP3K7CL", "PON3", "FOXA1", "CADPS2", "LRRK2", "CD47", "VSTM1", "SLN", "TRPC5OS", "RAB3C", "HOXB-AS1", "CHRNA3", "DIRAS2", "HOXA5", "SCN9A", "HOXC9", "PLD5", "C20orf85", "HOXA-AS3", "RS1", "MARCKSL1", "AK5", "NR4A3", "LOC100506314", "SCN3B", "ZIC2", "FOXA2", "SLC8A1-AS1","TFF3", "ABCD2", "SCN1A", "OR2H1", "FRMD3", "HOXA1", "STAG3")
  gene_set[["PF_EPN_A"]] <- c("CXorf67", "PLAG1", "EPOP", "ZNF556", "NSG1", "TKTL1", "MYO16", "ALDH1L1", "PPARG", "LDOC1", "IGF2BP3", "HCN1", "RYR3", "MAB21L2", "SERPINE2",  "SLC22A3", "IGSF1", "LAMA2", "STK26", "RELN", "PAX3", "GRIA4", "ANGPTL1", "DSC3", "EPHA7", "KCND2",  "VIPR2", "ADGRD1", "SSTR1", "CACNA2D1", "FAM69C", "MECOM", "VGLL3", "BMP5", "FXYD1", "CHODL", "SLC6A13", "ANKRD36BP2", "RORB", "FAP", "ARHGAP36", "ALX1", "CPNE6", "NKAIN4", "FAM110C", "EYA2", "CXCL5", "TNC", "SHISA6", "C16orf89", "PPP1R1A", "AIFM3")
  gene_set[["PF_SE"]] <- c("C2orf82", "LINC01539", "KAZALD1", "SLC38A1", "PRKAG2", "ACTC1", "ST3GAL6", "LOC100507477", "BHMT2", "DBH-AS1", "GRIN3A", "MSX2", "OSGIN2", "SCGB1D2", "METTL7B", "OTOGL", "RUNX3", "KIT", "THBS4", "C20orf197", "AKR1C1", "WSCD2", "PNPLA7", "SLC7A5", "TSHR", "MYOT", "C6orf141", "FSTL4", "LOC100128239", "HOXA-AS2", "KLHL34", "WNK4", "SFRP2", "SEMA3B", "PAH", "NRXN3", "PCDHB12", "UST", "CYP2J2", "LCAT", "ITIH1", "PRODH", "TMEM176B", "OLFML2B", "LOC90246", "ADARB2", "ITGBL1", "FERMT1", "JAKMIP1", "LINC01094", "GSTM3", "LPAR3", "TTR", "HAO1", "SLC22A4", "LRRN4CL", "HHATL", "ENPP6")
  gene_set[["ST_EPN_YAP1"]] <-c("DMRT3", "LGALSL", "SHISA9", "ARL4D", "TUNAR", "CLDN1", "COL26A1", "ANKRD1", "CPZ", "CXCL14", "KRT7", "LIPG", "GPR39", "ARHGAP22", "PBX3", "GRM7", "FNDC1", "EP300-AS1", "LOC101929122", "RBM20", "FABP7", "DOCK3", "HCG22", "FEZF2", "P3H2", "C1QTNF4", "DMRT1", "DIRAS3", "C21orf62", "LRRC3B", "NRXN2", "LINC01551", "DRD4", "SULF1", "SYNDIG1", "SHROOM3", "WBSCR17", "NLGN4X", "HEPH", "OGN", "LRRC2", "SCML1", "KIRREL3", "RAMP3", "ZNF204P", "FBLN5", "TNMD", "LINGO2", "RHOU", "ZIC1", "RHOBTB3", "HOPX", "PTN")
  gene_set[["SP_MPE"]] <- c("HOXB13", "PRAC1", "NEFL", "HOXA13", "KMO", "SLC39A2", "ARL15", "HOTAIR", "HNF1B", "C14orf105", "LOC101928731", "LOC101060400", "HOXC10", "CYTL1", "HSPB3", "SCGN", "HOXC13", "ACSM3", "LEFTY1", "TM4SF4", "DPP4", "ZNF385B", "NEFM", "DKK1", "MUC3A", "PRAC2", "HMGCS2", "SLC37A4", "SHC4", "CHL1", "DRAIC", "HOXB9", "MTCP1", "LOXL4", "ARMCX2", "CILP", "SLC24A2", "CHST7", "GNA14", "CFD", "MUC12", "LINC00643", "HOXD10", "SOCS2-AS1", "CACHD1", "IGFBP5", "ITGB3", "SLC35F4", "DLGAP1", "RBPMS-AS1", "LEFTY2", "PLA2R1", "CTTNBP2", "CALCRL", "COL4A3", "IL18")
  gene_set[["ST_SE"]] <- c( "SLC7A3", "GALNT13", "SPX", "CCK", "NXPH3", "HSD17B6", "UNC5D", "HTR2C", "NTSR2", "GJB6", "BHLHE22", "SORCS3", "SETD7", "MYH6", "SST", "NKX2-1", "GJB2", "ACVR1C", "CNGA3", "ARPP21", "STON1", "LINC00672", "PTER", "BRINP1", "KCNK2", "STAC", "RYR1", "ATRNL1", "MYZAP", "GABRA2", "HPCAL4", "CA4", "AK4", "MT1M", "PI15", "MT1G", "ETNPPL", "TSHZ2", "PRR16", "GABRB1", "MYBPC1", "LDLRAD3", "LRAT", "GLDN", "FGF7", "SIM2", "IP6K3", "NXPH1", "SIX3-AS1", "IL33", "SPECC1")


  #Gene set in bulk data
  gene_set <- lapply(gene_set, function(x){
    x[x %in% row.names(bulk)]
  })

  #Same gene set length
  gene_set <- lapply(gene_set,function(x){
    x[1:min(sapply(gene_set, length))]
  })

  cat(paste0(round(length(gene_set[[1]])/length(gene_set$ST_EPN_RELA),2)*100,"% of the EPN molecular subgroup marker genes are expressed in your transcriptomic data"))

  if (length(gene_set[[1]]) <= 5){
    message("\n5 or fewer EPN molecular subgroup marker genes are expressed in your data. This classification might be inaccurate")
  }
  if (length(gene_set[[1]]) == 0){
    message("\nNone of the EPN moleular subgroup marker genes are expressed in your data. Cannot classify data")
    return(NULL)
  }

  #Estimate the sampling distribution of enrichment scores
  null_dist <- Make_null(permutations, bulk[,sample(1:ncol(bulk), 1), drop=F], length(gene_set[[1]]))

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
