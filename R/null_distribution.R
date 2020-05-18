#'
#' Enrichment score for gene_set in single bulk RNA-seq sample
#'
#' @param num_bulk_genes number of genes the in bulk RNA-seq sample
#' @param num_genes number of genes in gene set
#' @param indx_gene_set numeric vector indexing the gene set in genes the list of genes expressed in the bulk sample
#'
#' @export
#'
null_EnrichmentScore <- function(num_bulk_genes, num_genes, indx_gene_set){
  ES <- 0
  best_ES <- 0
  #Number of hits
  Nh <- num_genes
  up_ES <- 1/Nh
  #Number of misses
  Nm <- num_bulk_genes - num_genes
  down_ES <- 1/Nm

  for (i in 2:length(indx_gene_set)){
    ES <- ES - down_ES*(indx_gene_set[i]-indx_gene_set[i-1]-1)
    if (abs(ES) > best_ES){
      best_sign <- sign(ES)
      best_ES <- abs(ES)
    }
    ES <- ES + up_ES
    if (abs(ES) > best_ES){
      best_sign <- ES/abs(ES)
      best_ES <- abs(ES)
    }
  }
  return(best_ES*best_sign)
}


#' Create null distribution of enrichment scores from random gene sets (estimate the enrichment score sampling distribution)
#'
#' @param permutations the number of enrichment scores that will be calculated in order to estimate the sample distribution of enrichment scores
#' @param bulk a matrix of a single bulk RNA-seq sample (genes by 1)
#' @param num_genes the number of genes in the gene set
#'
#' @importFrom plyr count
#'
#' @export
#'
Make_null <- function(permutations, bulk, num_genes){

  #Preallocate list of each subgroup
  null_es <- NULL
  null_es <- vector(mode="list", length=2)
  names(null_es) <- c("null values", "freq")

  #Order bulk RNA list of genes by expression
  bulk <- bulk[order(bulk, decreasing = TRUE),, drop = FALSE]

  for (i in 1:permutations){
    new_es <- NULL
    gene_set <- NULL
    gene_set <- sample(row.names(bulk), num_genes, replace = FALSE)
    indx_gene_set <- 0
    indx_gene_set <- c(indx_gene_set, which(row.names((bulk)) %in% gene_set))
    gene_set <- row.names(bulk)[indx_gene_set]
    new_es <- null_EnrichmentScore(nrow(bulk), num_genes, indx_gene_set)
    null_es[["null values"]] <- c(null_es[["null values"]], new_es)
  }

  null_es[["freq"]] <- count(round(null_es[["null values"]], 4))
  null_es[["freq"]][1, 'tot_area'] <- sum(abs(null_es[["freq"]]$x)*null_es[["freq"]]$freq)
  return(null_es)
}
