#'
#' Calculate the pvlaue associated with the enrichment score
#'
#' @param best_es list for each subgroup with its enrichment score located in "es"
#' @param null_dist list with an estimate of the sampling distribution of enrichment scores (null distribution)
#' @param gene_set character vector of genes
#'
#' @export
#'
FindPvalue <- function(best_es, null_dist, gene_set){

  #Preallocate matrix
  pvalue_es <- NULL
  pvalue_es <- matrix(0, ncol = length(gene_set), nrow = 1)
  colnames(pvalue_es) <- names(gene_set)

  for (i in 1:length(gene_set)){
    es <- round(best_es[[names(gene_set)[i]]]$es, 4)
    if (es > max(null_dist[["freq"]]$x)){
      pvalue_es[1,names(gene_set)[i]] <- 0
    } else {
      if ( es < min(null_dist[["freq"]]$x)){
        pvalue_es[1,names(gene_set)[i]] <- 1
      } else {
        val = sum(abs(null_dist[["freq"]][null_dist[["freq"]]$x > es,"x"])*null_dist[["freq"]][null_dist[["freq"]]$x > es, "freq"])
        pvalue_es[1,names(gene_set)[i]] <- (val / null_dist[["freq"]][1, 'tot_area'])
      }
    }
  }
  return(pvalue_es)
}
