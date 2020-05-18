#'
#' Calculate the enrichment score for a single bulk sample and gene set
#'
#' @param num_bulk_genes number of genes the in bulk RNA-seq sample
#' @param num_genes number of genes in gene set
#' @param indx_gene_set numeric vector indexing the gene set in the bulk sample
#'
#' @export
#'
internal_EnrichmentScore <- function(num_bulk_genes, num_genes, indx_gene_set){

  #Preallocate list of each subgroup
  final_es <- NULL
  final_es <- vector(mode="list", length=length(indx_gene_set))
  names(final_es) <- names(indx_gene_set)

  for (i in 1:length(indx_gene_set)){
    ES <- 0
    best_ES <- 0
    #Number of hits
    Nh <- num_genes
    up_ES <- 1/Nh
    #Number of misses
    Nm <- num_bulk_genes - num_genes
    down_ES <- 1/Nm

    for (j in 2:length(indx_gene_set[[i]])){
      ES <- ES - down_ES*(indx_gene_set[[i]][j]-indx_gene_set[[i]][j-1]-1)
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
    final_es[[names(indx_gene_set)[i]]] <- best_ES*best_sign
  }
  return(final_es)
}
