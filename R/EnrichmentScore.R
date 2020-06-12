#'
#' Calculate the enrichment score for a single bulk sample and gene set
#'
#' @param bulk a matrix of a single bulk RNA-seq sample (genes by 1)
#' @param gene_set List of character vector of genes
#'
#' @export
#'
EnrichmentScore <- function(bulk, gene_set){

  #Preallocate list of each subgroup
  final_es <- NULL
  final_es <- vector(mode="list", length=length(gene_set))
  names(final_es) <- names(gene_set)

  #Order bulk RNA genes by expression
  #bulk <- bulk[order(bulk, decreasing = TRUE),,drop = FALSE]

  for (i in 1:length(gene_set)){
    ES <- 0
    best_ES <- 0
    #Number of hits
    Nh <- length(gene_set[[names(gene_set)[i]]])
    up_ES <- 1/Nh
    #Number of misses
    Nm <- dim(bulk)[1] - length(gene_set[[names(gene_set)[i]]])
    down_ES <- 1 /Nm
    x <- NULL
    y <- NULL
    for (j in 1:nrow(bulk)){
      if (row.names(bulk)[j] %in% gene_set[[names(gene_set)[i]]]){
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

