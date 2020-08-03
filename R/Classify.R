#'
#' Assign classificatin
#'
#' @param classification a list of ist samples pvalues for each subgroup
#' @param min_pvalue the lowest p-value that will be used to classify a sample (numeric value)
#'
#' @export
#'
Classify <- function(classification, min_pvalue = 0.35){
  assign_class <- NULL
  assign_class <- lapply(classification, function(x){
    if (min(x$pvalue) < min_pvalue){
      assign_class <- c(assign_class, paste(colnames(x$pvalue)[x$pvalue %in% min(x$pvalue)], collapse = " & "))
    } else {
      assign_class <- c(assign_class,"NA")
    }
  })
  assign_class <- unlist(assign_class)
}
