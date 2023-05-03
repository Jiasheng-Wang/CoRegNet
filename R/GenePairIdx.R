

#' Gene pair index to gene names
#'
#' @param Idx The idx of gene pairs
#' @param CoMatrix The coregulation matrix
#'
#' @return The name of genes with corresponding index
#' @export
#'
#' @examples
#' GenePairIdx(Idx = fitting_summary$qualified_gene_idx, CoMatrx = CoMatrix)
GenePairIdx <- function(Idx, CoMatrix){
  temp <- which(upper.tri(CoMatrix, diag = FALSE), arr.ind = TRUE)
  temp <- temp[Idx,]
  temp[,1] <- rownames(CoMatrix)[as.numeric(temp[,1])]
  temp[,2] <- colnames(CoMatrix)[as.numeric(temp[,2])]
  temp <- as.data.frame(temp, stringsAsFactors = FALSE)
  temp <- cbind.data.frame(as.numeric(Idx), temp)

  colnames(temp) <- c("GenePairIdx","gene1","gene2")

  return(temp)
}
