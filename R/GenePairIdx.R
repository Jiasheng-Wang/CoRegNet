

#' Gene pair index to gene names
#'
#' @param Idx The idx of gene pairs
#' @param CoPMatrix The coperturbation matrix
#'
#' @return The name of genes with corresponding index
#' @export
#'
#' @examples
#' GenePairIdx(Idx = fitting_summary$qualified_gene_idx, CoPMatrx = CoPMatrix)
GenePairIdx <- function(Idx, CoPMatrix){
  temp <- which(upper.tri(CoPMatrix, diag = FALSE), arr.ind = TRUE)
  temp <- temp[Idx,]
  temp[,1] <- rownames(CoPMatrix)[as.numeric(temp[,1])]
  temp[,2] <- colnames(CoPMatrix)[as.numeric(temp[,2])]
  temp <- as.data.frame(temp, stringsAsFactors = FALSE)
  temp <- cbind.data.frame(as.numeric(Idx), temp)

  colnames(temp) <- c("GenePairIdx","gene1","gene2")

  return(temp)
}
