#' Example Data
#'
#' @return example DEG matrix
#' @export
#'
#' @examples
#' DEGMatrix <- DEGMatrix_dat()
create_DEGMatrix <- function(){
  DEGMatrix <- matrix(0, nrow = 20, ncol = 10)
  set.seed(2222)
  c <- sample(2:11,10,replace = TRUE)

  for (i in 1:10) {
    DEGMatrix[sample(1:20,c[i]),i] <- 1
  }


  rownames(DEGMatrix) <- paste0("exp_",seq(1:20))
  colnames(DEGMatrix) <- paste0("gene_",seq(1:10))

  return(DEGMatrix)
}
