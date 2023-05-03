#' Example Data
#'
#' @return example DEG matrix
#' @export
#'
#' @examples
#' DEGMatrix <- DEGMatrix_dat()
create_DEGMatrix <- function(){
  DEGMatrix <- matrix(0, nrow = 10, ncol = 20)
  set.seed(2222)
  c <- sample(2:11,20,replace = TRUE)

  for (i in 1:10) {
    DEGMatrix[i,sample(1:20,c[i])] <- 1
  }


  rownames(DEGMatrix) <- paste0("gene_",seq(1:10))
  colnames(DEGMatrix) <- paste0("exp_",seq(1:20))

  return(DEGMatrix)
}
