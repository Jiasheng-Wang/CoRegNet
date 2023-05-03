#' Sup function for matrix transformation
#'
#' This function generates a summary of transforming gene-gene co-regulation matrix to GenePairs-Simulations matrix
#'
#' @param RandCoM_directory The directory where simulated Co-regulation matrix saved
#' @param Ngenes Number of genes in the binary DEG matrix
#' @param NTransformedMatrix Number of GenePair vs. SimTimes matrix required
#'
#' @return A summary data frame helps to decide whether the transformation will be good
#' @export
#'
#' @examples
#' trans_summary <- Transform_CoM_sup(RandCoM_directory = "~/RandBM/", Ngenes = 10, NTransformedMatrix = 1)
Transform_CoM_sup <- function(RandCoM_directory, Ngenes, NTransformedMatrix = 1){
  NGenePairs <- (Ngenes^2 - Ngenes)/2
  files<-list.files(RandCoM_directory,full.names = T)
  vect_size <- sapply(files, file.size)
  size_files <- sum(vect_size)/1024/1024
  sup_dataframe <- data.frame(matrix(nrow = NTransformedMatrix, ncol = 4))
  colnames(sup_dataframe) <- c("idx","#GenePairs","MB","GB")
  sup_dataframe$idx <- seq(1,NTransformedMatrix)

  mean_genepairs <- NGenePairs %/% NTransformedMatrix
  mod_genepairs <- NGenePairs %% NTransformedMatrix

  if(mod_genepairs ==0){
    sup_dataframe$`#GenePairs` <- mean_genepairs
    sup_dataframe$MB <- size_files/NGenePairs*sup_dataframe$`#GenePairs`
    sup_dataframe$GB <- size_files/NGenePairs*sup_dataframe$`#GenePairs`/1024
  }else{
    sup_dataframe$`#GenePairs` <- mean_genepairs
    sup_dataframe$`#GenePairs`[1:mod_genepairs] <- sup_dataframe$`#GenePairs`[1:mod_genepairs]+1
    sup_dataframe$MB <- size_files/NGenePairs*sup_dataframe$`#GenePairs`
    sup_dataframe$GB <- size_files/NGenePairs*sup_dataframe$`#GenePairs`/1024
  }

  return(sup_dataframe)
}
