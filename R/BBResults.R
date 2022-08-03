#' Extract Results of BetaBinomial Fitting Results
#'
#' This function generates random co-perturbation Matrix. As each co-perturbation matrix is symmetric, it saves the upper_tir area of the matrix for saving the space
#'
#' @param Result_Dir The directory of BetaBinomial Fitting Results
#' @param GenePair_IdxTable The GenePair Index Table

#'
#' @return A combined Edge list of Co-perturbation Network
#' @export
#'
#' @examples
#' BBResults(Result_Dir = "~/FittingResults/")

BBResults <- function(Result_Dir,GenePair_IdxTable){
  if(dir.exists(Result_Dir)){
  }else{
    stop("Error, Result_Dir not found")
  }
  N_files <- list.files(Result_Dir)
  BBResult0 <- vector("list",length(N_files))
  for (i in 1:length(N_files)) {
    Temp <- read.table(paste0(Result_Dir,"BetaBinomial_Fitting_",i,".txt"),stringsAsFactors = FALSE,
                       header = T)
    Temp <- Temp[,c(1,2,5,6,7)]
    BBResult0[[i]] <- Temp
  }
  rm(Temp)
  BBResult0 <- dplyr::bind_rows(BBResult0)
  colnames(BBResult0)[2] <- "GenePairIdx"
  GenePair_IdxTable$gene1 <- as.character(GenePair_IdxTable$gene1)
  GenePair_IdxTable$gene2 <- as.character(GenePair_IdxTable$gene2)
  BBResult1 <- dplyr::left_join(x = BBResult0, y = GenePair_IdxTable, by = "GenePairIdx")


  return(BBResult1)
}
