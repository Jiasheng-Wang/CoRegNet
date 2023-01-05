

#' Prepare data for Betabinomial fitting
#'
#' This function prepares fitting files for the following step of betabinomial fitting.
#'
#' @param FittingFile_Dir The directory for betabinomial fitting
#' @param Transform_CoPM_Dir The directory of transformed simulation-based Co-perturbation matrix in gene by simulation format
#' @param CoPMatrix The real co-perturbation matrix
#' @param SigIdx Index of gene pairs with significant simluation-based P-values
#' @param P_cutoff The P-value threshold determining significant simulation-based P-values
#' @param FittingFiles Number of files for following betabinomial fitting
#'
#' @return Files ready for fitting in the FittingFile_Dir
#' @export
#'
#' @examples
#' FittingFile_Dir = "~/CopNet_Test/FittingData/", Transform_CoPM_Dir = "~/CopNet_Test/TransformedMatrix/", SigIdx = fitting_summary$qualified_gene_idx, P_cutoff = P_cutoff, FittingFiles = FittingFiles)

FittingDataPreparation <- function(FittingFile_Dir,Transform_CoPM_Dir,CoPMatrix,SigIdx,P_cutoff,FittingFiles){
  if(dir.exists(FittingFile_Dir)){
  }else{
    stop("Error, save_path not found")
  }

  ##Real Co-perturbation Times:
  RealCoTimes <- CoPMatrix[upper.tri(CoPMatrix,diag = FALSE)]
  RealCoTimes <- RealCoTimes[SigIdx]

  ##Max Co-perturbation Times:
  ColsumVector <- diag(CoPMatrix)
  Max_CoPTimes <-  which(upper.tri(CoPMatrix,diag = FALSE), arr.in=TRUE)
  Max_CoPTimes_sig <- Max_CoPTimes[SigIdx,]
  Max_CoPTimes <- cbind.data.frame(GenePair_Idx = SigIdx,
                                   N1 = ColsumVector[Max_CoPTimes_sig[,1]],
                                   N2 = ColsumVector[Max_CoPTimes_sig[,2]])
  Max_CoPTimes$diff <-  Max_CoPTimes$N1-Max_CoPTimes$N2
  Max_CoPTimes$MaxCoperturbationTimes <- ifelse(test = Max_CoPTimes$diff >0, yes = Max_CoPTimes$N2, no = Max_CoPTimes$N1)
  Max_CoperturbationTimes <- Max_CoPTimes[,c("GenePair_Idx","MaxCoperturbationTimes")]

  ##Combined the parameters:
  FITTINGPARAMETERS <- cbind.data.frame(Max_CoperturbationTimes,CoperturbationTimes = RealCoTimes)
  rm(Max_CoperturbationTimes,Max_CoPTimes,Max_CoPTimes_sig,RealCoTimes)
  ##SplitFiles:
  SplitContinuousNumber <- function(Number, nsplits){
    n <- Number
    l <- n %/% nsplits

    startend <- data.frame(matrix(nrow = nsplits, ncol = 3,0), stringsAsFactors = FALSE)
    colnames(startend) <- c("Start","End","Length")
    startend$Start <- seq(0,(nsplits - 1)*l, by = l)+1
    startend$End <- seq(l,nsplits*l, by = l)
    startend$End[nsplits] <- startend$End[nsplits] + n %% nsplits
    startend$Length <- startend$End - startend$Start + 1

    return(startend)
  }
  SplitMatrix <- SplitContinuousNumber(length(SigIdx),FittingFiles)
  for (i in 1:nrow(SplitMatrix)) {
    FittingParameters <- FITTINGPARAMETERS[SplitMatrix$Start[i]:SplitMatrix$End[i],]
    ##load the data of real co-perturbation times & max co-perturbation times:
    save(FittingParameters, file = paste0(FittingFile_Dir,"/FittingParameters_",i,".RData"))
  }

  ##Extract FittingData:
  Nfiles <- list.files(Transform_CoPM_Dir)
  Nfiles <- grep("TransformedMatrix_",Nfiles)
  Nfiles <- length(Nfiles)
  load(paste0(Transform_CoPM_Dir,"/TransformedMatrix_",1,".RData"))
  temp_TCoPM <- GenePairsVsSims[which(rownames(GenePairsVsSims) %in% SigIdx),]

  if (Nfiles >= 2) {
    k <- 2
    for (j in 1:nrow(SplitMatrix)) {
      if (j == nrow(SplitMatrix)) {
        while(nrow(temp_TCoPM) < SplitMatrix$Length[j]){
          load(paste0(Transform_CoPM_Dir,"/TransformedMatrix_",k,".RData"))
          GenePairsVsSims <- GenePairsVsSims[which(rownames(GenePairsVsSims) %in% SigIdx),]
          temp_TCoPM <- rbind.data.frame(temp_TCoPM,GenePairsVsSims)
        }
        FittingData <- temp_TCoPM
        save(x = FittingData, file = paste0(FittingFile_Dir,"/FittingFiles_",j,".RData"))
      }else{
        while(nrow(temp_TCoPM) <= SplitMatrix$Length[j]){
          load(paste0(Transform_CoPM_Dir,"/TransformedMatrix_",k,".RData"))
          GenePairsVsSims <- GenePairsVsSims[which(rownames(GenePairsVsSims) %in% SigIdx),]
          temp_TCoPM <- rbind.data.frame(temp_TCoPM,GenePairsVsSims)
        }
        FittingData <- temp_TCoPM[1:SplitMatrix$Length[j],]
        temp_TCoPM <- temp_TCoPM[-c(1:SplitMatrix$Length[j]),,drop = FALSE]
        save(x = FittingData, file = paste0(FittingFile_Dir,"/FittingFiles_",j,".RData"))
      }
    }
  }else{
    for (j in 1:nrow(SplitMatrix)) {
      FittingData <- temp_TCoPM[SplitMatrix$Start[j]:SplitMatrix$End[j],,drop = FALSE]
      save(x = FittingData, file = paste0(FittingFile_Dir,"/FittingFiles_",j,".RData"))
    }
  }

  return(paste0("Done, Please go and check the files in '",FittingFile_Dir,"'"))
}
