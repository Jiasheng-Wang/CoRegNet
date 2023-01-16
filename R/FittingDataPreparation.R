

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

FittingDataPreparation <- function(FittingFile_Dir,Transform_CoPM_Dir,CoPMatrix,SigIdx,P_cutoff,FittingFiles,trans_summary,Ncores){
  if (dir.exists(FittingFile_Dir)) {
  }
  else {
    stop("Error, save_path not found")
  }
  RealCoTimes <- CoPMatrix[upper.tri(CoPMatrix, diag = FALSE)]
  Nfiles <- nrow(trans_summary)
  Assignment_table <- data.frame(matrix(nrow = Nfiles,
                                        ncol = 4))
  colnames(Assignment_table) <- c("idx", "GenePairs", "istart",
                                  "iend")
  Assignment_table$idx <- seq(1, Nfiles)
  Assignment_table$GenePairs <- trans_summary$`#GenePairs`
  Assignment_table$istart[1] <- 1
  Assignment_table$iend[nrow(Assignment_table)] <- sum(Assignment_table$GenePairs)
  if (nrow(Assignment_table) >= 2) {
    for (ii in 1:(nrow(Assignment_table) - 1)) {
      Assignment_table$iend[ii] <- Assignment_table$istart[ii] +
        Assignment_table$GenePairs[ii] - 1
      Assignment_table$istart[ii + 1] <- Assignment_table$iend[ii] +
        1
    }
  }else {
  }
  RealCoTimes <- RealCoTimes[SigIdx]
  ColsumVector <- diag(CoPMatrix)
  Max_CoPTimes <- which(upper.tri(CoPMatrix, diag = FALSE),
                        arr.in = TRUE)
  Max_CoPTimes_sig <- Max_CoPTimes[SigIdx, ]
  Max_CoPTimes <- cbind.data.frame(GenePair_Idx = SigIdx,
                                   N1 = ColsumVector[Max_CoPTimes_sig[, 1]], N2 = ColsumVector[Max_CoPTimes_sig[,
                                                                                                                2]])
  Max_CoPTimes$diff <- Max_CoPTimes$N1 - Max_CoPTimes$N2
  Max_CoPTimes$MaxCoperturbationTimes <- ifelse(test = Max_CoPTimes$diff >
                                                  0, yes = Max_CoPTimes$N2, no = Max_CoPTimes$N1)
  Max_CoperturbationTimes <- Max_CoPTimes[, c("GenePair_Idx",
                                              "MaxCoperturbationTimes")]
  FITTINGPARAMETERS <- cbind.data.frame(Max_CoperturbationTimes,
                                        CoperturbationTimes = RealCoTimes)
  rm(Max_CoperturbationTimes, Max_CoPTimes, Max_CoPTimes_sig,
     RealCoTimes)
  gc()
  SplitContinuousNumber <- function(Number, nsplits) {
    n <- Number
    l <- n%/%nsplits
    k <- n%%nsplits
    startend <- data.frame(matrix(nrow = nsplits, ncol = 3,
                                  0), stringsAsFactors = FALSE)
    length_seq <- rep(l, nsplits)
    if (k == 0) {
    }
    else {
      length_seq[1:k] <- length_seq[1:k] + 1
    }
    colnames(startend) <- c("Start", "End", "Length")
    startend$Length <- length_seq
    startend$Start[1] <- 1
    startend$End[1] <- startend$Start[1] + startend$Length[1] -
      1
    if (nrow(startend) > 1) {
      for (ii in 2:nrow(startend)) {
        startend$Start[ii] <- startend$End[ii - 1] +
          1
        startend$End[ii] <- startend$Start[ii] + startend$Length[ii] -
          1
      }
    }
    else {
    }
    return(startend)
  }
  SplitMatrix <- SplitContinuousNumber(length(SigIdx), FittingFiles)
  SigIdx_list <- vector("list",length = FittingFiles)
  for (i in 1:nrow(SplitMatrix)) {
    FittingParameters <- FITTINGPARAMETERS[SplitMatrix$Start[i]:SplitMatrix$End[i], ]
    SigIdx_list[[i]] <- FittingParameters$GenePair_Idx
    save(FittingParameters, file = paste0(FittingFile_Dir,
                                          "/FittingParameters_", i, ".RData"))
  }

  LoadingList <- vector("list",length = FittingFiles)
  for(i in 1:length(SigIdx_list)){
    st_g <- SigIdx_list[[i]][1]
    ed_g <- SigIdx_list[[i]][length(SigIdx_list[[i]])]
    st_f <- max(which(Assignment_table$istart <= st_g))
    ed_f <- min(which(Assignment_table$iend >= ed_g))
    LoadingList[[i]] <- st_f:ed_f
  }

  if (Nfiles == 1) {
    load(paste0(Transform_CoPM_Dir, "/TransformedMatrix_", 1,".RData"))
    temp_TCoPM <- GenePairsVsSims[which(rownames(GenePairsVsSims) %in% SigIdx), , drop = FALSE]
    for (j in 1:nrow(SplitMatrix)) {
      FittingData <- temp_TCoPM[SplitMatrix$Start[j]:SplitMatrix$End[j], , drop = FALSE]
      save(x = FittingData, file = paste0(FittingFile_Dir, "/FittingFiles_", j, ".RData"))
    }
  }else{
    PP_FittingData <- function(FittingDataIdx){
      loading_k <- LoadingList[[FittingDataIdx]]
      temp_TCoPM <- data.frame()
      for (k in loading_k) {
        load(paste0(Transform_CoPM_Dir, "/TransformedMatrix_", k, ".RData"))
        GenePairsVsSims <- GenePairsVsSims[which(rownames(GenePairsVsSims) %in% SigIdx_list[[FittingDataIdx]]), , drop = FALSE]
        temp_TCoPM <- rbind.data.frame(temp_TCoPM, GenePairsVsSims)
      }
      FittingData <- temp_TCoPM
      save(x = FittingData, file = paste0(FittingFile_Dir, "/FittingFiles_", FittingDataIdx, ".RData"))
    }
    mc <- getOption("mc.cores", Ncores)
    res <- parallel::mclapply(seq(1:nrow(SplitMatrix)), PP_FittingData, mc.cores = mc)

  }

  return(paste0("Done, Please go and check the files in '",
                FittingFile_Dir, "'"))
}
