

#' Prepare data for Betabinomial fitting
#'
#' This function prepares fitting files for the following step of betabinomial fitting.
#'
#' @param FittingFile_Dir The directory for betabinomial fitting
#' @param Transform_CoM_Dir The directory of transformed simulation-based Co-perturbation matrix in gene by simulation format
#' @param CoMatrix The real co-perturbation matrix
#' @param SigIdx Index of gene pairs with significant simluation-based P-values
#' @param P_cutoff The P-value threshold determining significant simulation-based P-values
#' @param FittingFiles Number of files for following betabinomial fitting
#'
#' @return Files ready for fitting in the FittingFile_Dir
#' @export
#'
#' @examples
#' FittingFile_Dir = "~/FittingData/", Transform_CoM_Dir = "~/TransformedMatrix/", SigIdx = fitting_summary$qualified_gene_idx, P_cutoff = P_cutoff, FittingFiles = FittingFiles)

FittingDataPreparation <- function(FittingFile_Dir,Transform_CoM_Dir,CoMatrix,SigIdx,P_cutoff,FittingFiles,trans_summary,Ncores){
  if (dir.exists(FittingFile_Dir)) {
  }
  else {
    stop("Error, save_path not found")
  }
  RealCoTimes <- CoMatrix[upper.tri(CoMatrix, diag = FALSE)]
  Nfiles <- nrow(trans_summary)
  Assignment_table <- data.frame(matrix(nrow = Nfiles,
                                        ncol = 4))
  colnames(Assignment_table) <- c("idx", "GenePairs", "istart",
                                  "iend")
  Assignment_table$idx <- seq(1, Nfiles)
  Assignment_table$GenePairs <- as.integer(trans_summary$GenePairs)
  Assignment_table$istart[1] <- 1L
  Assignment_table$iend[nrow(Assignment_table)] <- as.integer(sum(Assignment_table$GenePairs))
  if (nrow(Assignment_table) >= 2) {
    for (ii in 1:(nrow(Assignment_table) - 1)) {
      Assignment_table$iend[ii] <- Assignment_table$istart[ii] + Assignment_table$GenePairs[ii] - 1L
      Assignment_table$istart[ii + 1] <- Assignment_table$iend[ii] + 1L
    }
  }else {
  }
  RealCoTimes <- RealCoTimes[SigIdx]
  ColsumVector <- diag(CoMatrix)
  Max_CoRegTimes <- which(upper.tri(CoMatrix, diag = FALSE),
                          arr.in = TRUE)
  Max_CoRegTimes_sig <- Max_CoRegTimes[SigIdx, ]
  Max_CoRegTimes <- cbind.data.frame(GenePair_Idx = SigIdx,
                                     N1 = ColsumVector[Max_CoRegTimes_sig[, 1]], N2 = ColsumVector[Max_CoRegTimes_sig[,
                                                                                                                      2]])
  Max_CoRegTimes$diff <- Max_CoRegTimes$N1 - Max_CoRegTimes$N2
  Max_CoRegTimes$MaxCoregulationTimes <- ifelse(test = Max_CoRegTimes$diff >
                                                  0, yes = Max_CoRegTimes$N2, no = Max_CoRegTimes$N1)
  Max_CoregulationTimes <- Max_CoRegTimes[, c("GenePair_Idx",
                                              "MaxCoregulationTimes")]
  FITTINGPARAMETERS <- cbind.data.frame(Max_CoregulationTimes,
                                        CoregulationTimes = RealCoTimes)
  FITTINGPARAMETERS$MaxCoregulationTimes <- as.integer(FITTINGPARAMETERS$MaxCoregulationTimes)
  FITTINGPARAMETERS$CoregulationTimes <- as.integer(FITTINGPARAMETERS$CoregulationTimes)
  rm(Max_CoregulationTimes, Max_CoRegTimes, Max_CoRegTimes_sig,
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
    load(paste0(Transform_CoM_Dir, "/TransformedMatrix_", 1,".RData"))
    temp_TCoM <- GenePairsVsSims[which(rownames(GenePairsVsSims) %in% SigIdx), , drop = FALSE]
    for (j in 1:nrow(SplitMatrix)) {
      FittingData <- temp_TCoM[SplitMatrix$Start[j]:SplitMatrix$End[j], , drop = FALSE]
      save(x = FittingData, file = paste0(FittingFile_Dir, "/FittingFiles_", j, ".RData"))
    }
  }else{
    PP_FittingData <- function(FittingDataIdx){
      loading_k <- LoadingList[[FittingDataIdx]]
      temp_TCoM <- data.frame()
      for (k in loading_k) {
        load(paste0(Transform_CoM_Dir, "/TransformedMatrix_", k, ".RData"))
        GenePairsVsSims <- GenePairsVsSims[which(rownames(GenePairsVsSims) %in% SigIdx_list[[FittingDataIdx]]), , drop = FALSE]
        temp_TCoM <- rbind.data.frame(temp_TCoM, GenePairsVsSims)
      }
      FittingData <- temp_TCoM
      save(x = FittingData, file = paste0(FittingFile_Dir, "/FittingFiles_", FittingDataIdx, ".RData"))
    }
    mc <- getOption("mc.cores", Ncores)
    res <- parallel::mclapply(seq(1:nrow(SplitMatrix)), PP_FittingData, mc.cores = mc)

  }

  return(paste0("Done, Please go and check the files in '",
                FittingFile_Dir, "'"))
}
