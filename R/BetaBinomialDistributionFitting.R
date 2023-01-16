#' BetaBinomialDistributionFitting
#'
#' This function enables Betabinomial fitting on the simulation results.
#'
#' @param FittingFile_Dir The directory contains fitting data and fitting parameters
#' @param Result_Dir The directory for saving the results
#' @param Ncores Number of cores to complete the task
#'
#' @return Result will be saved in .txt format
#' @export
#'
#' @examples
#' BetaBinomialDistributionFitting(FittingFile_Dir = "~/FittingData/",Result_Dir = "~/FittingResults/", Ncores = 2)

BetaBinomialDistributionFitting<- function(FittingFile_Dir,Result_Dir, Ncores = 1){
  if(dir.exists(Result_Dir)){
  }else{
    stop("Error, save_path not found")
  }

  NFittings <- list.files(FittingFile_Dir)
  NFittings <- grep("FittingFiles_",NFittings)
  NFittings <- length(NFittings)

  BetaBinomialFitting <- function(data_idx){
    ##load the data of real co-perturbation times & max co-perturbation times:
    load(paste0(FittingFile_Dir,"/FittingParameters_",data_idx,".RData"))
    ##load data for fitting:
    load(paste0(FittingFile_Dir,"/FittingFiles_",data_idx,".RData"))
    ##Check if IDX are matched:
    if(all(rownames(FittingParameters) == rownames(FittingData))){
      #MATCHED!Continue
    }else{
      message("Gene Pair Index in Fitting Data and Fitting Parameters are NOT MATCHED! Fitting Terminated.")
      return()
    }

    ##Generate the frame of output in "Result"
    Result <- data.frame(matrix(ncol = 7, nrow = nrow(FittingParameters), 0),stringsAsFactors = FALSE)
    colnames(Result) <- c("File_idx","GenePair_Idx","mu","sigma","CoperturbationTimes","MaxCoperturbationTimes","Fitting_Pvalue")

    ##Write part information into Result:
    Result$File_idx <- data_idx
    Result$GenePair_Idx <- FittingParameters$GenePair_Idx
    Result$CoperturbationTimes <- FittingParameters$CoperturbationTimes
    Result$MaxCoperturbationTimes <- FittingParameters$MaxCoperturbationTimes

    rm(FittingParameters)
    ##Start to do the fitting
    for (j in 1:nrow(Result)) {
      simdat <- as.numeric(FittingData[j,])
      N <- Result$MaxCoperturbationTimes[j]
      bdata <- cbind.data.frame("y" = simdat, "N" = N)
      invisible(capture.output(fit <- try(VGAM::vglm(cbind(y, N-y) ~ 1, family = VGAM::betabinomial, bdata, trace=TRUE), silent = TRUE)))

      if('try-error' %in% class(fit)){
        Result$mu[j] <- NA
        Result$sigma[j] <- NA
        Result$Fitting_Pvalue[j] <- NA
        next
      }else{
        mu = VGAM::Coef(fit)[1]
        sigma = VGAM::Coef(fit)[2]
        Result$mu[j] <- mu
        if (sigma == 0) {
          sigma <- 0.1e-9
          Result$sigma[j] <- sigma
        }else{
          Result$sigma[j] <- sigma
        }
        ##Calculate the corresponding Fitting_Pvalue:
        B <- seq(from = Result$CoperturbationTimes[j], to = Result$MaxCoperturbationTimes[j], by = 1)
        A <- sum(gamlss.dist::dBB(x = B,mu = mu,sigma = sigma,bd = Result$MaxCoperturbationTimes[j]))

        Result$Fitting_Pvalue[j] <- A
      }
    }
    write.table(Result,paste0(Result_Dir,"/BetaBinomial_Fitting_",data_idx,".txt"))
  }

  mc <- getOption("mc.cores", Ncores)
  res <- parallel::mclapply(seq(1:NFittings), BetaBinomialFitting, mc.cores = mc)

  return(paste0("Done, Please go and check the files in '",Result_Dir,"'"))
}
