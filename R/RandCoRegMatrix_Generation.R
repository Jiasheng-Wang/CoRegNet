#' Generate Random Binary DEG matrix
#'
#' This function generates random co-regulation Matrix. As each co-regulation matrix is symmetric, it saves the upper_tir area of the matrix for saving the space
#'
#' @param Ncores Number of CPU to simulate the binary matrix with given RowSum and ColSum, and calculate the co-regulation matrix.
#' @param EachCoreSims Number of matrix each core simulates
#' @param RowSumVector The fixed RowSum Vecotor of true binary DEG matrix - Each row is a gene
#' @param ColSumVector The fixed ColSum Vecotor of true binary DEG matrix - Each column is an experiment
#' @param save_path The directory to save the simulated co-regulation matrix
#' @param SaveSparse Save Results in Sparse Matrix/Vector

#' @return save Simulated co-regulation matrix with given rowsum and colsum of the true binary DEG matrix to the given directory
#' @export
#'
#' @examples
#' RandCoPMatrix_Generation(Ncores = 20,EachCoreSims = 50,RowSumVector = c(2,3,4,1,5,2,3,6,3),ColSumVector = c(4,3,2,2,2,5,3,2,3,3), save_path = "~/RandBM/")
RandCoRegMatrix_Generation <- function(Ncores = 1,EachCoreSims = 1000, RowSumVector,ColSumVector,save_path,SaveSparse = TRUE){
  ##Check save_path:
  if(dir.exists(save_path)){
  }else{
    stop("Error, save_path not found")
  }
  ##Check Simuation Times
  if(Ncores*EachCoreSims < 1000){
    warning("The total simulations is too small")
  }else{
  }

  RBMGeneration <- function(Ncores){
    iters <- EachCoreSims
    cc <- 0
    set.seed(Ncores*1000)
    while (cc < iters) {
      r <- RowSumVector
      s <- ColSumVector
      n <- length(r)
      m <- length(s)

      M0 <- matrix(0,nrow = n, ncol = m)
      ra <- r/m
      sa <- s/n

      while(sum(r)!=0 || sum(s)!=0){
        if(max(ra) >= max(sa)){
          Ridx <- which(max(ra) == ra)
          Ridx <- unname(Ridx[1])
          trysi <- try(which(s > 0)[sample(1:length(which(s > 0)),r[Ridx])], silent = TRUE)
          if('try-error' %in% class(trysi)){
            break
          }else{
            M0[Ridx,trysi] <- 1
            r[Ridx] <- 0
            s[trysi] <- s[trysi]-1
            n <- sum(r!=0)
            m <- sum(s!=0)
            ra <- r/m
            sa <- s/n
          }
        }else{
          Cidx <- which(max(sa) == sa)
          Cidx <- unname(Cidx[1])
          tryri <- try(which(r > 0)[sample(1:length(which(r > 0)),s[Cidx])], silent = TRUE)
          if('try-error' %in% class(tryri)){
            break
          }else{
            M0[tryri,Cidx] <- 1
            s[Cidx] <- 0
            r[tryri] <- r[tryri] - 1
            n <- sum(r!=0)
            m <- sum(s!=0)
            ra <- r/m
            sa <- s/n

          }
        }
      }

      if(all(colSums(M0) == ColSumVector) &&
         all(rowSums(M0) == RowSumVector)){
        cc <- cc + 1
        M1 <- M0 %*% t(M0)
        M2 <- M1[upper.tri(M1,diag = FALSE)]

        if(SaveSparse){
          M2  <- as(M2 , "sparseVector")
          save(x = M2,file = paste0(save_path,"/CRM_",(Ncores-1)*iters+cc,".RData"))
        }else{
          save(x = M2,file = paste0(save_path,"/CRM_",(Ncores-1)*iters+cc,".RData"))
        }
      
        rm(M0,M1,M2)
        gc()
      }else{
        rm(M0)
        gc()
      }
    }
  }
  mc <- getOption("mc.cores", Ncores)
  res <- parallel::mclapply(seq(1:Ncores), RBMGeneration, mc.cores = mc)
  return(paste0("Done, Please go and check the files in '",save_path,"'"))
}


