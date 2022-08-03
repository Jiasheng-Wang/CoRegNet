# CoPNet
 Co-perturbation Model: CopNet package Manual

##Load your DEGMatrix: Each row is a contrast; Each column is a gene
##Here, as an example, you can just use the function to create a qualified DEG Matrix
DEGMatrix <- create_DEGMatrix()
CoPMatrix <- t(DEGMatrix) %*% DEGMatrix
CoTimes <- CoPMatrix[upper.tri(CoPMatrix,diag = FALSE)]

#Generate simulation-based co-perturbation matrix with fixed rowsum and colsum of binary DEG matrix:
RandCoPMatrix_Generation(Ncores = 10, EachCoreSims = 50, RowSumVector = rowSums(DEGMatrix), ColSumVector = colSums(DEGMatrix),save_path = "~/CopNet_Test/RBM/")

##Transform the simulation-based co-perturbation matrix into gene by simulation format for p-value calculation:
n <- 1
trans_summary <- Transform_CoPM_sup(RandCoPM_directory = "~/CopNet_Test/RBM/", Ngenes = ncol(DEGMatrix), NTransformedMatrix = n)
Transform_CoPM(RandCoPM_directory = "~/CopNet_Test/RBM/", save_path = "~/CopNet_Test/TransformedMatrix/", NGenePairs = trans_summary[,2], NTransformedMatrix = n, Ncores = 1)

##Calculate the simulation-based P-values
SimPCalculation(transformed_CoPM_Dir = "~/CopNet_Test/TransformedMatrix/", save_path = "~/CopNet_Test/SimP/",CoTimes = CoTimes, trans_summary = trans_summary, Ncores = 2)

##Do beta binomial fitting only on those simulation p-value smaller than the cutoff to save time
## Decide the p-value cutoff and number of cores to do betabinomial fittings:
P_cutoff <- 0.05
FittingFiles <- 2

fitting_summary <- BBFitting_sup(SimP_Dir = "~/CopNet_Test/SimP/", P_cutoff = P_cutoff, FittingFiles = FittingFiles)
fitting_summary_1 <- fitting_summary$summary

##Prepare files for Betabinomial fitting:
FittingDataPreparation(FittingFile_Dir = "~/CopNet_Test/FittingData/",Transform_CoPM_Dir = "~/CopNet_Test/TransformedMatrix/",
                       CoPMatrix = CoPMatrix, SigIdx = fitting_summary$qualified_gene_idx, P_cutoff = P_cutoff, FittingFiles = FittingFiles)

##BetabinomialFitting:
BetaBinomialDistributionFitting(FittingFile_Dir = "~/CopNet_Test/FittingData/",Result_Dir = "~/CopNet_Test/FittingResults/", Ncores = 2)
GenePair_IdxTable <- GenePairIdx(Idx = fitting_summary$qualified_gene_idx, CoPMatrix = CoPMatrix)
