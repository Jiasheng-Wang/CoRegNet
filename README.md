# Please check our Recount3 CoRegNet via: http://CRN.liuzlab.org/

# Co-Regulation Model: CoRegNet package Manual

### Load your DEGMatrix: Each row is a gene; Each column is an experiment/contrast
### OR, here, as an example, you can use the following function to quickly create a qualified DEG Matrix
DEGMatrix <- create_DEGMatrix()

CoRegMatrix <- DEGMatrix %*% t(DEGMatrix)

CoTimes <- CoRegMatrix[upper.tri(CoRegMatrix,diag = FALSE)]


### Generate simulation-based co-regulation matrix with fixed rowsum and colsum of binary DEG matrix:
RandCoRegMatrix_Generation(Ncores = 10, EachCoreSims = 50, RowSumVector = rowSums(DEGMatrix), ColSumVector = colSums(DEGMatrix),save_path = "~/CoregNet_Test/RBM/")

### Transform the simulation-based co-regulation matrix into gene by simulation format for p-value calculation:
n <- 1

trans_summary <- Transform_CoM_sup(RandCoM_directory = "~/CoregNet_Test/RBM/", Ngenes = ncol(DEGMatrix), NTransformedMatrix = n)

Transform_CoM(RandCoM_directory = "\~/CoregNet_Test/RBM/",save_path = "~/CoregNet_Test/TransformedMatrix/", NGenePairs = trans_summary[,2], NTransformedMatrix = n, Ncores = 1)

### Calculate the simulation-based P-values
SimPCalculation(transformed_CoRegM_Dir = "\~/CoregNet_Test/TransformedMatrix/", save_path = "~/CoregNet_Test/SimP/",CoTimes = CoTimes,trans_summary = trans_summary,Ncores = 2)

### Do beta binomial fitting only on those simulation p-value smaller than the cutoff to save time
#Decide the p-value cutoff and number of cores to do betabinomial fittings:
P_cutoff <- 0.05

FittingFiles <- 2

fitting_summary <- BBFitting_sup(SimP_Dir = "~/CoregNet_Test/SimP/", P_cutoff = P_cutoff, FittingFiles = FittingFiles)

### Prepare files for Betabinomial fitting:
FittingDataPreparation(FittingFile_Dir = "\~/CoregNet_Test/FittingData/",Transform_CoM_Dir = "~/CoregNet_Test/TransformedMatrix/",
                       CoRegMatrix = CoRegMatrix, SigIdx = fitting_summary$qualified_gene_idx,
                       trans_summary = trans_summary, P_cutoff = P_cutoff, FittingFiles = FittingFiles)

### BetabinomialFitting:
BetaBinomialDistributionFitting(FittingFile_Dir = "\~/CoregNet_Test/FittingData/",Result_Dir = "~/CoregNet_Test/FittingResults/", Ncores = 2)

GenePair_IdxTable <- GenePairIdx(Idx = fitting_summary$qualified_gene_idx, CoRegMatrix = CoRegMatrix)
