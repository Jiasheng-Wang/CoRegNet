# Please check our Recount3 CoRegNet via:
http://CRN.liuzlab.org/

# Co-Regulation Model: CoRegNet package Manual

### Step 1: Load your DEGMatrix: Each row is a gene; Each column is an experiment/contrast
#### You can download *DEGMatrix_Recount3_Subset.csv* file, which contains a subset experiments from Recount3.<br>OR, here, as an example, you can use the following function to quickly create a qualified DEG Matrix
DEGMatrix <- create_DEGMatrix()

#### If you are using *DEGMatrix_Recount3_Subset.csv* file and you only take it as an example, it is recommended to run: 
DEGMatrix <- DEGMatrix[rowSums(DEGMatrix) > 40,]
#### This code keeps genes that are identified as DEGs in more than 40 experiments, and it will keep around 3000 genes in order to save the computational power;<br> If you would like to get more accurate results, it is recommended to set a lower threshold than 40 and keep more genes for the following analysis.

CoRegMatrix <- DEGMatrix %*% t(DEGMatrix)

CoTimes <- CoRegMatrix[upper.tri(CoRegMatrix,diag = FALSE)]


### Step 2: Generate simulation-based co-regulation matrix with fixed rowsum and colsum of binary DEG matrix:
RandCoRegMatrix_Generation(Ncores = 10, EachCoreSims = 50, RowSumVector = rowSums(DEGMatrix), ColSumVector = colSums(DEGMatrix),save_path = "~/CoregNet_Test/RBM/")

### Step 3: Transform the simulation-based co-regulation matrix into gene by simulation format for p-value calculation:
n <- 10

trans_summary <- Transform_CoM_sup(RandCoM_directory = "~/CoregNet_Test/RBM/", Ngenes = ncol(DEGMatrix), NTransformedMatrix = n)

Transform_CoM(RandCoM_directory = "\~/CoregNet_Test/RBM/",save_path = "~/CoregNet_Test/TransformedMatrix/", NGenePairs = trans_summary[,2], NTransformedMatrix = n, Ncores = 2)

### Step 4: Calculate the simulation-based P-values
SimPCalculation(transformed_CoRegM_Dir = "\~/CoregNet_Test/TransformedMatrix/", save_path = "~/CoregNet_Test/SimP/",CoTimes = CoTimes,trans_summary = trans_summary,Ncores = 2)

### Step 5: Do beta binomial fitting only on those simulation p-value smaller than the cutoff to save time
#Decide the p-value cutoff and number of cores to do betabinomial fittings:
P_cutoff <- 0.05

FittingFiles <- 2

fitting_summary <- BBFitting_sup(SimP_Dir = "~/CoregNet_Test/SimP/", P_cutoff = P_cutoff, FittingFiles = FittingFiles)

### Step 6: Prepare files for Betabinomial fitting:
FittingDataPreparation(FittingFile_Dir = "\~/CoregNet_Test/FittingData/",Transform_CoM_Dir = "~/CoregNet_Test/TransformedMatrix/",
                       CoRegMatrix = CoRegMatrix, SigIdx = fitting_summary$qualified_gene_idx,
                       trans_summary = trans_summary, P_cutoff = P_cutoff, FittingFiles = FittingFiles)

### Step 7: BetabinomialFitting:
BetaBinomialDistributionFitting(FittingFile_Dir = "\~/CoregNet_Test/FittingData/",Result_Dir = "~/CoregNet_Test/FittingResults/", Ncores = 2)

GenePair_IdxTable <- GenePairIdx(Idx = fitting_summary$qualified_gene_idx, CoRegMatrix = CoRegMatrix)

CoRegulationNetwork <- BBResults(Result_Dir = FittingResultPath,GenePair_IdxTable = GenePair_IdxTable)
