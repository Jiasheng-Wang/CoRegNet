# Please check our Recount3 CoRegNet via:
http://CRN.liuzlab.org/

# Co-Regulation Model: CoRegNet package Manual

### Step 1: Load your DEGMatrix: Each row is a gene; Each column is an experiment/contrast
#### You can download *DEGMatrix_Recount3_Subset.csv* file, which contains a subset experiments from Recount3.<br>OR, here, as an example, you can use the following function to quickly create a qualified DEG Matrix
<pre>
remotes::install_github("Jiasheng-Wang/CoRegNet")
library(CoRegNet)
DEGMatrix <- create_DEGMatrix()
</pre>
#### If you are using *DEGMatrix_Recount3_Subset.csv* file and you only take it as an example, it is recommended to run: 
<pre>
DEGMatrix <- DEGMatrix[rowSums(DEGMatrix) > 40,] #This code retains genes identified as DEGs in over 40 experiments,
#typically preserving approximately 3,000 genes to optimize computational efficiency.
#For a more precise analysis, consider setting a threshold lower than 40, thereby including a larger gene set in subsequent analyses.
</pre>
#### We advise users to exclude genes seldom identified as DEGs prior to analysis to enhance computational efficiency. It's suggested to retain fewer than 10,000 genes. However, for a more precise analysis, consider including as many genes as feasible.
<pre>
CoRegMatrix <- DEGMatrix %*% t(DEGMatrix)
CoTimes <- CoRegMatrix[upper.tri(CoRegMatrix,diag = FALSE)]
</pre>
### Step 2: Generate sampling-based co-regulation matrix with fixed rowsum and colsum of binary DEG matrix:
#### *Users must manually create all paths/directories used in the functions. As intermediate data/results are often substantial in size, automatic path creation isn't enabled to prevent inadvertent storage in unintended locations.* <br>dir.create() can be used to create a path/directory
<pre>
RandCoRegMatrix_Generation(Ncores = 10, EachCoreSims = 50, RowSumVector = rowSums(DEGMatrix), ColSumVector = colSums(DEGMatrix),save_path = "~/CoregNet_Test/RBM/")
</pre>
### Step 3: Transform the sampling-based co-regulation matrix into gene by simulation format for p-value calculation:
<pre>
n <- 10 # Set 'n' based on 'trans_summary' results. If file size from 'trans_summary' are too large, increase 'n'.
Ncores <- 2 #Adjust 'Ncores' based on your system's capabilities. Ensure 'Ncores' is no larger than 'n' and that 'Ncores' evenly divides 'n'
trans_summary <- Transform_CoM_sup(RandCoM_directory = "~/CoregNet_Test/RBM/",
                                   Ngenes = ncol(DEGMatrix),
                                   NTransformedMatrix = n)

Transform_CoM(RandCoM_directory = "~/CoregNet_Test/RBM/",
              save_path = "~/CoregNet_Test/TransformedMatrix/", 
              NGenePairs = trans_summary[,2], 
              NTransformedMatrix = n, 
              Ncores = Ncores)
</pre>
### Step 4: Calculate the sampling-based P-values
<pre>
SimPCalculation(transformed_CoRegM_Dir = "~/CoregNet_Test/TransformedMatrix/",
                save_path = "~/CoregNet_Test/SimP/",
                CoTimes = CoTimes,
                trans_summary = trans_summary,
                Ncores = 2)
</pre>
### Step 5: Do beta binomial fitting only on those sampling-based p-value smaller than the cutoff to save time
<pre>
#Decide the p-value cutoff and number of cores to do betabinomial fittings:
P_cutoff <- 0.05 #Gene pairs with a p-value exceeding 0.05 are generally not deemed statistically significant and
                 #pursuing further beta-binomial fitting might not be time-efficient.
                 #Nonetheless, users can choose to proceed if they deem it necessary.

FittingFiles <- 10 #Set 'FittingFiles' based on 'fitting_summary' results. If gene pairs for each file from 'fitting_summary' are too large, increase 'FittingFiles'.
fitting_summary <- BBFitting_sup(SimP_Dir = "~/CoregNet_Test/SimP/",
                                P_cutoff = P_cutoff,
                                FittingFiles = FittingFiles)
</pre>
### Step 6: Prepare files for Betabinomial fitting:
<pre>
FittingDataPreparation(FittingFile_Dir = "~/CoregNet_Test/FittingData/",
                       Transform_CoM_Dir = "~/CoregNet_Test/TransformedMatrix/",
                       CoRegMatrix = CoRegMatrix,
                       SigIdx = fitting_summary$qualified_gene_idx,
                       trans_summary = trans_summary,
                       P_cutoff = P_cutoff,
                       FittingFiles = FittingFiles,
                       Ncores = 10) #Set 'Ncores' based on the capabilities of your computer/server
</pre>
### Step 7: BetabinomialFitting:
<pre>
BetaBinomialDistributionFitting(FittingFile_Dir = "~/CoregNet_Test/FittingData/",
                                Result_Dir = "~/CoregNet_Test/FittingResults/",
                                Ncores = 10) #Set 'Ncores' based on the capabilities of your computer/server
  
GenePair_IdxTable <- GenePairIdx(Idx = fitting_summary$qualified_gene_idx,
                                 CoRegMatrix = CoRegMatrix)
CoRegulationNetwork <- BBResults(Result_Dir = FittingResultPath,
                                 GenePair_IdxTable = GenePair_IdxTable)
</pre>
