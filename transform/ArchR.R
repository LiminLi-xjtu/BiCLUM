
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

library(ArchR)
library(stringr)
library(Seurat)
library(SeuratObject)
library(sceasy)

inputFiles <- c("pbmc" = "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")

addArchRGenome("hg38") 
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  # filterTSS = 4, 
  # filterFrags = 1000,
  filterFrags = NULL,
  filterTSS = NULL,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

projHeme <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
projHeme
getAvailableMatrices(projHeme)

GeneScoreMatrix=getMatrixFromProject(
  ArchRProj = projHeme,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

genescore_name = GeneScoreMatrix@elementMetadata@listData[["name"]]
genescore_ex = GeneScoreMatrix@assays@data@listData[["GeneScoreMatrix"]]
rownames(genescore_ex) = genescore_name
coln = colnames(genescore_ex)
str_cell = str_split(coln,'#',simplify = T)
colns = str_cell[,2]
colnames(genescore_ex) = colns

scgeneactivity = CreateSeuratObject(counts = genescore_ex)
sceasy::convertFormat(scgeneactivity, from="seurat", to="anndata",
                      outFile='scGAM_ArchR.h5ad')
