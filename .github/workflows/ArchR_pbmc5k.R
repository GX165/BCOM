suppressPackageStartupMessages(library(ArchR))
addArchRGenome("hg19")
addArchRThreads(1)
#Get Input Fragment Files
inputFiles <- getInputFiles("PBMC5K")[1]
names(inputFiles) <- "PBMC5K"

#Create Arrow Files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,                             
  sampleNames = names(inputFiles),                                     
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  minTSS = 0,
  validBarcodes = validBC
)
#ArchRProject
proj <- ArchRProject(ArrowFiles)

########error#############
seRNA <- import10xFeatureMatrix(
  input = c("PBMC5k/atac_v1_pbmc_5k_filtered_peak_bc_matrix.h5"),
  names = c("PBMC5k")
)
#Add scRNA
proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)
#Filter Cells
proj <- proj[proj$TSSEnrichment > 6 & proj$nFrags > 2500 & !is.na(proj$Gex_nUMI)]


#Doublet Filtration. 
proj <- addDoubletScores(proj)
proj <- filterDoublets(proj)

#LSI-ATAC
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC"
)

#LSI-RNA
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)

#Combined Dims
proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
