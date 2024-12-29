#load matrix, gene, barcodes
process_data <- function(matrix_file, genes_file, barcodes_file) {
  mat <- readMM(matrix_file)
  gene_info <- read.delim(genes_file, stringsAsFactors = FALSE, sep = "\t", header = FALSE)
  barcodes <- read.delim(barcodes_file, stringsAsFactors = FALSE, sep = "\t", header = FALSE)
  
  rownames(mat) <- gene_info[, 1]
  colnames(mat) <- barcodes[, 1]
  mat <- as.matrix(mat)
  data_sample_1 <- as.matrix(mat)
  rownames(data_sample_1) <- gene_info[, 2]
  return(data_sample_1)
}
#lung, sample1-4
sample1 <- process_data("..\\matrix.mtx",
                        "..\\genes.tsv",
                        "..\\barcodes.tsv")
sample2 <- process_data("..\\matrix.mtx",
                        "..\\genes.tsv",
                        "..\\barcodes.tsv")
sample3 <- process_data("..\\matrix.mtx",
                        "..\\genes.tsv",
                        "..\\barcodes.tsv")
sample4 <- process_data("..\\matrix.mtx",
                        "..\\genes.tsv",
                        "..\\barcodes.tsv")
lung<-cbind(sample1,sample2,sample3,sample4)

# create seurat object
assay <- CreateAssayObject(lung[,names(label_factor)],analysis = "RNA")
lung<- CreateSeuratObject(assay, meta.data=metadata)

# normalization + PCA 
lung <- NormalizeData(lung, verbose = FALSE)
var_features <- lapply(split(row.names(lung@meta.data), lung@meta.data[["batch"]]), function(cells_use) {
  var_feats <- FindVariableFeatures(lung[, cells_use], selection.method = "vst", nfeatures = 2000)
  VariableFeatures(var_feats)
})
VariableFeatures(lung) <- unique(unlist(var_features))
lung <- ScaleData(lung, verbose = FALSE)
lung <- RunPCA(lung, features = VariableFeatures(lung), npcs = 20, verbose = FALSE)

#Harmony
lung<- RunHarmony(lung, "batch", plot_convergence = TRUE, 
                  nclust = 50, max_iter = 10, early_stop = TRUE)