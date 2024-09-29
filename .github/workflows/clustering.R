library(SNFtool) 
library(aricode)
library(mclust)
library(wordspace)
library(fastcluster)
library(BAMMSC) 
library(GSEABase) 
library(AUCell)
library(SingleCellExperiment)
library(corrplot)
library(Matrix)
library(data.table)
library(gridExtra)
library(Rtsne)
library(SCENIC)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(Seurat) 


# calculate auc
pathway_scoring <- function(gSet, mat_gene){
  cells_rankings <- AUCell_buildRankings(mat_gene, nCores=1, plotStats=TRUE)
  cells_AUC <- AUCell_calcAUC(gSet, cells_rankings, nCores=1,aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))
  aucMatrix <- getAUC(cells_AUC)
  aucMatrix = aucMatrix[rowSums(aucMatrix)>0.0,]
  return(aucMatrix)
}

#pathway selection
jaccard <- function(gSet, mat_gene) {
  pathway <- names(gSet)
  gene <- rownames(mat_gene)
  similarity <- data.frame(Jaccard = numeric(length(pathway)), pathway = pathway)
  for (i in 1:length(gSet)) {
    gset_i <- geneIds(gSet[[i]])
    intersection_size <- length(intersect(gset_i, gene))
    union_size <- length(union(gset_i, gene))
    jaccard_sim <- intersection_size / union_size
    similarity[i, "Jaccard"] <- jaccard_sim
  }
  similarity <- similarity[order(similarity$Jaccard, decreasing = TRUE), ]
  similarity <- similarity[1:200, , drop = FALSE]
  top_paths <- subset(gSet,names(gSet) %in% similarity$pathway)
  return(top_paths)
}

filter <- function(matrix, top_n = 100) {
  row_sums <- rowSums(matrix)
  sorted_matrix <- matrix[order(row_sums, decreasing = TRUE),]
  if (nrow(sorted_matrix) > top_n) {
    sorted_matrix <- sorted_matrix[1:top_n,]
  }
  return(sorted_matrix)
}

clean_sets <- function(gSet){
  min.size = 10; max.size = 500
  len_s = sapply(gSet, function(x) length(geneIds(x))) 
  idx = (len_s > min.size)&(len_s<max.size)
  gSet = gSet[idx]
  return(gSet)
}

#SNF fusion
integrating_pathway <- function(mat_gene, mat_path){
  K = 10; # number of neighbors, usually (10~30)
  alpha = 0.5; # hyperparameter, usually (0.3~0.8)
  T = 20; # Number of Iterations, usually (10~20)
  mat_gene = t(mat_gene)
  mat_gene = standardNormalization(mat_gene)
  mat_gene = (dist2(as.matrix(mat_gene),as.matrix(mat_gene)))^(1/2)
  mat_gene = affinityMatrix(mat_gene, K, alpha)
  mat_path = t(mat_path)
  mat_path = standardNormalization(mat_path)
  mat_path = (dist2(as.matrix(mat_path),as.matrix(mat_path)))^(1/2)
  mat_path = affinityMatrix(mat_path, K, alpha)
  W = SNF(list(mat_path, mat_gene), K, T)
  return(W)
}

# Hierarchical clustering
clustering<- function(mat_gene,k){
  dis_gene = dist.matrix(t(mat_gene),method = "euclidean",as.dist = TRUE)
  gcl <- fastcluster::hclust(dis_gene, method = 'ward.D');
  clust_results <- cutree(gcl, k)
  return(clust_results)
}

# ARI&NMI
evaluation <- function(result,label_int){
  ari <- ARI(result, label_int)
  nmi <- NMI(result, label_int)
  evaluation <- c(ARI = ari, NMI = nmi)
  return(evaluation)
}

#evaluation
evaluate_clustering <- function(count, norm, harmony, caa, mnn, gSet, label) {
  
  # Gene
  cluster_g <- clustering(count, k = length(unique(label)))
  evaluation_g <- evaluation(cluster_g, label)
  
  # Pathway
  KEGG <- getGmt("/home/223050029/BAMMC/KEGG_human.gmt")
  gset <- subsetGeneSets(KEGG, rownames(norm))
  gset <- jaccard(gset, norm)
  gset <- clean_sets(gset)
  mat_path <- pathway_scoring(gset, norm)
  W <- integrating_pathway(norm, mat_path)
  colnames(W) <- colnames(count)
  cluster_p <- clustering(W, k = length(unique(label)))
  evaluation_path <- evaluation(cluster_p, label)
  
  # AUC
  mat_auc <- pathway_scoring(gSet, norm)
  cluster_auc <- clustering(mat_auc, k = length(unique(label)))
  evaluation_auc <- evaluation(cluster_auc, label)
  
  # Harmony-gene
  cluster_h <- clustering(t(harmony), k = length(unique(label)))
  evaluation_h <- evaluation(cluster_h, label)
  
  # Harmony-pathway
  w <- integrating_pathway(t(harmony), mat_path)
  colnames(w) <- colnames(t(harmony))
  cluster_hp <- clustering(w, k = length(unique(label)))
  evaluation_hp <- evaluation(cluster_hp, label)
  
  # Harmony-AUC
  w1 <- integrating_pathway(t(harmony), mat_auc)
  colnames(w1) <- colnames(t(harmony))
  cluster_ha <- clustering(w1, k = length(unique(label)))
  evaluation_ha <- evaluation(cluster_ha, label)

  #CCA
  cluster_c<-clustering(t(cca),k=length(unique(label)))
  evaluation_c<-evaluation(cluster_c,label)

  #CCA-path
  w2<-integrating_pathway(t(cca),mat_path)
  colnames(w2) <- colnames(t(cca))
  cluster_cp<-clustering(w2,k=length(unique(label)))
  evaluation_cp<-evaluation(cluster_cp,label)

  #CCA_AUC
  w3<-integrating_pathway(t(cca),mat_auc)
  colnames(w3) <- colnames(t(cca))
  cluster_ca<-clustering(w3,k=length(unique(label)))
  evaluation_ca<-evaluation(cluster_ca,label)

  #MNN
  cluster_m<-clustering(t(mnn),k=length(unique(label)))
  evaluation_m<-evaluation(cluster_m,label)

  #MNN-path
  w4<-integrating_pathway(t(mnn),mat_path)
  colnames(w4) <- colnames(t(mnn))
  cluster_mp<-clustering(w4,k=length(unique(label)))
  evaluation_mp<-evaluation(cluster_mp,label)

  #MNN_AUC
  w5<-integrating_pathway(t(mnn),mat_auc)
  colnames(w5) <- colnames(t(mnn))
  cluster_ma<-clustering(w5,k=length(unique(label)))
  evaluation_ma<-evaluation(cluster_ma,label)
  
  # BAMMSC-gene 
  result1 <- DIMMSC(as.matrix(count), K = length(unique(label)))
  clust1 <- unlist(result1$mem)
  names(clust1) <- colnames(count)
  evaluation_bg <- evaluation(as.factor(clust1), label)
  
  # BAMMSC-pathway 
  mat_count <- round(mat_pathapply(count, 2, max))
  result2 <- DIMMSC(mat_count, K = length(unique(label)))
  clust2 <- unlist(result2$mem)
  names(clust2) <- colnames(count)
  evaluation_bp <- evaluation(as.factor(clust2), label)
  
  # BAMMSC-AUC 
  mat_count1 <- round(mat_aucapply(count, 2, max))
  result3 <- DIMMSC(mat_count1, K = length(unique(label)))
  clust3 <- unlist(result3$mem)
  names(clust3) <- colnames(count)
  evaluation_ba <- evaluation(as.factor(clust3), label)
  
  return(list(
    gene = evaluation_g,
    path = evaluation_path,
    auc = evaluation_auc,
    harmony = evaluation_h,
    harmony_path = evaluation_hp,
    harmony_auc = evaluation_ha,
    caa = evaluation_c,
    caa_path = evaluation_cp,
    caa_auc = evaluation_ca,
    mnn = evaluation_m,
    mnn_path = evaluation_mp,
    mnn_auc = evaluation_ma,
    bammsc_gene = evaluation_bg,
    bammsc_path = evaluation_bp,
    bammsc_auc = evaluation_ba
  ))
}

results <- evaluate_clustering(count, norm, harmony, caa, mnn, gSet, label)


#data preparation

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
sample1 <- process_data("E:\\bioinformatics\\MBI6013\\reference_0117\\GSE128066_RAW\\GSM3662224_K5_472\\matrix.mtx",
                        "E:\\bioinformatics\\MBI6013\\reference_0117\\GSE128066_RAW\\GSM3662224_K5_472\\genes.tsv",
                        "E:\\bioinformatics\\MBI6013\\reference_0117\\GSE128066_RAW\\GSM3662224_K5_472\\barcodes.tsv")
sample2 <- process_data("E:\\bioinformatics\\MBI6013\\reference_0117\\GSE128066_RAW\\GSM3662225_K6_473\\matrix.mtx",
                        "E:\\bioinformatics\\MBI6013\\reference_0117\\GSE128066_RAW\\GSM3662225_K6_473\\genes.tsv",
                        "E:\\bioinformatics\\MBI6013\\reference_0117\\GSE128066_RAW\\GSM3662225_K6_473\\barcodes.tsv")
sample3 <- process_data("E:\\bioinformatics\\MBI6013\\reference_0117\\GSE128066_RAW\\GSM3662226_K7_478\\matrix.mtx",
                        "E:\\bioinformatics\\MBI6013\\reference_0117\\GSE128066_RAW\\GSM3662226_K7_478\\genes.tsv",
                        "E:\\bioinformatics\\MBI6013\\reference_0117\\GSE128066_RAW\\GSM3662226_K7_478\\barcodes.tsv")
sample4 <- process_data("E:\\bioinformatics\\MBI6013\\reference_0117\\GSE128066_RAW\\GSM3662227_K8_477\\matrix.mtx",
                        "E:\\bioinformatics\\MBI6013\\reference_0117\\GSE128066_RAW\\GSM3662227_K8_477\\genes.tsv",
                        "E:\\bioinformatics\\MBI6013\\reference_0117\\GSE128066_RAW\\GSM3662227_K8_477\\barcodes.tsv")
lung<-cbind(sample1,sample2,sample3,sample4)

# create seurat object
assay <- CreateAssayObject(lung[,names(label_factor)],analysis = "RNA")
lung<- CreateSeuratObject(assay)

#meta.data
cellsTsne<-tsne_result[ApproxTruth[[1]],]
celltype<-ApproxTruth[[2]]
names(celltype)<-cellsTsne$Barcode
label<-celltype[names(celltype) %in% colnames(lung)]
label_factor <- factor(label)
names(label_factor) <- names(label)
vector1 <- rep(1, 577)
names(vector1)<-colnames(sample1)
vector2 <- rep(2, 724)
names(vector2)<-colnames(sample2)
vector3 <- rep(3, 684)
names(vector3)<-colnames(sample3)
vector4 <- rep(4, 649)
names(vector4)<-colnames(sample4)
result_vector <- c(vector1, vector2,vector3, vector4)
batch_factor <- factor(result_vector)
names(batch_factor) <- names(result_vector)
batch_info <- batch_factor[names(label_factor)]
lung@meta.data$batch <- batch_info
lung@meta.data$celltype <- label_factor

#seurat 
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

#CCA
lung.list <- SplitObject(lung, split.by = "batch")
lung.anchors <- FindIntegrationAnchors(object.list = lung.list, dims = 1:30)
lung.integrated <- IntegrateData(anchorset = lung.anchors, dims = 1:30)
DefaultAssay(lung.integrated) <- "integrated"
lung.integrated <- ScaleData(lung.integrated, verbose = FALSE)
lung.integrated <- RunPCA(lung.integrated, npcs = 30, verbose = FALSE)

#MNN
lung_mnn<- RunFastMNN(object.list = lung.list)

#input data
label<-lung@meta.data$celltype
norm<-lung@assays$RNA@data
count<-lung@assays$RNA@counts
harmony<- Embeddings(lung, 'harmony')
cca<-lung.integrated@reductions[["pca"]]@cell.embeddings
mnn<-lung@reductions[["mnn"]]@cell.embeddings

#regulon identification
motifAnnotations_mgi <- motifAnnotations
scenicOptions <- initializeScenic(org = "mgi", 
                                  dbDir = "/home/223050029/BAMMC/cisTarget_databases", 
                                  nCores = 20)
exprMat<-as.matrix(count)
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] 
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,coexMethod=c("top5perTarget")) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 

