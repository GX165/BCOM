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
integrating <- function(mat_list) {
  K = 10  
  alpha = 0.5  
  T = 20  
  processed_matrices <- lapply(mat_list, function(mat) {
    mat <- t(mat)
    mat <- standardNormalization(mat)
    mat <- (dist2(as.matrix(mat), as.matrix(mat)))^(1/2) 
    affinityMatrix(mat, K, alpha)
  })
  W <- SNF(processed_matrices, K, T)
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
  gset <- subsetGeneSets(KEGG, rownames(count))
  gset <- jaccard(gset, count)
  gset <- clean_sets(gset)
  mat_path <- pathway_scoring(gset, norm)
  cluster_p <- clustering(mat_path, k = length(unique(label)))
  evaluation_path <- evaluation(cluster_p, label)
  
  # Regulon
  mat_auc <- pathway_scoring(gSet, norm)
  cluster_auc <- clustering(mat_auc, k = length(unique(label)))
  evaluation_auc <- evaluation(cluster_auc, label)

	# Atac
	cluster_atac <- clustering(cellassign, k = length(unique(label)))
  evaluation_atac <- evaluation(cluster_atac, label)

	# Atac-pathway
	list1 <- list(cellassign,mat_path)
  wcp <- integrating(list1)
  colnames(wcp) <- colnames(cellassign)
  cluster_cp <- clustering(wcp, k = length(unique(label)))
  evaluation_cp <- evaluation(cluster_cp, label)
  
  # Atac-regulon
	list2 <- list(cellassign,mat_auc)
  wca <- integrating(list2)
  colnames(wca) <- colnames(cellassign)
  cluster_ca <- clustering(wca, k = length(unique(label)))
  evaluation_ca <- evaluation(cluster_ca, label)

	# Atac-pathway-regulon
	list3 <- list(cellassign,mat_path,mat_auc)
  whcp <- integrating(list3)
  colnames(whcp) <- colnames(cellassign)
  cluster_hcp <- clustering(whcp, k = length(unique(label)))
  evaluation_hcp <- evaluation(cluster_hcp, label)

	#Harmony
	cluster_h<-clustering(t(harmony),k=length(unique(label)))
	evaluation(cluster_h,label)
  
  # Harmony-pathway
	list4 <- list(t(harmony),mat_path)
  whp <- integrating(list4)
  colnames(whp) <- colnames(t(harmony))
  cluster_hp <- clustering(whp, k = length(unique(label)))
  evaluation_hp <- evaluation(cluster_hp, label)
  
  # Harmony-regulon
	list5 <- list(t(harmony), mat_auc)
  wha <- integrating(list5)
  colnames(wha) <- colnames(t(harmony))
  cluster_ha <- clustering(wha, k = length(unique(label)))
  evaluation_ha <- evaluation(cluster_ha, label)

  # Harmony-atac
	list6 <- list(t(harmony), cellassign)
  whc<-integrating(list6)
  colnames(whc) <- colnames(t(harmony))
  cluster_hc<-clustering(whc,k=length(unique(label)))
  evaluation_hc<-evaluation(cluster_hc,label)

  # Harmony-atac-pathway
  list7 <- list(t(harmony), cellassign, mat_path)
  whcp<-integrating(list7)
  colnames(whcp) <- colnames(t(harmony))
  cluster_hcp<-clustering(whcp,k=length(unique(label)))
  evaluation_hcp<-evaluation(cluster_hcp,label)
	
  # Harmony-atac-regulon
  list8 <- list(t(harmony), cellassign, mat_auc)
  whca<-integrating(list8)
  colnames(whca) <- colnames(t(harmony))
  cluster_hca<-clustering(whca,k=length(unique(label)))
  evaluation_hca<-evaluation(cluster_hca,label)

  # Harmony-atac-pathway-regulon
  list9 <- list(t(harmony), cellassign, mat_path, mat_auc)
  whcpa<-integrating(list9)
  colnames(whcpa) <- colnames(t(harmony))
  cluster_hcpa<-clustering(whcpa,k=length(unique(label)))
  evaluation_hcpa<-evaluation(cluster_hcpa,label)
  
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
		atac = evaluation_atac,
		atac_path = evaluation_cp,
		atac_regulon = evaluation_ca,
		atac_path_regulon = evaluation_cpa,
    harmony = evaluation_h,
    harmony_path = evaluation_hp,
    harmony_auc = evaluation_ha,
		harmony_atac = evaluation_hc,
    harmony_atac_auc = evaluation_hca,
		harmony_atac_pathway = evaluation_hcp,
    harmony_atac_auc_pathway = evaluation_hcpa,
    bammsc_gene = evaluation_bg,
    bammsc_path = evaluation_bp,
    bammsc_auc = evaluation_ba
  ))
}

results <- evaluate_clustering(count, harmony, cellassign, gSet, label)

