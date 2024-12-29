# calculate auc
pathway_scoring <- function(gSet, mat_gene){
  cells_rankings <- AUCell_buildRankings(mat_gene, nCores=1, plotStats=TRUE)
  cells_AUC <- AUCell_calcAUC(gSet, cells_rankings, nCores=1,aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))
  aucMatrix <- getAUC(cells_AUC)
  aucMatrix = aucMatrix[rowSums(aucMatrix)>0.0,]
  return(aucMatrix)
}

# pathway selection
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

# SNF fusion
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