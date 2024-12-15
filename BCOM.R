# Initialize the cisTopic object from the CellRanger ATAC count matrix
cisTopicObject<-createcisTopicObject(
  atac,
  project.name = "cisTopicProject",
  min.cells = 1,
  min.regions = 1,
  is.acc = 1,
  keepCountsMatrix = TRUE
)
# Build several models using LDA on the binary accessibility matrix 
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2,5, 10:25, 30,40), seed=987, nCores=20, iterations = 100, addModels=FALSE)
# Model selection: The second derivative measures the changed in the curvature from point to point.(i.e. the highest the second derivative means that the next model is not improving much more the log-likelihood) 
cisTopicObject <- selectModel(cisTopicObject, type='derivative')
# Retrieve the normalised topic assignments to the cells
cellassign <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
# Get the likelihood of each region in each cell
pred.matrix <- predictiveDistribution(cisTopicObject)
# Annotate the genomic regions
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
cisTopicObject <- annotateRegions(cisTopicObject, txdb=txdb, annoDb='org.Mm.eg.db')
# Get highly variable genes
top_genes <- getTopGenes(cisTopicObject, pred.matrix)
getTopGenes <- function(cisTopicObject, pred.matrix) {
  # Prepare the gene annotation dataframe with SYMBOL information
  name <- data.frame(
    RegionName = cisTopicObject@region.names,
    Symbol = cisTopicObject@region.data[["SYMBOL"]]
  )

  # Prepare the list of topics from the cisTopicObject
  topic_list <- lapply(1:length(cisTopicObject@binarized.cisTopics), function(i) {
    cisTopicObject@binarized.cisTopics[[paste0("Topic", i)]]
  })

  # Annotate regions with gene symbols
  getSymbols <- function(topic_list, name) {
    for (i in 1:length(topic_list)) {
      df <- topic_list[[i]]
      df$RegionName <- rownames(df)
      merged_df <- merge(df, name, by = "RegionName", all.x = TRUE)
      topic_list[[i]] <- merged_df
    }
    
    # Extract all gene symbols and filter those that are repeated across topics
    all_symbols <- unlist(lapply(topic_list, function(df) df$Symbol))
    symbol_counts <- table(all_symbols)
    repeated_symbols <- names(symbol_counts[symbol_counts = length(topic_list)])
    return(repeated_symbols)
  }

  # Get the common genes appearing in most topics
  gene <- getSymbols(topic_list, name)

  # Subset the matrix to include only the identified genes
  save <- subset(name, name$Symbol %in% gene)
  pred <- subset(pred.matrix, rownames(pred.matrix) %in% save$RegionName)
  region_to_symbol <- setNames(save$Symbol, save$RegionName)
  rownames(pred) <- region_to_symbol[rownames(pred)]
  
  # Aggregate the expression values by gene
  pred_avg <- aggregate(pred, by = list(rownames(pred)), FUN = mean)
  rownames(pred_avg) <- pred_avg$Group.1
  pred_avg <- pred_avg[, -1]
  matrix <- as.matrix(pred_avg)

  # Compute cell rankings
  rank <- AUCell_buildRankings(pred, plot=FALSE, verbose=FALSE)

  # Calculate gene variance (SD) across samples and select top 1000 genes
  gene_variance <- apply(rank, 1, sd)
  top_1000_genes <- order(gene_variance, decreasing = TRUE)[1:1000]
  top_genes <- rownames(rank)[top_1000_genes]

  # Return the top 1000 genes
  return(top_genes)
}
# Filter the raw gene expression matrix
filtered_count<-count[rownames(count) %in% top_gene,]
# Perform BAMMSC clustering
result<-BAMMSC(as.matrix(filtered_count),K=length(unique(label)))
cluster_bcom<-unlist(result$mem)
names(cluster_bcom)<-colnames(filtered_count)


