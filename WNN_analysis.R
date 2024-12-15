# Load the 10x hdf5 file 
inputdata.10x <- Read10X_h5("../filtered_feature_bc_matrix.h5")

# Extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
obj <- CreateSeuratObject(counts = rna_counts)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Add in the ATAC-seq data
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "../fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
obj[["ATAC"]] <- chrom_assay

# RNA analysis
DefaultAssay(obj) <- "RNA"
obj <- SCTransform(obj, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
# Apply Harmony for batch correction on PCA embeddings
obj <- RunHarmony(obj, group.by.vars = "Group", reduction = "pca", dims.use = 1:50)

# ATAC analysis
# The first dimension is excluded as this is typically correlated with sequencing depth
DefaultAssay(obj) <- "ATAC"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj)
obj <- RunUMAP(obj, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities
obj <- FindMultiModalNeighbors(obj, reduction.list = list("harmony", "lsi"), dims.list = list(1:50, 2:50))

# Perform hierarchical clustering
obj <- FindClusters(obj, resolution = 0.5, algorithm = 3, method = "ward.D2")

cluster_SeuratV4 <- obj$seurat_clusters
