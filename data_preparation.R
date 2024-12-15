# call consensus peak from a set of fragments files
get.peaks <- function(fragpaths, cutoff = 1000, annotation = NULL) {
  # Initialize an empty list to store peak data for each fragment path
  peaks_list <- list()
  
  # Loop through each fragment path in the provided list
  for (fragpath in fragpaths) {
    
    # Count the fragments in the current fragment path
    total_counts <- CountFragments(fragpath)
    
    # Filter out the barcodes where the frequency count is below the cutoff threshold
    barcodes <- total_counts[total_counts$frequency_count > cutoff, ]$CB
    
    # Create a fragment object for the selected barcodes from the current fragment path
    frags <- CreateFragmentObject(path = fragpath, cells = barcodes)
    
    # Call peaks from the fragment object
    peaks <- CallPeaks(frags)
    
    # Create a feature matrix
    counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)
    
    # Create a chromatin assay
    assay <- CreateChromatinAssay(
      counts,
      sep = c(":", "-"),  # Separator for fragment identifiers
      fragments = fragpath,  # Path to the fragment file
      annotation = annotation  # Optional annotation for the assay
    )
    
    # Extract the genomic ranges for the identified peaks from the chromatin assay
    peaks_granges <- granges(assay)
    
    # Store the peaks as genomic ranges in the peaks_list for the current fragment path
    peaks_list[[fragpath]] <- peaks_granges
  }
  
  # Combine all the peaks from different fragment paths into a single list
  combined.peaks <- Reduce(function(x, y) {
    c(x, y)  # Concatenate two genomic range objects
  }, peaks_list)
  
  # Calculate the width (size) of each peak in the combined list
  peakwidths <- width(combined.peaks)
  
  # Keep only those with widths between 20 and 10,000 base pairs
  combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]
  
  # Return the filtered list of combined peaks
  return(combined.peaks)
}


# Combine multiple scATAC-seq datasets into a single Seurat object
create.combined <- function(counts_paths, frag_paths, combined.peaks) {
  # Initialize a list of Seurat objects for each dataset
  seurat_objects <- lapply(1:length(counts_paths), function(i) {
    
    # Read in the gene expression data from the specified path (counts in 10X format)
    counts <- Read10X_h5(counts_paths[i])$Peaks  
    
    # Create a Fragment object based on the fragment file corresponding to the counts data
    frags <- CreateFragmentObject(
      path = frag_paths[i],  # Path to the fragment file
      cells = colnames(counts)  # Extract the cell names from the counts data
    )  
    
    # Create a feature matrix from the fragment object, using the given peaks
    matrix <- FeatureMatrix(
      fragments = frags,  # Fragment data
      features = combined.peaks,  # The set of combined peaks to extract
      cells = colnames(counts)  # List of cells
    )
    
    # Create a Chromatin Assay from the feature matrix and fragment data
    assay <- CreateChromatinAssay(matrix, fragments = frags)
    
    # Create a Seurat object with the Chromatin Assay
    obj <- CreateSeuratObject(assay, assay = "ATAC")
    
    # Return the Seurat object
    return(obj) 
  })
  
  # Merge all Seurat objects into one combined object
  combined <- merge(
    x = seurat_objects[[1]],
    y = seurat_objects[2:length(seurat_objects)],
    add.cell.ids = paste0("KTB", 1:length(seurat_objects))  # Add unique cell identifiers (prefix "KTB")
  )
  
  # Return the combined Seurat object containing data from all samples
  return(combined)
}

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)# EnsDb.Mmusculus.v79
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

peaks <- get.peaks(fragpath=frags, cutoff = 1000, annotation = annotation)
combined <- create.combined(counts, frags, peaks)

# Peak-cell counts
atac<- combined@assays[["ATAC"]]@counts
rownames(atac)<-sub("-", ":", rownames(atac))

