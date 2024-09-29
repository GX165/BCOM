library(SnapATAC)
library(GenomicRanges)
metadata <- read.table('metadata.tsv',
                       header = TRUE,
                       stringsAsFactors=FALSE,quote="",row.names=1)
metadata$label = as.character(metadata$label)

#Barcode selection
x.sp = createSnap(
  file="10xpbmc5k.snap",
  sample="10xpbmc5k",
  do.par = TRUE,
  num.cores=10)
x.sp = x.sp[which(x.sp@barcode %in% rownames(metadata)),]
#filter cells only using number of fragments and UMI with the following cutoffs
x.sp = filterCells(
 obj=x.sp, 
 subset.names=c("fragment.num", "UMI"),
 low.thresholds=c(1000,1000),
 high.thresholds=c(Inf, Inf)
 )
# show what bin sizes exist in file
showBinSizes("10xpbmc5k.snap");
x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=10)
# Matrix binarization
x.sp = makeBinary(x.sp, mat="bmat")
# Bin filtration
black_list = read.table('wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz')
black_list.gr = GRanges(black_list[,1], IRanges(black_list[,2], black_list[,3]))
idy1 = queryHits(findOverlaps(x.sp@feature, black_list.gr))
idy2 = grep("chrM|random", x.sp@feature)
idy = unique(c(idy1, idy2))
x.sp = x.sp[,-idy, mat="bmat"]
options(repr.plot.width=4, repr.plot.height=4)
plotBinCoverage(
  x.sp,
  pdf.file.name=NULL,
  col="grey",
  border="grey",
  breaks=10,
  xlim=c(-6,6)
)
x.sp = filterBins(
  x.sp,
  low.threshold=-2,
  high.threshold=2,
  mat="bmat"
)
############Jaccard matrix:error##############
x.sp = runJaccard(
  obj = x.sp,
  tmp.folder=tempdir(),
  mat = "bmat",
  max.var=2000,
  ncell.chunk=1000,
  do.par=FALSE,
  num.cores=10,
  seed.use=10
)
#normalization
x.sp = runNormJaccard(
    obj = x.sp,
    tmp.folder=tempdir(),
    ncell.chunk=1000,
    method="normOVE",
    row.center=TRUE,
    row.scale=TRUE,
    low.threshold=-5,
    high.threshold=5,
    do.par=FALSE,
    num.cores=10,
    seed.use=10
    )
#Linear Dimentionality Reduction
x.sp = runDimReduct(
    x.sp,
    pc.num=50,
    input.mat="jmat",
    method="svd",
    center=TRUE,
    scale=FALSE,
    seed.use=10
    )
#Determine statistically significant principal components
options(repr.plot.width=6, repr.plot.height=6)
plotDimReductElbow(
    obj=x.sp, 
    point.size=1.5,
    point.shape=19,
    point.color="red",
    point.alpha=1,
    pdf.file.name=NULL,
    pdf.height=7,
    pdf.width=7,
    labs.title="PCA Elbow plot",
    labs.subtitle=NULL
    )
plotDimReductPW(
    obj=x.sp, 
    pca.dims=1:50,
    point.size=0.3,
    point.color="grey",
    point.shape=19,
    point.alpha=0.6,
    down.sample=5000,
    pdf.file.name=NULL, 
    pdf.height=7, 
    pdf.width=7
    )