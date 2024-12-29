# create cisTopic object
library(cisTopic)
cisTopicObject<-createcisTopicObject(
  atac,
  project.name = "cisTopicProject",
  min.cells = 1,
  min.regions = 1,
  is.acc = 1,
  keepCountsMatrix = TRUE
)
cisTopicObject <- addCellMetadata(cisTopicObject, metadata)
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2, 5, 10:25, 30, 35, 40), seed=987, nCores=1, iterations = 500, addModels=FALSE)

par(mfrow=c(3,3))
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- selectModel(cisTopicObject, type='derivative')

cellassign <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
save(cellassign,file="cellassign.rda")