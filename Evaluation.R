evaluation <- function(result,label_int){
  ari <- ARI(result, label_int)
  nmi <- NMI(result, label_int)
  evaluation <- c(ARI = ari, NMI = nmi)
  return(evaluation)
}
