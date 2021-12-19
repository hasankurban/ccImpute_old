library(mclust)
library(Rtsne)
library(SummarizedExperiment)
library(SC3)
library(SingleCellExperiment)
# library(rlist)
library(stats)
library(Rcpp)
library(cluster)
sourceCpp("~/ccImpute/cpp/wCorr_m.cpp")
sourceCpp("~/ccImpute/cpp/solver.cpp")

#Compute ARI for each possibility
eval_alg <- function(X, X_log, labels, num_clusters,threshold, isFAST,dataset,filename) {

  xlog_t <- read.csv(file=paste("~/ccImpute/datasets/", dataset,"_", filename, ".csv", sep=""))
  rownames(xlog_t) <- xlog_t[,1]
  xlog_t <- xlog_t[,-1]
#  colnames(xlog_t) <- xlog_t[1,]
#  xlog_t <- xlog_t[-1,]
  # print(c(nrow(xlog_t), ncol(xlog_t)))
  # p <- 30
  # 
  # if (ncol(t(as.matrix(xlog_t))) <= p*2){
  #   p <- 9
  # }
  # 
  # cells <- ncol(X)
  # if(cells > 1000){
  #   print("Reducing rank")
  #   pca_red <- prcomp(as.matrix(xlog_t), rank. = 500)$x
  #   tsne_red <- Rtsne(pca_red, perplexity = p, check_duplicates = FALSE)$Y
  #   restarts <- 50
  # 
  # }
  # else{
    # pca_red <- prcomp(as.matrix(xlog_t))$x
  #   tsne_red <- Rtsne(as.matrix(xlog_t), perplexity = p, check_duplicates = FALSE)$Y
  #   restarts <- 1000
  # }
  # 
  # c1 = adjustedRandIndex(kmeans(
  #   pca_red,
  #   centers = num_clusters,
  #   iter.max = 1e+09,
  #   nstart = restarts
  # )$cluster,
  # labels)
  # 
  # print("PCA kmeans finished")
  # 
  # 
  # 
  # #tsne/kmeans
  # c2 = adjustedRandIndex(kmeans(
  #   tsne_red,
  #   centers = num_clusters,
  #   iter.max = 1e+09,
  #   nstart = restarts
  # )$cluster,
  # labels)
  # print("tsne kmeans finished")


  #c0 = adjustedRandIndex(eval(parse(text=paste("colData(sce)$sc3", toString(num_clusters), "clusters", sep="_"))),
#                         labels)


  prop_zeros_removed <- 1.00-(sum(xlog_t==0))/sum(X_log==0)

  pca_red <- prcomp(as.matrix(xlog_t))$x
  pca_dist <- as.matrix(stats::dist(pca_red, method = "euclidean", diag = TRUE, upper = TRUE))

  int_labels <- as.numeric((as.factor(labels)))

  silh_pca <- silhouette(int_labels, pca_dist)
  silh_pca_avr <- as.numeric(summary(silh_pca)['avg.width'])

  return(c(prop_zeros_removed, silh_pca_avr, threshold))
}


driver <- function(filename, repeats, threshold, isFAST,dataset){
  dataset_names <- list(dataset)
  # dataset_names <-list("blakeley", "deng", "pollen","darmanis", "segerstolpe")
  
  
  
  for(i in 1:length(dataset_names)){
    dataset = dataset_names[[i]]
    sce <- readRDS(file = paste("~/ccImpute/datasets/", dataset, ".rds", sep=""))
    
    X <- assays(sce)$counts
    X_log <- assays(sce)$logcounts
    
    print(paste(dataset, "Genes(rows):", nrow(X), "Cells(cols):", ncol(X), sep=" "))
    
    # colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "col")
    labels<-if(is.null(colData(sce)$cell_type2)) colData(sce)$cell_type1 else colData(sce)$cell_type2
    row_sums <- rowSums(X[,-1])
    X <- X[row_sums>0,] # remove genes that are not expressed at all
    X_log <- X_log[row_sums>0,] # remove genes that are not expressed at all
    
    num_clusters = length(unique(labels))
    
    data_aris <- replicate(repeats, eval_alg(X, X_log, labels, num_clusters,threshold, isFAST, dataset, filename))
    
    means <- rowMeans(data_aris)
    stdevs <- rowSds(data_aris)
    
    print(c("Clustering results: ", dataset))
    print(means)
    print(stdevs)
    fileConn<-eval(parse(text=paste('file("~/ccImpute/results/', "ccimpute__2_", dataset, '_', filename, '_', repeats, '_', threshold, 'txt")', sep="")))
    writeLines(c(paste(dataset, "Genes(rows):", nrow(X), "Cells(cols):", ncol(X), "clusters: ", num_clusters, sep=" "), means, stdevs), fileConn)
    close(fileConn)
  }
  print(sum)
}

args = commandArgs(trailingOnly=TRUE)

if(args[1] == "slow"){
  print(c("slow", args[2]))
  driver("slow-65", 1, .65, FALSE, args[2])
} else {driver("fast-95", 1, .95, TRUE, args[2])}
