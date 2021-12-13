if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("SC3", quietly = TRUE))
  BiocManager::install("SC3")
install.packages("./SC3", repos = NULL, type = "source")
if (!requireNamespace("mclust", quietly = TRUE))
  install.packages("mclust")
if (!requireNamespace("Rtsne", quietly = TRUE))
  install.packages("Rtsne")
if (!requireNamespace("mclust", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("mclust", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
  BiocManager::install("SummarizedExperiment")
if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
  BiocManager::install("SingleCellExperiment")
if (!requireNamespace("stats", quietly = TRUE))
  install.packages("stats")
if (!requireNamespace("Rcpp", quietly = TRUE))
  install.packages("Rcpp")
if (!requireNamespace("cluster", quietly = TRUE))
  install.packages("cluster")
if (!requireNamespace("Rmagic", quietly = TRUE))
  install.packages("Rmagic")

library(mclust)
library(Rtsne)
library(SummarizedExperiment)
library(SC3)
library(SingleCellExperiment)
library(wCorr)
library(rlist)
library(stats)
library(matrixStats)
library(Rcpp)
library(cluster)
library(Rmagic)

sourceCpp("/home/marcinmalec/Desktop/rnaseq_imputation/wCorr_m.cpp")
sourceCpp("/home/marcinmalec/Desktop/rnaseq_imputation/solver.cpp")

#Compute ARI for each possibility
eval_alg <- function(X, X_log, labels, num_clusters,threshold) {
  start_time <- Sys.time()
  xlog_t <- as.matrix(magic(t(X_log), genes="all_genes"))

  end_time <- Sys.time()
  
  print("Imputation finished")

  p <- 30
  
  if (ncol(t(as.matrix(xlog_t))) <= p*2){
    p <- 9
  }
  
  cells <- ncol(X)
  if(cells > 1000){
    print("Reducing rank")
    pca_red <- prcomp(as.matrix(xlog_t), rank. = 500)$x
    tsne_red <- Rtsne(pca_red, perplexity = p, check_duplicates = FALSE)$Y
    restarts <- 50
    
  }
  else{
    pca_red <- prcomp(as.matrix(xlog_t))$x
    tsne_red <- Rtsne(as.matrix(xlog_t), perplexity = p, check_duplicates = FALSE)$Y
    restarts <- 1000
  }
  
  c1 = adjustedRandIndex(kmeans(
    pca_red,
    centers = num_clusters,
    iter.max = 1e+09,
    nstart = restarts
  )$cluster,
  labels)
  
  print("PCA kmeans finished")
  
  
  
  #tsne/kmeans
  c2 = adjustedRandIndex(kmeans(
    tsne_red,
    centers = num_clusters,
    iter.max = 1e+09,
    nstart = restarts
  )$cluster,
  labels)
  print("tsne kmeans finished")
  
  
  temp <-t(X_log)
  temp[temp==0] <- xlog_t[temp==0]
  pca_red <- prcomp(as.matrix(temp))$x
  
  prop_zeros_removed <- 1.00-(sum(xlog_t==0))/sum(X_log==0)
  
  pca_dist <- as.matrix(stats::dist(pca_red, method = "euclidean", diag = TRUE, upper = TRUE))
  
  int_labels <- as.numeric((as.factor(labels)))
  
  silh_pca <- silhouette(int_labels, pca_dist)
  silh_pca_avr <- as.numeric(summary(silh_pca)['avg.width'])
  
  return(c(c1, c2, difftime(end_time, start_time, units="secs") , prop_zeros_removed, silh_pca_avr, threshold))
}


driver <- function(filename, repeats, threshold){
  dataset_names <- list("chen")
  # dataset_names <-list("blakeley", "deng", "pollen","darmanis", "segerstolpe")
  
  
  
  for(i in 1:length(dataset_names)){
    dataset = dataset_names[[i]]
    sce <- readRDS(file = paste("/home/marcinmalec/Desktop/rnaseq_imputation/datasets/", dataset, ".rds", sep=""))
    
    X <- assays(sce)$counts
    X_log <- assays(sce)$logcounts
    
    print(paste(dataset, "Genes(rows):", nrow(X), "Cells(cols):", ncol(X), sep=" "))
    
    # colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "col")
    labels<-if(is.null(colData(sce)$cell_type2)) colData(sce)$cell_type1 else colData(sce)$cell_type2
    row_sums <- rowSums(X[,-1])
    X <- X[row_sums>0,] # remove genes that are not expressed at all
    X_log <- X_log[row_sums>0,] # remove genes that are not expressed at all
    
    num_clusters = length(unique(labels))
    
    data_aris <- replicate(repeats, eval_alg(X, X_log, labels, num_clusters,threshold))
    
    means <- rowMeans(data_aris)
    stdevs <- rowSds(data_aris)
    
    print(c("Clustering results: ", dataset))
    print(means)
    print(stdevs)
    fileConn<-eval(parse(text=paste('file("/home/marcinmalec/results/', "magic_", dataset, '_', '_', repeats, '_', threshold, 'txt")', sep="")))
    writeLines(c(paste(dataset, "Genes(rows):", nrow(X), "Cells(cols):", ncol(X), "clusters: ", num_clusters, sep=" "), means, stdevs), fileConn)
    close(fileConn)
  }
  print(sum)
}

# driver("slow-65", 1, .50)
# driver("slow-65", 1, .55)
# driver("slow-65", 1, .60)

driver("slow-65", 1, .65)



# driver("slow-65", 1, .70)
# driver("fast-95", 10, 0, 2, 0, 0)

# driver("em_star", 1, 0, 2, 0, 0)


