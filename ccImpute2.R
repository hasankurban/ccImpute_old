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
eval_alg <- function(X, X_log, labels, num_clusters,threshold) {
  start_time <- Sys.time()
  distances <- list()
  names <- c()
  names <- c(names, "Spearman")
  distances <-list(w_cor_dist_m(X_log, rowVars(X_log)))
  names(distances) <- names

  sce <- SingleCellExperiment(
    assays = list(
      counts = as.matrix(X),
      logcounts = X_log
    ),
    colData = labels
  )
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sc3_prepare(sce, gene_filter = FALSE)

  metadata(sce)$sc3$distances <- distances
  sce <- sc3_calc_transfs(sce)
  sce <- sc3_kmeans(sce, num_clusters, FALSE)
  sce <- sc3_calc_consens(sce)

  # Get consensus matrix from the SC3
  cm <- eval(parse(text=paste("metadata(sce)$sc3$consensus$'", toString(num_clusters), "'$consensus", sep="")))

  # Remove diagonal entries
  cm <- cm - diag(nrow(cm))

  cm[cm < threshold] <- 0
  cm2 <- t(apply((cm), 2,  # Normalize the entries to get weighted average
                 function(i) i/sum(i)))

  # Replace NA values with 0
  cm2[is.na(cm2)] <- 0

  xlog_t = t(X_log)
  x_imp <- xlog_t
  # print(c(nrow(x_imp),ncol(x_imp)))

  t2 = x_imp == 0
  x_t_vote <- x_imp

  x_t_vote[t2] <- -1
  x_t_vote[x_imp > 0] <- 1

  #compute votes - majority wins - negative means an actual 0, otherwise it is some positive value
  votes <- matrix(0L, nrow = nrow(x_imp), ncol = ncol(x_imp))
  votes[t2] <- (cm2 %*% x_t_vote)[t2]

  t3 <- votes > 0
  x_imp <- solve_dropouts(cm2, x_imp, which(t3, arr.ind = TRUE))

  end_time <- Sys.time()

  print("Imputation finished")
  xlog_t[t3] <- x_imp[t3]

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


  c0 = adjustedRandIndex(eval(parse(text=paste("colData(sce)$sc3", toString(num_clusters), "clusters", sep="_"))),
                         labels)


  prop_zeros_removed <- 1.00-(sum(xlog_t==0))/sum(X_log==0)

  pca_dist <- as.matrix(stats::dist(pca_red, method = "euclidean", diag = TRUE, upper = TRUE))

  int_labels <- as.numeric((as.factor(labels)))

  silh_pca <- silhouette(int_labels, pca_dist)
  silh_pca_avr <- as.numeric(summary(silh_pca)['avg.width'])

  return(c(c0, c1, c2, difftime(end_time, start_time, units="secs") , prop_zeros_removed, silh_pca_avr, threshold))
}


driver <- function(filename, repeats, threshold){
  dataset_names <- list("chen")
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
    
    data_aris <- replicate(repeats, eval_alg(X, X_log, labels, num_clusters,threshold))
    
    means <- rowMeans(data_aris)
    stdevs <- rowSds(data_aris)
    
    print(c("Clustering results: ", dataset))
    print(means)
    print(stdevs)
    fileConn<-eval(parse(text=paste('file("~/ccImpute/results/', "ccimpute_", dataset, '_', filename, '_', repeats, '_', threshold, 'txt")', sep="")))
    writeLines(c(paste(dataset, "Genes(rows):", nrow(X), "Cells(cols):", ncol(X), "clusters: ", num_clusters, sep=" "), means, stdevs), fileConn)
    close(fileConn)
  }
  print(sum)
}

# driver("slow-65", 1, .50)
# driver("slow-65", 1, .55)
# driver("slow-65", 1, .60)
#driver("fast-95", 1, .95)
driver("slow-65", 1, .65)
