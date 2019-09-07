#' Fusion data impute by PPI and Pathway information, in order to identify cell types.
#'
#' @param ppi_data The data impute by PPI information.
#' @param pathway_data The data impute by pathwy information.
#' @param clust_dist_cutoff A numeric or a vector indicate the cutting parameter in hierarchical clustering.
#'
#' @return A list contaion predict labels in each height cut.
#' @export
#'
#' @examples
fusion_clustering <- function(ppi_data, pathway_data, clust_dist_cutoff){
  ppi_diss <- as.matrix(1 - cor(ppi_data, method = "pearson"))
  pathway_diss <- as.matrix(1 - cor(pathway_data, method = "pearson"))
  cons <- rbind(pathway_diss, ppi_diss)
  diss_1 <- as.matrix(1 - cor(cons, method = "pearson"))
  diss_2 <- as.matrix(1 - cor(cons, method = "spearman"))
  diss <- (diss_1 + diss_2) / 2
  fit <- stats::hclust(as.dist(diss))

  cluster_res <- list()
  if(is.null(clust_dist_cutoff)){
    clust_dist_cutoff <- seq(0, 1, 0.1)
  }
  for(i in 1:length(clust_dist_cutoff)){
    predict_lab <- cutree(fit, h = clust_dist_cutoff[i])
    cluster_res[[i]] <- predict_lab
  }
  names(cluster_res) <- clust_dist_cutoff
  return(cluster_res = cluster_res)
}
