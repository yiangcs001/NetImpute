#' Main function of NetImpute
#'
#' @param data The initial single-cell RNA-seq data.
#' @param is_log_data A logical value indicate whether the data is logarithm-form data.
#' @param fcm_cluster_num An integer usd in fuzzy c-means, recommend closer to the reue value.
#' @param filter_thre Zero number ratio in genes larger than filter_thre are filtered.
#' @param pca_dim An integer indicate how many number of dimentions used in fuzzy c-means.
#' @param ppi_impute A logical value indicate whether do imputation by PPI.
#' @param ppi PPI data used in imputation.
#' @param pathway_impute A logical value indicate whether do imputation by PPI.
#' @param pathway Pathway data used in imputation.
#' @param p_low_cut Low data imputated ratio, between 0 and 1.
#' @param p_low_cut_theta A parameter control low data imputated ratio, between 0 and 1.
#' @param p_high_cut High data imputated ratio, between 0 and 1.
#' @param p_high_cut_theta A parameter control high data imputated ratio, between 0 and 1.
#' @param lambda A parameter used in elastic network regularization.
#' @param alpha A parameter used in elastic network regularization.
#' @param do_fusion_clustering A logical value indicate whether to do the integrated step to identify cell types.
#' @param clust_dist_cutoff A numeric or a vector indicate the cutting parameter in hierarchical clustering.
#'
#' @return NetImpute returns imputated data use interactions between genes in PPI network/pathways, respectively.
#' @export
#'
#' @examples
netimpute <-function(data, is_log_data = TRUE, fcm_cluster_num, filter_thre = 0.1, pca_dim = NULL,
                     ppi_impute = TRUE, ppi = NULL, pathway_impute = TRUE, pathway = NULL,
                     p_low_cut = NULL, p_low_cut_theta = 0, p_high_cut = NULL, p_high_cut_theta = 0.5,
                     lambda = 10, alpha = 0.5, do_fusion_clustering = TRUE, clust_dist_cutoff = NULL){
  # Initialization
  # Log transformation data
  if(!is.logical(is_log_data)){
    stop("'is_log_data' must be a logical value!")
  }else{
    if(!is_log_data){
      data <-log2(data + 1)
    }
  }
  # Low data p-value cut and high data p-value
  if(is.null(p_low_cut)){
    p_low_cut <- ((1 - 1 / (sqrt(2) - p_low_cut_theta) ^ 2) + 0.5) / 2
  }else if(p_low_cut > 1 | p_low_cut < 0){
    stop("'p_low_cut' must be a numeric between 0 and 1!")
  }
  if(is.null(p_high_cut)){
    p_high_cut <- ((1 - 1 / (sqrt(2) - p_high_cut_theta) ^ 2) + 0.5) / 2
  }else if(p_high_cut > 1 | p_high_cut < 0){
    stop("'p_high_cut' must be a numeric between 0 and 1!")
  }
  # Elastic network parameter lambda1 for L1 and lambda2 for L2
  lambda1 <- lambda * alpha
  lambda2 <- lambda * (1 - alpha) / 2
  # Load PPI / Pathway data
  if(ppi_impute & is.null(ppi)){
    ppi <- NetImpute::ppi
  }
  if(pathway_impute & is.null(pathway)){
    pathway <- NetImpute::pathway
  }

  # Impute
  if(ppi_impute){
    ppi_impute_res <- ppi_impute(data, fcm_cluster_num, filter_thre, pca_dim, ppi,
                                 p_low_cut, p_high_cut, lambda1, lambda2)
  }else{
    ppi_impute_res <- NULL
  }
  if(pathway_impute){
    pathway_impute_res <- pathway_impute(data, fcm_cluster_num, filter_thre, pca_dim, pathway,
                                         p_low_cut, p_high_cut, lambda1, lambda2)
  }else{
    pathway_impute_res <- NULL
  }

  # Fusion clustering
  if(do_fusion_clustering){
    cluster_res <- fusion_clustering(ppi_impute_res$impute_data, pathway_impute_res$impute_data, clust_dist_cutoff)
    }else{
    cluster_res <- NULL
  }

  return(list(ppi = ppi_impute_res, pathway = pathway_impute_res, cluster_res = cluster_res))
}

