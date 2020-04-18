#' Impute data by pathway information
#'
#' @param data The initial single-cell RNA-seq data.
#' @param fcm_cluster_num An integer usd in fuzzy c-means, recommend closer to the reue value.
#' @param filter_thre Zero number ratio in genes larger than filter_thre are filtered.
#' @param pca_dim An integer indicate how many number of dimentions used in fuzzy c-means.
#' @param pathway Pathway data used in imputation.
#' @param p_low_cut Low data imputation ratio, between 0 and 1.
#' @param p_high_cut High data imputation ratio, between 0 and 1.
#' @param lambda1 L1 regularization parameter.
#' @param lambda2 L2 regularization parameter.
#'
#' @return The imputated data use interactions between genes in pathways.
#' @export
#'
#' @examples
pathway_impute <- function(data, fcm_cluster_num, filter_thre, pca_dim, pathway,
                           p_low_cut, p_high_cut, lambda1, lambda2){
  print("Impute through Pathway ...")
  # Filter data
  gene_fileter <- rowSums(data > 0) > (ncol(data) * filter_thre)
  filter_data <- data[gene_fileter,]
  # Normal the PPI network and data
  pathway_gene <- unique(unlist(pathway))
  pathway_data <- filter_data[rownames(filter_data) %in% pathway_gene, ]
  # Determin the PCA dimension number
  dist_res <- as.matrix(1 - cor(pathway_data, method = "pearson"))
  pca_res <- prcomp(dist_res, center = TRUE, scale. = FALSE)
  if(is.null(pca_dim)){
    pca_dim <- 0
    s_break <- 0
    variance <- (pca_res$sdev ^ 2) / sum((pca_res$sdev) ^ 2)
    signal <- (-diff(variance) / variance[2:length(variance)] > 0.6)
    for(s in signal){
      if(s_break < 2){
        pca_dim = pca_dim + 1
        if(s){
          s_break <- 0
        }else{
          s_break = s_break + 1
        }
      }
    }
    pca_dim <- pca_dim - 2
    if(pca_dim < 3){
      pca_dim <- 3
    }else if(pca_dim > 10){
      pca_dim <- 10
    }
  }

  # Fuzzy C-Means
  pca_tmp <- pca_res$rotation[, 1:pca_dim]
  fcm_res <- ppclust::fcm(pca_tmp, centers = fcm_cluster_num, nstart = 10, iter.max = 1e+4)
  fcm_member_degree <- fcm_res$u
  clust_res <- fcm_member_degree > 0.5
  unidentify_index <- which(rowSums(clust_res) != 1)
  clust_res[unidentify_index, ] <- (fcm_member_degree[unidentify_index, ] > (2 / fcm_cluster_num))
  unidentify_index <- which(rowSums(clust_res) > 1)
  outliers <- (rowSums(clust_res) == 0)

  # Imputation
  impute_data <- pathway_data
  unidentify_data <- vector("list", length(unidentify_index))
  names(unidentify_data) <- unidentify_index
  if(length(which(rowSums(clust_res) < 1)) > 0){
    unidentify_data <- unidentify_data[names(unidentify_data) %in% which(rowSums(clust_res) < 1) == FALSE]
  }
  for(clust_index in 1:fcm_cluster_num){
    print(paste("Imputation for cluster", clust_index, "..."))
    cell_index <- which(clust_res[, clust_index])
    if(length(cell_index) > 5){
      object <- impute_data[, cell_index]
      object <- object[rowSums(object > 0) > (ncol(object) * filter_thre),]
      # Calculate the dropout candidate
      m_data <- apply(object, 1, mean)
      v_data <- apply(object, 1, var)
      s_data <- t(apply(object, 1, function(x) (x - mean(x)) ^ 2))
      p_data <- 1 - (s_data / ((sqrt(2) ^ 2) * v_data))
      dropout_candidate <- ((object < m_data) & (p_data < p_low_cut)) | ((object > m_data) & (p_data < p_high_cut))
      # Imputation data in one cluster by network information
      temp <- object * (!dropout_candidate)
      for(row_index in 1:nrow(temp)){
        temp_yimpute <- list()
        d_flag <- dropout_candidate[row_index, ]
        if(sum(d_flag) > 0){
          g_name <- rownames(temp)[row_index]
          temp_pathway <- pathway[sapply(pathway, function(x) all(g_name %in% x))]
          for(temp_pathway_gene in temp_pathway){
            nb_gene <- setdiff(intersect(temp_pathway_gene, rownames(temp)), g_name)
            if(length(nb_gene) > 1){
              xx <- t(temp[nb_gene, !d_flag])
              row_name <- rownames(xx)[which(rowSums(xx) != 0)]
              if(length(row_name) > 1){
                yy <- temp[row_index, !d_flag]
                ximpute <- t(temp[nb_gene, d_flag])
                nnls_res <- penalized::penalized(response = yy, penalized = xx, unpenalized = ~0,
                                                 positive = TRUE, lambda1 = lambda1, lambda2 = lambda2,
                                                 maxiter = 1e+3, trace = FALSE)
                if(length(penalized::coefficients(nnls_res)) == 0){
                  yimpute <- rep(0, nrow(ximpute))
                  names(yimpute) <- rownames(ximpute)
                }else{
                  yimpute_res <- penalized::predict(nnls_res, penalized = ximpute, unpenalized = ~0)
                  if(class(yimpute_res) == "matrix"){
                    yimpute <- yimpute_res[,1]
                  }else if(nrow(ximpute) > 1){
                    yimpute <- yimpute_res
                  }else{
                    yimpute <- yimpute_res[1]
                  }
                }
                temp_yimpute <- c(temp_yimpute, list(yimpute))
              }

            }
          }
          if(length(temp_yimpute) > 0){
            impute_data[g_name,][match(names(which(d_flag==TRUE)), colnames(impute_data))] <- do.call(pmax, temp_yimpute)
          }
        }
      }
      for(index in as.integer(names(unidentify_data))){
        if(clust_res[index, clust_index]){
          unidentify_data[[which(names(unidentify_data) == index)]][as.character(clust_index)] <- list(impute_data[, index])
        }
      }
    }
  }

  if(length(unidentify_data) == 0){
    pathway_impute_data <- impute_data
  }else{
    pathway_impute_data <- impute_data
    for(temp_index in 1:length(unidentify_data)){
      cell_index <- as.integer(names(unidentify_data[temp_index]))
      temp_data <- data.frame(unidentify_data[[temp_index]])
      pathway_impute_data[, cell_index] <- apply(temp_data, 1, max)
    }
  }
  return(list(raw_data = pathway_data, impute_data = pathway_impute_data,
              pca_dim = pca_dim, outliers = outliers))
}
