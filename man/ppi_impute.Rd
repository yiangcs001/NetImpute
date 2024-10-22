% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppi_impute.R
\name{ppi_impute}
\alias{ppi_impute}
\title{Impute data by PPI information}
\usage{
ppi_impute(data, fcm_cluster_num, filter_thre, pca_dim, ppi, p_low_cut,
  p_high_cut, lambda1, lambda2)
}
\arguments{
\item{data}{The initial single-cell RNA-seq data.}

\item{fcm_cluster_num}{An integer usd in fuzzy c-means, recommend closer to the reue value.}

\item{filter_thre}{Zero number ratio in genes larger than filter_thre are filtered.}

\item{pca_dim}{An integer indicate how many number of dimentions used in fuzzy c-means.}

\item{ppi}{PPI data used in imputation.}

\item{p_low_cut}{Low data imputation ratio, between 0 and 1.}

\item{p_high_cut}{High data imputation ratio, between 0 and 1.}

\item{lambda1}{L1 regularization parameter.}

\item{lambda2}{L2 regularization parameter.}
}
\value{
The imputated data use interactions between genes in PPI network.
}
\description{
Impute data by PPI information
}
