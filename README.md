# NetImpute

Introduction
----------
NetImpute is an approach towards the identification of cell types from scRNA-seq data by integrating multiple types of biological networks. It includes two main models: the imputation model and the integration model. 

- **NetImpute-imputation model**: Use a Chebyshev inequality based statistic method to detect the noise data items in scRNA-seq data, and take into account the gene associations from the PPI network and gene pathways to impute these data noise. 
- **NetImpute-integration model**: Fuse the similarity information across cells from both the PPI-based and pathway-based imputation data to identify cell types using the hierarchical clustering algorithm.

Installation
----------
You can install the latest version of NetImpute from the GitHub repository:

	install.packages("devtools")
	devtools::install_github("yiangcs001/NetImpute")

Example
----------
	library(NetImpute)
	res <- netimpute(data = camp, fcm_cluster_num = 5, clust_dist_cutoff = 0.3)
If you do not input the PPI or pathway information, NetImpute will use the default PPI and pathway as biological networks. Otherwise, you can input your own specific PPI and pathway information as follows:

	res <- netimpute(data = camp, fcm_cluster_num = 5, ppi = ppi, pahtway = pathway, clust_dist_cutoff = 0.3)

Copyright and License Information
----------
Copyright (C) 2019 Northwestern Polytechnical University, Xi’an, China.
