### SETUP ###

# DATA FRAME MANIPULATION
library(tidyverse)
library(magrittr)
library(matrixStats)
library(plyr)
library(dplyr)
# library(reshape2)
# library(useful) # For 'corner' function

# zmad <- function(x) {
#   med_x <- median(x, na.rm = T)
#   mad_x <- mad(x, constant = 0.6744897501960817, na.rm = T)
#   
#   return((x - med_x)/(mad_x))
# }

# Plotting
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(grid)
library(gridExtra)
# library(ggridges)
library(ggpubr)
# library(GGally)
library(plotly)
# library(VennDiagram)
library(RColorBrewer)
library(extrafont)
# library(cowplot)
# library(networkD3)
suppressMessages(loadfonts())

# For tables
library(data.table)
library(DT)
library(htmltools)

# Reading in data
library(readr)
library(readxl)

# GSEA
# library(GSEABase)
library(clusterProfiler)

# CDS
# library(taigr)
# library(cdsr)
# library(celllinemapr)

# Other
library(pbapply)
library(parallel)
# library(sva) # For ComBat
library(limma)
# library(stringr)
# library(irr)
# library(zoo)
library(scales)
# library(rlang)
# library(rmarkdown)
# library(knitr)
# library(lsa)

# Clustering
# library(clues)
# library(mclust)
library(pheatmap)
library(Rtsne)
# library(NMF)
library(Seurat)
library(umap)
# FIt-SNE implementation (functon to use is fftRtsne. Use fast_tsne_path indicated below)
source('FIt-SNE/fast_tsne.R')
fast_tsne_path <- 'FIt-SNE/bin/fast_tsne'

# Define repeated kmeans function here
# ds: sample x feature matrix
# pca: whether or not to run PCA first
# num_pcs: number of PCs to keep
# iterations: number of iterations to run (number of kmeans)
# n_clusters: number of clusters to separate into to generate the similarity score
# plus_minus: if you want to vary the n_clusters to be within a random range (e.g, in range 190-210, then set n_clusers to 200 and plus_minus to 10)
repeated_k_means <- function(ds, pca=TRUE, num_pcs=200, iterations=1000, n_clusters=2, plus_minus=0) {
  # Run PCA
  for_clustering <- ds
  if (pca) {
    # Run PCA
    num_pcs <- min(num_pcs, ncol(ds))
    pca_res <- prcomp(ds)$x
    
    for_clustering <- pca_res[,paste0('PC', seq(1, num_pcs))]
  }
  
  # Run the clustering
  keep_count <- NULL
  pb <- txtProgressBar(min = 0, max = iterations, style = 3)
  
  for (r in seq(1, iterations)) {
    increment_clust <- sample(seq((n_clusters-plus_minus), n_clusters+plus_minus), 1)
    clusters <- suppressWarnings(kmeans(as.matrix(for_clustering), increment_clust)$cluster)
    
    if (is.null(keep_count)) {
      keep_count <- data.frame(run_1=clusters)
    } else {
      keep_count[names(clusters),paste0('run_', r)] <- clusters[names(clusters)]
    }
    setTxtProgressBar(pb, r)
  }
  print('Completed clustering iterations. Aggregating results...')
  
  # Now calculate the similarity score
  repeated_kmeans_similarity <- apply(keep_count, 1, FUN = function(x) apply(keep_count, 1, FUN = function(y) length(which(x==y))))
  repeated_kmeans_similarity <- repeated_kmeans_similarity/iterations

  return(repeated_kmeans_similarity)
}

## Publication theme ggplot
theme_Publication <- function(base_size=14, base_family="Helvetica") {
  theme_foundation(base_size=base_size, base_family=base_family) +
    theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          # panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.justification = 'top',
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold")
    )
}

# cols_to_use_for_groups <- NULL
# 
# scale_fill_Publication <- function(...){
#   discrete_scale("fill","Publication",manual_pal(values = cols_to_use_for_groups), ...)
# }
# 
# scale_colour_Publication <- function(...){
#   discrete_scale("colour","Publication",manual_pal(values = cols_to_use_for_groups), ...)
# }
# 
# scale_color_Publication <- function(...){
#   discrete_scale("colour","Publication",manual_pal(values = cols_to_use_for_groups), ...)
# }

options(stringsAsFactors = FALSE)

colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color, 
            x)
  } else x
}

as.num = function(x, na.strings = c("NA", "")) {
  stopifnot(is.character(x))
  na = (x %in% na.strings) | (grepl("[A-Ba-b]", x))
  x[na] = 0
  x = as.numeric(x)
  x[na] = NA_real_
  x
}
