#' batch.table.freq
#'
#' @param x factor vector: batch vector
#' @param batch_names character vector: batch names
#' @description
#' Internal helper function.
#' @export
batch.table.freq <- function(x, batch_names) {
  env.table <- table(x)
  res <- rep(0, length(batch_names))
  names(res) <- batch_names
  res[names(env.table)] <- env.table
  res <- res/sum(res)
  return(res)
}

#' KLD
#'
#' @param P probability distribution P
#' @param Q probability distribution Q
#' @description
#' Computes the Kullback-Leibler-Divergence for distributions P and Q.
#' @export
KLD <- function(P,Q) {
  KL <- 0
  for(i in 1:length(P)) {
    if(P[i] != 0) {
      KL <- KL + P[i]*log2(P[i]/Q[i])
    }
  }
  return(KL)
}

#' NESTscore
#'
#' @param obj Seurat object
#' @param group_by name of metadata-column with batch assignment
#' @param k_nn number of nearest neighbors to be considered. If set to "cellwise", a cell-wise k is computed based on a the connectivity of a fuzzy umap graph (computed with `uwot::umap` default parameters on the considered dimensionality reduction).
#' @param reduction_in name of dimensionality reduction to be considered for nearest neighbor search
#' @param ndims number of dimensions in dimensionality reduction to be considered for nearest neighbor search
#' @param assay_in name of seurat object assay
#' @param return_pairwise_eval return pairwise batch-to-batch evaluation as matrix
#' @param show_heatmap plots pairwise batch-to-batch matrix as heatmap (only considered if `return_pairwise_eval = TRUE`)#' @param name_suffix additional suffix for new NEST-Score metadata column name
#' @param heatmap_col color for the optimal batch frequency in the pairwise batch-to-batch heatmap (only considered if `return_pairwise_eval = TRUE` and `show_heatmap = TRUE`)
#' @param name_suffix additional suffix for new NEST-Score metadata column name
#' @returns Returns a list including the Seurat object `seuratobj` with cell-wise NEST-Score stored in an additional metadata column and the theoretical upper and lower bound of the NEST-Score `NESTscore_limits` for the considered batch assignment.
#'          If `return_pairwise_eval = TRUE` the list also includes a pairwise batch-to-batch evaluation matrix `pairwise_matrix` and helpful plotting parameters `proposed_params_pheatmap` for using `pheatmap()`.
#' @description
#' Computes the cell-wise NEST-Score for a given batch variable. For the NEST-Score the frequency vector of batches is computed for each cell considering their k-nearest neighborhood (knn).
#' This vector is divided by the global frequency vector of batches and scaled to add up to 1. To compare the resulting distribution vector against a uniform distribution, the negative Kullback-Leibler-Divergence is computed.
#' A high NEST-Score (close to 0) means, that the distribution of samples in the cell's knn matches the global distribution. A low NEST-Score (negative) means that a cell's knn is a less homogeneous, i.e. consists of fewer batches and/or different batch frequencies than the global distribution.
#' If all considered batches would mix randomly throughout the data, we would expect a NEST-Score close to 0 for all cells.
#'
#' @export
NESTscore <- function(obj,
                      group_by,
                      k_nn = 100,
                      reduction_in = "pca",
                      ndims = 50,
                      assay_in = Seurat::DefaultAssay(obj),
                      return_pairwise_eval = FALSE,
                      show_heatmap = FALSE,
                      heatmap_col = "darkred",
                      name_suffix = "") {

  # catch missing reductions, metadata or wrong assays
  if(methods::is(obj@reductions[[reduction_in]])[1] != "DimReduc") {
    stop("Perform ", reduction_in, " dimensionality reduction on ",assay_in," assay or adjust assay.")
  } else {
    if(obj@reductions[[reduction_in]]@assay.used != assay_in) {
      stop("Perform ", reduction_in, " dimensionality reduction on ",assay_in," assay or adjust assay.")
    }
  }
  if(!(group_by %in% colnames(obj@meta.data))) {
    stop(group_by, " not found in object's metadata.")
  }

  batch <- factor(obj@meta.data[,group_by], levels = unique(obj@meta.data[,group_by]))
  global_frequency <- table(batch)/sum(table(batch))
  uniform_freq <- rep(1/length(levels(batch)),length(levels(batch)))

  # compute lower bound of nest score
  worst_freq <- rep(0,length(levels(batch)))
  worst_freq[1] <- 1
  worst_value <- -KLD(worst_freq, uniform_freq)

  if(k_nn == "cellwise") {
    umap_new <- uwot::umap(X = obj@reductions[[reduction_in]]@cell.embeddings,
                           ret_extra = "fgraph")
    cellwise_knn <- rowSums(as.matrix(umap_new$fgraph)>0)
    k_nn <- max(cellwise_knn)

    obj <- Seurat::FindNeighbors(obj,
                                 reduction = reduction_in,
                                 k.param = k_nn+1,
                                 return.neighbor = T,
                                 dims = 1:ndims,
                                 assay = assay_in)

    env.total <- obj@neighbors[[paste0(assay_in,".nn")]]@nn.idx
    NA_indx <- as.matrix(dplyr::bind_rows(apply(cbind(1:length(cellwise_knn),cellwise_knn), 1,
                                function(x) data.frame("x" = x[1],"y" = (x[2]+1):(k_nn+1)))))
    env.total[NA_indx] <- NA

  } else {
    # find k-nearest neighbors for all cells
    obj <- Seurat::FindNeighbors(obj,
                                 reduction = reduction_in,
                                 k.param = k_nn+1,
                                 return.neighbor = T,
                                 dims = 1:ndims,
                                 assay = assay_in)

    env.total <- obj@neighbors[[paste0(assay_in,".nn")]]@nn.idx
  }

  env <- env.total[,2:(k_nn+1)]
  dim.store <- dim(env)
  env <- c(env)
  env.batch <- batch[env]
  env.batch <- matrix(env.batch, dim.store[1],dim.store[2])

  # frequency of batch samples per cell
  freq.mat <- t(apply(env.batch, 1, batch.table.freq, batch_names = names(table(batch))))
  # scale with global batch sample frequency
  freq.mat.glob <- t(apply(freq.mat, 1, function(x) x/global_frequency))
  # scale to sum to 1
  freq.mat.glob.scaled <- t(apply(freq.mat, 1, function(x) x/sum(x)))

  KLD_vec <- t(apply(freq.mat.glob.scaled[,colnames(freq.mat.glob.scaled) %in% levels(batch)], 1, function(x) KLD(x, uniform_freq)))[1,]

  obj <- Seurat::AddMetaData(obj, as.numeric(-KLD_vec),
                             col.name = paste0("NESTscore_",group_by,name_suffix))
  # add command log
  command <- Seurat::LogSeuratCommand(object = obj, return.command = TRUE)
  methods::slot(object = command, name = "assay.used") <- assay_in
  command.name <- methods::slot(object = command, name = "name")
  obj[[command.name]] <- command

  if(return_pairwise_eval) {
    heatmap.matrix.ratio <- matrix(NA, nrow = length(levels(batch)), ncol = length(levels(batch)))
    colnames(heatmap.matrix.ratio) <- rownames(heatmap.matrix.ratio) <- levels(batch)

    # average scaled frequencies for each sample
    for(b in levels(batch)) {
      freq.mat.sub <- freq.mat.glob.scaled[which(batch == b),]
      heatmap.matrix.ratio[b,] <- colMeans(freq.mat.sub)[colnames(heatmap.matrix.ratio)]
    }

    heatmap.matrix.ratio <- apply(heatmap.matrix.ratio, 2, function(x) x/sum(x, na.rm = T))
    heatmap.matrix.ratio <- apply(heatmap.matrix.ratio, 1, function(x) x/sum(x, na.rm = T))

    color_vec = c(grDevices::colorRampPalette(c("grey90",heatmap_col))(which.min(abs(seq(0,1,0.005)-1/length(levels(batch))))),
                  grDevices::colorRampPalette(c(heatmap_col,"grey20"))((length(seq(0,1,0.005)) - 1 - which.min(abs(seq(0,1,0.005)-1/length(levels(batch)))))))
    breaks_vec = seq(0,1,0.005)

    if(show_heatmap) {
      pheatmap::pheatmap(heatmap.matrix.ratio,
                         cluster_rows = F,
                         cluster_cols = F,
                         color = grDevices::adjustcolor(color_vec,.8),
                         breaks = breaks_vec,
                         na_col = "transparent",
                         display_numbers = round(heatmap.matrix.ratio,2),
                         fontsize_number=20)
    }

    return(list("seuratobj" = obj,
                "NESTscore_limits" = c(worst_value,0),
                "pairwise_matrix" = heatmap.matrix.ratio,
                "proposed_params_pheatmap" = list("breaks" = breaks_vec,
                                                  "color" = color_vec)))
  } else {
    return(list("seuratobj" = obj,
                "NESTscore_limits" = c(worst_value,0)))
  }
}

#' NEST_threshold_proposal
#'
#' @param obj Seurat object
#' @param group_by name of metadata-column with batch assignment
#' @param min_nsamples minimum number of samples that should be uniformly mixed for lower NEST-Score bound
#' @param thresh_distribution (only considered if `min_nsamples == NULL`) limiting frequency distribution for lower NEST-Score bound. The length of the vector has to match the number of batches.
#' @returns proposed NEST-Score threshold to binarily assign cells to have a homogeneous or heterogeneous neighborhood regarding the considered batch assignment.
#' @description
#' If `thresh_distribution = NULL` the lower bound distribution is a vector of zeros, where only `min_nsamples` entries have value `1/min_nsamples`. The output threshold
#'
#' @export
NEST_threshold_proposal <- function(obj,
                                    group_by,
                                    min_nsamples = 2,
                                    thresh_distribution = NULL) {

  # catch missing metadata or too high min_nsamples
  if(!(group_by %in% colnames(obj@meta.data))) {
    stop(group_by, " not found in object's metadata.")
  } else {
    if(!(paste0("NESTscore_",group_by) %in% colnames(obj@meta.data))) {
      warning(paste0("No matching NESTscore found for assignment ",group_by,". `NESTscore()` is called."))
      obj <- NESTscore(obj = obj, group_by = group_by)$seuratobj
    }
  }

  batch <- factor(obj@meta.data[,group_by], levels = unique(obj@meta.data[,group_by]))

  if(min_nsamples > length(levels(batch))) {
    warning("`min_nsamples` can not be higher than number  of unique batches.")
    warning("`min_nsamples` is set to the number of unique batches.")
    min_nsamples == length(levels(batch))
  }

  uniform_freq <- rep(1/length(levels(batch)),
                      length(levels(batch)))

  if(is.null(thresh_distribution)) {
    thresh_freq <- rep(0,length(levels(batch)))
    thresh_freq[1:min_nsamples] <- 1/min_nsamples
  } else {
    thresh_freq <- thresh_distribution
    thresh_freq <- thresh_freq/sum(thresh_freq)
  }
  
  thresh <- -KLD(thresh_freq, uniform_freq)
  
  return(thresh)
}

#' NEST_search_k
#'
#' @param obj Seurat object
#' @param group_by name of metadata-column with batch assignment
#' @param k_search_vector minimum and maximum number of nearest neighbors to be considered, maximal number has to be #cells -1.
#' @param color_plot_initial output plot is colored by cell-wise binary assignment based on `k_initial` and `thresh_initial`
#' @param k_initial number of nearest neighbors used for initial binary assignment for plotting
#' @param thresh_initial NEST-score threshold used for initial binary assignment for plotting. If `NULL` (default), it will be computed by `NEST_threshold_proposal()` with default parameters.
#' @param reduction_in name of dimensionality reduction to be considered for nearest neighbor search
#' @param ndims number of dimensions in dimensionality reduction to be considered for nearest neighbor search
#' @param assay_in name of seurat object assay
#' @description
#' Plots the NEST-Score for all cells across multiple choices of `k` nearest neighbors.
#'
#' @export
NEST_search_k <- function(obj,
                          group_by,
                          k_search_vector = seq(2, 100, length.out = 50),
                          color_plot_initial = T,
                          k_initial = k_search_vector[round(length(k_search_vector)/2)],
                          thresh_initial = NULL,
                          reduction_in = "pca",
                          ndims = 50,
                          assay_in = Seurat::DefaultAssay(obj)) {

  # catch missing reductions, metadata or wrong assays
  if(methods::is(obj@reductions[[reduction_in]])[1] != "DimReduc") {
    stop("Perform ", reduction_in, " dimensionality reduction on ",assay_in," assay or adjust assay.")
  } else {
    if(obj@reductions[[reduction_in]]@assay.used != assay_in) {
      stop("Perform ", reduction_in, " dimensionality reduction on ",assay_in," assay or adjust assay.")
    }
  }
  if(!(group_by %in% colnames(obj@meta.data))) {
    stop(group_by, " not found in object's metadata.")
  }
  if(!(k_initial %in% k_search_vector) & color_plot_initial) {
    stop("`k_initial` has to be in `k_search_vector`.")
  }
  if(!is.integer(k_search_vector)) {
    warning("`k_search_vector` contains non-integer values and hence is rounded.")
    k_search_vector <- round(k_search_vector)
  }
  if(max(k_search_vector) > (dim(obj)[2]-1)) {
    warning("Some value of `k_search_vector` are bigger than #cells-1. Those will be ignored.")
    k_search_vector <- k_search_vector[(k_search_vector <= (dim(obj)[2]-1))]
  }

  if(is.null(thresh_initial)) {
    thresh_initial <- NEST_threshold_proposal(obj = obj, group_by = group_by)
  }

  batch <- factor(obj@meta.data[,group_by], levels = unique(obj@meta.data[,group_by]))

  global_frequency <- table(batch)/sum(table(batch))
  uniform_freq <- rep(1/length(levels(batch)),length(levels(batch)))

  obj <- Seurat::FindNeighbors(obj,
                               reduction = reduction_in,
                               k.param = max(k_search_vector)+1,
                               return.neighbor = T,
                               dims = 1:ndims,
                               assay = assay_in)

  env.total <- obj@neighbors[[paste0(assay_in,".nn")]]@nn.idx

  mix_score_mat <- matrix(NA, nrow = dim(obj)[2], ncol = length(k_search_vector))
  count <- 0

  for(nnn in k_search_vector) {

    print(paste0(round(count/length(k_search_vector),4)*100,"% of iterations"))
    count <- count + 1

    # NESTScore Code
    env <- env.total[,2:(nnn+1)]
    dim.store <- dim(env)
    env <- c(env)
    env.batch <- batch[env]
    env.batch <- matrix(env.batch, dim.store[1],dim.store[2])

    freq.mat <- t(apply(env.batch, 1, batch.table.freq, batch_names = names(table(batch))))
    freq.mat.glob <- t(apply(freq.mat, 1, function(x) x/global_frequency))
    freq.mat.glob.scaled <- t(apply(freq.mat, 1, function(x) x/sum(x)))

    KLD_vec <- t(apply(freq.mat.glob.scaled[,colnames(freq.mat.glob.scaled) %in% levels(batch)], 1, function(x) KLD(x, uniform_freq)))[1,]
    mix_score_mat[,count] <- -as.numeric(KLD_vec)
  }

  df_plot <- reshape2::melt(data.frame(t(mix_score_mat)))
  df_plot$rep <- rep(k_search_vector, dim(mix_score_mat)[1])

  if(color_plot_initial) {
    df_plot$id <- rep(mix_score_mat[,which(k_search_vector == k_initial)] > thresh_initial,
                      each = dim(mix_score_mat)[2])

    df_plot$NEST_bin <- "heterogeneous"
    df_plot$NEST_bin[df_plot$id] <- "homogeneous"
    df_plot$NEST_bin <- factor(df_plot$NEST_bin,
                               levels = c("homogeneous","heterogeneous"))

    p <- ggplot2::ggplot(df_plot,ggplot2::aes(x = rep,
                                              y = value,
                                              group.by = variable,
                                              col = NEST_bin)) +
      ggplot2::geom_line(alpha = 0.01, lwd = 0.2) +
      ggplot2::guides(col = ggplot2::guide_legend(title="")) +
      ggplot2::scale_color_manual(values = c("heterogeneous" = "#ffc425",
                                             "homogeneous" = "#00aedb")) +
      ggplot2::geom_line(alpha = 0.01, lwd = 0.2) +
      ggplot2::xlim(0,max(k_search_vector)) +
      ggplot2::geom_hline(yintercept = thresh_initial, lty = 2) +
      ggplot2::geom_vline(xintercept = k_initial, lty = 2) +
      ggplot2::ylab("NEST-Score") +
      ggplot2::xlab("number of nearest neighbors") +
      ggplot2::theme_minimal() +
      #Seurat::NoLegend() +
      ggplot2::theme(text = ggplot2::element_text(size = 8, colour = "black"),
                     axis.text = ggplot2::element_text(size = 6, colour = "black"),
                     axis.title = ggplot2::element_text(size = 8, colour = "black"),
                     title = ggplot2::element_text(size = 8, colour = "black"),
                     legend.text = ggplot2::element_text(size = 8, colour = "black"),
                     strip.text = ggplot2::element_text(size = 8, colour = "black"))

  } else {
    p <- ggplot2::ggplot(df_plot,ggplot2::aes(x = rep,
                                              y = value,
                                              group.by = variable)) +
      ggplot2::geom_line(alpha = 0.01, lwd = 0.2) +
      ggplot2::guides(col = ggplot2::guide_legend(title="")) +
      ggplot2::geom_line(alpha = 0.01, lwd = 0.2) +
      ggplot2::xlim(min(k_search_vector),max(k_search_vector)) +
      ggplot2::geom_hline(yintercept = thresh_initial, lty = 2) +
      ggplot2::geom_vline(xintercept = k_initial, lty = 2) +
      ggplot2::ylab("NEST-Score") +
      ggplot2::xlab("number of nearest neighbors") +
      ggplot2::theme_minimal() +
      #Seurat::NoLegend() +
      ggplot2::theme(text = ggplot2::element_text(size = 8, colour = "black"),
                     axis.text = ggplot2::element_text(size = 6, colour = "black"),
                     axis.title = ggplot2::element_text(size = 8, colour = "black"),
                     title = ggplot2::element_text(size = 8, colour = "black"),
                     legend.text = ggplot2::element_text(size = 8, colour = "black"),
                     strip.text = ggplot2::element_text(size = 8, colour = "black"))

  }

  return(p)
}
