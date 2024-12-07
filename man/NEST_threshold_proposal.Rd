% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NESTscore.R
\name{NEST_threshold_proposal}
\alias{NEST_threshold_proposal}
\title{NEST_threshold_proposal}
\usage{
NEST_threshold_proposal(
  obj,
  group_by,
  min_nsamples = 2,
  thresh_distribution = NULL
)
}
\arguments{
\item{obj}{Seurat object}

\item{group_by}{name of metadata-column with batch assignment}

\item{min_nsamples}{minimum number of samples that should be uniformly mixed for lower NEST-Score bound}

\item{thresh_distribution}{(only considered if \code{min_nsamples == NULL}) limiting frequency distribution for lower NEST-Score bound. The length of the vector has to match the number of batches.}
}
\value{
proposed NEST-Score threshold to binarily assign cells to have a homogeneous or heterogeneous neighborhood regarding the considered batch assignment.
}
\description{
If \code{thresh_distribution = NULL} the lower bound distribution is a vector of zeros, where only \code{min_nsamples} entries have value \code{1/min_nsamples}. The output threshold
}
