% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pi_estimation.R
\name{pi_estimation}
\alias{pi_estimation}
\title{A function to get robust cell-type composition estimation based on the prior estimates and transcriptomic data.}
\usage{
pi_estimation(expression_matrix, n_iteration = 5, prior)
}
\arguments{
\item{expression_matrix}{G by N gene expression matrix for estimating cell type composition of the tissue.
Each row is a gene and each column is tissue sample. Gene name and sample ID can be included as row and column names,
but not in the expression matrix.}

\item{n_iteration}{Number of bootstrap samples to estimate the pi.}

\item{prior}{Prior estimation of cell-type composition.}
}
\value{
Updated cell-type composition estimates. A ribble of 2 columns: The first columns shows the sample id, the second column shows the robust cell-type composition estimation
}
\description{
A function to get robust cell-type composition estimation based on the prior estimates and transcriptomic data.
}
