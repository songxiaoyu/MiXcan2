% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MiXcan2_Ensemble.R
\name{MiXcan2_ensemble}
\alias{MiXcan2_ensemble}
\title{MiXcan2 Ensemble}
\usage{
MiXcan2_ensemble(
  y,
  x,
  cov,
  pi,
  yName = NULL,
  xNameMatrix = NULL,
  B = 10,
  seed = NULL
)
}
\arguments{
\item{y}{The pre-cleaned expression level data for a single gene in N samples.}

\item{x}{A N by P matrix for all the genetic predictors used to predict the genetically regulated expression  of the gene.}

\item{cov}{A N by Q matrix for the covariates adjusted in the model (e.g. age, population stratification).}

\item{pi}{An estimation of cell-type faction of the cell type of interest.}

\item{yName}{Default is NULL. A row vector to save the expression information, such as gene ID, gene name.}

\item{xNameMatrix}{Default is NULL. A matrix to save the X matrix information,
such as variable ID, position, rsid, ref_allele, eff_allele.}

\item{B}{No. of ensemble models}

\item{seed}{seed}
}
\value{
list with 9 elements. It contains

\item{ensemble_summary}{Summarized model descriptions from all B ensemble models.}
\item{ensemble_summary_by_type}{Summarized model descriptions for Cell-Type-Specific and Non-Specific
models separately.}
\item{ensemble_weight}{Prediction weights from B ensemble models summarized separately for
Cell-Type-Specific and Non-Specific models.}
\item{ensemble_intecept}{Prediction intercepts from B ensemble models summarized separately for
Cell-Type-Specific and Non-Specific models.}
\item{CTS_weight}{Prediction weights summarized from all
Cell-Type-Specific models.}
\item{NS_weight}{Prediction weights summarized from all
Non-Specific models.}
\item{all_summary}{All B model descriptions w/o summarization (average).}
\item{all_weights}{All B weights w/o summarization (average).}
}
\description{
MiXcan2 Ensemble
}
