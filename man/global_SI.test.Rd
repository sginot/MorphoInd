% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global_SI.test.R
\name{global_SI.test}
\alias{global_SI.test}
\title{Global Size Index (LSI, VSI, VSI*) test}
\usage{
global_SI.test(dat, x, ref = 1, compsp = 1:length(dat), iter = 1000)
}
\arguments{
\item{dat}{Input data should be a list of matrices/data frames (same as for morpho_indices function).}

\item{x}{Output from the morpho_indices function, containing the dataframes with SI values, as well as bootstraped values.}

\item{ref}{Index of the reference sample (numeric or character), default is 1, so the first sample in 'dat' is taken to be the reference. Must be the same as defined for the morpho_indices function.}

\item{compsp}{Vector of index(es) of the samples for which the difference is to be tested (numeric or character). Default is 1:length(dat), i.e. all samples, including the reference.}

\item{iter}{Number of values generated under H0, i.e. the size of the H0 distribution (default = 1000).}
}
\value{
mat.res A matrix containing observed values of the statistic, 0.95 quantile values for the statistic under H0 (i.e. critical statistic values for significance with alpha=0.05), and P-values (i.e. P(SI >= SI_obs)).
}
\description{
This function will test for the difference in SI between sample(s) and the reference sample of the dataset. This test is global, i.e. across all morphological descriptors (variables).
}
\details{
The test is based on H0: no significant difference between sample and reference. The H0 distribution of the statistic is computed based on the number of iterations (iter). The statistic is the average of absolute mean SI values for each sample (which represent the global absolute difference between sample and reference). The expected input is 1) the original dataset, as formatted for the morpho_indices function (argument dat) and 2) the output from the morpho_indices function (argument x).
}
\examples{
data("equus", package="MorphoInd")
o <- morpho_indices(dat=equus, ref=1, data.type="summary", bootstrap="p", plot=F)
test <- global_SI.test(dat=equus, x=o, ref=1, compsp=2:10, iter=1000)
}
