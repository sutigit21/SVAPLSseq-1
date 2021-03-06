% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svplsSurr.R
\name{svplsSurr}
\alias{svplsSurr}
\title{svplsSurr}
\usage{
svplsSurr(dat, group, controls = NULL, phi = function(x) log(x + const),
  const = 1, pls.method = "oscorespls", max.surrs = 3, cutoff = 10^-7,
  parallel = FALSE, num.cores = NULL, plot = FALSE)
}
\arguments{
\item{dat}{A gene expression count matrix or a 'SummarizedExperiment' object 
or a 'DGEList' object.}

\item{group}{a factor representing the sample indices belonging to the two 
different groups.}

\item{controls}{The set of control probes with no differential expression 
between the two groups (set to NULL by default).}

\item{phi}{The transforming function to be applied on the original gene 
expression count data (set to be log function with an offset \code{const}).}

\item{const}{The offset parameter for the transforming function \code{phi} 
(set to 1 by default).}

\item{pls.method}{The non-linear partial least squares method to be used. 
The different options available are: the classical orthogonal scores 
algorithm ("oscorespls", default), the kernel algorithm ("kernelpls") and 
wide kernel algorithm ("widekernelpls"). Using the "oscorespls" option is 
recommended for producing mutually orthogonal surrogate variables.}

\item{max.surrs}{The maximum number if surrogate variables to be extracted 
from the NPLS algorithm (set to 3 by default).}

\item{cutoff}{The user-specified pvalue cutoff for testing the significance
of the extracted surrogate variables (set to 1e-07 by default).}

\item{parallel}{Logical, indicating if the computations should be 
parallelized or not (set to \code{FALSE} by default).}

\item{num.cores}{The requested number of cores to be used in the parallel 
computations inside the function (used only when \code{parallel} is 
\code{TRUE}, \code{NULL} by default).}

\item{plot}{Logical, if \code{TRUE} a barplot of the variance proportions 
explained by the significant surrogate variables is returned (set to 
\code{FALSE} by default).}
}
\value{
An \code{\link{svplsSurr}} object.
}
\description{
This function extracts the surrogated estimates of the hidden variables in the 
data by using the partial least squares (PLS) algorithm 
on two multivariate random matrices. It provides the user with two options: 

(1) \bold{Unsupervised SVAPLS}: Here a standard linear regression model is first used on
a transformed version of the expression count matrix to estimate the primary signals 
of differential expression for all the genes. The fitted model residuals and the 
transformed count matrix are then organized respectively into two multivariate matrices 
\code{E} and \code{Y}, in such a way that each column corresponds to a certain gene.  
\code{E} is then regressed on \code{Y} using a Non-linear partial least squares (NPLS) 
algorithm and the extracted scores in the column-space of \code{Y} are deemed as the surrogate 
variables.

(2) \bold{Supervised SVAPLS}: In case information on a set of control genes (probes) is provided,
this function uses a Non-linear partial least squares (NPLS) algorithm to regress \code{Y} on a 
submatrix of \code{Y} (\code{Y.sub}) corresponding to the set of controls and scores in the column-
space of \code{Y.sub} are considered as the surrogate variables.

The function then regresses the first eigenvector of the residual matrix \code{E} (for Unsupervised
SVAPLS or the control matrix \code{Y.sub} for Supervised SVAPLS) on these surrogate variables and 
tests them for statistical significance with a certain user-specified pvalue cutoff. The variables 
yielding a pvalue below the cutoff are returned.
}
\examples{
##Loading the simulated dataset
data(sim.dat)

##Extracting the significant surrogate variables
group = as.factor(c(rep(1, 10), rep(-1, 10)))
sv <- svplsSurr(dat = sim.dat, group = group)
print(sv)

}

