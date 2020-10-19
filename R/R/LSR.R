#' Log-Shape Ratios (from Rmorph, M. Baylac)
#'
#' @description This function computes log-shape ratios and size from a matrix of morphological measurements (columns) from individual specimens (rows).
#' @details The code was copied and modified from the Rmorph package from M. Baylac. By default, the function applies a natural log transform to the data, so if your data matrix is already log-transformed, please set transf=F.
#'
#' @param mat The data matrix with specimens as rows and measurements as columns. Missing values (NAs) will cause problems.
#' @param transf Logical. If TRUE (default) a natural log is applied to te data matrix.
#'
#' @return size. A vector containing sizes (geometric mean of measurements).
#' @return lsr. A matrix containing log-shape ratio values.
#' @return nobj. The number of specimens (objects).
#' @return nvar. The number of measurements (variables).
#'
#' @export
#' @examples
#' data("equus", package="MorphoInd")
#' m <- matrix(NA, 22, 15) #Make empty matrix
#' for (i in 1:15) {m[,i] <- rnorm(n=22,mean=equus[[1]]["Mean",i], sd=equus[[1]]["St-Dev",i])} #Generate normally distributed morphological data based on summary statistics of the first sample in the equus dataset.
#' LSR <-lsr(m) #Compute size and LSR values (warning: NAs can be produced if the generated data has negative values, which is not possible for real data).


lsr <- function(mat, transf=TRUE) {

	if(ncol(mat) < 2) stop("lsr: mat should have at least two columns")
	mat <- as.matrix(mat)
	nobj <- nrow(mat)
	nvar <-ncol(mat)
	if (transf == TRUE) mat <- log(mat)
	size <- apply(mat,1,mean)
	lsr <- scale(t(mat), center = TRUE, scale = FALSE)
	lsr <- t(lsr)
	attr(lsr, "scaled:center") <- NULL
	structure(list(size= size, lsr= lsr, nobj= nobj, nvar= nvar, transf= transf))
}

