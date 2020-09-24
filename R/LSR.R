#' Mosimann's (1970) Log-Shape Ratio (LSR)
#'
#' @description The function computes Log-Shape ratios from a linear measurements matrix.
#' @param mat Input data in the form of a matrix with individuals as lines and measurements (variables) as columns.
#' @param transf Should the matrix be log-transformed or not. If TRUE (default), the input matrix should be raw measurements, if FALSE, the input matrix should contain values that are already log-transformed.

#' @return size Size value(s) computed as the geometric mean of all measurements for each individual.
#' @lsr Log-Shape ratio values.
#' @nobj Number of 'objects', i.e. number of individuals.
#' @nvar Number of 'variables, i.e. number of measurements.
#'
#' @export


lsr<-function(mat, transf=TRUE) {
	print("This code was copied and modified from M. Baylac's Rmorph package.")
	print("NB : transf = T (default) means mat will be converted to log(mat). If mat already contains log values, use transf=F.")
	if(ncol(mat) < 2) stop("mat should have at least two columns")

	mat <- as.matrix(mat)
	nobj <- nrow(mat)
	nvar <-ncol(mat)
	if (transf == TRUE) mat <- log(mat)
	size <- apply(mat,1,mean)
	lsr <- scale(t(mat), center = TRUE, scale = FALSE)
	lsr <- t(lsr)

	structure(list(size= size, lsr= lsr, nobj= nobj, nvar= nvar))
}
