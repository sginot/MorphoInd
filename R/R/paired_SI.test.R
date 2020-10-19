#' Paired Size Index (LSI, VSI, VSI*) test
#'
#' @description This function will test for the difference in SI between sample(s) and the reference sample of the dataset. This tests descriptors individually, i.e. each morphological descriptor will have its own test.
#' @details NB: This test should be used only when the global test (global_SI.test function) showed significant global differences between sample(s) and reference. The test is based on H0: no significant difference between sample and reference. The statistic follows a Welch's T distribution. The expected input is 1) the original dataset, as formatted for the morpho_indices function (argument dat) and 2) the output from the morpho_indices function (argument x).
#'
#' @param dat Input data should be a list of matrices/data frames (same as for morpho_indices function).
#' @param x Output from the morpho_indices function, containing the dataframes with SI values, as well as bootstraped values.
#' @param ref Index of the reference sample (numeric or character), default is 1, so the first sample in 'dat' is taken to be the reference. Must be the same as defined for the morpho_indices function.
#' @param compsp Vector of index(es) of the samples for which the difference is to be tested (numeric or character). Default is 1:length(dat), i.e. all samples, including the reference.
#'
#' @return mat.res A matrix containing observed values of the statistic, 0.95 quantile values for the statistic under H0 (i.e. critical statistic values for significance with alpha=0.05), and P-values (i.e. P(SI >= SI_obs)). 
#'
#' @export
#' @examples 
#' data("equus", package="MorphoInd")
#' o <- morpho_indices(dat=equus, ref=1, data.type="summary", bootstrap="p", plot=F)
#' test <- paired_SI.test(dat=equus, x=o, ref=1, compsp=2:10)


paired_SI.test <- function(dat, x, ref=1, compsp=1:length(dat)) {

N <- ncol(x$LSI.df)
k <- length(c(ref,compsp))

Ri <- unlist(dat[[ref]]["Mean",])
s_Ri <- unlist(dat[[ref]]["St-Dev",])
n_Ri <- unlist(dat[[ref]]["N",])

mat.t.LSI <- matrix(NA, ncol=N, nrow=k)
mat.df.LSI <- matrix(NA, ncol=N, nrow=k)
mat.P.LSI <- matrix(NA, ncol=N, nrow=k)

mat.t.VSI <- matrix(NA, ncol=N, nrow=k)
mat.df.VSI <- matrix(NA, ncol=N, nrow=k)
mat.P.VSI <- matrix(NA, ncol=N, nrow=k)

mat.t.VSI2 <- matrix(NA, ncol=N, nrow=k)
mat.df.VSI2 <- matrix(NA, ncol=N, nrow=k)
mat.P.VSI2 <- matrix(NA, ncol=N, nrow=k)

   for (j in compsp) { #Start of loop for specified comparisons (compsp)
   n_Xi <- unlist(dat[[j]]["N",])
   s_Xi <- unlist(dat[[j]]["St-Dev",])

   Xib_h0 <- Rib_h0 <- matrix(NA, ncol=length(Ri), nrow=iter)
      for (i in 1:length(Ri)) {
      Xib_h0[,i] <- rnorm(n=iter, mean=Ri[i], sd=s_Ri[i]/sqrt(n_Xi[i]))
      Rib_h0[,i] <- rnorm(n=iter, mean=Ri[i], sd=s_Ri[i]/sqrt(n_Ri[i]))
      }

   s_Xib <- s_Rib <- matrix(NA, ncol=length(Ri), nrow=iter)
      for (i in 1:length(Ri)) {
      s_Xib[,i] <- sqrt(((n_Xi[i]-1)*s_Ri[i]^2)/rchisq(n=iter, df=n_Xi[i]-1))
      s_Rib[,i] <- sqrt(((n_Ri[i]-1)*s_Ri[i]^2)/rchisq(n=iter, df=n_Ri[i]-1))
      }

   s_Xb <- sqrt(((n_Xi-1)*(s_Xib^2)+(n_Ri-1)*(s_Rib^2))/(n_Xi+n_Ri-2)) #Not entirely sure this is correct

   LSI_h0 <- log(Xib_h0/Rib_h0)
   VSI_h0 <- 25*(Xib_h0 - Rib_h0)/s_Rib
   VSI2_h0 <- (Xib_h0 - Rib_h0)/s_Xb #See calculation of s_Xb
   
   meanLSIi_b <- apply(x$LSI.bootstrap[[j]], 2, mean) 
   meanVSIi_b <- apply(x$VSI.bootstrap[[j]], 2, mean) 
   meanVSI2i_b <- apply(x$VSI2.bootstrap[[j]], 2, mean) 

   s_LSIi_b <- apply(x$LSI.bootstrap[[j]], 2, sd)
   s_VSIi_b <- apply(x$VSI.bootstrap[[j]], 2, sd)
   s_VSI2i_b <- apply(x$VSI2.bootstrap[[j]], 2, sd)

   s_LSIi_h0 <- apply(LSI_h0, 2, sd)
   s_VSIi_h0 <- apply(VSI_h0, 2, sd)
   s_VSI2i_h0 <- apply(VSI2_h0, 2, sd)

   ti_LSI <- abs(meanLSIi_b - meanLSIi_h0) / sqrt((s_LSIi_b^2 + s_LSIi_h0^2)/n_Xi)
   ti_VSI <- abs(meanVSIi_b - meanVSIi_h0) / sqrt((s_VSIi_b^2 + s_VSIi_h0^2)/n_Xi)
   ti_VSI2 <- abs(meanVSI2i_b - meanVSI2i_h0) / sqrt((s_VSI2i_b^2 + s_VSI2i_h0^2)/n_Xi)

   vi_LSI <- ((s_LSIi_b^2 + s_LSIi_h0^2)/n_Xi)^2/((s_LSIi_b^4 + s_LSIi_h0^4)/(n_Xi^3-n_Xi^2))
   vi_VSI <- ((s_VSIi_b^2 + s_VSIi_h0^2)/n_Xi)^2/((s_VSIi_b^4 + s_VSIi_h0^4)/(n_Xi^3-n_Xi^2))
   vi_VSI2 <- ((s_VSI2i_b^2 + s_VSI2i_h0^2)/n_Xi)^2/((s_VSI2i_b^4 + s_VSI2i_h0^4)/(n_Xi^3-n_Xi^2))

   Pi_LSI <- 1 - pt(ti_LSI, df=vi_LSI, lower.tail=T) + pt(ti_LSI, df=vi_LSI, lower.tail=F)
   Pi_VSI <- 1 - pt(ti_VSI, df=vi_VSI, lower.tail=T) + pt(ti_VSI, df=vi_VSI, lower.tail=F)
   Pi_VSI2 <- 1 - pt(ti_VSI2, df=vi_VSI2, lower.tail=T) + pt(ti_VSI2, df=vi_VSI2, lower.tail=F)

   mat.t.LSI[j,] <- ti_LSI
   mat.t.VSI[j,] <- ti_VSI
   mat.t.VSI2[j,] <- ti_VSI2
 
   mat.df.LSI[j,] <- vi_LSI
   mat.df.VSI[j,] <- vi_VSI
   mat.df.VSI2[j,] <- vi_VSI2

   mat.P.LSI[j,] <- Pi_LSI
   mat.P.VSI[j,] <- Pi_VSI
   mat.P.VSI2[j,] <- Pi_VSI2
   }

colnames(mat.t.LSI) <- colnames(mat.t.VSI) <- colnames(mat.t.VSI2) <- colnames(mat.df.LSI) <- colnames(mat.df.VSI) <- colnames(mat.df.VSI2) <- colnames(mat.P.LSI) <- colnames(mat.P.VSI) <- colnames(mat.P.VSI2) <- colnames(x$LSI.df)

rownames(mat.t.LSI) <- rownames(mat.t.VSI) <- rownames(mat.t.VSI2) <- rownames(mat.df.LSI) <- rownames(mat.df.VSI) <- rownames(mat.df.VSI2) <- rownames(mat.P.LSI) <- rownames(mat.P.VSI) <- rownames(mat.P.VSI2) <- rownames(x$LSI.df)
  
message("P-values for LSI")
print(data.frame(mat.P.LSI))

message("P-values for VSI")
print(data.frame(mat.P.VSI))

message("P-values for VSI2")
print(data.frame(mat.P.VSI2))


return(list(t.LSI=mat.t.LSI, df.LSI=mat.df.LSI, P.LSI=mat.P.LSI, t.VSI=mat.t.VSI, df.VSI=mat.df.VSI, P.VSI=mat.P.VSI, t.VSI2=mat.t.VSI2, df.VSI2=mat.df.VSI2, P.VSI2=mat.P.VSI2))
}

