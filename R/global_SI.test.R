#' Global Size Index (LSI, VSI, VSI*) test
#'
#' @description This function will test for the difference in SI between sample(s) and the reference sample of the dataset. This test is global, i.e. across all morphological descriptors (variables).
#' @details The test is based on H0: no significant difference between sample and reference. The H0 distribution of the statistic is computed based on the number of iterations (iter). The statistic is the average of absolute mean SI values for each sample (which represent the global absolute difference between sample and reference). The expected input is 1) the original dataset, as formatted for the morpho_indices function (argument dat) and 2) the output from the morpho_indices function (argument x).
#'
#' @param dat Input data should be a list of matrices/data frames (same as for morpho_indices function).
#' @param x Output from the morpho_indices function, containing the dataframes with SI values, as well as bootstraped values.
#' @param ref Index of the reference sample (numeric or character), default is 1, so the first sample in 'dat' is taken to be the reference. Must be the same as defined for the morpho_indices function.
#' @param compsp Vector of index(es) of the samples for which the difference is to be tested (numeric or character). Default is 1:length(dat), i.e. all samples, including the reference.
#' @param iter Number of values generated under H0, i.e. the size of the H0 distribution (default = 1000).
#'
#' @return mat.res A matrix containing observed values of the statistic, 0.95 quantile values for the statistic under H0 (i.e. critical statistic values for significance with alpha=0.05), and P-values (i.e. P(SI >= SI_obs)). 
#'
#' @export
#' @examples 
#' data("equus", package="MorphoInd")
#' o <- morpho_indices(dat=equus, ref=1, data.type="summary", bootstrap="p", plot=F)
#' test <- global_SI.test(dat=equus, x=o, ref=1, compsp=2:10, iter=1000)

global_SI.test <- function(dat, x, ref=1, compsp=1:length(dat), iter=1000) { #Start of function

N <- ncol(x$LSI.df)
k <- nrow(x$LSI.df)

Ri <- unlist(dat[[ref]]["Mean",])
s_Ri <- unlist(dat[[ref]]["St-Dev",])
n_Ri <- unlist(dat[[ref]]["N",])

mat.res <- matrix(NA, ncol=9, nrow=k)

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
   
   meanLSIi_h0 <- apply(LSI_h0, 2, mean)
   meanVSIi_h0 <- apply(VSI_h0, 2, mean)
   meanVSI2i_h0 <- apply(VSI2_h0, 2, mean)

   meanLSI_h0 <- mean(meanLSIi_h0)
   meanVSI_h0 <- mean(meanVSIi_h0)
   meanVSI2_h0 <- mean(meanVSI2i_h0)

   abs_LSI_h0 <- abs(LSI_h0)
   abs_VSI_h0 <- abs(VSI_h0)
   abs_VSI2_h0 <- abs(VSI2_h0)

   distri_abs_LSI <- apply(abs_LSI_h0, 1, mean)
   distri_abs_VSI <- apply(abs_VSI_h0, 1, mean)
   distri_abs_VSI2 <- apply(abs_VSI2_h0, 1, mean) #These are the distributions for our statistic under H0

   absmeanLSIi_h0 <- apply(abs_LSI_h0, 2, mean)
   absmeanVSIi_h0 <- apply(abs_VSI_h0, 2, mean)
   absmeanVSI2i_h0 <- apply(abs_VSI2_h0, 2, mean)

   absmeanLSI_h0 <- mean(absmeanLSIi_h0)
   absmeanVSI_h0 <- mean(absmeanVSIi_h0)
   absmeanVSI2_h0 <- mean(absmeanVSI2i_h0)

   abs_real_LSI <- mean(abs(unlist(x$LSI.df[j,])))
   abs_real_VSI <- mean(abs(unlist(x$VSI.df[j,])))
   abs_real_VSI2 <- mean(abs(unlist(x$VSI2.df[j,]))) #These are the observed values of our statistic (mean across all descriptors of absolute SI value=>difference between reference and compared sample). 
   crit_LSI <- sort(distri_abs_LSI)[0.95*iter]
   crit_VSI <- sort(distri_abs_VSI)[0.95*iter]
   crit_VSI2 <- sort(distri_abs_VSI2)[0.95*iter]

   P_LSI <- sum(distri_abs_LSI >= abs_real_LSI)/iter
   P_VSI <- sum(distri_abs_VSI >= abs_real_VSI)/iter
   P_VSI2 <- sum(distri_abs_VSI2 >= abs_real_VSI2)/iter #P-values. These are unilateral because our statistic is strictly positive.

   mat.res[j,] <- c(abs_real_LSI, abs_real_VSI, abs_real_VSI2, crit_LSI, crit_VSI, crit_VSI2, P_LSI, P_VSI, P_VSI2)
   } #End of big loop for compsp

colnames(mat.res) <- c("Obs. mean |LSI|", "Obs. mean |VSI|", "Obs. mean |VSI2|", "LSI 0.95 quantile H0", "VSI 0.95 quantile H0", "VSI2 0.95 quantile H0", "P (LSI)","P (VSI)","P (VSI2)")
rownames(mat.res) <- rownames(x$LSI.df)

return(mat.res)
}

