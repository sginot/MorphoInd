#' Combining outputs from separate runs of function 'morpho.ind'
#'
#' This is to be used ONLY if both runs had the same reference sample. Generally it wwill be used when mixed data is used and one run of morpho.ind will have parametric bootstraps while the other has non-parametric bootstrap.

combin.mo.ind <- function(out1=out1, out2=out2) {

out.combin <- vector("list", length=length(out1))

names(out.combin) <- names(testindi)

out.combin$LSI.df <- rbind(out1$LSI.df, out2$LSI.df)
out.combin$VSI.df <- rbind(out1$VSI.df, out2$VSI.df)
out.combin$VSI2.df <- rbind(out1$VSI2.df, out2$VSI2.df)

out.combin$LSI.meanboot <- rbind(out1$LSI.meanboot, out2$LSI.meanboot)
out.combin$VSI.meanboot <- rbind(out1$VSI.meanboot, out2$VSI.meanboot)
out.combin$VSI2.meanboot <- rbind(out1$VSI2.meanboot, out2$VSI2.meanboot)

out.combin$LSI.medianboot <- rbind(out1$LSI.medianboot, out2$LSI.medianboot)
out.combin$VSI.medianboot <- rbind(out1$VSI.medianboot, out2$VSI.medianboot)
out.combin$VSI2.medianboot <- rbind(out1$VSI2.medianboot, out2$VSI2.medianboot)

out.combin$LSI.sdboot <- rbind(out1$LSI.sdboot, out2$LSI.sdboot)
out.combin$VSI.sdboot <- rbind(out1$VSI.sdboot, out2$VSI.sdboot)
out.combin$VSI2.sdboot <- rbind(out1$VSI2.sdboot, out2$VSI2.sdboot)

out.combin$LSI.bootstrap <- c(out1$LSI.bootstrap, out2$LSI.bootstrap)
out.combin$VSI.bootstrap <- c(out1$VSI.bootstrap, out2$VSI.bootstrap)
out.combin$VSI2.bootstrap <- c(out1$VSI2.bootstrap, out2$VSI2.bootstrap)

out.combin$LSI.CI <- c(out1$LSI.CI, out2$LSI.CI)
out.combin$VSI.CI <- c(out1$VSI.CI, out2$VSI.CI)
out.combin$VSI2.CI <- c(out1$VSI2.CI, out2$VSI2.CI)

out.combin$LSI.mISD <- c(out1$LSI.mISD, out2$LSI.mISD)
out.combin$VSI.mISD <- c(out1$VSI.mISD, out2$VSI.mISD)
out.combin$VSI2.mISD <- c(out1$VSI2.mISD, out2$VSI2.mISD)

return(out.combin)
}
