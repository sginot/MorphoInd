#' Simpson's Log-Size Index (LSI) and Uerpmann's Variability Size Index, with bootstrap
#'
#' @description The function computes Simpson's Log-Size Index (LSI) and Uerpmann's Variability Size Index (and bootstrap). 
#' @details The function computes Simpson's Log-Size Index (LSI) and Uerpmann's Variability Size Index. Works with several samples, one of which is the reference. Samples are constituted by: Several measurements for several individuals, or summary statistics for each descriptor (N, mean, sd).
#' The function can be used with: List of matrices (samples x individuals x measurements), or list of matrices (samples x summary stats x measurements), or mixed list of matrices (some with individual measurements, other with summary stats). In the latter case, matrices with individual measurements will be converted to matrices with summary stats.
#' The function includes parametric and non-parametric bootstraps, this allows to test differences between samples. Tests can be global first then check which measurement matters.
#'
#' @param dat Input data should be a list of matrices/data frames
#' @param ref Index of the reference sample, default is 1, so the first sample in 'dat' is taken to be the reference.
#' @param k Value for k (arbitrary constant for VSI calculation).
#' @param data.type Type of the data. Can be 'individuals', if all samples in 'dat' contain detailed measurements by individual. Can be 'summary' (default), if all samples in 'dat' have N, Mean, and St-Dev (summary statistics) as lines. Can be 'both' if both types of data are mixed, in that case, argument 'which.ind' must be specified. When data.type="both", all samples ith individuals are converted to summary statistics before computing anyhting.
#' @param which.ind Numeric vector containing index of samples which have individual measurements (including reference if needed).
#' @param which.desc Which morphological descriptors (traits) are to be used. By default all are used. Note: Descriptors can also be filtered when plotting with the 'plot.mo.ind' function.
#' @param N.sample Number of samples (including reference).
#' @param bootstrap Type of bootstrap to be used. Default is 'off' = no bootstrap. Can be 'p' = parametric bootstrap. Can be 'np' = non-parametric bootstrap, only possible when using 'individuals' as data type. If data.type="both", only the parametric bootstrap is used (samples are converted to summary statistics). To combine 'np' and 'p' bootstrap when using mixed data, run the function separately for individual samples and summary samples (both with the SAME reference sample, which should be individuals). Outputs can then be combined using the 'combin.mo.ind' function.
#' @param iter Number of iterations for the resampling (default = 1000).
#' @param VCV Maintain variance-covariance structure when resampling or not (default = F). VCV = T only available for bootstrap="np" so far. This means that individuals can be resampled 'as a whole' or resampled fo each descriptor. Differences in confidence intervals should be minor between both.
#' @param center.mISD Should the data be centered with mISD (isometric size)? Default=F.
#' @param plots Should basic plots be plotted (default = F). More flexible plotting can be done with 'plot.mo.ind' function.
#' @param col Vector of colors for the plot (length=number of samples).
#'
#' @return LSI.df Dataframe of values of Simpson's Log-Size Index (Simpson 1941).
#' @return VSI.df Dataframe of values of Uerpmann's Variability Size Index (Upermann 1982).
#' @return VSI2.df Dataframe of values of modified Uerpmann's Variability Size Index (Escarguel 2008).
#' @return LSI.bootstrap List of matrices of bootstrapped values of LSI. Each matrix represents boostraps for one sample.
#' @return VSI.bootstrap Same as previous for VSI.
#' @return VSI2.bootstrap Same as previous for VSI2.
#' @return LSI.mISD Vector of values of mean isometric size difference (mISD) for LSI.
#' @return VSI.mISD Same as previous for VSI.
#' @return VSI2.mISD Same as previous for VSI2
#' @return LSI.CI List of matrices of 95% confidence intervals for LSI. Line 1 is the lower interval values, line 2 is the upper interval values.
#' @return VSI.CI Same as previous for VSI.
#' @return VSI2.CI Same as previous for VSI2.
#' @return LSI.meanboot Matrix of LSI bootstrapped average values.
#' @return VSI.meanboot Same as previous for VSI.
#' @return VSI2.meanboot Same as previous for VSI2.
#' @return LSI.medianboot Matrix of LSI bootstrapped median values.
#' @return VSI.medianboot Same as previous for VSI.
#' @return VSI2.medianboot Same as previous for VSI2.
#' @return LSI.sdboot Matrix of LSI standard deviation values.
#' @return VSI.sdboot Same as previous for VSI.
#' @return VSI2.sdboot Same as previous for VSI2.
#'
#' @examples
#' o <- morpho.indices(dat=dat, ref=1, data.type="both", which.ind=c(1,2,8,10), bootstrap="p", plot=T) #This will compute LSI, VSI and VSI2 value, parametric bootstrap values, CIs, and plots.

morpho.indices <- function(
	dat = data.list, 
	ref = 1, 
	k = 25, 
	data.type = "summary", 
	which.ind = NULL, 
	which.desc = 1:dim(dat[[1]])[2], 
	N.sample = length(dat), 
	bootstrap = "off", 
        iter = 1000, 
	VCV = F, 
	center.mISD = F, 
	plots = F, 
	col=c("black", rainbow(N.sample-1)) 
	) {

   if (!is.list(dat)) {stop("Dataset must be a list")}
   if (!data.type %in% c('summary', 'individuals', 'both')) {stop("Data type specified incorrectly: argument 'data.type' must be one of 'summary', 'individuals', 'both'.")}

   N.descriptors <- length(which.desc)
   for (i in 1:N.sample) {dat[[i]] <- dat[[i]][, which.desc]}

   LSI.df <- data.frame(matrix(NA, ncol = N.descriptors, nrow = N.sample), row.names = names(dat))
   VSI.df <- data.frame(matrix(NA, ncol = N.descriptors, nrow = N.sample), row.names = names(dat))
   VSI2.df <- data.frame(matrix(NA, ncol = N.descriptors, nrow = N.sample), row.names = names(dat))

   ### This part is for data with individual measurements for ALL the samples####
   if (data.type == "individuals") { ### Basic calculation of indices
   
   ref.samp <- dat[[ref]]

   R <- apply(ref.samp, 2, mean, na.rm=T)
   Sr <- apply(ref.samp, 2, sd, na.rm=T)
   Nr <- dim(ref.samp)[1] - apply(is.na(ref.samp), 2, sum)
      
      for (i in 1:N.sample) { ### Start of for loop to go through all samples
      samp <- dat[[i]]
      Xj <- apply(samp, 2, mean, na.rm=T)
      Sj <- apply(samp, 2, sd, na.rm=T)
      Nj <- dim(samp)[1] - apply(is.na(samp), 2, sum)

      Sx <- sqrt((((Nj - 1) * Sj^2) + ((Nr - 1) * Sr^2)) / (Nj + Nr -2))

      VSI2.df[i,] <- (Xj - R) / Sx
      LSI.df[i,] <- log(Xj) - log(R)
      VSI.df[i,] <- k * (Xj - R) / Sr
      message("Sample ", i, ": Indices computed")
      } ### End of for loop (basic calculation of indices)
   
      if (bootstrap == "np") { ### Non-parametric bootstrap estimates, ie resampling within samples with replacement. Should the reference values also be bootstraped?
      message("Non-parametric bootstrap")

      LSI.bootstrap <- VSI.bootstrap <- VSI2.bootstrap <- list()
      LSI.CI <- VSI.CI <- VSI2.CI <- list()

         if (VCV) { ### Start of bootstrap estimates with VCV structure preserved (resampling done by row within samples)
         message("Variance-covariance structure PRESERVED")

            for (i in 1:N.sample) {
            message("Samples bootstrap: ", i)
            mat.boot.VSI2 <- mat.boot.VSI <- mat.boot.LSI <- matrix(NA, nrow=iter, ncol=N.descriptors)

               for (j in 1:iter) {
               samp <- dat[[i]]

                  repeat { ### Repeat loop : if by resampling one of the descriptors does not have at least 2 measurements, it will restart.
                  bref <- ref.samp[sample.int(dim(ref.samp)[1], replace=T),]
                  bNr <- dim(bref)[1] - apply(is.na(bref), 2, sum)
                  if (sum(bNr == 0 | bNr == 1) == 0) {break}
                  }

                  repeat { ### Repeat loop: if when resampling one of the descriptors does not have at least 2 measurements, it will restart.
                  bsamp <- samp[sample.int(dim(samp)[1], replace=T),]
                  bNj <- dim(bsamp)[1] - apply(is.na(bsamp), 2, sum)
                  if (sum(bNj == 0 | bNj == 1) == 0) {break}
                  }

               bR <- apply(bref, 2, mean, na.rm=T)
               bSr <- apply(bref, 2, sd, na.rm=T)
               bXj <- apply(bsamp, 2, mean, na.rm=T)
               bSj <- apply(bsamp, 2, sd, na.rm=T)
               bSx <- sqrt((((bNj - 1) * bSj^2) + ((bNr - 1) * bSr^2)) / (bNj + bNr -2))

               mat.boot.VSI2[j,] <- (bXj - bR) / bSx
               mat.boot.LSI[j,] <- log(bXj) - log(bR)
               mat.boot.VSI[j,] <- k * (bXj - bR) / bSr
               }

            LSI.bootstrap[[i]] <- mat.boot.LSI
            VSI.bootstrap[[i]] <- mat.boot.VSI
            VSI2.bootstrap[[i]] <- mat.boot.VSI2

            VSI.CI[[i]] <- apply(mat.boot.VSI, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            VSI2.CI[[i]] <- apply(mat.boot.VSI2, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            LSI.CI[[i]] <- apply(mat.boot.LSI, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            }

         } ### End of bootstrap estimates with VCV structure preserved

         if (!VCV) { ### Start of bootstrap estimates with VCV structure removed (resampling done for each descriptors separately within samples)
         message("Variance-covariance structure REMOVED")

            for (i in 1:N.sample) {
            message("Samples bootstrap: ", i)

            mat.boot.LSI <- matrix(NA, nrow=iter, ncol=N.descriptors)
            mat.boot.VSI <- matrix(NA, nrow=iter, ncol=N.descriptors)
            mat.boot.VSI2 <- matrix(NA, nrow=iter, ncol=N.descriptors)

               for (j in 1:iter) {
               samp <- dat[[i]]

                  repeat { ### Repeat loop : if by resampling one of the descriptors does not have at least 2 measurements, it will restart.
                  bref <- apply(ref.samp, 2, sample,  dim(ref.samp)[1], replace=T)
                  bNr <- dim(bref)[1] - apply(is.na(bref), 2, sum)
                  if (sum(bNr == 0 | bNr == 1) == 0) {break}
                  }

                  repeat { ### Repeat loop : if by resampling one of the descriptors does not have at least 2 measurements, it will restart.
                  bsamp <- apply(samp, 2, sample,  dim(samp)[1], replace=T)
                  bNj <- dim(bsamp)[1] - apply(is.na(bsamp), 2, sum)
                  if (sum(bNj == 0 | bNj == 1) == 0) {break}
                  }

               bR <- apply(bref, 2, mean, na.rm=T)
               bSr <- apply(bref, 2, sd, na.rm=T)
               bXj <- apply(bsamp, 2, mean, na.rm=T)
               bSj <- apply(bsamp, 2, sd, na.rm=T)
               bSx <- sqrt((((bNj - 1) * bSj^2) + ((bNr - 1) * bSr^2)) / (bNj + bNr -2))

               mat.boot.VSI2[j,] <- (bXj - bR) / bSx
               mat.boot.LSI[j,] <- log(bXj) - log(bR)
               mat.boot.VSI[j,] <- k * (bXj - bR) / bSr
               }

            LSI.bootstrap[[i]] <- mat.boot.LSI
            VSI.bootstrap[[i]] <- mat.boot.VSI
            VSI2.bootstrap[[i]] <- mat.boot.VSI2

            VSI.CI[[i]] <- apply(mat.boot.VSI, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            VSI2.CI[[i]] <- apply(mat.boot.VSI2, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            LSI.CI[[i]] <- apply(mat.boot.LSI, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            }

         }  ### End of bootstrap estimates with VCV structure removed
      } ### End of non-parametric bootstrap computation
 

      if (bootstrap == "p") { ### Start of parametric bootstrap part
      message("Parametric bootstrap (normal distribution is assumed)")

      LSI.bootstrap <- VSI.bootstrap <- VSI2.bootstrap <- list()
      LSI.CI <- VSI.CI <- VSI2.CI <- list()

         if (VCV) {message("Parametric bootstrap with variance-covariance structure based on Cholesky factorization")

         #mat.means <- matrix(NA, nrow=length(dat), ncol=N.descriptors)

            #for (i in 1:length(dat)) {
            #mat.means[i,] <- apply(dat[[i]], 2, mean, na.rm=T)
            #}

         #GC <- cor(mat.means)
         #chol.GC <- chol(GC)
         }

         if (!VCV) { ### Start of bootstrap estimates with VCV structure removed (resampling done for each descriptors separately within samples)
         message("No variance-covariance structure between descriptors")
         
         ref.samp <- rbind(dim(ref.samp)[1] - apply(is.na(ref.samp), 2, sum), apply(ref.samp, 2, mean, na.rm=T), apply(ref.samp, 2, sd, na.rm=T))

            for (i in 1:N.sample) {
            message("Samples bootstrap: ", i)

            mat.boot.LSI <- matrix(NA, nrow=iter, ncol=N.descriptors)
            mat.boot.VSI <- matrix(NA, nrow=iter, ncol=N.descriptors)
            mat.boot.VSI2 <- matrix(NA, nrow=iter, ncol=N.descriptors)

            Ndat <- dim(dat[[i]])[1] - apply(is.na(dat[[i]]), 2, sum)
            Mdat <- apply(dat[[i]], 2, mean, na.rm=T)
            SDdat <- apply(dat[[i]], 2, sd, na.rm=T)

            samp <- rbind(Ndat, Mdat, SDdat)

               for (j in 1:iter) {
               #Parametric pseudo values generation and computation of indices.

               bNj <- unlist(samp[1,])
               bSr <- bR <- bXj <- bSj <- bSx <- rep(NA, N.descriptors)

                  for (h in 1:N.descriptors) {
                     repeat { #Repeat: if a <= 0, it will restart.
                     a <- rnorm(1, mean=samp[2,h], sd=sqrt(samp[3,h]^2/samp[1,h]))
                     b <- rnorm(1, mean=ref.samp[2,h], sd=sqrt(ref.samp[3,h]^2/ref.samp[1,h]))
                     if (a > 0 & b > 0) {bXj[h] <- a ; bR[h] <- b ; break}
                     } #WARNING: This avoids NAs from log() BUT produces a truncated normal distribution when SD >= mean.

                  bSj[h] <- sqrt(((bNj[h]-1)*samp[3,h]^2)/rchisq(1, df=unlist(bNj[h])-1))
                  bSr[h] <- sqrt(((ref.samp[1,h]-1)*ref.samp[3,h]^2)/rchisq(1, df=unlist(ref.samp[1,h]-1)))
                  } 

               bSx <- sqrt((((bNj - 1) * bSj^2) + ((ref.samp[1,] - 1) * bSr^2)) / (bNj + ref.samp[1,] -2))

               mat.boot.VSI2[j,] <- (bXj - bR) / bSx
               mat.boot.LSI[j,] <- log(bXj) - log(bR)
               mat.boot.VSI[j,] <- k * (bXj - bR) / bSr
               }

            LSI.bootstrap[[i]] <- mat.boot.LSI
            VSI.bootstrap[[i]] <- mat.boot.VSI
            VSI2.bootstrap[[i]] <- mat.boot.VSI2

            VSI.CI[[i]] <- apply(mat.boot.VSI, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            VSI2.CI[[i]] <- apply(mat.boot.VSI2, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            LSI.CI[[i]] <- apply(mat.boot.LSI, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            }
         }

      }  ### End of parametric bootstrap part

      if (bootstrap == "off") { ### Start = If no bootstrap computation required list results for indices
      message("No bootstrap")
      } ### End of no bootstrap

   } #### End of computation function for cases with individual data for ALL samples ####


   ### This part is for data with mixed samples, some with individual measurements, others with summaries ####
   if (data.type == "both" & is.null(which.ind)) {
   stop("For mixed data the argument which.ind must be specified.")
   }

   if (data.type == "both" & !is.null(which.ind)) {

      for (i in which.ind) { #Loop to convert individual data samples into summaries.
      m <- dat[[i]]
      mm <- apply(m, 2, mean, na.rm=T)
      msd <- apply(m, 2, sd, na.rm=T)
      mn <- dim(m)[1] - apply(is.na(m), 2, sum)
      dat[[i]] <- data.frame(rbind(mn, mm, msd), row.names=c("N", "Mean", "St-Dev"))
      }

   ref.samp <- as.matrix(dat[[ref]])

   R <- ref.samp["Mean",]
   Sr <- ref.samp["St-Dev",]
   Nr <- ref.samp["N",]

      for (i in 1:N.sample) {
      samp <- dat[[i]]
      Xj <- samp["Mean",]
      Sj <- samp["St-Dev",]
      Nj <- samp["N",]

      Sx <- sqrt((((Nj - 1) * Sj^2) + ((Nr - 1) * Sr^2)) / (Nj + Nr -2))

      VSI2.df[i,] <- (Xj - R) / Sx
      LSI.df[i,] <- log(Xj) - log(R)
      VSI.df[i,] <- k * (Xj - R) / Sr
      message("Sample ", i, ": Indices computed")
      }
      
      if (bootstrap == "off") { ### Start = If no bootstrap computation required list results for indices
      message("No bootstrap")
      } ### End of no bootstrap

      if (bootstrap == "np") { ### Start = non-parametric bootstrap
      warning("Non-parametric bootstrap is only possible for data.type=individuals")
      LSI.bootstrap <- VSI.bootstrap <- VSI2.bootstrap <- list()
      LSI.CI <- VSI.CI <- VSI2.CI <- list()
      } ### End of non-param bootstrap

      if (bootstrap == "p") { ### Start = parametric bootstrap
      message("Parametric bootstrap (normal distribution assumed)")
      LSI.bootstrap <- VSI.bootstrap <- VSI2.bootstrap <- list()
      LSI.CI <- VSI.CI <- VSI2.CI <- list()

         if (VCV) {warning("Parametric bootstrap with variance-covariance structure in development")}

         if (!VCV) { ### Start of bootstrap estimates with VCV structure removed (resampling done for each descriptors separately within samples)
         message("No variance-covariance structure between descriptors")

            for (i in 1:N.sample) {
            message("Samples bootstrap: ", i)

            mat.boot.LSI <- matrix(NA, nrow=iter, ncol=N.descriptors)
            mat.boot.VSI <- matrix(NA, nrow=iter, ncol=N.descriptors)
            mat.boot.VSI2 <- matrix(NA, nrow=iter, ncol=N.descriptors)

            samp <- as.matrix(dat[[i]])
            
               for (j in 1:iter) {
               #Parametric pseudo values generation and computation of indices.

               bNj <- unlist(samp[1,])
               bSr <- bR <- bXj <- bSj <- bSx <- rep(NA, N.descriptors)

                  for (h in 1:N.descriptors) {
                     repeat { #Repeat: if a <= 0, it will restart.
                     a <- rnorm(1, mean=samp[2,h], sd=sqrt(samp[3,h]^2/samp[1,h]))
                     b <- rnorm(1, mean=ref.samp[2,h], sd=sqrt(ref.samp[3,h]^2/ref.samp[1,h]))
                     if (a > 0 & b > 0) {bXj[h] <- a ; bR[h] <- b ; break}
                     } #WARNING: This avoids NAs from log() BUT produces a truncated normal distribution when SD >= mean.

                  bSj[h] <- sqrt(((bNj[h]-1)*samp[3,h]^2)/rchisq(1, df=unlist(bNj[h])-1))
                  bSr[h] <- sqrt(((ref.samp[1,h]-1)*ref.samp[3,h]^2)/rchisq(1, df=unlist(ref.samp[1,h]-1)))
                  } 

               bSx <- sqrt((((bNj - 1) * bSj^2) + ((ref.samp[1,] - 1) * bSr^2)) / (bNj + ref.samp[1,] -2))

               mat.boot.VSI2[j,] <- (bXj - bR) / bSx
               mat.boot.LSI[j,] <- log(bXj) - log(bR)
               mat.boot.VSI[j,] <- k * (bXj - bR) / bSr
               }

            LSI.bootstrap[[i]] <- mat.boot.LSI
            VSI.bootstrap[[i]] <- mat.boot.VSI
            VSI2.bootstrap[[i]] <- mat.boot.VSI2

            VSI.CI[[i]] <- apply(mat.boot.VSI, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            VSI2.CI[[i]] <- apply(mat.boot.VSI2, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            LSI.CI[[i]] <- apply(mat.boot.LSI, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            } #End of for loop (i)

         } #End of bootstrap estimate with VCV. 

      } ### End of parametric bootstrap

   } #### End of computation for cases with mixed data ####

   ### This part is for data with only summary statistics for ALL samples ####
   if (data.type == "summary") {
   ref.samp <- as.matrix(dat[[ref]])

   R <- ref.samp["Mean",]
   Sr <- ref.samp["St-Dev",]
   Nr <- ref.samp["N",]

      for (i in 1:N.sample) {
      samp <- as.matrix(dat[[i]])
      Xj <- samp["Mean",]
      Sj <- samp["St-Dev",]
      Nj <- samp["N",]

      Sx <- sqrt((((Nj - 1) * Sj^2) + ((Nr - 1) * Sr^2)) / (Nj + Nr -2))

      VSI2.df[i,] <- (Xj - R) / Sx
      LSI.df[i,] <- log(Xj) - log(R)
      VSI.df[i,] <- k * (Xj - R) / Sr
      message("Sample ", i, ": Indices computed")
      }

      if (bootstrap == "np") { ### Start = non-parametric bootstrap
      warning("Non-parametric bootstrap is only possible for data.type=individuals")
      LSI.bootstrap <- VSI.bootstrap <- VSI2.bootstrap <- list()
      LSI.CI <- VSI.CI <- VSI2.CI <- list()
      } ### End of non-param bootstrap

      if (bootstrap == "p") { ### Start = parametric bootstrap
      message("Parametric bootstrap (normal distribution assumed)")
      LSI.bootstrap <- VSI.bootstrap <- VSI2.bootstrap <- list()
      LSI.CI <- VSI.CI <- VSI2.CI <- list()
          
         if (VCV) {warning("Parametric bootstrap with variance-covariance structure in development")}

         if (!VCV) { ### Start of bootstrap estimates with VCV structure removed (resampling done for each descriptors separately within samples)
         message("No variance-covariance structure between descriptors")

            for (i in 1:N.sample) {
            message("Samples bootstrap: ", i)

            mat.boot.LSI <- matrix(NA, nrow=iter, ncol=N.descriptors)
            mat.boot.VSI <- matrix(NA, nrow=iter, ncol=N.descriptors)
            mat.boot.VSI2 <- matrix(NA, nrow=iter, ncol=N.descriptors)

            samp <- as.matrix(dat[[i]])

               for (j in 1:iter) {
               #Parametric pseudo values generation and computation of indices.

               bNj <- unlist(samp[1,])
               bSr <- bR <- bXj <- bSj <- bSx <- rep(NA, N.descriptors)

                  for (h in 1:N.descriptors) {
                     repeat { #Repeat: if a <= 0, it will restart.
                     a <- rnorm(1, mean=samp[2,h], sd=sqrt(samp[3,h]^2/samp[1,h]))
                     b <- rnorm(1, mean=ref.samp[2,h], sd=sqrt(ref.samp[3,h]^2/ref.samp[1,h]))
                     if (a > 0 & b > 0) {bXj[h] <- a ; bR[h] <- b ; break}
                     } #WARNING: This avoids NAs from log() BUT produces a truncated normal distribution when SD >= mean.

                  bSj[h] <- sqrt(((bNj[h]-1)*samp[3,h]^2)/rchisq(1, df=unlist(bNj[h])-1))
                  bSr[h] <- sqrt(((ref.samp[1,h]-1)*ref.samp[3,h]^2)/rchisq(1, df=unlist(ref.samp[1,h]-1)))
                  } 

               bSx <- sqrt((((bNj - 1) * bSj^2) + ((ref.samp[1,] - 1) * bSr^2)) / (bNj + ref.samp[1,] -2))

               mat.boot.VSI2[j,] <- (bXj - bR) / bSx
               mat.boot.LSI[j,] <- log(bXj) - log(bR)
               mat.boot.VSI[j,] <- k * (bXj - bR) / bSr
               }

            LSI.bootstrap[[i]] <- mat.boot.LSI
            VSI.bootstrap[[i]] <- mat.boot.VSI
            VSI2.bootstrap[[i]] <- mat.boot.VSI2

            VSI.CI[[i]] <- apply(mat.boot.VSI, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            VSI2.CI[[i]] <- apply(mat.boot.VSI2, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            LSI.CI[[i]] <- apply(mat.boot.LSI, 2, sort)[c(round(0.025*iter),round(0.975*iter)),]
            } #End of for loop (i)

         } #End of bootstrap estimate with VCV. 

      } ### End of parametric bootstrap
      
   } #### End of computation for cases with summary data ####

   if (bootstrap=="p" | bootstrap =="np") {
   LSI.meanboot <- LSI.medianboot <- LSI.sdboot <- VSI.meanboot <- VSI.medianboot <- VSI.sdboot <-VSI2.meanboot <- VSI2.medianboot <- VSI2.sdboot <- matrix(NA, N.sample, N.descriptors)

   LSI.mISD <- VSI.mISD <- VSI2.mISD <- rep(NA, N.sample)

      for (i in 1:N.sample) {
      LSI.meanboot[i,] <- apply(LSI.bootstrap[[i]], 2, mean)
      LSI.medianboot[i,] <- apply(LSI.bootstrap[[i]], 2, quantile, 0.5)
      LSI.sdboot[i,] <- apply(LSI.bootstrap[[i]], 2, sd)
      VSI.meanboot[i,] <- apply(VSI.bootstrap[[i]], 2, mean)
      VSI.medianboot[i,] <- apply(VSI.bootstrap[[i]], 2, quantile, 0.5)
      VSI.sdboot[i,] <- apply(VSI.bootstrap[[i]], 2, sd)
      VSI2.meanboot[i,] <- apply(VSI2.bootstrap[[i]], 2, mean)
      VSI2.medianboot[i,] <- apply(VSI2.bootstrap[[i]], 2, quantile, 0.5)
      VSI2.sdboot[i,] <- apply(VSI2.bootstrap[[i]], 2, sd)

      LSI.mISD[i] <- sum(1/LSI.sdboot[i,]*LSI.meanboot[i,]) / sum(1/LSI.sdboot[i,])
      VSI.mISD[i] <- sum(1/VSI.sdboot[i,]*VSI.meanboot[i,]) / sum(1/VSI.sdboot[i,])
      VSI2.mISD[i] <- sum(1/VSI2.sdboot[i,]*VSI2.meanboot[i,]) / sum(1/VSI2.sdboot[i,])
      }
   }

   if (plots) { ### Plotting part, if required
   layout(matrix(1:3, nrow=1))

   plot(unlist(LSI.df[1,]), type="o", lwd = 1.5, pch=20, ylim = c(min(LSI.df), max(LSI.df)), xlab="Descriptor", ylab="LSI", col=col[1])

   abline(v=c(1:N.descriptors), col="gray", lty=2)

   legend("bottomleft", legend=names(dat), col=col, lty=1, lwd=1.5, bty="n")
   for (i in 1:(N.sample-1)) {points(unlist(LSI.df[i+1,]), type="o", lwd = 1.5, pch=20, col = col[i+1])}

   plot(unlist(VSI.df[1,]), type="o", lwd = 1.5, pch=20, ylim = c(min(VSI.df), max(VSI.df)), xlab="Descriptor", ylab="VSI", col=col[1])

   abline(v=c(1:N.descriptors), col="gray", lty=2)

   for (i in 1:(N.sample-1)) {points(unlist(VSI.df[i+1,]), type="o", pch=20, lwd = 1.5, col = col[i+1])}

   plot(unlist(VSI2.df[1,]), type="o", lwd = 1.5, pch=20, ylim = c(min(VSI2.df), max(VSI2.df)), xlab="Descriptor", ylab="VSI*", col=col[1])

   abline(v=c(1:N.descriptors), col="gray", lty=2)

   for (i in 1:(N.sample-1)) {points(unlist(VSI2.df[i+1,]), type="o", pch=20, lwd = 1.5, col = col[i+1])}
   } ### End of plotting part

if (bootstrap == "off") {
	return(list(
	LSI.df=LSI.df,
	VSI.df=VSI.df,
	VSI2.df=VSI2.df
	))}

if (!center.mISD & bootstrap == "np" | !center.mISD & bootstrap == "p") {
	return(list(
	LSI.df=LSI.df,
	VSI.df=VSI.df,
	VSI2.df=VSI2.df,
	LSI.bootstrap=LSI.bootstrap,
	VSI.bootstrap=VSI.bootstrap,
	VSI2.bootstrap=VSI2.bootstrap, 
	LSI.CI=LSI.CI,
	VSI.CI=VSI.CI,
	VSI2.CI=VSI2.CI,
	LSI.meanboot=LSI.meanboot,
	LSI.medianboot=LSI.medianboot,
	VSI.meanboot=VSI.meanboot,
	VSI.medianboot=VSI.medianboot,
	VSI2.meanboot=VSI2.meanboot,
	VSI2.medianboot=VSI2.medianboot,
	LSI.sdboot=LSI.sdboot,
	VSI.sdboot=VSI.sdboot,
	VSI2.sdboot=VSI2.sdboot,
	LSI.mISD=LSI.mISD,
	VSI.mISD=VSI.mISD,
	VSI2.mISD=VSI2.mISD
	))}

if (center.mISD & bootstrap == "np" | center.mISD & bootstrap == "p") {

   for (i in 1:N.sample) {
   LSI.CI[[i]] <- LSI.CI[[i]] - LSI.mISD[i]
   VSI.CI[[i]] <- VSI.CI[[i]] - VSI.mISD[i]
   VSI2.CI[[i]] <- VSI2.CI[[i]] - VSI2.mISD[i]
   LSI.df[i,] <- LSI.df[i,] - LSI.mISD[i]
   VSI.df[i,] <- VSI.df[i,] - VSI.mISD[i]
   VSI2.df[i,] <- VSI2.df[i,] - VSI2.mISD[i]
   LSI.meanboot[i,] <- LSI.meanboot[i,] - LSI.mISD[i]
   VSI.meanboot[i,] <- VSI.meanboot[i,] - VSI.mISD[i]
   VSI2.meanboot[i,] <- VSI2.meanboot[i,] - VSI2.mISD[i]
   LSI.medianboot[i,] <- LSI.medianboot[i,] - LSI.mISD[i]
   VSI.medianboot[i,] <- VSI.medianboot[i,] - VSI.mISD[i]
   VSI2.medianboot[i,] <- VSI2.medianboot[i,] - VSI2.mISD[i]
   }
	return(list(
	LSI.df=LSI.df,
	VSI.df=VSI.df,
	VSI2.df=VSI2.df,
	LSI.bootstrap=LSI.bootstrap,
	VSI.bootstrap=VSI.bootstrap,
	VSI2.bootstrap=VSI2.bootstrap, 
	LSI.mISD=LSI.mISD,
	VSI.mISD=VSI.mISD,
	VSI2.mISD=VSI2.mISD,
	LSI.CI=LSI.CI,
	VSI.CI=VSI.CI,
	VSI2.CI=VSI2.CI,
	LSI.meanboot=LSI.meanboot,
	LSI.medianboot=LSI.medianboot,
	VSI.meanboot=VSI.meanboot,
	VSI.medianboot=VSI.medianboot,
	VSI2.meanboot=VSI2.meanboot,
	VSI2.medianboot=VSI2.medianboot,
	LSI.sdboot=LSI.sdboot,
	VSI.sdboot=VSI.sdboot,
	VSI2.sdboot=VSI2.sdboot
	))}

} #### End of function ####



