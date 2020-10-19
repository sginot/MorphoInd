#' Plotting function for 'morpho_indices' outputs
#'
#' @param dat Data is the output from function morpho.indices().
#' @param name Character vector of names for ALL samples INCLUDING reference.
#' @param ref Index number (or name as character string) of reference sample. Default = 1, so the first sample is assumed to be the reference.
#' @param samp Numeric (or character) vector containing the index numbers (or names) of samples to be plotted. Defaults to all except reference sample.
#' @param var.names Character vector containing the names of ALL morphological descriptors. By default the column names of data frame from the "dat" argument.
#' @param which.var Numeric (or character) vector specifying the morphological descriptors to be plotted. By default all of them.
#' @param index Which morphological index will be plotted: "LSI" (default), "VSI" or "VSI2".
#' @param plot.ref Plot reference sample? Default = TRUE
#' @param plot.index Plot original values of LSI, VSI or VSI2. Default = TRUE.
#' @param bootmean Plot average bootstraped value of LSI, VSI or VSI2. Default = FALSE.
#' @param bootmedian Plot median bootstraped value of LSI, VSI or VSI2. Default = FALSE.
#' @param plot.mISD Should mISD values be plotted?
#' @param CI Plot 95 confidence intervals? Default = FALSE.
#' @param col.ref Color for reference sample. Default is "black".
#' @param col Vector of colors for samples to be plotted. Length must be equal to length(samp).
#' @param alpha.lev Opacity level (0 = transparent, 1 = opaque). Requires packages 'scales'.
#' @param lty.index, lty.bootmean, lty.bootmed Line types for index, bootstrap mean and bootstrap median, respectively.
#' @param lgnd Plot legend? Default=TRUE.
#' @param lgnd.pos Position of legend, can be specified as for function legend(). Default = 'topleft'.
#' @param cex.lgnd Size of legend. Default = 0.5.
#' @param xlab, ylab As for standard 'plot'.
#' @param tick.labels Labels for tickmarks.
#'
#' @export

plot_mo_ind <- function(dat=dat, 
                        name=name, 
                        ref=1, 
                        samp=name[-ref], 
                        var.names=colnames(dat$LSI.df), 
                        which.var=1:length(var.names), 
                        index="LSI", 
                        plot.ref=T, 
                        plot.index=T, 
                        bootmean=F, 
                        bootmedian=F, 
                        plot.mISD=T,
                        CI=F, 
                        col.ref="black", 
                        col=rainbow(length(samp)),
                        alpha.lev = 0.2,
                        lty.index=1,
                        lty.bootmean=2,
                        lty.bootmed=3, 
                        lgnd=T, 
                        lgnd.pos="topleft",
                        cex.lgnd=0.5, 
                        xlab="Descriptors", ylab=index, tick.labels=colnames(dfs))
                        {

dfs <- dat[[paste(index, "df", sep=".")]]
CIs <- dat[[paste(index, "CI", sep=".")]]
boot <- dat[[paste(index, "bootstrap", sep=".")]]
mISD <- dat[[paste(index, "mISD", sep=".")]]

rownames(dfs) <- names(CIs) <- names(boot) <- names(mISD) <- name
colnames(dfs) <- var.names

ref.name <- rownames(dfs[ref,])
samp.name <- rownames(dfs[samp,])

dfs <- as.matrix(dfs[, which.var])

N.s <- dim(dfs)[1]
N.d <- dim(dfs)[2]

   for (i in 1:length(CIs)) {
   colnames(CIs[[i]]) <- var.names
   colnames(boot[[i]]) <- var.names
   CIs[[i]] <- CIs[[i]][,which.var]
   boot[[i]] <- boot[[i]][,which.var]
   }

   if (bootmean) { # Can be deleted if added in morpho.indices function
   bootmeans <- matrix(NA, N.s, N.d)
      for (i in 1:N.s) {
      bootmeans[i,] <- apply(boot[[i]], 2, mean)
      }
   }

   if (bootmedian) { # Can be deleted if added in previous function
   bootmedians <- matrix(NA, N.s, N.d)
      for (i in 1:N.s) {
      bootmedians[i,] <- apply(boot[[i]], 2, quantile, 0.5)
      }
   }

   if (!CI) {

   plot(dfs[ref,], ylim = c(min(dfs[samp,]), max(dfs[samp,])), type="n", xlab=xlab, ylab=ylab, xaxt="n", las=2)
   axis(side=1, labels=tick.labels, at=1:N.d, las=2)

      for (i in 1:length(samp)) {

         if (plot.mISD) {abline(h=mISD[samp[i]], col=col[i], lty=4)}
         if (plot.index) {
         points(dfs[samp[i],], col=col[i], type="o", pch=20, lty=lty.index)
            if (plot.ref) {points(dfs[ref,], type="o", pch=20, col=col.ref, lty=lty.index)}
         }
         if (bootmean) {
         points(bootmeans[samp[i],], col=col[i], type="o", pch=20, lty=lty.bootmean)
            if (plot.ref) {points(bootmeans[ref,], type="o", pch=20, col=col.ref, lty=lty.bootmean)}
         }
         if (bootmedian) {
         points(bootmedians[samp[i],], col=col[i], type="o", pch=20, lty=lty.bootmed)
            if (plot.ref) {points(bootmedians[ref,], type="o", pch=20, col=col.ref, lty=lty.bootmed)}
         }
      }
   }

   if (CI) {

   plot(dfs[ref,], ylim = c(min(unlist(CIs[samp])), max(unlist(CIs[samp]))), type="n", xlab=xlab, ylab=ylab, xaxt="n", las=2)
   axis(side=1, labels=tick.labels, at=1:N.d, las=2)
   
      for (i in 1:length(samp)) {

      if (plot.mISD) {abline(h=mISD[samp[i]], col=col[i], lty=4)}
      if (plot.ref) {
      polygon(c(1:N.d, N.d:1), c(CIs[[ref]][1,], rev(CIs[[ref]][2,])), col=alpha(col.ref, alpha.lev), border=NA)
      }

      polygon(c(1:N.d, N.d:1), c(CIs[[samp[i]]][1,], rev(CIs[[samp[i]]][2,])), col=alpha(col[i], alpha.lev), border=NA)

         if (plot.index) {
         points(dfs[samp[i],], col=col[i], type="o", pch=20, lty=lty.index)
            if (plot.ref) {points(dfs[ref,], type="o", pch=20, col=col.ref, lty=lty.index)}
         }
         if (bootmean) {
         points(bootmeans[samp[i],], col=col[i], type="o", pch=20, lty=lty.bootmean)
            if (plot.ref) {points(bootmeans[ref,], type="o", pch=20, col=col.ref, lty=lty.bootmean)}
         }
         if (bootmedian) {
         points(bootmedians[samp[i],], col=col[i], type="o", pch=20, lty=lty.bootmed)
            if (plot.ref) {points(bootmedians[ref,], type="o", pch=20, col=col.ref, lty=lty.bootmed)}
         }
      }
   }

   if (plot.mISD) {
   leg <- c(samp.name, "mISD")
   lty.leg <- c(rep(1, length(samp.name)), 4)
   col <- c(col,"black")
   pchm <- c(rep(20, length(samp.name)), NA)
   } else {
   leg <- samp.name
   lty.leg <- rep(1, length(samp.name))
   pchm <- c(rep(20, length(samp.name)))
   }

   if (lgnd & plot.ref) {
   legend(lgnd.pos,
   legend=c(paste(ref.name, "(reference sample)"), leg),
   col=c(col.ref, col),
   pch=c(20,pchm),
   lty=c(1,lty.leg),
   bty="n",
   cex=cex.lgnd)
   }

   if (lgnd & !plot.ref) {
   legend(lgnd.pos,
   legend=leg,
   col=col,
   pch=pchm,
   lty=lty.leg,
   bty="n",
   cex=cex.lgnd)
   }


}

