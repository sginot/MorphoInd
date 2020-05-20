#' Statistical tests for differences in LSI, VSI and VSI*
#'
#' This function is designed to globally test for differences between samples based on output from function morpho_indices. If this global test is significant, specific test can be used to find out which morphological variable are responsible.
#' @arguments x The output from the morpho_indices function.
#' @arugments ref Numeric string specifying which sample is the reference.


global_SI.test <- function(x, ref=1) { #Start of function

N <- ncol(x$LSI.df)
k <- nrow(x$LSI.df)

m.LSI <- apply(x$LSI.df, 1, mean)
m.VSI <- apply(x$VSI.df, 1, mean)
m.VSI2 <- apply(x$VSI2.df, 1, mean)

v.LSI <- apply(x$LSI.df, 1, var)
v.VSI <- apply(x$VSI.df, 1, var)
v.VSI2 <- apply(x$VSI2.df, 1, var)

t.LSI <- abs(m.LSI-m.LSI[ref]) / sqrt((v.LSI+v.LSI[ref])/N)
t.VSI <- abs(m.VSI-m.VSI[ref]) / sqrt((v.VSI+v.VSI[ref])/N)
t.VSI2 <- abs(m.VSI2-m.VSI2[ref]) / sqrt((v.VSI2+v.VSI2[ref])/N)

df.LSI <- (((v.LSI+v.LSI[ref])/N)^2) / ((v.LSI^2+v.LSI[ref]^2)/(N^3-N^2))
df.VSI <- (((v.VSI+v.VSI[ref])/N)^2) / ((v.VSI^2+v.VSI[ref]^2)/(N^3-N^2))
df.VSI2 <- (((v.VSI2+v.VSI2[ref])/N)^2) / ((v.VSI2^2+v.VSI2[ref]^2)/(N^3-N^2))

p.LSI <- 2*pt(-abs(t.LSI), df.LSI)
p.VSI <- 2*pt(-abs(t.VSI), df.VSI)
p.VSI2 <- 2*pt(-abs(t.VSI2), df.VSI2)

s.LSI <- s.VSI <- s.VSI2 <- rep("", k)
s.LSI[which(p.LSI <= 0.05)] <- "*"
s.VSI[which(p.VSI <= 0.05)] <- "*"
s.VSI2[which(p.VSI2 <= 0.05)] <- "*"

end.LSI <- data.frame(t=t.LSI[-ref], df=df.LSI[-ref], p=round(p.LSI[-ref],5), sig=s.LSI[-ref])

end.VSI <- data.frame(t=t.VSI[-ref], df=df.VSI[-ref], p=round(p.VSI[-ref],5), sig=s.VSI[-ref])

end.VSI2 <- data.frame(t=t.VSI2[-ref], df=df.VSI2[-ref], p=round(p.VSI2[-ref],5), sig=s.VSI2[-ref])

message("T-test statistics for differences in LSI, VSI and VSI2")
return(list(LSI.test=end.LSI, VSI.test=end.VSI, VSI2.test=end.VSI2))
} #End of function
