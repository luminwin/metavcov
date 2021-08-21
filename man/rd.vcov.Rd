\name{rd.vcov}
\alias{rd.vcov}
\title{Covariance matrix for risk differences}
\author{
Min Lu
}
\usage{
rd.vcov(r, nt, nc, st, sc, n_rt = 0, n_rc = 0)
}
\description{
 Compute variance-covariance matrix for multivariate meta-analysis when effect size is risk difference.
}
\arguments{
\item{r }{A list of correlation coefficient matrices of the outcomes from the studies. r[[k]][i,j] is the correlation coefficient between outcome i and outcome j from study k.}
\item{nt }{ Sample sizes in the treatment group reporting each of the outcome. nt[i,j] means the sample size from study i reporting outcome j.}
\item{nc }{ Defined in a similar way as nt for control group.}
\item{st }{ Number of participants with event for all outcomes (dichotomous) in treatment group. st[i,j] is number of participants with event for outcome j in treatment group.}
\item{sc }{ Defined in a similar way as st for control group.}
\item{n_rt }{A list of matrices storing sample sizes in the treatment group reporting pairwised outcomes in the off diagonal elements. n_rt[[k]][i,j] means the sample size reporting both outcome i and outcome j from study k. Diagonal elements of these matrices are not used. The default value is zero, which means the smaller sample size reporting the corresponding two outcomes: i.e. n_rt[[k]][i,j]=min(nt[k,i],nt[k,j]). }
\item{n_rc }{ Defined in a similar way as n_rt for control group.}
}

\value{
    \tabular{ll}{
  rd \tab Computed risk difference from input.\cr
  list.lgOR.cov \tab A list of computed variance-covariance matrices.\cr
  lgOR.cov \tab A matrix whose rows are computed variance-covariance vectors.\cr}}


\references{
Ahn, S., Lu, M., Lefevor, G.T., Fedewa, A. & Celimli, S. (2016). Application of meta-analysis in sport and exercise science. In N. Ntoumanis, & N. Myers (Eds.), \emph{An Introduction to Intermediate and Advanced Statistical Analyses for Sport and Exercise Scientists} (pp.233-253). Hoboken, NJ: John Wiley and Sons, Ltd.

Wei, Y., &  Higgins, J. (2013). Estimating within study covariances in multivariate meta-analysis with multiple outcomes. \emph{Statistics in Medicine, 32}(7), 119-1205.
}

\examples{
###########################################################################
# Example: Geeganage2010 data
# Preparing risk differences and covariances for multivariate meta-analysis
###########################################################################
data(Geeganage2010)
## set the correlation coefficients list r
r12 <- 0.71
r.Gee <- lapply(1:nrow(Geeganage2010), function(i){matrix(c(1, r12, r12, 1), 2, 2)})


computvocv <- rd.vcov(nt = subset(Geeganage2010, select = c(nt_DD, nt_D)),
                      nc = subset(Geeganage2010, select = c(nc_DD, nc_D)),
                      st = subset(Geeganage2010, select = c(st_DD, st_D)),
                      sc = subset(Geeganage2010, select = c(sc_DD, sc_D)),
                      r = r.Gee)
# name computed relative risk as an input
Input <- computvocv$rd
colnames(Input) <- c("rd.DD", "rd.D")
# name variance-covariance matrix of trnasformed z scores as covars
covars <- computvocv$rd.cov
###############################################################
# Running random-effects model using package "mvmeta"
################################################################
#library(mvmeta)
#mvmeta_RE <- summary(mvmeta(cbind(rd.DD, rd.D),
#                          S = covars, data = as.data.frame(Input),
#                          method = "reml"))
#mvmeta_RE
}
