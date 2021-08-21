\name{metavcov-package}
\alias{metavcov-package}
\alias{metavcov}
\docType{package}
\title{Variances and Covariances for Multivariate Meta-Analysis}
\description{
R package metavov computes variances and covariances for multivariate meta-analysis. Effect sizes include correlation (r), mean difference (MD), standardized mean difference (SMD), log odds ratio (logOR), log risk ratio (logRR), and risk difference (RD).}
\author{
Min Lu (Maintainer,<m.lu6@umiami.edu>)
}

\references{
Ahn, S., Lu, M., Lefevor, G.T., Fedewa, A. & Celimli, S. (2016). Application of meta-analysis in sport and exercise science. In N. Ntoumanis, & N. Myers (Eds.), \emph{An Introduction to Intermediate and Advanced Statistical Analyses for Sport and Exercise Scientists} (pp.233-253). Hoboken, NJ: John Wiley and Sons, Ltd.

Olkin, I., & Ishii, G. (1976). Asymptotic distribution of functions of a correlation matrix. In S. Ikeda (Ed.), \emph{Essays in probability and statistics: A volume in honor of Professor Junjiro Ogawa} (pp.5-51). Tokyo, Japan: Shinko Tsusho.

Wei, Y., &  Higgins, J. (2013). Estimating within study covariances in multivariate meta-analysis with multiple outcomes. \emph{Statistics in Medicine, 32}(7), 119-1205.
}

\examples{
###############################################
# Effect size : correlation coefficients
###############################################
data(Craft2003)
computvocv <- r.vcov(n = Craft2003$N,
                  corflat = subset(Craft2003, select = C1:C6),
                  method = "average")
# name transformed z scores as an input
Input <- computvocv$zr
# name variance covariance matrix of trnasformed z scores as covars
covars <- computvocv$zcov
# Next step: Overall analysis: Running random effects model
#                              using package "mvmeta"
#library(mvmeta)
#mvmeta_RE <- summary(mvmeta(cbind(C1, C2, C3, C4, C5, C6),
#                           S = covars, data = Input, method = "reml"))
#mvmeta_RE
}
