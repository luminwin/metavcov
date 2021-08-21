\name{smd.vcov}
\alias{smd.vcov}
\title{Covariance matrix for standardized mean differences}
\author{
Min Lu
}
\usage{
smd.vcov(nt, nc, d, r, n_rt = 0, n_rc = 0)
}
\description{
 Compute variance-covariance matrix for multivariate meta-analysis when effect size is standardized mean difference.
}
\arguments{
\item{nt }{A matrix with sample sizes in the treatment group reporting each of the outcome. nt[i,j] is the sample size from study i reporting outcome j.}
\item{nc }{ Defined in a similar way as nt for control group.}
\item{d }{ A matrix with standardized mean differences from each of the outcome. d[i,j] is the value from study i for outcome j.}
\item{r }{ A list of correlation coefficient matrices of the outcomes from the studies. r[[k]][i,j] is the correlation coefficient between outcome i and outcome j from study k.}
\item{n_rt }{A list of matrices storing sample sizes in the treatment group reporting pairwised outcomes in the off diagonal elements. n_rt[[k]][i,j] is the sample size reporting both outcome i and outcome j from study k. Diagonal elements of these matrices are not used. The default value is zero, which means the smaller sample size reporting the corresponding two outcomes: i.e. n_rt[[k]][i,j]=min(nt[k,i],nt[k,j]). }
\item{n_rc }{ Defined in a similar way as n_rt for control group.}
}

\value{
    \tabular{ll}{
  list.mix.cov \tab A list of computed variance-covariance matrices.\cr
  mix.cov \tab A matrix whose rows are computed variance-covariance vectors.\cr}}


\references{
Ahn, S., Lu, M., Lefevor, G.T., Fedewa, A. & Celimli, S. (2016). Application of meta-analysis in sport and exercise science. In N. Ntoumanis, & N. Myers (Eds.), \emph{An Introduction to Intermediate and Advanced Statistical Analyses for Sport and Exercise Scientists} (pp.233-253). Hoboken, NJ: John Wiley and Sons, Ltd.

Wei, Y., &  Higgins, J. (2013). Estimating within study covariances in multivariate meta-analysis with multiple outcomes. \emph{Statistics in Medicine, 32}(7), 119-1205.
}

\examples{
######################################################
# Example: Geeganage2010 data
# Preparing covarianceS for multivariate meta-analysis
######################################################
data(Geeganage2010)
## set the correlation coefficients list r
r12 <- 0.71
r.Gee <- lapply(1:nrow(Geeganage2010), function(i){matrix(c(1, r12, r12, 1), 2, 2)})

computvocv <- smd.vcov(nt = subset(Geeganage2010, select = c(nt_SBP, nt_DBP)),
               nc = subset(Geeganage2010, select = c(nc_SBP, nc_DBP)),
               d = subset(Geeganage2010, select = c(SMD_SBP, SMD_DBP)), r = r.Gee)
# name variance-covariance matrix as covars
covars <- computvocv$smd.cov
#####################################################
# Running random-effects model using package "mvmeta"
#####################################################
#library(mvmeta)
#mvmeta_RE <- summary(mvmeta(cbind(SMD_SBP, SMD_DBP),
#                     S = covars,
#                     data = subset(Geeganage2010, select = c(SMD_SBP, SMD_DBP)),
#                     method = "reml"))
#mvmeta_RE
}
