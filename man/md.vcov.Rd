\name{md.vcov}
\alias{md.vcov}
\title{Covariance matrix for mean differences}
\author{
Min Lu
}
\usage{
md.vcov(r, nt, nc, n_rt = 0, n_rc = 0, sdt, sdc)
}
\description{
Compute variance-covariance matrix for multivariate meta-analysis when effect size is mean difference.
}
\arguments{
\item{r }{A list of correlation coefficient matrices of the outcomes from the studies. r[[k]][i,j] is the correlation coefficient between outcome i and outcome j from study k.}
\item{nt }{A matrix with sample sizes in the treatment group reporting each of the outcome. nt[i,j] is the sample size from study i reporting the outcome j.}
\item{nc }{ Defined in a similar way as nt for control group.}
\item{n_rt }{A list of matrices storing sample sizes in the treatment group reporting pairwised outcomes in the off diagonal elements. n_rt[[k]][i,j] is the sample size reporting both outcome i and outcome j from study k. Diagonal elements of these matrices are not used. The default value is zero, which means the smaller sample size reporting the corresponding two outcomes: i.e. n_rt[[k]][i,j]=min(nt[k,i],nt[k,j]). }
\item{n_rc }{ Defined in a similar way as n_rt for control group.}
\item{sdt }{ Sample standard deviation from each of the outcome. sdt[i,j] is the sample standard deviation from study i for outcome j.}
\item{sdc }{ Defined in a similar way as sdt for control group.}
}

\value{
    \tabular{ll}{
  list.md.cov \tab A list of computed variance-covariance matrices.\cr
  md.cov \tab A matrix whose rows are computed variance-covariance vectors.\cr}}


\references{
Ahn, S., Lu, M., Lefevor, G.T., Fedewa, A. & Celimli, S. (2016). Application of meta-analysis in sport and exercise science. In N. Ntoumanis, & N. Myers (Eds.), \emph{An Introduction to Intermediate and Advanced Statistical Analyses for Sport and Exercise Scientists} (pp.233-253). Hoboken, NJ: John Wiley and Sons, Ltd.

Wei, Y., &  Higgins, J. (2013). Estimating within study covariances in multivariate meta-analysis with multiple outcomes. \emph{Statistics in Medicine, 32}(7), 119-1205.
}

\examples{
######################################################
# Example: Geeganage2010 data
# Preparing covariances for multivariate meta-analysis
######################################################
## set the correlation coefficients list r
r12 <- 0.71
r.Gee <- lapply(1:nrow(Geeganage2010), function(i){matrix(c(1, r12, r12, 1), 2, 2)})

computvocv <- md.vcov(nt = subset(Geeganage2010, select = c(nt_SBP, nt_DBP)),
                    nc = subset(Geeganage2010, select = c(nc_SBP, nc_DBP)),
                    sdt = subset(Geeganage2010, select=c(sdt_SBP, sdt_DBP)),
                    sdc = subset(Geeganage2010, select=c(sdc_SBP, sdc_DBP)),
                    r = r.Gee)
# name variance-covariance matrix as covars
covars <- computvocv$md.cov

#####################################################
# Running random-effects model using package "mvmeta"
#####################################################
#library(mvmeta)
#mvmeta_RE <- summary(mvmeta(cbind(MD_SBP, MD_DBP), S = covars,
#                         data = subset(Geeganage2010, select = c(MD_SBP, MD_DBP)),
#                         method = "reml"))
#mvmeta_RE
}
