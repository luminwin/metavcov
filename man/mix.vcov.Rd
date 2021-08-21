\name{mix.vcov}
\alias{mix.vcov}
\title{Covariance matrix for mixed effect sizes }
\author{
Min Lu
}
\usage{
mix.vcov(d, r, nt, nc, st, sc, n_rt = 0, n_rc = 0, sdt, sdc, type)
}
\description{
Compute variance-covariance matrices between different effect sizes. Effect sizes include mean difference (MD), standardized mean difference (SMD), log odds ratio (logOR), log risk ratio (logRR), and risk difference (RD). Formulas are in Table I of Wei et al.'s paper.
}
\arguments{
\item{d }{ A matrix with standard mean differences from each of the outcome. d[i,j] is the value from study i for outcome j. If outcome j is not mean difference, NA has to be imputed in column j. }
\item{r }{ A list of correlation coefficient matrices of the outcomes from the studies. r[[k]][i,j] is the correlation coefficient between outcome i and outcome j from study k.}
\item{nt }{ A matrix with sample sizes in the treatment group reporting each of the outcome. nt[i,j] is the sample size from study i reporting outcome j.}
\item{nc }{ Defined in a similar way as nt for control group.}
\item{st }{ Number of participants with event for all outcomes (dichotomous) in treatment group. st[i,j] reports number of participants with event for outcome j in treatment group. If outcome j is not dichotomous, NA has to be imputed in column j.}
\item{sc }{ Defined in a similar way as st for control group.}
\item{n_rt }{A list of matrices storing sample sizes in the treatment group reporting pairwised outcomes in the off diagonal elements. n_rt[[k]][i,j] is the sample size reporting both outcome i and outcome j from study k. Diagonal elements of these matrices are not used. The default value is zero, which means the smaller sample size reporting the corresponding two outcomes: i.e. n_rt[[k]][i,j]=min(nt[k,i],nt[k,j]). }
\item{n_rc }{ Defined in a similar way as n_rt for control group.}
\item{sdt }{ Sample standard deviation from each of the outcome. sdt[i,j] is the sample standard deviation from study i for outcome j. If outcome j is not mean difference, NA has to be imputed in the jth column.}
\item{sdc }{ Defined in a similar way as sdt for control group.}
\item{type }{ A vector indicating types of effect sizes. "MD" stands for mean difference, "SMD" stands for standardized mean difference,"logOR" stands for log odds ratio, "logRR" stands for log risk ratio, and "RD" stands for risk difference. }
}
\value{
   \tabular{ll}{
  list.mix.cov \tab A list of computed variance-covariance matrices.\cr
  mix.cov \tab A matrix whose rows are computed variance-covariance vectors.\cr}
}

\references{
Wei, Y., &  Higgins, J. (2013). Estimating within study covariances in multivariate meta-analysis with multiple outcomes. \emph{Statistics in Medicine, 32}(7), 119-1205.
}

\examples{
#########################################################################
# Example: Geeganage2010 data
# Preparing covariances for multivariate meta analysis
# Choose variable SBP, DBP, DD, D with effect sizes "MD","MD","RD","lgOR"
#########################################################################
data(Geeganage2010)
## set the correlation coefficients list r
r12 <- 0.71
r13 <- 0.5
r14 <- 0.25
r23 <- 0.6
r24 <- 0.16
r34 <- 0.16
r <- vec2sm(c(r12, r13, r14, r23, r24, r34), diag = FALSE)
diag(r) <- 1
mix.r <- lapply(1:nrow(Geeganage2010), function(i){r})
attach(Geeganage2010)
## compute variance co-variances
computvocv <- mix.vcov(type = c("MD", "MD", "RD", "lgOR"),
                     d = cbind(MD_SBP, MD_DBP, NA, NA),
                     sdt = cbind(sdt_SBP, sdt_DBP, NA, NA),
                     sdc = cbind(sdc_SBP, sdc_DBP, NA, NA),
                     nt = cbind(nt_SBP, nt_DBP, nt_DD, nt_D),
                     nc = cbind(nc_SBP, nc_DBP, nc_DD, nc_D),
                     st = cbind(NA, NA, st_DD, st_D),
                     sc = cbind(NA, NA, sc_DD, sc_D),
                     r = mix.r)
# name different effect sizes as an input
Input <- subset(Geeganage2010, select = c(MD_SBP, MD_DBP))
Input$RD_DD <- st_DD/nt_DD - sc_DD/nc_DD
Input$lgOR_D <- log((st_D/(nt_D - st_D))/(sc_D/(nc_D - sc_D)))
# name variance-covariance matrix as covars
covars <- computvocv$mix.cov
###############################################
# Running random-effects model using package "mvmeta"
###############################################
#library(mvmeta)
#mvmeta_RE <- summary(mvmeta(cbind(MD_SBP, MD_DBP, RD_DD, lgOR_D),
#                           S = covars, data = Input, method = "reml"))
#mvmeta_RE
}
