\name{metavcov-package}
\alias{metavcov-package}
\alias{metavcov}
\docType{package}
\title{Computing Variances and Covariances, Visualization and Missing Data Solution for Multivariate Meta Analysis}
\description{
The package \code{metavcov} computes variances and covariances for effect sizes as preparations for multivariate meta-analysis. Effect sizes include correlation (r), mean difference (MD), standardized mean difference (SMD), log odds ratio (logOR), log risk ratio (logRR), and risk difference (RD). Functions for plotting confidence intervals and multiple imputation for missing data are offered. It can fit a fixed-effect model for demonstration purposes; users are highly encouraged to use packages \code{mixmeta} and \code{metaSEM} for random-effect models. }
\author{
Min Lu (Maintainer,<m.lu6@umiami.edu>)
}

\references{
Ahn, S., Lu, M., Lefevor, G.T., Fedewa, A. & Celimli, S. (2016). Application of meta-analysis in sport and exercise science. In N. Ntoumanis, & N. Myers (Eds.), \emph{An Introduction to Intermediate and Advanced Statistical Analyses for Sport and Exercise Scientists} (pp.233-253). Hoboken, NJ: John Wiley and Sons, Ltd.

Olkin, I., & Ishii, G. (1976). Asymptotic distribution of functions of a correlation matrix. In S. Ikeda (Ed.), \emph{Essays in probability and statistics: A volume in honor of Professor Junjiro Ogawa} (pp.5-51). Tokyo, Japan: Shinko Tsusho.

B. J. Becker. (2009) Model-based meta-analysis. In H. Cooper, L. V. Hedges, and J. C. Valentine, (Ed.), \emph{The handbook of research synthesis and meta-analysis}, chapter 20, pages 377-395. Russell Sage Foundation.

Wei, Y., &  Higgins, J. (2013). Estimating within study covariances in multivariate meta-analysis with multiple outcomes. \emph{Statistics in Medicine, 32}(7), 119-1205.
}

\examples{
\donttest{
###############################################
# Effect size: correlation coefficients
###############################################
data(Craft2003)
# extract correlation from the dataset (craft)
corflat <- subset(Craft2003, select=C1:C6)
# transform correlations to z and compute variance-covariance matrix.
computvcov <- r.vcov(n = Craft2003$N, corflat = corflat, method = "average")
# name transformed z scores as y
y <- computvcov$ef
# name variance-covariance matrix of trnasformed z scores as S
S <- computvcov$matrix.vcov
S[1, ]
## fixed-effect model
MMA_FE <- summary(metafixed(y = y, Slist = computvcov$list.vcov))
MMA_FE
# Restricted maximum likelihood (REML) estimator from the mixmeta package
#library(mixmeta)
#mvmeta_RE <- summary(mixmeta(cbind(C1, C2, C3, C4, C5, C6)~1,
#                           S = S, data = y, method = "reml"))
#mvmeta_RE

# maximum likelihood estimators from the metaSEM package
# library(metaSEM)
# metaSEM_RE <- summary(meta(y = y, v = S))
# metaSEM_RE
# Plotting the result:
obj <- MMA_FE
# obj <- mvmeta_RE
# obj <- metaSEM_RE
# pdf("CI.pdf", width = 4, height = 7)
plotCI(y = computvcov$ef, v = computvcov$list.vcov,
        name.y = NULL, name.study = Craft2003$ID,
        y.all = obj$coefficients[,1],
        y.all.se = obj$coefficients[,2])
# dev.off()
#########################################################################
# Other effect sizes of the same or different type
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
r <- vecTosm(c(r12, r13, r14, r23, r24, r34))
diag(r) <- 1
mix.r <- lapply(1:nrow(Geeganage2010), function(i){r})
attach(Geeganage2010)
## compute variance co-variances
computvcov <- mix.vcov(type = c("MD", "MD", "RD", "lgOR"),
                       d = cbind(MD_SBP, MD_DBP, NA, NA),
                       sdt = cbind(sdt_SBP, sdt_DBP, NA, NA),
                       sdc = cbind(sdc_SBP, sdc_DBP, NA, NA),
                       nt = cbind(nt_SBP, nt_DBP, nt_DD, nt_D),
                       nc = cbind(nc_SBP, nc_DBP, nc_DD, nc_D),
                       st = cbind(NA, NA, st_DD, st_D),
                       sc = cbind(NA, NA, sc_DD, sc_D),
                       r = mix.r,
                       name = c("MD.SBP", "MD.DBP", "RD.DD", "lgOR.D"))
# save different effect sizes in y
y <- computvcov$ef
head(y)
# save variances and covariances of all the effect sizes in a matrix S
S <- computvcov$matrix.vcov
S[1, ]
## fixed-effect model
MMA_FE <- summary(metafixed(y = y, Slist = computvcov$list.vcov))
# Restricted maximum likelihood (REML) estimator from the mixmeta package
# library(mixmeta)
# mvmeta_RE <- summary(mixmeta(cbind(MD.SBP, MD.DBP, RD.DD, lgOR.D) ~1,
#                            S = S, data = y, method = "reml"))
# mvmeta_RE

# maximum likelihood estimators from the metaSEM package
# library(metaSEM)
# metaSEM_RE <- summary(meta(y = y, v = S))
# metaSEM_RE
# Plotting the result:
obj <- MMA_FE
# obj <- mvmeta_RE
# obj <- metaSEM_RE
# pdf("CI.pdf", width = 4, height = 7)
plotCI(y = computvcov$ef, v = computvcov$list.vcov,
        name.y = NULL, name.study = Geeganage2010$studyID,
        y.all = obj$coefficients[,1],
        y.all.se = obj$coefficients[,2],
        hline = c(0, 0, 0, 1))
# dev.off()
#####################################################################################
# Multiple Imputation for missing data
#####################################################################################
# prepare a dataset with missing values and input arguments for meta.mi
Craft2003.mnar <- Craft2003[, c(2, 4:10)]
Craft2003.mnar[sample(which(Craft2003$C4 < 0), 6), "C4"] <- NA
dat <- Craft2003.mnar
n.name <- "N"
ef.name <- c("C1", "C2", "C3", "C4", "C5", "C6")
# fixed-effect model
obj <- metami(dat, M = 2, vcov = "r.vcov",
             n.name, ef.name,
             func = "metafixed")
# Plotting the result
computvcov <- r.vcov(n = Craft2003$N,
                     corflat = subset(Craft2003.mnar, select = C1:C6),
                     method = "average")
plotCI(y = computvcov$ef, v = computvcov$list.vcov,
        name.y = NULL, name.study = Craft2003$ID,
        y.all = obj$coefficients[,1],
        y.all.se = obj$coefficients[,2])
}}
