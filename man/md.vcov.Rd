\name{md.vcov}
\alias{md.vcov}
\title{Computing Variance-Covariance Matrices for Mean Differences}
\author{
Min Lu
}
\usage{
md.vcov(r, nt, nc, n_rt = NA, n_rc = NA, sdt, sdc)
}
\description{
The function \code{md.vcov} computes effect sizes and variance-covariance matrix for multivariate meta-analysis when the effect sizes of interest are all measured by mean difference. See \code{\link{mix.vcov}} for effect sizes of the same or different types.
}
\arguments{
\item{r }{ A \eqn{N}-dimensional list of \eqn{p \times p}{p x p} correlation matrices for the \eqn{p} outcomes from the \eqn{N} studies. \code{r[[k]][i,j]} is the correlation coefficient between outcome \eqn{i} and outcome \eqn{j} from study \eqn{k}.}
\item{nt }{ A \eqn{N \times p}{N x p} matrix storing sample sizes in the treatment group reporting the  \eqn{p} outcomes. \code{nt[i,j]} is the sample size from study \eqn{i} reporting outcome \eqn{j}.}
\item{nc }{ A matrix defined in a similar way as \code{nt} for the control group.}
\item{n_rt }{A \eqn{N}-dimensional list of \eqn{p \times p}{p x p} matrices storing sample sizes in the treatment group reporting pairwise outcomes in the off-diagonal elements. \code{n_rt[[k]][i,j]} is the sample size reporting both outcome \eqn{i} and outcome \eqn{j} from study \eqn{k}. Diagonal elements of these matrices are discarded. The default value is \code{NA}, which means that the smaller sample size reporting the corresponding two outcomes is imputed: i.e. \code{n_rt[[k]][i,j]=min(nt[k,i],nt[k,j])}. }
\item{n_rc }{ A list defined in a similar way as \code{n_rt} for the control group.}
\item{sdt }{ A \eqn{N \times p}{N x p} matrix storing sample standard deviations for each outcome from treatment group. \code{sdt[i,j]} is the sample standard deviation from study \eqn{i} for outcome \eqn{j}. If outcome \eqn{j} is not continuous such as MD or SMD, NA has to be imputed in the \eqn{j}th column.}
\item{sdc }{ A matrix defined in a similar way as \code{sdt} for the control group.}
}

\value{
\item{list.vcov }{ A \eqn{N}-dimensional list of \eqn{p(p+1)/2 \times p(p+1)/2}{p(p+1)/2 x p(p+1)/2} matrices of computed variance-covariance matrices.}
\item{matrix.vcov }{ A \eqn{N \times p(p+1)/2}{N x p(p+1)/2} matrix whose rows are computed variance-covariance vectors.}}


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

computvcov <- md.vcov(nt = subset(Geeganage2010, select = c(nt_SBP, nt_DBP)),
                    nc = subset(Geeganage2010, select = c(nc_SBP, nc_DBP)),
                    sdt = subset(Geeganage2010, select=c(sdt_SBP, sdt_DBP)),
                    sdc = subset(Geeganage2010, select=c(sdc_SBP, sdc_DBP)),
                    r = r.Gee)
# name variance-covariance matrix as S
S <- computvcov$matrix.vcov
## fixed-effect model
y <- as.data.frame(subset(Geeganage2010, select = c(MD_SBP, MD_DBP)))
MMA_FE <- summary(metafixed(y = y, Slist = computvcov$list.vcov))
MMA_FE
#######################################################################
# Running random-effects model using package "mixmeta" or "metaSEM"
#######################################################################
# Restricted maximum likelihood (REML) estimator from the mixmeta package
#library(mixmeta)
#mvmeta_RE <- summary(mixmeta(cbind(MD_SBP, MD_DBP)~1, S = S,
#                         data = subset(Geeganage2010, select = c(MD_SBP, MD_DBP)),
#                         method = "reml"))
#mvmeta_RE

# maximum likelihood estimators from the metaSEM package
# library(metaSEM)
# metaSEM_RE <- summary(meta(y = y, v = S))
# metaSEM_RE
##############################################################
# Plotting the result:
##############################################################
# obj <- MMA_FE
# obj <- mvmeta_RE
# obj <- metaSEM_RE

# plotCI(y = y, v = computvcov$list.vcov,
#         name.y = c("MD_SBP", "MD_DBP"), name.study = Geeganage2010$studyID,
#         y.all = obj$coefficients[,1],
#         y.all.se = obj$coefficients[,2])
}
