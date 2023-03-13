\name{smd.vcov}
\alias{smd.vcov}
\title{Computing Variance-Covariance Matrices for Standardized Mean Differences}
\author{
Min Lu
}
\usage{
smd.vcov(nt, nc, d, r, n_rt = NA, n_rc = NA, name = NULL)
}
\description{
The function \code{lgOR.vcov} computes effect sizes and variance-covariance matrix for multivariate meta-analysis when the effect sizes of interest are all measured by standardized mean difference. See \code{\link{mix.vcov}} for effect sizes of the same or different types.
}
\arguments{
\item{nt }{ A \eqn{N \times p}{N x p} matrix storing sample sizes in the treatment group reporting the  \eqn{p} outcomes. \code{nt[i,j]} is the sample size from study \eqn{i} reporting outcome \eqn{j}.}
\item{nc }{ A matrix defined in a similar way as \code{nt} for the control group.}
\item{d }{ A \eqn{N \times p}{N x p} matrix or data frame with standard mean differences (SMD) from the \eqn{N} studies. \code{d[i,j]} is the value from study \eqn{i} for outcome \eqn{j}. }
\item{r }{ A \eqn{N}-dimensional list of \eqn{p \times p}{p x p} correlation matrices for the \eqn{p} outcomes from the \eqn{N} studies. \code{r[[k]][i,j]} is the correlation coefficient between outcome \eqn{i} and outcome \eqn{j} from study \eqn{k}.}
\item{n_rt }{A \eqn{N}-dimensional list of \eqn{p \times p}{p x p} matrices storing sample sizes in the treatment group reporting pairwise outcomes in the off-diagonal elements. \code{n_rt[[k]][i,j]} is the sample size reporting both outcome \eqn{i} and outcome \eqn{j} from study \eqn{k}. Diagonal elements of these matrices are discarded. The default value is \code{NA}, which means that the smaller sample size reporting the corresponding two outcomes is imputed: i.e. \code{n_rt[[k]][i,j]=min(nt[k,i],nt[k,j])}. }
\item{n_rc }{ A list defined in a similar way as \code{n_rt} for the control group.}
\item{name }{ Names for the outcomes.}
}

\value{
\item{ef}{ A \eqn{N \times p}{N x p} data frame that transforms the input argument \code{d} into Hedges's g (Wei and Higgins, 2013).}
\item{list.vcov}{ A \eqn{N}-dimensional list of \eqn{p(p+1)/2 \times p(p+1)/2}{p(p+1)/2 x p(p+1)/2} variance-covariance matrices for Hedges's g (Wei and Higgins, 2013).}
\item{matrix.vcov}{A \eqn{N \times p(p+1)/2}{N x p(p+1)/2} whose rows are computed variance-covariance vectors for Hedges's g (Wei and Higgins, 2013).}
\item{list.dvcov}{A \eqn{N}-dimensional list of \eqn{p(p+1)/2 \times p(p+1)/2}{p(p+1)/2 x p(p+1)/2} variance-covariance matrices for SMD (Olkin and Gleser, 2009).}
\item{matrix.dvcov}{A \eqn{N \times p(p+1)/2}{N x p(p+1)/2} matrix whose rows are computed variance-covariance vectors for SMD (Olkin and Gleser, 2009).}
  }


\references{
Ahn, S., Lu, M., Lefevor, G.T., Fedewa, A. & Celimli, S. (2016). Application of meta-analysis in sport and exercise science. In N. Ntoumanis, & N. Myers (Eds.), \emph{An Introduction to Intermediate and Advanced Statistical Analyses for Sport and Exercise Scientists} (pp.233-253). Hoboken, NJ: John Wiley and Sons, Ltd.

Wei, Y. &  Higgins, J. (2013). Estimating within study covariances in multivariate meta-analysis with multiple outcomes. \emph{Statistics in Medicine, 32}(7), 119-1205.

Olkin, I. & Gleser, L. (2009). Stochastically dependent effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (pp. 357-376). New York: Russel Sage Foundation.
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

computvcov <- smd.vcov(nt = subset(Geeganage2010, select = c(nt_SBP, nt_DBP)),
               nc = subset(Geeganage2010, select = c(nc_SBP, nc_DBP)),
               d = subset(Geeganage2010, select = c(SMD_SBP, SMD_DBP)), r = r.Gee,
               name = c("SMD_SBP", "SMD_DBP"))
# name variance-covariance matrix as S
S <- computvcov$matrix.vcov
## fixed-effect model
y <- computvcov$ef
MMA_FE <- summary(metafixed(y = y, Slist = computvcov$list.vcov))
#######################################################################
# Running random-effects model using package "mixmeta" or "metaSEM"
#######################################################################
# Restricted maximum likelihood (REML) estimator from the mixmeta package
#library(mixmeta)
#mvmeta_RE <- summary(mixmeta(cbind(SMD_SBP, SMD_DBP)~1,
#                     S = S,
#                     data = y,
#                     method = "reml"))
#mvmeta_RE

# maximum likelihood estimators from the metaSEM package
# library(metaSEM)
# metaSEM_RE <- summary(meta(y = y, v = S))
# metaSEM_RE
##############################################################
# Plotting the result:
##############################################################
obj <- MMA_FE
# obj <- mvmeta_RE
# obj <- metaSEM_RE
# pdf("CI.pdf", width = 4, height = 7)
plotCI(y = computvcov$ef, v = computvcov$list.vcov,
        name.y = NULL, name.study = Geeganage2010$studyID,
        y.all = obj$coefficients[,1],
        y.all.se = obj$coefficients[,2])
# dev.off()
}
