\name{lgRR.vcov}
\alias{lgRR.vcov}
\title{Computing Variance-Covariance Matrices for Log Risk Ratios}
\author{
Min Lu
}
\usage{
lgRR.vcov(r, nt, nc, st, sc, n_rt = NA, n_rc = NA)
}
\description{
The function \code{lgOR.vcov} computes effect sizes and variance-covariance matrix for multivariate meta-analysis when the effect sizes of interest are all measured by log risk ratio (or log relative risk). See \code{\link{mix.vcov}} for effect sizes of the same or different types.
}
\arguments{
\item{r}{ A \eqn{N}-dimensional list of \eqn{p \times p}{p x p} correlation matrices for the \eqn{p} outcomes from the \eqn{N} studies. \code{r[[k]][i,j]} is the correlation coefficient between outcome \eqn{i} and outcome \eqn{j} from study \eqn{k}.}
\item{nt }{ A \eqn{N \times p}{N x p} matrix storing sample sizes in the treatment group reporting the  \eqn{p} outcomes. \code{nt[i,j]} is the sample size from study \eqn{i} reporting outcome \eqn{j}.}
\item{nc }{ A matrix defined in a similar way as \code{nt} for the control group.}
\item{st }{ A \eqn{N \times p}{N x p} matrix recording number of participants with event for all outcomes (dichotomous) in treatment group. \code{st[i,j]} reports number of participants with event for outcome \eqn{j} in treatment group for study \eqn{i}. If outcome \eqn{j} is not dichotomous, NA has to be imputed in column \eqn{j}.}
\item{sc }{ Defined in a similar way as \code{st} for the control group.}
\item{n_rt }{A \eqn{N}-dimensional list of \eqn{p \times p}{p x p} matrices storing sample sizes in the treatment group reporting pairwise outcomes in the off-diagonal elements. \code{n_rt[[k]][i,j]} is the sample size reporting both outcome \eqn{i} and outcome \eqn{j} from study \eqn{k}. Diagonal elements of these matrices are discarded. The default value is \code{NA}, which means that the smaller sample size reporting the corresponding two outcomes is imputed: i.e. \code{n_rt[[k]][i,j]=min(nt[k,i],nt[k,j])}. }
\item{n_rc }{ A list defined in a similar way as \code{n_rt} for the control group.}
}

\value{
\item{ ef}{A \eqn{N \times p}{N x p} data frame whose columns are computed log risk ratios.}
\item{list.vcov }{ A \eqn{N}-dimensional list of \eqn{p(p+1)/2 \times p(p+1)/2}{p(p+1)/2 x p(p+1)/2} matrices of computed variance-covariance matrices.}
\item{matrix.vcov }{ A \eqn{N \times p(p+1)/2}{N x p(p+1)/2} matrix whose rows are computed variance-covariance vectors.}
  }


\references{
Ahn, S., Lu, M., Lefevor, G.T., Fedewa, A. & Celimli, S. (2016). Application of meta-analysis in sport and exercise science. In N. Ntoumanis, & N. Myers (Eds.), \emph{An Introduction to Intermediate and Advanced Statistical Analyses for Sport and Exercise Scientists} (pp.233-253). Hoboken, NJ: John Wiley and Sons, Ltd.

Wei, Y., &  Higgins, J. (2013). Estimating within study covariances in multivariate meta-analysis with multiple outcomes. \emph{Statistics in Medicine, 32}(7), 119-1205.
}

\examples{
##########################################################################
# Example: Geeganage2010 data
# Preparing log risk ratios and covariances for multivariate meta-analysis
##########################################################################
data(Geeganage2010)
## set the correlation coefficients list r
r12 <- 0.71
r.Gee <- lapply(1:nrow(Geeganage2010), function(i){matrix(c(1, r12, r12, 1), 2, 2)})


computvcov <- lgRR.vcov(nt = subset(Geeganage2010, select = c(nt_DD, nt_D)),
                      nc = subset(Geeganage2010, select = c(nc_DD, nc_D)),
                      st = subset(Geeganage2010, select=c(st_DD, st_D)),
                      sc = subset(Geeganage2010, select=c(sc_DD, sc_D)),
                      r = r.Gee)
# name computed log risk ratio as y
y <- computvcov$ef
colnames(y) = c("lgRR.DD", "lgRR.D")
# name variance-covariance matrix of trnasformed z scores as covars
S <- computvcov$matrix.vcov
## fixed-effect model
MMA_FE <- summary(metafixed(y = y, Slist = computvcov$list.vcov))
MMA_FE
#######################################################################
# Running random-effects model using package "mixmeta" or "metaSEM"
#######################################################################
#library(mixmeta)
#mvmeta_RE = summary(mixmeta(cbind(lgRR.DD, lgRR.D)~1,
#                            S = S, data = as.data.frame(y),
#                            method = "reml"))
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
        name.y = c("lgRR.DD", "lgRR.D"),
        name.study = Geeganage2010$studyID,
        y.all = obj$coefficients[,1],
        y.all.se = obj$coefficients[,2],
        hline = 1)
# dev.off()
}
