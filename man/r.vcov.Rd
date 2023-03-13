\name{r.vcov}
\alias{r.vcov}
\alias{rToz}
\alias{smTovec}
\title{Computing Variance-Covariance Matrices for Correlation Coefficients}
\author{
Min Lu
}
\usage{
r.vcov(n, corflat, zscore = FALSE, name = NULL, method = "average", na.impute = NA)
}
\description{
The function \code{r.vcov} computes variance-covariance matrix for multivariate meta-analysis when the effect size is measured by correlation coefficient.
}
\arguments{
\item{n }{ A \eqn{N}-dimensional vector containing sample sizes from \eqn{N} studies.}
\item{corflat }{ A \eqn{N \times p}{N x p} matrix or data frame storing \eqn{p} correlation coefficients from each of the \eqn{N} studies.}
\item{zscore }{ Whether the correlation coefficients in \code{corflat} are already transformed into Fisher's \eqn{z} scores.}
\item{name }{ A \eqn{p}-dimensional vector containing names for the \eqn{p} correlation coefficients.}
\item{method }{ Method \code{"average"} computes variance and covariances with mean correlation coefficients that are sample-size weighted from all the \eqn{N} studies (missing values are omitted); method \code{"each"} computes variance and covariances with each of the corresponding correlation coefficients.}
\item{na.impute }{ Missing values can be imputed by a numeric value, such as zero by setting \code{na.impute = 0}. With the default setting \code{na.impute = NA}, missing values are not imputed. If specifying \code{na.impute = "average"}, missing values are imputed by mean correlation coefficients that are sample-size weighted from the complete records. }
}

\details{
How to arrange correlation coefficients of each study from matrix to vector is in Cooper et al book page 385 to 386. Details for average method are in book of Cooper et al page 388. Let \eqn{r_{ist}} denote the sample correlation coefficient that describes the relationship between variables \eqn{s} and \eqn{t} in study \eqn{i}. We can calculate its variance as
\eqn{var({{r}_{ist}})={{(1-\rho_{ist}^{2})}^{2}}/{{n}_{i}},} and the covariance between two correlation coefficients is
\eqn{
	cov({{r}_{ist}},{{r}_{iuv}})=[.5{{\rho }_{ist}}{{\rho }_{iuv}}({{\rho }}_{isu}^{2}+{{\rho }}_{isv}^{2}+{{\rho }}_{itu}^{2}+{{\rho }}_{itv}^{2})+{{\rho }_{isu}}{{\rho }_{itv}}+{{\rho }_{isv}}{{\rho }_{itu}} \nonumber \\
	 -({{\rho }_{ist}}{{\rho }_{isu}}{{\rho }_{isv}}+{{\rho }_{its}}{{\rho }_{itu}}{{\rho }_{itv}}+{{\rho }_{ius}}{{\rho }_{iut}}{{\rho }_{iuv}}+{{\rho }_{ivs}}{{\rho }_{ivt}}{{\rho }_{ivu}})]/{{n}_{i}},
}
where \eqn{\rho_{i..}} represents the population value. In practice, \eqn{\rho_{i..}} can be substituted by the observed sample correlation or sample-size weighted mean correlation coefficients from all studies.
}
\value{
\item{r }{ A \eqn{N \times p}{N x p} data frame that contains the input argument \code{corflat} with column names and imputed values according to the input argument \code{na.impute}. If the input argument \code{zscore=TRUE}, \code{r} is transformed from Fisher's \eqn{z} score in \code{corflat}. }
\item{list.rvcov }{ A \eqn{N}-dimensional list of \eqn{p(p+1)/2 \times p(p+1)/2}{p(p+1)/2 x p(p+1)/2} matrices of computed variance-covariance matrices.}
\item{matrix.rvcov }{ A \eqn{N \times p(p+1)/2}{N x p(p+1)/2} matrix whose rows are computed variance-covariance vectors.}
\item{ef }{ A \eqn{N \times p}{N x p} data frame that contains Fisher's \eqn{z} transformed correlation coefficients from the input argument \code{corflat} with column names and imputed values according to the input argument \code{na.impute}. }
\item{list.vcov }{ A list in the same format of \code{list.rvcov} for Fisher's \eqn{z} transformed correlation coefficients.}
\item{matrix.vcov }{ A matrix \code{matrix.rvcov} for Fisher's \eqn{z} transformed correlation coefficients.}}

\references{
Ahn, S., Lu, M., Lefevor, G.T., Fedewa, A. & Celimli, S. (2016). Application of meta-analysis in sport and exercise science. In N. Ntoumanis, & N. Myers (Eds.), \emph{An Introduction to Intermediate and Advanced Statistical Analyses for Sport and Exercise Scientists} (pp.233-253). Hoboken, NJ: John Wiley and Sons, Ltd.

Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.) (2009). \emph{The handbook of research synthesis and
meta-analysis.} New York: Russell Sage Foundation.

Olkin, I., & Ishii, G. (1976). Asymptotic distribution of functions of a correlation matrix. In S. Ikeda (Ed.),
\emph{Essays in probability and statistics: A volume in honor of Professor Junjiro Ogawa} (pp.5-51). Tokyo,
Japan: Shinko Tsusho.
}

\examples{
##############################################################
# A simple example:
# Checking the example in Harris Cooper et al.'s book page 388
##############################################################
r <- matrix(c(-0.074, -0.127, 0.324, 0.523, -0.416, -0.414), 1)
n <- 142
computvcov <- r.vcov(n = n, corflat = r,
                     name = paste("C", c("st", "su", "sv", "tu", "tv", "uv"), sep = ""),
                     method = "each")
round(computvcov$list.rvcov[[1]], 4)
round(computvcov$ef, 4)
round(computvcov$list.vcov[[1]], 4)
round(computvcov$matrix.vcov, 4)
######################################################
# Example: Craft2003 data
# Preparing covariances for multivariate meta-analysis
######################################################
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
#######################################################################
# Running random-effects model using package "mixmeta" or "metaSEM"
#######################################################################
# Restricted maximum likelihood (REML) estimator from the mixmeta package
#library(mixmeta)
#mvmeta_RE <- summary(mixmeta(cbind(C1, C2, C3, C4, C5, C6) ~ 1,
#                           S = S, data = y, method = "reml"))
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
        name.y = NULL, name.study = Craft2003$ID,
        y.all = obj$coefficients[,1],
        y.all.se = obj$coefficients[,2])
# dev.off()
}
