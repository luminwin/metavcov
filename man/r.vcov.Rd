\name{r.vcov}
\alias{r.vcov}
\title{Covariance matrix for correlation coefficients}
\author{
Min Lu
}
\usage{
r.vcov(n, corflat, method = "average")
}
\description{
 Compute variance-covariance matrix for multivariate meta-analysis when effect size is correlation coefficient.
}
\arguments{
\item{n }{ Sample sizes from studies.}
\item{corflat }{ Correlation coefficients from studies.}
\item{method }{ Method "average" computes variance covariances with sample-size weighted mean correlation coefficients from all studies; method "each" computes variance covariances with each of the corresponding correlation coefficients.}
}

\details{
  How to arrange correlation coefficients of each study from matrix to vector is in Cooper et al book page 385 to 386. Details for average method are in book of Cooper et al page 388.
}
\value{
\item{list.rcov }{ A list of computed Variance-covariance matrices.}
\item{rcov }{ A matrix whose rows are computed Variance-covariance vectors.}
\item{zr }{ Z transformed correlation coefficients from input "corflat".}
\item{list.rcov }{ A list of computed Variance-covariance matrices from z transformed correlation coefficients.}
\item{zcov }{ A matrix whose rows are computed Variance-covariance vectors from z transformed correlation coefficients.}}

\references{
Ahn, S., Lu, M., Lefevor, G.T., Fedewa, A. & Celimli, S. (2016). Application of meta-analysis in sport and exercise science. In N. Ntoumanis, & N. Myers (Eds.), \emph{An Introduction to Intermediate and Advanced Statistical Analyses for Sport and Exercise Scientists} (pp.233-253). Hoboken, NJ: John Wiley and Sons, Ltd.

Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.) (2009). \emph{The handbook of research synthesis and 
meta-analysis.} New York: Russell Sage Foundation.

Olkin, I., & Ishii, G. (1976). Asymptotic distribution of functions of a correlation matrix. In S. Ikeda (Ed.), 
\emph{Essays in probability and statistics: A volume in honor of Professor Junjiro Ogawa} (pp.5-51). Tokyo, 
Japan: Shinko Tsusho.
}

\examples{
  ######################################################
  # Example: Craft2003 data
  # Preparing covariances for multivariate meta-analysis
  ######################################################
  data(Craft2003)
  #extract correlation from the dataset (craft)
  corflat <- subset(Craft2003, select=C1:C6)
  # transform correlations to z and compute variance-covariance matrix.
  computvocv <- r.vcov(n = Craft2003$N, corflat = corflat, method = "average")
  # name transformed z scores as an input
  Input <- computvocv$zr
  # name variance-covariance matrix of trnasformed z scores as covars
  covars <- computvocv$zcov
  #####################################################
  # Running random-effects model using package "mvmeta"
  #####################################################
  #library(mvmeta)
  #mvmeta_RE <- summary(mvmeta(cbind(C1, C2, C3, C4, C5, C6),
  #                           S = covars, data = Input, method = "reml"))
  #mvmeta_RE
  ##############################################################
  # Another example:
  # Checking the example in Harris Cooper et al.'s book page 388
  ##############################################################
r1 <- c(-0.074, -0.127, 0.324, 0.523, -0.416, -0.414)
r <- rbind(r1, r1) ### the r.vcov is to handle at least two studies
n <- c(142, 142)
computvcov <- r.vcov(n = n, corflat = r, method = "average")
round(computvcov$list.rcov[[1]], 4)
round(computvcov$list.zcov[[1]], 4)
}
