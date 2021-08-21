\name{md_smd}
\alias{md_smd}
\title{Covariance between mean difference and standardized mean difference}
\author{
Min Lu
}
\usage{
md_smd(r, n1c, n2c, n1t, n2t,
       n12c = min(n1c, n2c), n12t = min(n1t, n2t),
       sd1t, sd2t, sd1c, sd2c)
}
\description{
Compute covariance between mean difference and standardized mean difference, when effect sizes are different.
}
\arguments{
\item{r }{ Correlation coefficient of the two outcomes.}
\item{n1c }{ Number of participants reporting outcome 1 in control group.}
\item{n2c }{ Number of participants reporting outcome 2 in control group.}
\item{n1t }{ Number of participants reporting outcome 1 in treatment group.}
\item{n2t }{ Number of participants reporting outcome 2 in treatment group.}
\item{n12c }{ Number of participants reporting both outcome 1 and outcome 2 in control group. By default, it is equal to the smaller number between n1c and n2c.}
\item{n12t }{ Number defined in a similar way as n12c for treatment group. }
\item{sd1t }{ Sample standard deviation of outcome 1.}
\item{sd2t }{ Sample standard deviation of outcome 2.}
\item{sd1c }{ Defined in a similar way as sd1t for control group.}
\item{sd2c }{ Defined in a similar way as sd2t for control group.}
}

\value{
Return the computed covariance.
}


\references{
Ahn, S., Lu, M., Lefevor, G.T., Fedewa, A. & Celimli, S. (2016). Application of meta-analysis in sport and exercise science. In N. Ntoumanis, & N. Myers (Eds.), \emph{An Introduction to Intermediate and Advanced Statistical Analyses for Sport and Exercise Scientists} (pp.233-253). Hoboken, NJ: John Wiley and Sons, Ltd.

Wei, Y., &  Higgins, J. (2013). Estimating within study covariances in multivariate meta-analysis with multiple outcomes. \emph{Statistics in Medicine, 32}(7), 119-1205.
}

\examples{
## simple example
md_smd(r = 0.71, n1c = 34, n2c = 35, n1t = 25, n2t = 32,
       sd1t = 0.6, sd2t = 0.4, sd1c = 8, sd2c = 0.9)
## calculate covariances for variable SBP and DBP in Geeganage2010 data
attach(Geeganage2010)
SBP_DBP <- unlist(lapply(1:nrow(Geeganage2010), function(i){md_smd(r = 0.71,
                n1c = nc_SBP[i], n2c = nc_DBP[i], n1t = nt_SBP[i], n2t = nt_DBP[i],
                sd1t = sdt_SBP[i], sd2t = sdt_DBP[i],
                sd1c = sdc_SBP[i], sd2c = sdc_SBP[i])}))
SBP_DBP
}
