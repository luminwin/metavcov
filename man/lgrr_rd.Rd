\name{lgrr_rd}
\alias{lgrr_rd}
\title{Covariance between log risk ratio and risk difference}
\author{
Min Lu
}
\usage{
lgrr_rd(r, n1c, n2c, n1t, n2t,
        n12c = min(n1c, n2c),n12t = min(n1t, n2t),
}
\description{
Compute covariance between log risk ratio and risk difference, when the two outcomes are binary.
}
\arguments{
\item{r }{ Correlation coefficient of the two outcomes.}
\item{n1c }{ Number of participants reporting outcome 1 in control group.}
\item{n2c }{ Number of participants reporting outcome 2 in control group.}
\item{n1t }{ Number of participants reporting outcome 1 in treatment group.}
\item{n2t }{ Number of participants reporting outcome 2 in treatment group.}
\item{n12c }{ Number of participants reporting both outcome 1 and outcome 2 in control group. By default, it is equal to the smaller number between n1c and n2c.}
\item{n12t }{ Number defined in a similar way as n12c for treatment group. }
\item{s2c }{ Number of participants with event for outcome 2 (dichotomous) in control group.}
\item{s2t }{ Defined in a similar way as s2c for treatment group. }
\item{f2c }{ Number of participants without event for outcome 2 (dichotomous) in control group.}
\item{f2t }{ Defined in a similar way as f2c for treatment group. }
\item{s1c }{ Number of participants with event for outcome 1 (dichotomous) in control group.}
\item{s1t }{ Defined in a similar way as s1c for treatment group. }
\item{f1c }{ Number of participants without event for outcome 1 (dichotomous) in control group.}
\item{f1t }{ Defined in a similar way as f1c for treatment group. }
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
lgrr_rd(r = 0.71, n1c = 30, n2c = 35, n1t = 28, n2t = 32,
        s2c = 5, s2t = 8, f2c = 30, f2t = 24,
        s1c = 5, s1t = 8, f1c = 25, f1t = 20)
## calculate covariances for variable D and DD in Geeganage2010 data
attach(Geeganage2010)
D_DD <- unlist(lapply(1:nrow(Geeganage2010), function(i){lgrr_rd(r = 0.71,
        n1c = nc_SBP[i], n2c = nc_DD[i],
        n1t = nt_SBP[i], n2t = nt_DD[i], s2t = st_DD[i], s2c = sc_DD[i],
        f2c = nc_DD[i] - sc_DD[i], f2t = nt_DD[i] - st_DD[i],
        s1t = st_D[i], s1c = sc_D[i], f1c = nc_D[i] - sc_D[i], f1t = nt_D[i] - st_D[i])}))
D_DD
}
