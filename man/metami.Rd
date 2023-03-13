\name{metami}
\alias{metami}
\alias{print.summary.metami}
\alias{summary.metami}
\alias{rubinpool}
\title{Multiple Imputation for Missing Data in Meta-Analysis}
\author{
Min Lu
}
\usage{
metami(data, M = 20, vcov = "r.vcov",
        r.n.name, ef.name, x.name = NULL,
        rvcov.method = "average", rvcov.zscore = TRUE,
        type = NULL,
        d = NULL, sdt = NULL, sdc = NULL,
        nt = NULL, nc = NULL,
        st = NULL, sc = NULL,
        n_rt = NA, n_rc = NA,
        r = NULL,
        func = "mixmeta",
        formula = NULL,
        method = "fixed",
        pool.seq = NULL,
        return.mi = FALSE,
        ci.level = 0.95)
}
\description{
Multiple imputation allows for the uncertainty about the missing data by generating several different plausible imputed data sets and appropriately combining results obtained from each of them. Let  \eqn{\hat{\theta}_{*m}} be the estimated coefficient from the \eqn{m}th imputed dataset for one of the \eqn{p} dimensions in the multivariate outcome, where \eqn{m=1,\dots,M}. The coefficient from MI \eqn{\bar{\theta}} is  simply just an arithmetic mean of the individual coefficients estimated from each of the \eqn{M} meta-analysis. We have
\deqn{\bar{\theta}=\frac{\sum_{m=1}^{M}\hat{\theta}_{*m}}{M}.}
Estimation of the standard error for each variable is little more complicated. Let \eqn{V_W} be the within imputation variance, which is the average of the variance of the estimated coefficient from each imputed dateset:
  \deqn{V_W=\frac{\sum_{m=1}^{M}V ({\hat{\theta}_{*m}})}{M},}
where \eqn{V ({\hat{\theta}_{*m}})} is the variance of the estimator calculated from generalized least squares methods using the imputed dataset. Let \eqn{V_B} be the between imputation variance, which is calculated as
\deqn{V_B=\frac{\sum_{m=1}^{M}({\hat{\theta}_{*m}}-\bar{\theta})^2}{M-1}.}
From \eqn{V_W} and \eqn{V_B}, the variance of the pooled coefficients is calculated as
\deqn{V(\bar{\theta})=V_W+V_B+\frac{V_B}{M}}
The above variance is statistically principled since \eqn{V_W} reflects the sampling variance and \eqn{V_B} reflects the extra variance due to the missing data.
}
\arguments{
\item{data }{ A  \eqn{N \times p}{N x p} data frame that contains effect sizes and predictors for meta-regression, if any.}
\item{M}{ Number of imputed data sets.}
\item{vcov }{ Method for computing effect sizes; options including \code{vcov = "r.vcov"} for correlation coefficients and \code{vcov = "mix.vcov"} for other types of effect sizes. See \code{\link{r.vcov}} and \code{\link{mix.vcov}} for details.}
\item{r.n.name }{ A string defining the column name for sample sizes in \code{data} when the effect sizes are correlation coefficients (\code{vcov = "r.vcov"}).}
\item{ef.name }{ A \eqn{p}-dimensional vector that stores the column names for sample sizes in \code{data} when the effect sizes are correlation coefficients (\code{vcov = "r.vcov"}). }
\item{x.name }{ A vector that stores the column names in \code{data} for predictors for meta-regression. }
\item{rvcov.method}{Method used for \code{\link{r.vcov}}; options including \code{"average"} and \code{"each"}. }
\item{rvcov.zscore }{ Whether the correlation coefficients in \code{data} are already transformed into Fisher's \eqn{z} scores.}
\item{type }{ A \eqn{p}-dimensional vector indicating types of effect sizes for the argument \code{vcov = "mix.vcov"}. \code{"MD"} stands for mean difference, \code{"SMD"} stands for standardized mean difference, \code{"logOR"} stands for log odds ratio, \code{"logRR"} stands for log risk ratio, and \code{"RD"} stands for risk difference. }
\item{d }{ A \eqn{p}-dimensional vector that stores the column names in \code{data} for continuous effect sizes such as MD or SMD. If outcome \eqn{j} is dichotomous, NA has to be imputed in for \eqn{d[j]}. }
\item{sdt }{ A \eqn{p}-dimensional vector that stores the column names in \code{data} for the sample standard deviations of each outcome from the treatment group. If outcome \eqn{j} is dichotomous, NA has to be imputed in for \eqn{d[j]}. }
\item{sdc }{ A vector defined in a similar way as \code{sdt} for the control group.}
\item{nt }{ A \eqn{p}-dimensional vector that stores the column names in \code{data} for sample sizes of \eqn{p} outcomes from treatment group. }
\item{nc }{ A vector defined in a similar way as \code{nt} for the control group.}
\item{st }{ A \eqn{p}-dimensional vector that stores the column names in \code{data} for the number of participants with event for all outcomes (dichotomous) in the treatment group. If outcome \eqn{j} is dichotomous, NA has to be imputed in for \eqn{st[j]}. }
\item{sc }{ A vector defined in a similar way as \code{st} for the control group.}
\item{n_rt}{A \eqn{N}-dimensional list of \eqn{p \times p}{p x p} correlation matrices storing sample sizes in the treatment group reporting pairwise outcomes in the off-diagonal elements. See \code{\link{mix.vcov}} for details.}
\item{n_rc}{ A list defined in a similar way as \code{n_rt} for the control group.}
\item{r}{A \eqn{N}-dimensional list of \eqn{p \times p}{p x p} correlation matrices for the \eqn{p} outcomes from the \eqn{N} studies. See \code{\link{mix.vcov}} for details.}
\item{func}{ A string defining the function to be used for fitting the meta-analysis. Options include \code{func = "metafixed"} for fixed-effect meta-analysis (see \code{\link{metafixed}} for details). \code{func = "mixmeta"}, for which the \code{mixmeta} package must be installed beforehand, and \code{func = "meta"}, for which the \code{metaSEM} package must be installed beforehand. }
\item{formula}{ Formula used for the function  \code{func = "mixmeta"} from the \code{mixmeta} package when \code{func = "mixmeta"}. }
\item{method}{ Method used for the function  \code{func = "mixmeta"} from the \code{mixmeta} package when \code{func = "mixmeta"}.}
\item{pool.seq}{ A numeric vector indicating if the results are pooled from subsets of the \code{M} data sets. By default, the results are only pooled from all \code{M} data sets. }
\item{return.mi}{ Should the \code{M} imputed data sets be returned?}
\item{ci.level}{ Significant level for the pooled confidence intervals. The default is 0.05. }
}
\details{
For the imputation phase, this function imports the \code{mice} package that imputes incomplete multivariate data by chained equations. The pooling phase is performed via the Rubin's rules.
}
\value{
\item{coefficients}{ A data.frame that contains the pooled results from the \code{M} imputed data sets.}
\item{results.mi}{ A \eqn{M}-dimensional list of results from each imputed data set.}
\item{data.mi }{ A \eqn{M}-dimensional list of imputed data sets if the argument \code{return.mi = TRUE}.}
\item{result.seq}{ A list of results from the pooled results from the subsets of the \code{M} imputed data sets if the argument \code{pool.seq = TRUE}.}
}

\references{
Ahn, S., Lu, M., Lefevor, G.T., Fedewa, A. & Celimli, S. (2016). Application of meta-analysis in sport and exercise science. In N. Ntoumanis, & N. Myers (Eds.), \emph{An Introduction to Intermediate and Advanced Statistical Analyses for Sport and Exercise Scientists} (pp.233-253). Hoboken, NJ: John Wiley and Sons, Ltd.

Cooper, H., Hedges, L.V., & Valentine, J.C. (Eds.) (2009). \emph{The handbook of research synthesis and
meta-analysis.} New York: Russell Sage Foundation.

Olkin, I., & Ishii, G. (1976). Asymptotic distribution of functions of a correlation matrix. In S. Ikeda (Ed.),
\emph{Essays in probability and statistics: A volume in honor of Professor Junjiro Ogawa} (pp.5-51). Tokyo,
Japan: Shinko Tsusho.

Van Buuren, S. and Groothuis-Oudshoorn, K., 2011. mice: Multivariate imputation by chained equations in R. \emph{Journal of statistical software}, 45(1), pp.1-67.

Gasparrini A., Armstrong, B., Kenward M. G. (2012). Multivariate meta-analysis for non-linear and other multi-parameter associations. \emph{Statistics in Medicine}. 31(29):3821-3839.

Cheung, M.W.L. (2015). {metaSEM}: An R Package for Meta-Analysis using Structural Equation Modeling. \emph{Frontiers in Psychology} 5, 1521.

Rubin, D.B., 2004. \emph{Multiple imputation for nonresponse in surveys} (Vol. 81). John Wiley & Sons.
}

\examples{
\donttest{
#####################################################################################
# Example: Craft2003 data
# Preparing input arguments for meta.mi() and fixed-effect model
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

########################
# Plotting the result
########################
computvcov <- r.vcov(n = Craft2003$N,
                     corflat = subset(Craft2003.mnar, select = C1:C6),
                     method = "average")
plotCI(y = computvcov$ef, v = computvcov$list.vcov,
        name.y = NULL, name.study = Craft2003$ID,
        y.all = obj$coefficients[,1],
        y.all.se = obj$coefficients[,2])

########################
# Pooling from subsets
########################
# o1 <- metami(dat, M = 10, vcov = "r.vcov",
#               n.name, ef.name,
#               func = "metafixed",
#               pool.seq = c(5, 10))
# pooled results from M = 5 imputed data sets
# o1$result.seq$M5$coefficients
# pooled results from M = 10 imputed data sets
# o1$result.seq$M10$coefficients
#########################################################################################
# Running random-effects and meta-regression model using packages "mixmeta" or "metaSEM"
#########################################################################################
# Restricted maximum likelihood (REML) estimator from the mixmeta package
# library(mixmeta)
# o2 <- metami(dat, M = 10, vcov = "r.vcov",
#                     n.name, ef.name,
#                     formula = as.formula(cbind(C1, C2, C3, C4, C5, C6) ~ 1),
#                     func = "mixmeta",
#                     method = "reml")
# maximum likelihood estimators from the metaSEM package
# library(metaSEM)
# o3 <- metami(dat, M = 10, vcov = "r.vcov",
#               n.name, ef.name,
#               func = "meta")

# meta-regression
# library(metaSEM)
# o4 <- metami(dat, M = 10, vcov = "r.vcov",
#                     n.name, ef.name, x.name = "p_male",
#                     func = "meta")
# library(mixmeta)
# o5 <- metami(dat, M = 20, vcov = "r.vcov",
#               n.name, ef.name, x.name = "p_male",
#               formula = as.formula(cbind(C1, C2, C3, C4, C5, C6) ~ p_male ),
#               func = "mixmeta",
#               method = "reml")
#####################################################################################
# Example: Geeganage2010 data
# Preparing input arguments for meta.mi() and fixed-effect model
#####################################################################################
# Geeganage2010.mnar <- Geeganage2010
# Geeganage2010.mnar$MD_SBP[sample(1:nrow(Geeganage2010),7)] <- NA
# r12 <- 0.71
# r13 <- 0.5
# r14 <- 0.25
# r23 <- 0.6
# r24 <- 0.16
# r34 <- 0.16
# r <- vecTosm(c(r12, r13, r14, r23, r24, r34))
# diag(r) <- 1
# mix.r <- lapply(1:nrow(Geeganage2010), function(i){r})
# o <- metami(data = Geeganage2010.mnar, M = 10, vcov = "mix.vcov",
#               ef.name = c("MD_SBP", "MD_DBP", "RD_DD", "lgOR_D"),
#               type = c("MD", "MD", "RD", "lgOR"),
#               d = c("MD_SBP", "MD_DBP", NA, NA),
#               sdt = c("sdt_SBP", "sdt_DBP", NA, NA),
#               sdc = c("sdc_SBP", "sdc_DBP", NA, NA),
#               nt = c("nt_SBP", "nt_DBP", "nt_DD", "nt_D"),
#               nc = c("nc_SBP", "nc_DBP", "nc_DD", "nc_D"),
#               st = c(NA, NA, "st_DD", "st_D"),
#               sc = c(NA, NA, "sc_DD", "sc_D"),
#               r = mix.r,
#               func = "metafixed")
}}
