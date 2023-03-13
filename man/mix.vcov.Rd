\name{mix.vcov}
\alias{mix.vcov}
\alias{vecTosm}
\title{Computing Variance-Covariance Matrices for Effect Sizes of the Same or Different Types}
\author{
Min Lu
}
\usage{
mix.vcov(d, r, nt, nc,
         st, sc, n_rt = NA, n_rc = NA,
         sdt, sdc, type,
         name = NULL, na.impute = NA)
}
\description{
The function \code{r.vcov} computes effect sizes and variance-covariance matrices between effect sizes of the same of different types. Effect sizes include mean difference (MD), standardized mean difference (SMD), log odds ratio (logOR), log risk ratio (logRR), and risk difference (RD). Formulas are in Table I of Wei et al.'s paper (2013).
}
\arguments{
\item{d }{ A \eqn{N \times p}{N x p} matrix or data frame with mean differences (MD) and/or standard mean differences (SMD) from the \eqn{N} studies. \code{d[i,j]} is the value from study \eqn{i} for outcome \eqn{j}. If outcome \eqn{j} is not MD or SMD, NA has to be imputed in column \eqn{j}. }
\item{r }{ A \eqn{N}-dimensional list of \eqn{p \times p}{p x p} correlation matrices for the \eqn{p} outcomes from the \eqn{N} studies. \code{r[[k]][i,j]} is the correlation coefficient between outcome \eqn{i} and outcome \eqn{j} from study \eqn{k}.}
\item{nt }{ A \eqn{N \times p}{N x p} matrix storing sample sizes in the treatment group reporting the  \eqn{p} outcomes. \code{nt[i,j]} is the sample size from study \eqn{i} reporting outcome \eqn{j}.}
\item{nc }{ A matrix defined in a similar way as \code{nt} for the control group.}
\item{st }{ A \eqn{N \times p}{N x p} matrix recording number of participants with event for all outcomes (dichotomous) in treatment group. \code{st[i,j]} reports number of participants with event for outcome \eqn{j} in treatment group for study \eqn{i}. If outcome \eqn{j} is not dichotomous, NA has to be imputed in column \eqn{j}.}
\item{sc }{ Defined in a similar way as \code{st} for the control group.}
\item{n_rt }{A \eqn{N}-dimensional list of \eqn{p \times p}{p x p} matrices storing sample sizes in the treatment group reporting pairwise outcomes in the off-diagonal elements. \code{n_rt[[k]][i,j]} is the sample size reporting both outcome \eqn{i} and outcome \eqn{j} from study \eqn{k}. Diagonal elements of these matrices are discarded. The default value is \code{NA}, which means that the smaller sample size reporting the corresponding two outcomes is imputed: i.e. \code{n_rt[[k]][i,j]=min(nt[k,i],nt[k,j])}. }
\item{n_rc }{ A list defined in a similar way as \code{n_rt} for the control group.}
\item{sdt }{ A \eqn{N \times p}{N x p} matrix storing sample standard deviations for each outcome from treatment group. \code{sdt[i,j]} is the sample standard deviation from study \eqn{i} for outcome \eqn{j}. If outcome \eqn{j} is not continuous such as MD or SMD, NA has to be imputed in the \eqn{j}th column.}
\item{sdc }{ A matrix defined in a similar way as \code{sdt} for the control group.}
\item{type }{ A \eqn{p}-dimensional vector indicating types of effect sizes. \code{"MD"} stands for mean difference, \code{"SMD"} stands for standardized mean difference, \code{"logOR"} stands for log odds ratio, \code{"logRR"} stands for log risk ratio, and \code{"RD"} stands for risk difference. }
\item{name }{ A \eqn{p}-dimensional vector storing names for the effect sizes.}
\item{na.impute }{ Missing values in \code{d} can be imputed by a numeric value, such as zero by setting \code{na.impute = 0}. With the default setting \code{na.impute = NA}, missing values are not imputed. If specifying \code{na.impute = "average"}, missing values are imputed by the mean of \code{d} that is sample-size weighted from the complete records. }
}
\value{
\item{ ef}{A \eqn{N \times p}{N x p} data frame whose columns are computed effect sizes according to the input argument \code{"type"}. SMD will be converted to Hedges's g (Wei and Higgins, 2013).}
\item{list.vcov }{ A \eqn{N}-dimensional list of \eqn{p(p+1)/2 \times p(p+1)/2}{p(p+1)/2 x p(p+1)/2} matrices of computed variance-covariance matrices.}
\item{matrix.vcov }{ A \eqn{N \times p(p+1)/2}{N x p(p+1)/2} matrix whose rows are computed variance-covariance vectors.}
}

\references{
Wei, Y., &  Higgins, J. (2013). Estimating within study covariances in multivariate meta-analysis with multiple outcomes. \emph{Statistics in Medicine, 32}(7), 119-1205.
}

\examples{
#########################################################################
# Example: Geeganage2010 data
# Preparing covariances for a multivariate meta-analysis
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
#######################################################################
# Running random-effects model using package "mixmeta" or "metaSEM"
#######################################################################
# Restricted maximum likelihood (REML) estimator from the mixmeta package
# library(mixmeta)
# mvmeta_RE <- summary(mixmeta(cbind(MD.SBP, MD.DBP, RD.DD, lgOR.D) ~ 1,
#                            S = S, data = y, method = "reml"))
# mvmeta_RE

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
        y.all.se = obj$coefficients[,2],
        hline = c(0, 0, 0, 1))
# dev.off()
}
