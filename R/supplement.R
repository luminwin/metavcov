#######################################################################
#  Step 1: Data preparation with the package metavcov
#  Example 1: Correlation Coefficients with r.vcov()
#
#######################################################################
library(metavcov)
r <- matrix(c(-0.074, -0.127, 0.324, 0.523, -0.416, -0.414), 1)
n <- 142
computvcov <- r.vcov(n = n, corflat = r,
                     name = paste("C", c("st", "su", "sv", "tu", "tv", "uv"), sep = ""),
                     method = "each")
round(computvcov$list.rvcov[[1]], 4)
round(computvcov$ef, 4)
round(computvcov$list.vcov[[1]], 4)
round(computvcov$matrix.vcov, 4)

computvcov <- r.vcov(n = 142,
                     corflat = matrix(c(-0.074, -0.127, 0.324, 0.523, -0.416, NA), 1),
                     na.impute = 0)
computvcov$r



data(Craft2003)
computvcov <- r.vcov(n = Craft2003$N,
                     corflat = subset(Craft2003, select = C1:C6),
                     method = "average")
y <- computvcov$ef
S <- computvcov$matrix.vcov
S[1, ]

#######################################################################
#  Step 1: Data preparation with the package metavcov
#  Example 2: Standardized Mean Difference with smd.vcov()
#
#######################################################################
data(Geeganage2010)
## correlation coefficients between outcomes are missing in the data
## impute the correlation coefficient list based on expert knowledge
r12 <- 0.71
r.Gee <- lapply(1:nrow(Geeganage2010), function(i){matrix(c(1, r12, r12, 1), 2, 2)})

computvcov <- smd.vcov(nt = Geeganage2010[ ,c("nt_SBP", "nt_DBP")],
                       nc = Geeganage2010[ ,c("nc_SBP", "nc_DBP")],
                       d = Geeganage2010[ ,c("SMD_SBP", "SMD_DBP")],
                       r = r.Gee,
                       name = c("SBP", "DBP"))
head(computvcov$ef)   ## Hedge's g
head(computvcov$matrix.vcov)   ##  variances/covariances for Hedge's g
head(computvcov$matrix.dvcov)   ##  variances/covariances for SMD

#################################################################################
#  Step 1: Data preparation with the package metavcov
#  Example 3: Covariance between Mean Difference and Log Odds Ratio with md_lgor()
#
#################################################################################

md_lgor(r = 0.71, sd1t = 0.4, sd1c = 8,
        n1c = 34, n2c = 35,
        n1t = 25, n2t = 32,
        s2c = 5, s2t = 8,
        f2c = 30, f2t = 24)
#######################################################################
#  Step 1: Data preparation with the package metavcov
#  Example 4: Combination of All Effect Sizes with mix.vcov()
#
#######################################################################
data(Geeganage2010)
## correlation coefficients between outcomes are missing in the data
## impute the correlation coefficient list based on expert knowledge
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
## compute variances and  covariances
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
## save different effect sizes in y
y <- computvcov$ef
head(y)
## save variances and covariances of all the effect sizes in a matrix S
S <- computvcov$matrix.vcov
computvcov$list.vcov[[1]]
S[1, ]
#######################################################################
#  Step 2: Fixed Effect Meta-Analysis with the package metavcov or mixmeta
#  Example 1: Correlation Coefficients with metafixed() or mixmeta()
#
#######################################################################
data(Craft2003)
computvcov <- r.vcov(n = Craft2003$N,
                     corflat = subset(Craft2003, select = C1:C6),
                     method = "average")
y <- computvcov$ef
Slist <- computvcov$list.vcov
MMA_FE <- summary(metafixed(y = y, Slist = Slist))
MMA_FE

## An alternative way 
library(mixmeta)
S <- computvcov$matrix.vcov
MMA_FE <- summary(mixmeta(cbind(C1, C2, C3, C4, C5, C6)~1,
                            S = S, data = y, method = "fixed"))
## In most cases, a random effect model is more appropriate.

#######################################################################
#  Step 2: Random Effect Meta-Analysis with the package mixmeta or metaSEM 
#  Example 1: Correlation Coefficients with mixmeta(), meta(), or rml()
#
#######################################################################

library(mixmeta)
S <- computvcov$matrix.vcov
MMA_RE <- summary(mixmeta(cbind(C1, C2, C3, C4, C5, C6)~1,
                         S = S, data = y, method = "reml"))
## meta-regression
summary(mixmeta(cbind(C1, C2, C3, C4, C5, C6)~ p_male,
                            S = S, data = data.frame(y, p_male = Craft2003$p_male),
                      method = "reml"))

library(metaSEM)
# For the maximum likelihood (ML) estimation method
summary(meta(y = y, v = S, data = data.frame(y,S))	)
# For the restricted maximum likelihood (REML) estimation method
# summary(reml(y = y, v = S, data = data.frame(y,S))	)

dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

## meta-regression
summary(meta(y = y, v = S, data = data.frame(y, S, p_male = Craft2003$p_male)))
# summary(reml(y = y, v = S, data = data.frame(y, S, p_male = Craft2003$p_male)))
#######################################################################
#  Step 3: Visualization with the package metavcov
#  Example 1: Visualizing Correlation Coefficients with plotCI()
#
#######################################################################
obj <- MMA_RE
#pdf("CI.pdf", width = 4.5, height = 7)

plotCI(y = computvcov$ef, v = computvcov$list.vcov,
        name.y = c("Correlation: cognitive anxiety & somatic anxiety \n (C1)",
                   "Correlation: cognitive anxiety & self concept \n (C2)",
                   "Correlation: cognitive anxiety & athletic performance \n (C3)",
                   "Correlation: somatic anxiety & self concept \n (C4)",
                   "Correlation: somatic anxiety & athletic performance \n (C5)",
                   "Correlation: self concept & athletic performance \n (C6)"), 
        name.study = Craft2003$ID,
        y.all = obj$coefficients[,1],
        y.all.se = obj$coefficients[,2],
        up.bound = Inf, low.bound = -Inf)

#dev.off()

#######################################################################
#
#  Multiple Imputation for Missing Data
#
#######################################################################
# prepare a dataset with missing values and input arguments for meta.mi
Craft2003.mnar <- Craft2003[, c(2, 4:10)]
Craft2003.mnar[sample(which(Craft2003$C4 < 0), 6), "C4"] <- NA
dat <- Craft2003.mnar
n.name <- "N"
ef.name <- c("C1", "C2", "C3", "C4", "C5", "C6")

# fixed effect model
####################
o1 <- metami(dat, M = 2, vcov = "r.vcov",        ## M = 2 for fast demonstration; it should be M = 20
             n.name, ef.name,
             func = "metafixed")

computvcov <- r.vcov(n = Craft2003$N,
                     corflat = subset(Craft2003.mnar, select = C1:C6),
                     method = "average")
obj <- o1
plotCI(y = computvcov$ef, v = computvcov$list.vcov,
        name.y = NULL, name.study = Craft2003$ID,
        y.all = obj$coefficients[,1],
        y.all.se = obj$coefficients[,2])

# random effect model
####################

# restricted maximum likelihood (REML) estimator from the mixmeta package
library(mixmeta)
o2 <- metami(dat, M = 2, vcov = "r.vcov",        ## M = 2 for fast demonstration; it should be M = 20  
             n.name, ef.name,
             formula = as.formula(cbind(C1, C2, C3, C4, C5, C6) ~ 1),
             func = "mixmeta",
             method = "reml")
# maximum likelihood estimators from the metaSEM package
library(metaSEM)
o2 <- metami(dat, M = 2, vcov = "r.vcov",        ## M = 2 for fast demonstration; it should be M = 20
              n.name, ef.name,
              func = "meta")

# meta-regression
library(metaSEM)
o3 <- metami(dat, M = 2, vcov = "r.vcov",        ## M = 2 for fast demonstration; it should be M = 20
                    n.name, ef.name, x.name = "p_male",
                    func = "meta")
library(mixmeta)
o3 <- metami(dat, M = 2, vcov = "r.vcov",          ## M = 2 for fast demonstration; it should be M = 20
              n.name, ef.name, x.name = "p_male",
              formula = as.formula(cbind(C1, C2, C3, C4, C5, C6) ~ p_male ),
              func = "mixmeta",
              method = "reml")


#######################################################################
#
#  Simulation Study for Missing Data
#
#######################################################################

# prepare for true parameters
##############################
computvcov <- r.vcov(n = Craft2003$N,
                     corflat = subset(Craft2003, select = C1:C6),
                     method = "average")
mixmeta_RE <- mixmeta(cbind(C1, C2, C3, C4, C5, C6) ~ 1,
                    S = computvcov$matrix.vcov, data = computvcov$ef, method = "reml")

Truepar <- summary(mixmeta_RE)$coefficients[, 1]


# Simulation
##############################
Dat <- Craft2003
m <- nrow(Craft2003)
n.name <- "N"
ef.name <- c("C1", "C2", "C3", "C4", "C5", "C6")
f <- as.formula(cbind(C1, C2, C3, C4, C5, C6) ~ 1 )

perc <- 0.33 ## missing percentage
M <- 100 ## the number of imputations
rep <- 4 ## replications for the simulation. 4 for fast demonstration; it should be 100

resl.dl <- resl.avg <- resl.imp10 <- resl.imp20 <- resl.imp50 <- resl.imp100 <- NULL
for (i in 1:rep){
temp  <- Dat[, c(n.name, ef.name)]
temp[, ef.name] <- rToz(temp[, ef.name])
dat.mnar <- temp

set.seed(i)

dat.mnar[sample(which(Dat$C4 < 0), round(m*perc)), "C4"] <- NA

temp <- metami(dat.mnar, M = M, vcov = "r.vcov",
              n.name, ef.name,
              rvcov.method = "average",
              formula = f,
              func = "mixmeta",
              method = "reml",
              pool.seq = c(10, 20, 50, 100),
              return.mi = FALSE)
resl.imp10 <- rbind(resl.imp10, (temp$result.seq$M10$coefficients[, 1] - Truepar))
resl.imp20 <- rbind(resl.imp20, (temp$result.seq$M20$coefficients[, 1] - Truepar))
resl.imp50 <- rbind(resl.imp50, (temp$result.seq$M50$coefficients[, 1] - Truepar))
resl.imp100 <- rbind(resl.imp100, (temp$result.seq$M100$coefficients[, 1] - Truepar))

computvcov <- r.vcov(n = Dat[, n.name],
                      zscore = TRUE,
                      corflat = subset(dat.mnar, select = ef.name),
                      method = "average", na.impute = "average")
 y <- subset(dat.mnar, select = ef.name)
 S <- computvcov$matrix.vcov

 temp <- mixmeta(f, S = S, data = y, method = "reml")
 result <- summary(temp)$coefficients[, 1] - Truepar
 resl.dl <- rbind(resl.dl,result)

 y <- computvcov$ef

 temp <- mixmeta(f, S = S, data = y, method = "reml")
 result <- summary(temp)$coefficients[, 1] - Truepar
 resl.avg <- rbind(resl.avg,result)

 print(paste("replication", i))
}

method <- c("Omission with mean imputed covariances", "Mean imputation",
            "Multiple imputation, M = 10", "Multiple imputation, M = 20",
            "Multiple imputation, M = 50", "Multiple imputation, M = 100")
raw.r <- list(resl.dl, resl.avg, resl.imp10, resl.imp20, resl.imp50, resl.imp100)

datbox <- list()
for (i in 1:rep){
  datbox[[i]] <- list()
  for (j in 1:length(method)){
  temp <-   data.frame(bias = raw.r[[j]][i, ], MSE = (raw.r[[j]][i, ])^2,
               ef = ef.name, method = method[j])
  datbox[[i]] <- rbind(datbox[[i]], temp)
  }
}

Datbox <- do.call(rbind, datbox)
Datbox$ef <- factor(Datbox$ef, levels = ef.name)
Datbox$method <- factor(Datbox$method, levels = method)

boxplot(bias ~  method + ef, data = Datbox)
boxplot(MSE ~  method + ef, data = Datbox)


#######################################################################
#  Another example of effect size type:
#  Standardized Mean Difference 
#
#######################################################################
#######################################################################
#  Step 1: Data preparation with the package metavcov
#######################################################################
library("metavcov")
# load and print the data
Geeganage2010 <- as.data.frame(Geeganage2010) 
subset(Geeganage2010, select = c(nt_SBP, nt_DBP,nc_SBP, nc_DBP,SMD_SBP, SMD_DBP))
## correlation coefficients between outcomes are missing in the data
## impute the correlation coefficient list based on expert knowledge
r12 <- 0.71
r.Gee <- lapply(1:nrow(Geeganage2010), function(i){matrix(c(1, r12, r12, 1), 2, 2)})
## compute the variance-covariance matrix 
computvcov <- smd.vcov(nt = subset(Geeganage2010, select = c(nt_SBP, nt_DBP)),
                       nc = subset(Geeganage2010, select = c(nc_SBP, nc_DBP)),
                       d = subset(Geeganage2010, select = c(SMD_SBP, SMD_DBP)),
                       r = r.Gee,
                       name = c("SBP", "DBP"))
head(computvcov$ef)            ## Hedge's g
head(computvcov$matrix.vcov)   ## variance-covariances for Hedge's g (matrix)
computvcov$list.vcov           ## variance-covariances for Hedge's g (list)



#######################################################################
#  Step 2: Fixed Effect Meta-Analysis with the package metavcov or mixmeta
#######################################################################
## from the package "metavcov"
y <- computvcov$ef
Slist <- computvcov$list.vcov
MMA_FE <- summary(metafixed(y = y, Slist = Slist))
MMA_FE

## An alternative way
library(mixmeta)
S <- computvcov$matrix.vcov
MMA_FE <- summary(mixmeta(cbind(SBP,DBP)~1,
                          S = S, data = y, method = "fixed"))
MMA_FE

#######################################################################
#  Step 3: Visualization with the package metavcov
#######################################################################
obj <- MMA_FE
## from the package "metavcov"

# pdf("CI.pdf", width = 4, height = 7)

plotCI(y = computvcov$ef, v = computvcov$list.vcov,
       name.y = NULL, name.study = Geeganage2010$studyID,
       y.all = obj$coefficients[,1],
       y.all.se = obj$coefficients[,2],
       up.bound = Inf, low.bound = -Inf)

# dev.off()