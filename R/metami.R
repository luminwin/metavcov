metami <- function(data, M = 20, vcov = "r.vcov",
                    r.n.name, ef.name, x.name = NULL,
                    rvcov.method = "average",
                    rvcov.zscore = TRUE,
                    type = NULL,
                    d = NULL,
                    sdt = NULL,
                    sdc = NULL,
                    nt = NULL,
                    nc = NULL,
                    st = NULL,
                    sc = NULL,
                    n_rt = NA,
                    n_rc = NA,
                    r = NULL,
                    func = "mvmeta",
                    formula = NULL,
                    method = "fixed",
                    pool.seq = NULL,
                    return.mi = FALSE,
                    ci.level = 0.95){
# if (!exists(func)) {stop(paste("Please load the package for the function", func))}
pool <- c("coefficients")
dat <- data; rm(data)
p <- ncol(dat)
N <- nrow(dat)

predMatrix <- mice::make.predictorMatrix(dat)

cmplt <- colnames(dat)[unlist(lapply(1:p, function(i){length(which(is.na(dat[,i]) == TRUE)) == 0}))]

if (length(cmplt) == p) stop('There is no missing values in your data')
predMatrix[cmplt, ] <- 0  # don't impute these variables since they have no NAs

imp <- mice::mice(dat, print = FALSE, m = M, predictorMatrix = predMatrix, method = mice::make.method(dat))

mis <- lapply(1:M, function(i){ unlist(lapply(1:p, function(j){imp$imp[[j]][, i]}) )})

if (return.mi) dat.mi <- list() else dat.mi <- NULL
pp <- 2 # length(pool)
out.l <- list()
for (i in 1:M){
dat.imp <- as.data.frame(dat)
dat.imp[is.na(dat.imp)] <- mis[[i]]

vcov <- c("r.vcov", "mix.vcov")[match(vcov, c("r.vcov", "mix.vcov"))]
if (vcov == "r.vcov"){
obj <- r.vcov(n = dat.imp[, r.n.name],
                     corflat = subset(dat.imp, select = ef.name),
                     zscore = TRUE,
                     method = rvcov.method, name = ef.name) } else if (vcov == "mix.vcov"){
                       if (!is.na(n_rt)) {n_rt <- dat.imp[, n_rt]}
                       if (!is.na(n_rc)) {n_rc <- dat.imp[, n_rc]}
                       eval.d <- as.data.frame(matrix(NA, N, length(d)))
                       colnames(eval.d) <- ef.name
                       eval.d[, !is.na(d)] <- dat.imp[, d[!is.na(d)]]

                       eval.sdt <- as.data.frame(matrix(NA, N, length(sdt)))
                       eval.sdt[, !is.na(sdt)] <- dat.imp[, sdt[!is.na(sdt)]]

                       eval.sdc <- as.data.frame(matrix(NA, N, length(sdc)))
                       eval.sdc[, !is.na(sdc)] <- dat.imp[, sdc[!is.na(sdc)]]

                       eval.nt <- as.data.frame(matrix(NA, N, length(nt)))
                       eval.nt[, !is.na(nt)] <- dat.imp[, nt[!is.na(nt)]]

                       eval.nc <- as.data.frame(matrix(NA, N, length(nc)))
                       eval.nc[, !is.na(nc)] <- dat.imp[, nc[!is.na(nc)]]

                       eval.st <- as.data.frame(matrix(NA, N, length(st)))
                       eval.st[, !is.na(st)] <- dat.imp[, st[!is.na(st)]]

                       eval.sc <- as.data.frame(matrix(NA, N, length(sc)))
                       eval.sc[, !is.na(sc)] <- dat.imp[, sc[!is.na(sc)]]

                       obj <- mix.vcov(type = type,
                                        d = eval.d,
                                        sdt = eval.sdt,
                                        sdc = eval.sdc,
                                        nt = eval.nt,
                                        nc = eval.nc,
                                        st = eval.st,
                                        sc = eval.sc,
                                        n_rt = n_rt,
                                        n_rc = n_rc,
                                        r = r,
                                        name = ef.name)
                     }
if (rvcov.zscore == FALSE) {
  if (vcov == "mix.vcov") {stop("rvcov.zscore == FALSE only makes sense if argument vcov is r.vcov")}
  y.name <- "r"
   if (func == "metafixed") {y.v.name <- "list.rvcov"} else {
   y.v.name <- "rvcov"}} else {
  y.name <- "ef"
   if (func == "metafixed") {y.v.name <- "list.vcov"} else {
     y.v.name <- "matrix.vcov"}}

ef <- obj[[y.name]]
ef.v <- obj[[y.v.name]]

if (return.mi) dat.mi[[i]] <- list(dat.imp = dat.imp, ef = ef, ef.v = ef.v)

if (func == "mvmeta") {
  pool <- c("coefficients", "qstat")
if("mvmeta" %in% rownames(installed.packages()) == FALSE) {install.packages("mvmeta")}
  if (is.null(formula)) {
  stop("Formula must be specified for mvmeta") } else {
    if (is.null(x.name)) {
      o <- mvmeta::mvmeta(formula = formula, S = ef.v, data = ef, method = method) } else {
        xdat <- subset(dat.imp, select = x.name)
      o <- mvmeta::mvmeta(formula = formula, S = ef.v, method = method,
                  data = data.frame(ef, xdat))
      }

  } }
if (func == "metafixed") { o <- metafixed(y = ef, Slist = ef.v); pool <- c("coefficients", "qstat")}

if (func == "meta") {
  pool <- c("coefficients", "Q.stat")
  if("metaSEM" %in% rownames(installed.packages()) == FALSE) {install.packages("metaSEM")}
if (is.null(x.name)) {
  o <- metaSEM::meta(y = ef, v = ef.v, data = data.frame(ef,ef.v)) } else {
  xdat <- subset(dat.imp, select = x.name)
  o <- metaSEM::meta(y = ef, v = ef.v, x = xdat,
            data = data.frame(ef, ef.v, xdat))
    }}

oo <- summary(o)

output <- vector(mode = "list", length = pp)
names(output) <- pool

for (j in 1:pp){
  output[[pool[j]]]  <- as.data.frame(oo[[pool[j]]])
}
out.l[[i]] <- output
}

pp <- 1
out <- vector(mode = "list", length = pp)
names(out) <- "coefficients" # pool

for (j in 1:pp){
  temp <- lapply(1:M, function(i) { out.l[[i]][[j]] })
  out[[j]] <- Reduce("+", temp) / M
}

rnames <- rownames(out[["coefficients"]])
out[["coefficients"]] <- rubinpool(out.l, ci.level, rnames)

result <- out
result$data.mi <- dat.mi
result$results.mi <- out.l

if (!is.null(pool.seq)){

  temp <- vector(mode = "list", length = length(pool.seq))
  names(temp) <- paste("M", pool.seq, sep ="")
  for (i in 1:length(pool.seq)){
    for (j in 1:pp){
      tmp <- lapply(1:pool.seq[i], function(i) { out.l[[i]] })
      temp[[i]][[j]] <- rubinpool(tmp, ci.level, rnames)
    }
    names(temp[[i]]) <- "coefficients" # pool
  }
  result$result.seq <- temp
}


class(out) <- class(result) <- "metami"
cat(paste("pooled results from", M, "imputations for missing values in", paste(setdiff(colnames(dat), cmplt), collapse = ","), "\n"))
print(summary(out))
result
}

print.summary.metami <- function(x, ...){
  digits = 4
  cat("Fixed-effects coefficients","\n",sep="")
  signif <- symnum(x$coefficients[,"Pr(>|z|)"],corr=FALSE,na=FALSE,
                   cutpoints=c(0, 0.001,0.01,0.05,0.1,1),symbols=c("***","**","*","."," "))
  ###########################################################################
  # Summary table
  #
  tabletot <- formatC(x$coefficients,digits=digits,format="f")
  tabletot <- cbind(tabletot,signif)
  colnames(tabletot)[7] <- ""

  print(tabletot,quote=FALSE,right=TRUE,print.gap=2)
  cat("---\nSignif. codes: ",attr(signif,"legend"),"\n\n")
}

summary.metami <- function(object, ...){
  fit = object
  ci.level = 0.95
  x <- list(coefficients = fit$coefficients[!is.na(fit$coefficients[,2]),])
  class(x) <- "summary.metami"
  x
}

maketable <- function(fit, ci.level = 0.95, names){
  coef <- as.numeric(fit$coef)
  coef.se <- as.numeric(fit$vcov)
  zval <- coef/coef.se
  zvalci <- qnorm((1 - ci.level)/2,lower.tail = FALSE)
  pvalue <- 2*(1-pnorm(abs(zval)))
  ci.lb <- coef-zvalci*coef.se
  ci.ub <- coef+zvalci*coef.se
  cilab <- paste(signif(ci.level,2)*100,"%ci.",c("lb","ub"),sep = "")

  tab <- cbind(coef, coef.se, zval, pvalue, ci.lb, ci.ub)
  dimnames(tab) <- list(names,
                        c("Estimate","Std. Error","z","Pr(>|z|)",cilab))
  tab
}
rubinpool <- function(o.list, ci.level, names){
  M <- length(o.list)
  theta <- do.call(rbind, lapply(1:M, function(i){ o.list[[i]]$coefficients[,1]}))
  vw <- do.call(rbind, lapply(1:M, function(i){ o.list[[i]]$coefficients[,2]}))
  Vw <- colMeans(vw^2)
  thetabar <- colMeans(theta)
  Vb <- colSums(theta - matrix(rep(thetabar, M), nrow = M, byrow = TRUE))^2/(M-1)
  Vtotal <- Vw + Vb + Vb/M
  fito <- list(coef = thetabar, vcov = sqrt(Vtotal))
  maketable(fito, ci.level = ci.level, names)
}

