md_smd <- function(smd, r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), sd1t, sd2t, sd1c, sd2c){
  #  s1p: pooled standard deviation for outcome 1
  s1p <- sqrt(((n1t - 1)*sd1t^2+(n1c - 1)*sd1c^2)/(n1t + n1c - 2))
  J <- function(x = 1000){1 - 3/(4*x - 1)}
  v1 <- n1c + n1t - 2
  jv1 <- J(v1)
  rslt <- list()
  rslt$g <- smd/jv1
  rslt$v <- r*(jv1*n12c*sd1c*sd2c/n1c/n2c + jv1*n12t*sd1t*sd2t/n1t/n2t)/s1p
  rslt}

md_lgor <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, sd1c, sd1t){
  rslt <- list()
  rslt$lgor <- log((s2t/f2t)/(s2c/f2c))
  rslt$v <- r*sd1c*n12c*sqrt(n2c)*sqrt(1/s2c + 1/f2c)/n1c/n2c + r*sd1t*n12t*sqrt(n2t)*sqrt(1/s2t + 1/f2t)/n1t/n2t
  rslt}

md_lgrr <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, sd1c, sd1t){
  rslt <- list()
  rslt$lgrr <- log((s2t/n2t)/(s2c/n2c))
  rslt$v <- r*sd1c*n12c*sqrt(f2c/s2c)/n1c/n2c+r*sd1t*n12t**sqrt(f2t/s2t)/n1t/n2t
  rslt}

md_rd <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, sd1c, sd1t){
  rslt <- list()
  rslt$rd <- (s2t/n2t) - (s2c/n2c)
  rslt$v <- r*sd1c*n12c*sqrt(f2c*s2c)/n1c/(n2c^2) + r*sd1t*n12t*sqrt(f2t*s2t)/n1t/(n2t^2)
  rslt}

smd_lgor <- function(d, r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, sd1c, sd1t){
  #  s1p: pooled standard deviation for outcome 1
  s1p <- sqrt(((n1t - 1)*sd1t^2 + (n1c - 1)*sd1c^2)/(n1t + n1c - 2))
  J <- function(x = 1000){1 - 3/(4*x - 1)}
  v1 <- n1c + n1t - 2
  jv1 <- J(v1)
  rslt <- list()
  rslt$g <- d/jv1
  rslt$lgor <- log((s2t/f2t)/(s2c/f2c))
  rslt$v <- jv1*n12c*sqrt(n2c)*r*sd1c*sqrt(1/s2c + 1/f2c)/s1p/n1c/n2c + jv1*n12t*sqrt(n2t)*r*sd1t*sqrt(1/s2t+1/f2t)/s1p/n1t/n2t
  rslt}

smd_lgrr <- function(d, r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, sd1c, sd1t){
  #  s1p: pooled standard deviation for outcome 1
  s1p <- sqrt(((n1t - 1)*sd1t^2 + (n1c - 1)*sd1c^2)/(n1t + n1c - 2))
  J <- function(x = 1000){1 - 3/(4*x - 1)}
  v1 <- n1c + n1t - 2
  jv1 <- J(v1)
  rslt <- list()
  rslt$g <- d/jv1
  rslt$lgrr <- log((s2t/n2t)/(s2c/n2c))
  rslt$v <- jv1*n12c*r*sd1c*sqrt(f2c/s2c)/s1p/n1c/n2c + jv1*n12t*r*sd1t*sqrt(f2t/s2t)/s1p/n1t/n2t
  rslt}

smd_rd <- function(d, r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, sd1c, sd1t){
  #  s1p: pooled standard deviation for outcome 1
  s1p <- sqrt(((n1t - 1)*sd1t^2 + (n1c - 1)*sd1c^2)/(n1t + n1c - 2))
  J <- function(x = 1000){1 - 3/(4*x - 1)}
  v1 <- n1c + n1t - 2
  jv1 <- J(v1)
  rslt <- list()
  rslt$g <- d/jv1
  rslt$rd <- (s2t/n2t) - (s2c/n2c)
  rslt$v <- jv1*n12c*r*sd1c*sqrt(f2c*s2c)/s1p/n1c/(n2c^2) + jv1*n12t*r*sd1t*sqrt(f2t*s2t)/s1p/n1t/(n2t^2)
  rslt}

lgor_lgrr <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c,n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, s1c, s1t, f1t, f1c)  {
  rslt <- list()
  rslt$lgor <- log((s1t/f1t)/(s1c/f1c))
  rslt$lgrr <- log((s2t/n2t)/(s2c/n2c))
  rslt$v <- r*n12c*sqrt((1/s1c + 1/f1c)*f2c/s2c)/sqrt(n1c)/n2c + r*n12t*sqrt((1/s1t + 1/f1t)*f2t/s2t)/sqrt(n1t)/n2t
  rslt}

lgor_rd <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, s1c, s1t, f1t, f1c) {
  rslt <- list()
  rslt$lgor <- log((s1t/f1t)/(s1c/f1c))
  rslt$rd <- (s2t/n2t) - (s2c/n2c)
  rslt$v <- r*n12c*sqrt((s2c*f2c/s1c/f1c)*f2c/s2c)/n1c/n2c + r*n12t*sqrt((s2t*f2t/s1t/f1t)*f2t/s2t)/n1t/n2t
  rslt}

lgrr_rd <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, s1c, s1t, f1c, f1t) {
  rslt <- list()
  rslt$lgrr <- log((s1t/n1t)/(s1c/n1c))
  rslt$rd <- (s2t/n2t) - (s2c/n2c)
  rslt$v <- r*n12c*sqrt((f1c*s2c*f2c/s1c)*f2c/s2c)/n1c/(n2c^2) + r*n12t*sqrt((f1t*s2t*f2t/s1t)*f2t/s2t)/n1t/(n2t^2)
  rslt}


mix.vcov <- function(d, r, nt, nc, st, sc, n_rt = NA, n_rc = NA, sdt, sdc, type, name = NULL, na.impute = NA){
  colum.number <- length(type)
  if (is.null(name)){
    if (is.null(colnames(d))) {
      name <- paste("V", 1:colum.number, sep = "")
    } else {
      name <- colnames(d)
    }
  } else {
    name <- name[1:colum.number]
  }
  colnames(d) <- name

  cov.name <- c()
  for (i in 1:colum.number){
    if (i == colum.number){ temp <- NULL} else {
      temp <-  paste("cov",name[i], name, sep = "_")[(i+1):colum.number]
    }
    cov.name <- c(cov.name, paste("var", name[i], sep = "_"), temp)
  }


  ft <- nt - st
  fc <- nc - sc
  K <- nrow(nt)
  col.vac.number <- (colum.number + 1)*colum.number/2



  if (is.na(na.impute)) { d[is.na(d)] <- na.impute }
  else if (na.impute == "average") {
    n <- rowSums(nt + nc, na.rm = TRUE)
    temp <- unlist(lapply(1:colum.number, function(i){weighted.mean(d[, i], n, na.rm = TRUE)}))
    corflat.m <- matrix(rep(temp, K), K, colum.number, byrow = "TRUE")

    d[is.na(d)] <- corflat.m[is.na(d)]} else {
    d[is.na(d)] <- na.impute
  }


  type.sm <- matrix("SMD", ncol(nt), ncol(nt))
  for (i in 1:ncol(nt)){
    for (j in 1:ncol(nt)){
      type.sm[i,j] <- paste(type[i], type[j])
    }}
  g <- d
  result <- list()
  result$list.vcov <- list()

  if (is.na(n_rt)&(length(n_rt) == 1)){ n_rt <- rep(list(matrix(NA, colum.number, colum.number)), K) }

  for (k in 1:K) {
    for (i in 1:colum.number){
      for (j in 1:colum.number){
        if (is.na(n_rt[[k]][i, j]))
          n_rt[[k]][i, j] <- min(nt[k, i], nt[k, j])
      }
    }
  }

  if (is.na(n_rc)&(length(n_rc) == 1)){ n_rc <- rep(list(matrix(NA, colum.number, colum.number)), K) }

  for (k in 1:K) {
    for (i in 1:colum.number){
      for (j in 1:colum.number){
        if (is.na(n_rc[[k]][i, j]))
          n_rc[[k]][i, j] <- min(nc[k, i], nc[k, j])
      }
    }
  }

  for (k in 1:K){
    result$list.vcov[[k]] <- matrix(NA, colum.number, colum.number)
    colnames( result$list.vcov[[k]]) <- rownames( result$list.vcov[[k]]) <- name

    for (i in 1:colum.number){
      for (j in 1:colum.number)
      {
        if ((type[i] == "MD")&&(type[j] == "MD"))
        {
          tempr <- diag(2)
          tempr[1, 2] <- tempr[2, 1] <- r[[k]][i, j]
          tempnrt <- list(diag(2))
          diag(tempnrt[[1]]) <- c(n_rt[[k]][i, i], n_rt[[k]][j, j])
          tempnrt[[1]][1, 2] <- tempnrt[[1]][2, 1] <- n_rt[[k]][i, j]
          tempnrc <- list(diag(2))
          diag(tempnrc[[1]]) <- c(n_rc[[k]][i, i], n_rc[[k]][j, j])
          tempnrc[[1]][1,2] <- tempnrc[[1]][2, 1] <- n_rc[[k]][i, j]
          temp <- md.vcov(r = list(tempr), nt = cbind(nt[k, i], nt[k, j]), nc = cbind(nc[k, i], nc[k, j]), n_rt = tempnrt, n_rc = tempnrc, sdt = cbind(sdt[k, i], sdt[k, j]), sdc = cbind(sdc[k, i], sdc[k, j]))
          result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$list.vcov[[1]][1, 2]
        }

        if ((type[i] == "SMD")&&(type[j] == "SMD"))
        {          tempr <- diag(2)
        tempr[1, 2] <- tempr[2, 1] <- r[[k]][i, j]
        tempnrt <- list(diag(2))
        diag(tempnrt[[1]]) <- c(n_rt[[k]][i, i], n_rt[[k]][j, j])
        tempnrt[[1]][1, 2] <- tempnrt[[1]][2, 1] <- n_rt[[k]][i, j]
        tempnrc <- list(diag(2))
        diag(tempnrc[[1]]) <- c(n_rc[[k]][i, i], n_rc[[k]][j, j])
        tempnrc[[1]][1, 2] <- tempnrc[[1]][2, 1] <- n_rc[[k]][i, j]
        temp <- smd.vcov(r = list(tempr), d = cbind(d[k, i], d[k, j]), nt = cbind(nt[k, i], nt[k, j]), nc = cbind(nc[k, i], nc[k, j]), n_rt = tempnrt, n_rc = tempnrc)
        g[k, i] <- temp$ef[1, 1]
        g[k, j] <- temp$ef[1, 2]
        result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$list.vcov[[1]][1, 2]
        }

        if ((type[i] == "lgOR")&&(type[j] == "lgOR"))
        { tempr <- diag(2)
        tempr[1, 2] <- tempr[2, 1] <- r[[k]][i, j]
        tempnrt <- list(diag(2))
        diag(tempnrt[[1]]) <- c(n_rt[[k]][i, i], n_rt[[k]][j, j])
        tempnrt[[1]][1, 2] <- tempnrt[[1]][2, 1] <- n_rt[[k]][i, j]
        tempnrc <- list(diag(2))
        diag(tempnrc[[1]]) <- c(n_rc[[k]][i, i], n_rc[[k]][j, j])
        tempnrc[[1]][1, 2] <- tempnrc[[1]][2, 1] <- n_rc[[k]][i, j]
        temp <- lgOR.vcov(r = list(tempr), nt = cbind(nt[k, i], nt[k, j]), nc = cbind(nc[k, i], nc[k, j]),
                        n_rt = tempnrt, n_rc = tempnrc,
                        st = cbind(st[k, i], st[k, j]),
                        sc = cbind(sc[k, i], sc[k, j]))
        result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$list.vcov[[1]][1, 2]
        g[k, i] <- temp$ef[1, 1]
        g[k, j] <- temp$ef[1, 2]
        }

        if ((type[i] == "lgRR")&&(type[j] == "lgRR"))
        { tempr <- diag(2)
        tempr[1, 2] <- tempr[2, 1] <- r[[k]][i, j]
        tempnrt <- list(diag(2))
        diag(tempnrt[[1]]) <- c(n_rt[[k]][i, i], n_rt[[k]][j, j])
        tempnrt[[1]][1, 2] <- tempnrt[[1]][2, 1] <- n_rt[[k]][i, j]
        tempnrc <- list(diag(2))
        diag(tempnrc[[1]]) <- c(n_rc[[k]][i, i], n_rc[[k]][j, j])
        tempnrc[[1]][1, 2] <- tempnrc[[1]][2, 1] <- n_rc[[k]][i, j]
        temp <- lgRR.vcov(r = list(tempr), nt = cbind(nt[k, i], nt[k, j]), nc = cbind(nc[k, i], nc[k, j]),
                        n_rt = tempnrt, n_rc = tempnrc,
                        st = cbind(st[k, i], st[k, j]),
                        sc = cbind(sc[k, i], sc[k, j]))
        result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$list.vcov[[1]][1, 2]
        g[k, i] <- temp$ef[1, 1]
        g[k, j] <- temp$ef[1, 2]
        }

        if ((type[i] == "RD")&&(type[j] == "RD"))
        { tempr <- diag(2)
        tempr[1, 2] <- tempr[2, 1] <- r[[k]][i, j]
        tempnrt <- list(diag(2))
        diag(tempnrt[[1]]) <- c(n_rt[[k]][i, i], n_rt[[k]][j, j])
        tempnrt[[1]][1, 2] <- tempnrt[[1]][2, 1] <- n_rt[[k]][i, j]
        tempnrc <- list(diag(2))
        diag(tempnrc[[1]]) <- c(n_rc[[k]][i, i], n_rc[[k]][j, j])
        tempnrc[[1]][1, 2] <- tempnrc[[1]][2, 1] <- n_rc[[k]][i, j]
        temp <- rd.vcov(r = list(tempr), nt = cbind(nt[k, i], nt[k, j]), nc = cbind(nc[k, i], nc[k, j]),
                      n_rt = tempnrt, n_rc = tempnrc,
                      st = cbind(st[k, i], st[k, j]),
                      sc = cbind(sc[k, i], sc[k, j]))
        result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$list.vcov[[1]][1, 2]
        g[k, i] <- temp$ef[1, 1]
        g[k, j] <- temp$ef[1, 2]
        }

        ##############################
        if ((type[i] == "MD")&&(type[j] == "SMD"))
        { temp <- md_smd(smd = d[k, j], r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j], sd2t = sdt[k, j], sd2c = sdc[k, j],
                 sd1t = sdt[k, i], sd1c = sdc[k, i])
          result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$v
          g[k, j] <- temp$g
        }
        if ((type[i] == "MD")&&(type[j] == "lgOR"))
        { temp <- md_lgor(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                  sd1t = sdt[k, i], s2t = st[k, j], sd1c = sdc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j])
          result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$v
          g[k, j] <- temp$lgor
        }

        if ((type[i] == "MD")&&(type[j] == "lgRR"))
        { temp <- md_lgrr(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                          sd1t = sdt[k, i], s2t = st[k, j], sd1c = sdc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j])
          result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$v
          g[k, j] <- temp$lgrr
        }

        if ((type[i] == "MD")&&(type[j] == "RD"))
        { temp <- md_rd(r = r[[k]][i, j], n1c=nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                        sd1t = sdt[k, i], s2t = st[k, j], sd1c = sdc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j])
          result$list.vcov[[k]][i,j]<-result$list.vcov[[k]][j, i] <- temp$v
          g[k, j] <- temp$rd
        }

        if ((type[i] == "SMD")&&(type[j] == "lgOR"))
        { temp <- smd_lgor(d = d[k, i], r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k,j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                           sd1t = sdt[k, i], s2t = st[k, j], sd1c = sdc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j])
          result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$v
          g[k, i] <- temp$g
          g[k, j] <- temp$lgor
        }

        if ((type[i] == "SMD")&&(type[j] == "lgRR"))
        { temp <- smd_lgrr(d = d[k, i], r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                           sd1t = sdt[k, i], s2t = st[k, j], sd1c = sdc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j])
          result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$v
          g[k, i] <- temp$g
          g[k, j] <- temp$lgrr
        }

        if ((type[i] == "SMD")&&(type[j] == "RD"))
        { temp <- smd_rd(d = d[k, i], r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                         sd1t = sdt[k,i], s2t = st[k, j], sd1c = sdc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j])
          result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$v
          g[k, i] <- temp$g
          g[k, j] <- temp$rd
        }

        if ((type[i] == "lgOR")&&(type[j] == "lgRR"))
        { temp <- lgor_lgrr(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                    s1t = st[k, i], s2t = st[k, j], s1c = sc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j], f1c = fc[k, i], f1t = ft[k, i])
          result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$v
          g[k, i] <- temp$lgor
          g[k, j] <- temp$lgrr
        }

        if ((type[i] == "lgOR")&&(type[j] == "RD"))
        { temp <- lgor_rd(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                          s1t = st[k,i],s2t=st[k,j], s1c = sc[k, i],s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j], f1c = fc[k, i], f1t = ft[k, i])
          result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$v
          g[k, i] <- temp$lgor
          g[k, j] <- temp$rd
        }

        if ((type[i] == "lgRR")&&(type[j] == "RD"))
        { temp <- lgrr_rd(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                          s1t = st[k, i], s2t = st[k, j], s1c = sc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j], f1c = fc[k, i], f1t = ft[k, i])
          result$list.vcov[[k]][i, j] <- result$list.vcov[[k]][j, i] <- temp$v
          g[k, i] <- temp$lgrr
          g[k, j] <- temp$rd
        }
      }
    }
  }
  result$ef <- as.data.frame(g)
  temp <- matrix(unlist(lapply(1:K, function(k){smTovec(result$list.vcov[[k]])})), K, col.vac.number, byrow = TRUE)
  temp <- as.data.frame(temp)
  colnames(temp) <- cov.name
  result$matrix.vcov <- temp
  result
}
