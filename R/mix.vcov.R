md_smd <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), sd1t, sd2t, sd1c, sd2c){
  #  s1p: pooled standard deviation for outcome 1
  s1p <- sqrt(((n1t - 1)*sd1t^2+(n1c - 1)*sd1c^2)/(n1t + n1c - 2))
  J <- function(x = 1000){1 - 3/(4*x - 1)}
  v1 <- n1c + n1t - 2
  jv1 <- J(v1)
  r*(jv1*n12c*sd1c*sd2c/n1c/n2c + jv1*n12t*sd1t*sd2t/n1t/n2t)/s1p}

md_lgor <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, sd1c, sd1t){
  r*sd1c*n12c*sqrt(n2c)*sqrt(1/s2c + 1/f2c)/n1c/n2c + r*sd1t*n12t*sqrt(n2t)*sqrt(1/s2t + 1/f2t)/n1t/n2t }

md_lgrr <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, sd1c, sd1t){
  r*sd1c*n12c*sqrt(f2c/s2c)/n1c/n2c+r*sd1t*n12t**sqrt(f2t/s2t)/n1t/n2t  }

md_rd <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, sd1c, sd1t){
  r*sd1c*n12c*sqrt(f2c*s2c)/n1c/(n2c^2) + r*sd1t*n12t*sqrt(f2t*s2t)/n1t/(n2t^2) }

smd_lgor <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, sd1c, sd1t){
  #  s1p: pooled standard deviation for outcome 1
  s1p <- sqrt(((n1t - 1)*sd1t^2 + (n1c - 1)*sd1c^2)/(n1t + n1c - 2))
  J <- function(x = 1000){1 - 3/(4*x - 1)}
  v1 <- n1c + n1t - 2
  jv1 <- J(v1)
  jv1*n12c*sqrt(n2c)*r*sd1c*sqrt(1/s2c + 1/f2c)/s1p/n1c/n2c + jv1*n12t*sqrt(n2t)*r*sd1t*sqrt(1/s2t+1/f2t)/s1p/n1t/n2t }

smd_lgrr <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, sd1c, sd1t){
  #  s1p: pooled standard deviation for outcome 1
  s1p <- sqrt(((n1t - 1)*sd1t^2 + (n1c - 1)*sd1c^2)/(n1t + n1c - 2))
  J <- function(x = 1000){1 - 3/(4*x - 1)}
  v1 <- n1c + n1t - 2
  jv1 <- J(v1)
  jv1*n12c*r*sd1c*sqrt(f2c/s2c)/s1p/n1c/n2c + jv1*n12t*r*sd1t*sqrt(f2t/s2t)/s1p/n1t/n2t }

smd_rd <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, sd1c, sd1t){
  #  s1p: pooled standard deviation for outcome 1
  s1p <- sqrt(((n1t - 1)*sd1t^2 + (n1c - 1)*sd1c^2)/(n1t + n1c - 2))
  J <- function(x = 1000){1 - 3/(4*x - 1)}
  v1 <- n1c + n1t - 2
  jv1 <- J(v1)
  jv1*n12c*r*sd1c*sqrt(f2c*s2c)/s1p/n1c/(n2c^2) + jv1*n12t*r*sd1t*sqrt(f2t*s2t)/s1p/n1t/(n2t^2)  }

lgor_lgrr <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c,n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, s1c, s1t, f1t, f1c)  {
  r*n12c*sqrt((1/s1c + 1/f1c)*f2c/s2c)/sqrt(n1c)/n2c + r*n12t*sqrt((1/s1t + 1/f1t)*f2t/s2t)/sqrt(n1t)/n2t  }

lgor_rd <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, s1c, s1t, f1t, f1c) {
  r*n12c*sqrt((s2c*f2c/s1c/f1c)*f2c/s2c)/n1c/n2c + r*n12t*sqrt((s2t*f2t/s1t/f1t)*f2t/s2t)/n1t/n2t }

lgrr_rd <- function(r, n1c, n2c, n1t, n2t, n12c = min(n1c, n2c), n12t = min(n1t, n2t), s2c, s2t, f2c, f2t, s1c, s1t, f1t, f1c) {
  r*n12c*sqrt((f1c*s2c*f2c/s1c)*f2c/s2c)/n1c/(n2c^2) + r*n12t*sqrt((f1t*s2t*f2t/s1t)*f2t/s2t)/n1t/(n2t^2)  }


mix.vcov <- function(d, r, nt, nc, st, sc, n_rt = 0, n_rc = 0, sdt, sdc, type){
  ft <- nt - st
  fc <- nc - sc
  colum.number <- length(type)
  K <- nrow(nt)
  col.vac.number <- (colum.number + 1)*colum.number/2
  type.sm <- matrix("SMD", ncol(nt), ncol(nt))
  for (i in 1:ncol(nt)){
    for (j in 1:ncol(nt)){
      type.sm[i,j] <- paste(type[i], type[j])
    }}
  print(type.sm)
  list.corr.st.varcovar <- list()

  if (n_rt == 0){
    n_rt <- list()
    for (k in 1:K)
    {      temp <- diag(colum.number)
    for (i in 1:colum.number){
      for (j in 1:colum.number){
        temp[i,j] <- min(nt[k, i], nt[k, j])
      }
    }
    n_rt[[k]] <- temp
    }
  }
  if (n_rc == 0){
    n_rc <- list()
    for (k in 1:K)
    {temp <- diag(colum.number)
    for (i in 1:colum.number){
      for (j in 1:colum.number){
        temp[i, j] <- min(nc[k, i], nc[k, j])}
    }
    n_rc[[k]] <- temp
    }
  }

  for (k in 1:K){
    list.corr.st.varcovar[[k]] <- matrix(NA, colum.number, colum.number)

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
          list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- temp$list.md.cov[[1]][1, 2]
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
        list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- temp$list.smd.cov[[1]][1, 2]
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
        list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- temp$list.lgOR.cov[[1]][1, 2]
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
        list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- temp$list.lgRR.cov[[1]][1, 2]
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
        list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- temp$list.rd.cov[[1]][1, 2]
        }

        ##############################
        if ((type[i] == "MD")&&(type[j] == "SMD"))
        {list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- md_smd(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j], sd2t = sdt[k, j], sd2c = sdc[k, j],
                                                                                  sd1t = sdt[k, i], sd1c = sdc[k, i])}
        if ((type[i] == "MD")&&(type[j] == "lgOR"))
        {list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- md_lgor(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                                                                                   sd1t = sdt[k, i], s2t = st[k, j], sd1c = sdc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j])}

        if ((type[i] == "MD")&&(type[j] == "lgRR"))
        {list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- md_lgrr(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                                                                                   sd1t = sdt[k, i], s2t = st[k, j], sd1c = sdc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j])}

        if ((type[i] == "MD")&&(type[j] == "RD"))
        {list.corr.st.varcovar[[k]][i,j]<-list.corr.st.varcovar[[k]][j, i] <- md_rd(r = r[[k]][i, j], n1c=nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                                                                                 sd1t = sdt[k, i], s2t = st[k, j], sd1c = sdc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j])}

        if ((type[i] == "SMD")&&(type[j] == "lgOR"))
        {list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- smd_lgor(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k,j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                                                                                    sd1t = sdt[k, i], s2t = st[k, j], sd1c = sdc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j])}

        if ((type[i] == "SMD")&&(type[j] == "lgRR"))
        {list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- smd_lgrr(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                                                                                    sd1t = sdt[k, i], s2t = st[k, j], sd1c = sdc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j])}
        if ((type[i] == "SMD")&&(type[j] == "RD"))
        {list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- smd_rd(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                                                                                  sd1t = sdt[k,i], s2t = st[k, j], sd1c = sdc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j])}

        if ((type[i] == "lgOR")&&(type[j] == "lgRR"))
        {list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- lgor_lgrr(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                                                                                     s1t = st[k, i], s2t = st[k, j], s1c = sc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j], f1c = fc[k, i], f1t = ft[k, i])}

        if ((type[i] == "lgOR")&&(type[j] == "RD"))
        {list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- lgor_rd(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                                                                                   s1t = st[k,i],s2t=st[k,j], s1c = sc[k, i],s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j], f1c = fc[k, i], f1t = ft[k, i])}

        if ((type[i] == "lgRR")&&(type[j] == "RD"))
        {list.corr.st.varcovar[[k]][i, j] <- list.corr.st.varcovar[[k]][j, i] <- lgrr_rd(r = r[[k]][i, j], n1c = nc[k, i], n2c = nc[k, j], n1t = nt[k, i], n2t = nt[k, j], n12c = n_rc[[k]][i, j], n12t = n_rt[[k]][i, j],
                                                                                   s1t = st[k, i], s2t = st[k, j], s1c = sc[k, i], s2c = sc[k, j], f2c = fc[k, j], f2t = ft[k, j], f1c = fc[k, i], f1t = ft[k, i])}
      }
    }
  }

  corr.st.varcovar <- matrix(unlist(lapply(1:K, function(k, list.corr.st.varcovar){sm2vec(list.corr.st.varcovar[[k]], diag = TRUE)},list.corr.st.varcovar = list.corr.st.varcovar)), K, col.vac.number, byrow = TRUE)
  list(list.mix.cov = list.corr.st.varcovar, mix.cov = corr.st.varcovar)
}
