lgOR.vcov <- function(r, nt, nc, st, sc, n_rt = 0, n_rc = 0){
  ft <- nt - st
  fc <- nc - sc
  if (length(as.vector(ft)) == length(as.matrix(ft)[, 1]))    {colum.number <- 1} else    {     colum.number <- ncol(ft)}

  if (length(as.vector(ft)) == length(as.matrix(ft)[, 1]))   {K <- length(ft)}   else   {     K <- nrow(ft)}
  col.vac.number <- (colum.number + 1)*colum.number/2
  if (sum(unlist(n_rt)) == 0){
    n_rt <- list()
    for (k in 1:K)
    {      temp <- diag(colum.number)
    for (i in 1:colum.number){
      for (j in 1:colum.number){
        temp[i, j] <- min(nt[k, i], nt[k, j])
      }
    }
    n_rt[[k]] <- temp
    }
  }
  if (sum(unlist(n_rc)) == 0){
    n_rc <- list()
    for (k in 1:K)
    {temp <- diag(colum.number)
    for (i in 1:colum.number){
      for (j in 1:colum.number){
        temp[i,j] <- min(nc[k,i], nc[k,j])}
    }
    n_rc[[k]] <- temp
    }
  }
  list.corr.st.varcovar <- list()
  for (k in 1:K){
    list.corr.st.varcovar[[k]] <- matrix(NA, colum.number, colum.number)
    for (i in 1:colum.number){
      for (j in 1:colum.number)
      {list.corr.st.varcovar[[k]][i, j] <- r[[k]][i, j]*n_rc[[k]][i, j]*sqrt((1/sc[k, i] + 1/fc[k, i])*(1/sc[k, j] + 1/fc[k, j]))/sqrt(nc[k, i]*nc[k, j]) + r[[k]][i, j]*n_rt[[k]][i, j]*sqrt((1/st[k, i] + 1/ft[k, i])*(1/st[k, j] + 1/ft[k, j]))/sqrt(nt[k, i]*nt[k, j])
      }
    }
  }
  corr.st.varcovar <- matrix(unlist(lapply(1:K,function(k, list.corr.st.varcovar){sm2vec(list.corr.st.varcovar[[k]], diag=TRUE)}, list.corr.st.varcovar = list.corr.st.varcovar)), K, col.vac.number, byrow = TRUE)
  lgOR <- matrix(NA, K, colum.number)
  for (k in 1:K) {
    for (i in 1:colum.number){
      lgOR[k, i] <- log((st[k, i]/ft[k, i])/(sc[k, i]/fc[k, i]))
    }}
  list(list.lgOR.cov = list.corr.st.varcovar, lgOR.cov = corr.st.varcovar, lgOR = lgOR)
}