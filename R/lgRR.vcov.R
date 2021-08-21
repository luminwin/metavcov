lgRR.vcov <- function(r, nt, nc, st, sc, n_rt = 0, n_rc = 0)
{
  ft <- nt - st
  fc <- nc - sc
  if (length(as.vector(ft)) == length(as.matrix(ft)[, 1]))    {colum.number <- 1}    else   {    colum.number <- ncol(ft)}

  if (length(as.vector(ft)) == length(as.matrix(ft)[, 1]))    {K <- length(ft)}   else   {     K <- nrow(ft)}
  col.vac.number <- (colum.number + 1)*colum.number/2
  if (sum(unlist(n_rt)) == 0){
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
  if (sum(unlist(n_rc)) == 0){
    n_rc <- list()
    for (k in 1:K)
    {temp <- diag(colum.number)
    for (i in 1:colum.number){
      for (j in 1:colum.number){
        temp[i,j] <- min(nc[k, i],nc[k, j])}
    }
    n_rc[[k]] <- temp
    }
  }
  list.corr.st.varcovar <- list()
  for (k in 1:K){
    list.corr.st.varcovar[[k]] <- matrix(NA, colum.number, colum.number)
    for (i in 1:colum.number){
      for (j in 1:colum.number)
      {list.corr.st.varcovar[[k]][i, j]<-r[[k]][i, j]*n_rc[[k]][i, j]*sqrt(fc[k, i]*fc[k, j]/sc[k, i]/sc[k, j])/sqrt(nc[k, i]*nc[k, j]) + r[[k]][i, j]*n_rt[[k]][i, j]*sqrt(ft[k, i]*ft[k, j]/st[k, i]/st[k, j])/sqrt(nt[k, i]*nt[k, j])
      }
    }
  }

  lgRR <- matrix(NA, K, colum.number)
  for (k in 1:K) {
    for (i in 1:colum.number){
      lgRR[k, i] <- log((st[k, i]/nt[k, i])/(sc[k, i]/nc[k, i]))
    }}

  corr.st.varcovar <- matrix(unlist(lapply(1:K, function(k, list.corr.st.varcovar){sm2vec(list.corr.st.varcovar[[k]], diag = TRUE)}, list.corr.st.varcovar = list.corr.st.varcovar)), K, col.vac.number, byrow = TRUE)
  list(list.lgRR.cov = list.corr.st.varcovar, lgRR.cov = corr.st.varcovar, lgRR = lgRR)
}
