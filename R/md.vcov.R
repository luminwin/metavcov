md.vcov <- function(r, nt, nc, n_rt = 0, n_rc = 0, sdt, sdc)
{
  if (length(as.vector(nt)) == length(as.matrix(nt)[, 1]))    {colum.number <- 1}    else    {      colum.number <- ncol(nt)}
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
        temp[i, j] <- min(nc[k, i], nc[k, j])}
    }
    n_rc[[k]] <- temp
    }
  }
  list.corr.st.varcovar <- list()
  for (k in 1:K){
    list.corr.st.varcovar[[k]] <- matrix(NA, colum.number, colum.number)
    for (i in 1:colum.number){
      for (j in 1:colum.number)
      {list.corr.st.varcovar[[k]][i, j] <- r[[k]][i, j]*n_rc[[k]][i, j]*sdc[k, i]*sdc[k, j]/(nc[k, i]*nc[k, j]) + r[[k]][i, j]*n_rt[[k]][i, j]*sdt[k, i]*sdt[k, j]/(nt[k, i]*nt[k, j])
      }
    }
  }
  corr.st.varcovar <- matrix(unlist(lapply(1:K, function(k, list.corr.st.varcovar){sm2vec(list.corr.st.varcovar[[k]], diag = TRUE)}, list.corr.st.varcovar = list.corr.st.varcovar)), K, col.vac.number, byrow = TRUE)
  list(list.md.cov = list.corr.st.varcovar, md.cov = corr.st.varcovar)
}
