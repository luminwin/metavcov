smd.vcov <- function(nt, nc, d, r, n_rt = 0, n_rc = 0){

  if (length(as.vector(d)) == length(as.matrix(d)[, 1]))   {colum.number <- 1}  else {colum.number <- ncol(d)}
  if (length(r) == 1)   {K <- 1}    else   {    K <- nrow(d)}
  col.vac.number <- (colum.number + 1)*colum.number/2

  list.corr.st.varcovar <- list()

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
          temp[i,j] <- min(nc[k, i], nc[k, j])}
      }
      n_rc[[k]] <- temp
      }
    }
    J <- function(x = 1000){1 - 3/(4*x - 1)}
    for (k in 1:K){
      list.corr.st.varcovar[[k]] <- matrix(NA, colum.number, colum.number)
      #  i=j=2
      for (i in 1:colum.number){
        for (j in 1:colum.number)
        {
          vi <- nc[k,i] + nt[k,i] - 2
          jvi <- J(vi)
          vj <- nc[k,j] + nt[k,j] - 2
          jvj <- J(vj)
          ki <- (2*nt[k, i] - 2)/((nc[k, i] + nt[k, i] - 2)^2) + (2*nc[k, i] - 2)/((nc[k, i] + nt[k, i] - 2)^2)
          kj <- (2*nt[k, j] - 2)/((nc[k, j] + nt[k, j] - 2)^2) + (2*nc[k, j] - 2)/((nc[k, j] + nt[k, j] - 2)^2)
          kij <- 2*(nt[k, i]*nt[k, j]/(nt[k, i] + nt[k, j] - 1) + nc[k, i]*nc[k, j]/(nc[k, i] + nc[k,j] - 1) - 2)/(nc[k, i] + nt[k, i] - 2)/(nc[k, j] + nt[k, j] - 2)
          list.corr.st.varcovar[[k]][i, j] <- (n_rc[[k]][i, j]/nc[k, i]/nc[k, j] + n_rt[[k]][i, j]/nt[k, i]/nt[k, j])*r[[k]][i, j] + kij*((r[[k]][i, j])^2)*d[k, i]*d[k, j]*jvi*jvj*sqrt((vi/(vi - 2) - 1/(jvi^2))*(vj/(vj - 2) - 1/(jvj^2)))/sqrt(ki*kj)
        }
      }
    }
    corr.st.varcovar <- matrix(unlist(lapply(1:K, function(k, list.corr.st.varcovar){sm2vec(list.corr.st.varcovar[[k]], diag = TRUE)}, list.corr.st.varcovar = list.corr.st.varcovar)), K, col.vac.number, byrow = TRUE)
    list(list.smd.cov = list.corr.st.varcovar, smd.cov = corr.st.varcovar)
  }
