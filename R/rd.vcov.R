rd.vcov <- function(r, nt, nc, st, sc, n_rt = NA, n_rc = NA)
{
  ft <- nt - st
  fc <- nc - sc
  if (length(as.vector(ft)) == length(as.matrix(ft)[, 1]))   {colum.number <- 1}   else   {     colum.number <- ncol(ft)}

  if (length(as.vector(ft)) == length(as.matrix(ft)[, 1]))  {K <- length(ft)}  else  {    K <- nrow(ft)}
  col.vac.number <- (colum.number + 1)*colum.number/2

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

  list.corr.st.varcovar <- list()
  for (k in 1:K){
    list.corr.st.varcovar[[k]] <- matrix(NA, colum.number, colum.number)
    for (i in 1:colum.number){
      for (j in 1:colum.number)
      {list.corr.st.varcovar[[k]][i, j] <- unlist(r[[k]][i, j]*n_rc[[k]][i, j]*sqrt(fc[k, i]*fc[k, j]*sc[k, i]*sc[k, j])/((nc[k, i]*nc[k, j])^2)+r[[k]][i, j]*n_rt[[k]][i, j]*sqrt(ft[k, i]*ft[k, j]*st[k, i]*st[k, j])/((nt[k, i]*nt[k, j])^2))
      }
    }
  }
  rd <- matrix(NA, K, colum.number)
  for (k in 1:K) {
    for (i in 1:colum.number){
      rd[k, i] <- unlist((st[k, i]/nt[k, i]) - (sc[k, i]/nc[k, i]))
    }}
  corr.st.varcovar <- matrix(unlist(lapply(1:K, function(k, list.corr.st.varcovar){smTovec(list.corr.st.varcovar[[k]])}, list.corr.st.varcovar = list.corr.st.varcovar)), K, col.vac.number, byrow = TRUE)
  list(list.vcov = list.corr.st.varcovar, matrix.vcov = corr.st.varcovar, ef = as.data.frame(rd))
}
