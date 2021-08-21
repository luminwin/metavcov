r.vcov <- function(n, corflat, method = "average")
{ # requireNamespace("corpcor")
 #library(corpcor)
  vec2sm <- function (vec, diag = FALSE, order = NULL)
  {
    n = (sqrt(1 + 8 * length(vec)) + 1)/2
    if (diag == TRUE)
      n <- n - 1
    if (ceiling(n) != floor(n))
      stop("Length of vector incompatible with symmetric matrix")
    m <- matrix(NA, nrow = n, ncol = n)
    lo <- lower.tri(m, diag)
    if (is.null(order)) {
      m[lo] <- vec
    }
    else {
      vec.in.order <- rep(NA, length(order))
      vec.in.order[order] <- vec
      m[lo] <- vec.in.order
    }
    for (i in 1:(n - 1)) for (j in (i + 1):n) m[i, j] <- m[j,
                                                          i]
    return(m)
  }

  function (m, diag = FALSE)
  {
    return(as.vector(m[lower.tri(m, diag)]))
  }

  colum.number <- ncol(corflat)
  V <- round(sqrt(2*colum.number + 1/4) + 0.5,0)

  if (length(n) == 1)
  {K <- 1}
  else
  { K <- nrow(corflat)}


  col.vac.number <- (colum.number+1)*colum.number/2
  corr.st.varcovar <- matrix(NA, K, col.vac.number)
  list.corr.st.varcovar <- list()

  z.corr.st.varcovar <- matrix(NA, K, col.vac.number)
  list.z.corr.st.varcovar <- list()
  zz.corflat <- corflat
  for (i in 1:K){
    for (j in 1:colum.number)
    {
      zz.corflat[i, j] <- 0.5*log((1 + corflat[i, j])/(1 - corflat[i, j]))
    }
  }

  if (method == "each")
  {corflat <- corflat
   z.corflat <- corflat
   for (i in 1:K){
     for (j in 1:colum.number)
     {
       z.corflat[i, j] <- 0.5*log((1 + corflat[i, j])/(1 - corflat[i, j]))
     }
   }
  }

  if (method == "average")
  {temp <- unlist(lapply(1:colum.number, function(i, corflat,n){weighted.mean(corflat[, i], n)}, corflat = corflat, n = n))
   corflat <- matrix(rep(temp, K), K, colum.number, byrow = "TRUE")

   temp <- unlist(lapply(1:colum.number, function(i, corflat,n){weighted.mean(corflat[,i], n)}, corflat = zz.corflat, n = n))
   z.corflat <- matrix(rep(temp, K), K, colum.number, byrow = "TRUE")}

  sum <- 0
  iii <- 1
  foot <- matrix(NA, colum.number, 2)
  f <- NA
  for (ii in 1:(V-1))
  {
    f[ii] <- (V-ii)
    foot[((sum(f) - f[ii] + 1):(sum(f))), 1] <- c(((ii + 1):V))
    foot[((sum(f) - f[ii] + 1):(sum(f))), 2] <- ii

  }

  sub <- as.data.frame(matrix(NA, col.vac.number, 4))
  d <- colum.number
  colnames(sub) <- c("s","t","u","v")
  for (t in 1:(V-1))
  {
    for (s in (t+1):V)
    {
      sub[((sum + 1):(sum + d)), 1] <- s
      sub[((sum + 1):(sum + d)), 2] <- t
      sub[((sum + 1):(sum + d)), 3] <- foot[((colum.number - d + 1):colum.number), 1]
      sub[((sum + 1):(sum + d)), 4] <- foot[((colum.number - d + 1):colum.number), 2]
      sum <- sum+d
      d <- d - 1

    }
  }

  for (k in 1:K)
  {
    r <- vec2sm(t(corflat[k, ]), diag = FALSE)
    r[is.na(r)] <- 1

    for  (ii in 1:col.vac.number)
    {

      corr.st.varcovar[k, ii] <- (0.5*r[sub[ii,]$s,sub[ii,]$t]*r[sub[ii,]$u, sub[ii,]$v]*((r[sub[ii,]$s,sub[ii,]$u])^2
                                                                                      + (r[sub[ii,]$s,sub[ii,]$v])^2
                                                                                      + (r[sub[ii,]$t,sub[ii,]$u])^2+(r[sub[ii,]$t,sub[ii,]$v])^2) + r[sub[ii,]$s,sub[ii,]$u]*r[sub[ii,]$t,sub[ii,]$v]
                               + r[sub[ii,]$s, sub[ii,]$v]*r[sub[ii,]$t, sub[ii,]$u]
                               - (r[sub[ii,]$s, sub[ii,]$t]*r[sub[ii,]$s, sub[ii,]$u]*r[sub[ii,]$s, sub[ii,]$v]
                                 + r[sub[ii,]$t, sub[ii,]$s]*r[sub[ii,]$t, sub[ii,]$u]*r[sub[ii,]$t, sub[ii,]$v]
                                 + r[sub[ii,]$u, sub[ii,]$s]*r[sub[ii,]$u, sub[ii,]$t]*r[sub[ii,]$u, sub[ii,]$v]
                                 + r[sub[ii,]$v, sub[ii,]$s]*r[sub[ii,]$v, sub[ii,]$t]*r[sub[ii,]$v, sub[ii,]$u]))/n[k]
    }

    r <- vec2sm(t(z.corflat[k,]), diag = FALSE)
    r[is.na(r)] <- 1

    for  (ii in 1:col.vac.number)
    {
      z.corr.st.varcovar[k, ii] <- (0.5*r[sub[ii,]$s, sub[ii,]$t]*r[sub[ii,]$u,sub[ii,]$v]*((r[sub[ii,]$s,sub[ii,]$u])^2
                                                                                        + (r[sub[ii,]$s,sub[ii,]$v])^2
                                                                                        + (r[sub[ii,]$t,sub[ii,]$u])^2 + (r[sub[ii,]$t, sub[ii,]$v])^2) + r[sub[ii,]$s, sub[ii,]$u]*r[sub[ii,]$t, sub[ii,]$v]
                                 + r[sub[ii,]$s, sub[ii,]$v]*r[sub[ii,]$t, sub[ii,]$u]
                                 - (r[sub[ii,]$s, sub[ii,]$t]*r[sub[ii,]$s, sub[ii,]$u]*r[sub[ii,]$s, sub[ii,]$v]
                                   + r[sub[ii,]$t, sub[ii,]$s]*r[sub[ii,]$t, sub[ii,]$u]*r[sub[ii,]$t, sub[ii,]$v]
                                   + r[sub[ii,]$u, sub[ii,]$s]*r[sub[ii,]$u, sub[ii,]$t]*r[sub[ii,]$u, sub[ii,]$v]
                                   + r[sub[ii,]$v, sub[ii,]$s]*r[sub[ii,]$v, sub[ii,]$t]*r[sub[ii,]$v, sub[ii,]$u]))/n[k]/(1 - (r[sub[ii,]$s,sub[ii,]$t])^2)/(1 - (r[sub[ii,]$u,sub[ii,]$v])^2)
    }

    list.corr.st.varcovar[[k]] <- vec2sm(corr.st.varcovar[k,], diag = TRUE)


    list.z.corr.st.varcovar[[k]] <- vec2sm(z.corr.st.varcovar[k,], diag = TRUE)
  }

  temp <- list.z.corr.st.varcovar
  list.z.corr.st.varcovar <- lapply(1:K, function(k, temp, n){tempp <- diag(temp[[k]])
                                                         diag(temp[[k]]) <- tempp*(n[k]/(n[k] - 3))
                                                         temp[[k]]}, temp <- temp, n = n)

  z.corr.st.varcovar <- matrix(unlist(lapply(1:K, function(k, list.z.corr.st.varcovar){sm2vec(list.z.corr.st.varcovar[[k]], diag = TRUE)}, list.z.corr.st.varcovar = list.z.corr.st.varcovar)), K, col.vac.number, byrow = TRUE)

  return(list(list.rcov = list.corr.st.varcovar, rcov = corr.st.varcovar, zr = zz.corflat, zcov = z.corr.st.varcovar, list.zcov = list.z.corr.st.varcovar))

}
