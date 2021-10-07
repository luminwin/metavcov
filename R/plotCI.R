ciplotgg <- function(d, yinter = 0, eff.nam){
  # d is a data frame with 4 columns
  # d$x gives variable names
  # d$y gives center point
  # d$ylo gives lower limits
  # d$yhi gives upper limits
 # require(ggplot2)
  p <- ggplot2::ggplot(d, ggplot2::aes(x=x, y=y, ymin=ylo, ymax=yhi))+
    ggplot2::geom_pointrange()+
    ggplot2::geom_hline(yintercept = yinter, linetype=2)+
    ggplot2::coord_flip()+
    ggplot2::xlab(' ') +
    ggplot2::ylab(eff.nam) +
    ggplot2::theme_classic() +
    ggplot2::geom_point(ggplot2::aes(x=d[which(d$x == "Overall"),]$x,
                                     y=d[which(d$x == "Overall"),]$y), colour="blue", size = 3, shape =15)
  return(p)
}

plotCI <- function(y, v, name.y = NULL, name.study = NULL, y.all, y.all.se, hline = 0, up.bound = Inf, low.bound = -Inf, return.data = FALSE) {
  K <- nrow(y)
  p <- ncol(y)
  if (length(hline) == 1) {hline = rep(hline, p)}
  if (is.null(name.study)) { if (is.null(rownames(y))) { name.study <- paste("Study", 1:K)} else { name.study <- rownames(y) } }
  if (is.null(name.y)) { if (is.null(colnames(y))) { name.y <- paste("Outcome", 1:p)} else { name.y <- colnames(y) } }
  y.var <- do.call('rbind', lapply(1:K, function(i) {
    diag(v[[i]])
  }) )

  z <- qnorm(0.975)
  obj <- plots <- list()

  for (i in 1:p){
    d <- data.frame(x = name.study,
                    y = y[ ,i],
                    ylo = y[ ,i] - z*sqrt(y.var[ ,i]),
                    yhi= y[ ,i] + z*sqrt(y.var[ ,i])
                    )
    d <- rbind(data.frame( x = "Overall",
                               y = y.all[i],
                               ylo = y.all[i] - z*(y.all.se[i]),
                               yhi= y.all[i] + z*(y.all.se[i])
                              ),d)
    d$x <- factor(d$x, levels = d$x)
    d[which(d$ylo < low.bound),] <- low.bound
    d[which(d$yhi > up.bound),] <- up.bound
    obj[[i]] <- d
    colnames(obj[[i]]) <- c("name.study", "Estimate", "95%ci.lb",  "95%ci.ub" )
    (plots[[i]] <- (ciplotgg(d, yinter = hline[i], eff.nam <- name.y[i])))
  }
  names(plots) <- paste("Plotting", name.y)
  if (return.data == FALSE) { return(plots)} else {
return(list(data = obj,plots = plots))}
}
