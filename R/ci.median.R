ci.median <- function(x, conf=0.95){ # from R/asbio
  n <- nrow(as.matrix(x))
  if(qbinom((1 - conf)/2, n, 0.5) == 0){ stop("CI not calculable") }
  L <- qbinom((1 - conf)/2, n, 0.5)
  U <- n - L + 1
  if (L >= U){ stop("CI not calculable") }
  order.x <- sort(x)
  ci <- c(lower = order.x[L], upper = order.x[n - L + 1])
  return(ci)
}