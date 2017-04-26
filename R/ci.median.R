ci.median <- function(x, conf=0.95){ # from R/asbio
  n <- nrow(as.matrix(x))
  if(qbinom((1 - conf)/2, n, 0.5) == 0){ stop("CI not calculable") }
  L <- qbinom((1 - conf)/2, n, 0.5)
  U <- n - L + 1
  if (L >= U){ stop("CI not calculable") }
  order.x <- sort(x)
  results <- list()
  results$head <- paste(paste(as.character(conf * 100), "%", sep = ""), c("Confidence interval for population median"))
  results$ci <- c(median = median(x), lower = order.x[L], upper = order.x[n - L + 1])
  results$ends <- c("Estimate", paste(as.character(c((1 - conf)/2, 1 - ((1 - conf)/2)) * 100), "%", sep = ""))
  results$coverage <- 1 - (2 * pbinom(q = L - 1, n, 0.5))
  class(results) <- "ci"
  return(results)
}