make.EM.alt.formula <- function(formula, X){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- "y ~"
  this.formula <- as.formula(paste(this.formula.string, paste(gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE), collapse=" + ")))
  return(this.formula)
}
