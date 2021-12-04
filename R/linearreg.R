#'Linearreg
#'
#'Fits a linear regression model and does inference
#'
#'@param y vector of outcome values
#'
#'@param X matrix whose columns are predictor values, ordered like y. If missing fits intercept-only
#'
#'@param prednames optional names for predictors, in order of columns of X
#'
#'@examples
#'y <- c(1,2,3,4,5,6,7,8,9,10)
#'x1 <- c(10,8,9,6,5,9,3,2,1,3)
#'x2 <- c(6,4,8,6,2,9,0,1,4,6)
#'prednames <- c("Speed", "Weight")
#'linearreg(y, cbind(x1, x2))
#'
#'@export
#'
linearreg <- function(y, X, prednames) {
  #If no X given, fit intercept model
  if (missing(X)) {
    X <- rep(1,length(y))
    ncol <- 0
  }
  else if (is.null(ncol(X))) {
    ncol <- 1
  }
  else {
    ncol <- ncol(X)
  }
  #define design matrix
  X <- cbind(rep(1,length(y)), X)
  #assign generic variable names if not given
  if (missing(prednames)) {
    prednames <- c("Intercept")
    if (ncol > 0) {
      for (i in 1:ncol) {
        prednames <- c(prednames, paste0("x",i))
      }
    }
  }
  else {
    prednames <- c("Intercept", prednames)
  }

  #must be invertible or cannot solve
  if (det(t(X) %*% X) == 0) {
    stop("Error: XTX matrix not invertible")
  }

  beta <- solve(t(X) %*% X) %*% t(X) %*% y

  estmeans <- X %*% beta

  #model variance estimate
  varest <- sum((y - estmeans)^2)/(length(y)-ncol-1)
  sesmat <- solve(t(X) %*% X)
  #standard errors of each predictor
  ses <- sqrt(varest * diag(sesmat))
  tvals <- beta/ses
  pvals <- 2 * pt(abs(tvals), df = length(y) - ncol - 1, lower.tail = FALSE)

  beta <- round(beta, 5)
  ses <- round(ses, 5)
  tvals <- round(tvals, 5)
  pvals <- round(pvals, 5)

  #makes dataframe of coefficients and corresponding values
  coefs <- data.frame(prednames, beta, ses, tvals, pvals)
  colnames(coefs) <- c("Predictor", "Estimate", "S.E.", "t value", "p")

  R2 <- 1 - sum((y - estmeans)^2)/sum((y - mean(y))^2)
  adjR2 <- 1 - (1 - R2)*(length(y)-1)/(length(y)-ncol-1)
  R2 <- round(R2,5)
  adjR2 <- round(adjR2, 5)

  dfp <- ncol
  dfn <- length(y)-ncol-1
  f <- (sum((estmeans - mean(y))^2)/dfp) / (sum((y - estmeans)^2)/dfn)
  fp <- pf(f, dfp, dfn, lower.tail = FALSE)
  f <- round(f, 5)
  fp <- round(fp, 5)
  freturn <- paste(f, "on")
  freturn <- paste(freturn, dfp)
  freturn <- paste(freturn, "and")
  freturn <- paste(freturn, dfn)
  freturn <- paste(freturn, "degrees of freedom")
  freturn <- paste(freturn, ", p-value:")
  freturn <- paste(freturn, fp)

  #only include adjusted R-squared in output if more than one predictor
  if(ncol > 1) {
    toreturn <- list("Coefficients" = coefs, "R-squared" = R2, "Adjusted R-squared" = adjR2,
                     "F-test" = freturn)
  }
  else {
    toreturn <- list("Coefficients" = coefs, "R-squared" = R2, "F-test" = freturn)
  }

  return(toreturn)
}
