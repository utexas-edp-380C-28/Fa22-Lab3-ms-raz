#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 3: An Overview of Linear Regression 
#        and Data Generation
#
#' Discuss what this function does...
#'
#' @param n number of observations requested
#' @param mu mean parameters
#' @param Sigma covariance matrix
#' @return the data

rmvnorm <- function(n, mu, Sigma) {
    try(if
        (length(mu) != ncol(Sigma)) stop("Your input dimensions don't match"))
    # create random variable to start with that is n x p
    Z = matrix(rnorm(n*length(mu)), n, length(mu))
    # create a variable that is mvn that has specified mus and sigma
    X = matrix(1, n, 1) %*% t(mu) + Z %*% chol(Sigma) 
    # create column names here so you dont have to in lm
    colnames(X) <- c(sprintf("X%1d", 1:length(mu)))
    return(X)
}