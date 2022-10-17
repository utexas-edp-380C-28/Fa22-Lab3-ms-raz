#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 3: An Overview of Linear Regression 
#        and Data Generation
#
# Name:  RAZ: Rebecca A. Zarate
#---------------------------------------------------#

## SETUP -------------------------------------------
# Source rmvnorm()
source("scripts/rmvnorm.R")

## Begin Code --------------------------------------

n <- 100
reps <- 500

#########-----------------------1.a ------------------------##########
p_x <- list(
    mu = -2,
    sigma = 3
)

p_y <- list(
    sigma = 7, 
    b0 = 12, 
    b1= 4
)

# Create a function that generates X and Y based on Equation (9)
genFun <- function(n, p_x, p_y){
    X <- rnorm(n, p_x$mu, p_x$sigma)
    Y <- rnorm(n, (p_y$b0 + p_y$b1*X), p_y$sigma) #eq 9
    dat <- cbind(X = X, Y = Y)
    return(dat)
}

genDat <- genFun(n, p_x, p_y)

#3 Create an analyze function that takes one data set and outputs a one-dimensional vector (i.e., c()) with the following labels and ordering: 
analysisFun <- function(genDat){
    m_x = mean(genDat[,"X"])
    m_y = mean(genDat[,"Y"])
    s_x = sd(genDat[,"X"])
    s_y = sd(genDat[,"Y"])
    cor = cor(genDat[,"X"], genDat[,"Y"])
    mod = lm(genDat[,"Y"] ~ genDat[,"X"])
    output = c(m_x = m_x, 
               m_y = m_y, 
               s_x = s_x, 
               s_y = s_y, 
               cor = cor, 
               b0 = mod$coefficients[[1]], 
               b1 = mod$coefficients[[2]], 
               s_e = sigma(mod)) # residual sd
    return(output)
}

set.seed(17290)

# Set up and run a simulation to generate 500 reps of n = 100 using the replicate function
smallSim <- replicate(reps, genFun(n, p_x, p_y))

# Use the analyze function to analyze and compute the empirical mean and standard deviation of the parameter estimates
results <- apply(smallSim, 3, analysisFun)

rowMeans(results)
apply(results, 1, sd)

# Create a list(mean, sd) that holds the mean and standard deviation of the parameter estimates from the simulation. Make sure to label each item of the list. Print the list and include its output as a comment in your code.
num1 <- list(
    mean = rowMeans(results),
    sd = apply(results, 1, sd)
)

# m_x        m_y        s_x        s_y        cor         b0         b1        s_e 
# -2.0061317  3.9810216  3.0034995 13.8596003  0.8631812 11.9780298  3.9873856  6.9749258 
# > apply(results, 1, sd)
# m_x        m_y        s_x        s_y        cor         b0         b1        s_e 
# 0.31115637 1.41924855 0.20911389 1.03767609 0.02432744 0.86119359 0.24140973 0.49070350 

#########-----------------------2.a ------------------------##########

# 2.4 Construct a function that produces a solution to OLS regression. The function should map onto the following inputs: function(y, X)

myLM <- function(Y, X){
    betas <- solve((t(X) %*% X)) %*% t(X) %*% Y # eq 11
    RSE <- sqrt((t( Y - X %*% betas) %*% (Y - X %*% betas))/(length(Y) - ncol(X))) # eq 13 
    SEs <- t(sqrt(RSE^2 %*% diag(solve(t(X) %*% X)))) # eq 17 but read on
    R2 <- (t(betas) %*% cov(X) %*% betas) /(t(betas) %*% cov(X) %*% betas + RSE^2)
    tval <- betas/SEs
    PrT <- 2 * pt(abs(tval), length(Y)-ncol(X), lower.tail = FALSE)
    
    tmp <- cbind(
        Estimate = betas[,1], 
        SE = SEs[, 1], 
        't value' = tval[, 1], 
        'p value' = PrT[, 1]
    )
    
    rownames(tmp) <- c(sprintf("b%1d", 0:(length(betas)-1)))
    
    results <- rbind(tmp, 
                     'SD(e)' = c(RSE, rep(NA, ncol(tmp)-1)), 
                     'R2' = c(R2, rep(NA, ncol(tmp)-1))
    )
    
    return(results)
}

X <- as.matrix(mtcars[, c("wt", "cyl", "gear")])
X <- cbind(matrix(1, NROW(mtcars)), X) # add the row of 1s for the intercept 

myLM(Y = mtcars$mpg, X)

# Estimate        SE    t value      p value
# b0    42.3863641 4.3789952  9.6794726 1.965929e-10
# b1    -3.3920819 0.8208025 -4.1326409 2.941570e-04
# b2    -1.5280010 0.4197533 -3.6402363 1.092609e-03
# b3    -0.5228629 0.7788703 -0.6713094 5.075244e-01
# SD(e)  2.5921848        NA         NA           NA
# R2     0.8182681        NA         NA           NA

# check
# summary(lm(mpg ~ wt + cyl + gear, data = mtcars))

#########-----------------------3------------------------##########

#-------------------------3.1a------------------------------------------#
set.seed(21389)

# It will be advantageous to make a new function to translate the parameter list to the inputs expected by rmvnorm().
# translate parameters into form rvnorm works with
translateFun <- function(mu, sigma, rho){
    translated = list(
        mu = mu, 
        # create covariance matrix
        Sigma = matrix(diag(sigma), length(mu)) %*% 
            matrix(c(1, rho, rho, 1), nrow = 2) %*% 
            matrix(diag(sigma), length(mu)),
        rho = rho 
    )
}

# translate 
params <- translateFun(mu = c(5, 10), sigma = c(1,2), rho = .3)

n <- 100000

# generate correlated Xs from mus and covar matrix after the parameters have been translated using the translate function
X <- rmvnorm(n, params$mu, params$Sigma)
X <- cbind(matrix(1, NROW(X)), X) # add the row of 1s for the intercept 

#-------------------------3.2b------------------------------------------#
# Method 1: Setting the Total Variance Explained
set.seed(23921)

# new translate function that gets all params and saves them in a form rmvnorm can deal with and generates the covariance matrix
translateFun <- function(mu, sigma, rho, betas = NULL, R2 = NULL, muY = NULL, sigmaY = NULL){
    translated = list(
        mu = mu, 
        Sigma = matrix(diag(sigma), length(mu)) %*% 
            matrix(c(1, rho, rho, 1), nrow = 2) %*% 
            matrix(diag(sigma), length(mu)),
        rho = rho, 
        sigma = sigma,
        betas = betas,
        R2 = R2, 
        muY = muY, 
        sigmaY = sigmaY
    )
}

# translate parameters into form rmvnorm takes
params <- translateFun(mu = c(5, 10), 
                       sigma = c(1,2), 
                       rho = .3, 
                       betas = c(1, 1), 
                       R2 = 0.6, 
                       muY = 10, 
                       sigmaY = 5)

# matrix way
varY <- t(params$betas) %*% params$Sigma %*% params$betas # eq 22

# determine variance of error SD(e)
varE <- varY * (1/params$R2 - 1) # eq 23

# expectation of Y
Ey <- params$muY - t(params$mu) %*% params$betas

# [-1] bc of intercept, need to remove when multiplying by betas
Y <- Ey[1] + X[,-1] %*% params$betas + rnorm(0, sqrt(varE), n = nrow(X))

# analyze OLS using myLM function
myMod <- myLM(Y = Y, X = X)
myMod

# 
#           Estimate          SE   t value p value
# b0    -4.9370600 0.040463429 -122.0129       0
# b1     0.9891621 0.006771543  146.0763       0
# b2     0.9997359 0.003373688  296.3332       0
# SD(e)  2.0361731          NA        NA      NA
# R2     0.5978875          NA        NA      NA

# create a difference bw the population parameters for Y and the results from you generated data

# myMod
#           Estimate          SE      t value   p value
# b0    11.0629399882 0.040463429 273.40589430 0.0000000
# b1     0.9891621375 0.006771543 146.07632891 0.0000000
# b2    -0.0002641061 0.003373688  -0.07828409 0.9376022
# SD(e)  2.0361730796          NA           NA        NA
# R2     0.1900023581          NA           NA        NA

pop <- c(Ey, params$betas, varE, params$R2)
estimate <- myMod[,1]
difference <- pop - estimate

difTab <- cbind(pop, estimate, difference)

#           pop   estimate    difference
# b0    -5.000000 -4.9370600 -0.0629399882
# b1     1.000000  0.9891621  0.0108378625
# b2     1.000000  0.9997359  0.0002641061
# SD(e)  4.133333  2.0361731  2.0971602537
# R2     0.600000  0.5978875  0.0021125328

#-------------------------3.2c------------------------------------------#
# Method 2: Setting the Marginal Correlations
set.seed(123782)

newParams <- list(
    muX = c(5, 10),
    sigmaX1 = 1, 
    sigmaX2 = 2, 
    rhoY1 = .3, 
    rhoY2 = -.4, 
    muY = 10, 
    sigmaY = 5
)

# corrleation matrix of y, x1, x2
Rm <- matrix(c(1, newParams$rhoY1, newParams$rhoY2, 
               newParams$rhoY1, 1, params$rho, 
               newParams$rhoY2, params$rho, 1), nrow = 3)

# full covariance matrix
covariancematirx <- diag(c(newParams$sigmaY, 
                           newParams$sigmaX1, 
                           newParams$sigmaX2)) %*% 
    Rm %*% 
    diag(c(newParams$sigmaY,
           newParams$sigmaX1,
           newParams$sigmaX2))

# solve for betas
betaZ <- solve(covariancematirx[-1, -1]) %*% covariancematirx[1, -1]

# ereor variance
sigmaE <- covariancematirx[1,1] - t(betaZ) %*%  covariancematirx[1, -1]

# expectation of Y
Ey <- newParams$muY - t(newParams$muX) %*% betaZ

# generate y
yz <- Ey[1] + X[,-1] %*% betaZ + rnorm(n, 0, sqrt(sigmaE))

# analyze
myMod2 <- myLM(yz,  X)

#           Estimate          SE   t value p value
# b0    11.9179656 0.079797134  149.3533       0
# b1     2.3128462 0.013354027  173.1947       0
# b2    -1.3470146 0.006653184 -202.4616       0
# SD(e)  4.0154969          NA        NA      NA
# R2     0.3542002          NA        NA      NA

# get R2
R2 <- Rm[-1,1] %*% solve(Rm[-1, -1]) %*% Rm[-1, 1]

pop <- c(Ey, betaZ[1], betaZ[2], sqrt(sigmaE), R2)
estimate <- myMod2[,1]
difference <- pop - estimate

difTab <- cbind(pop, estimate, difference)
#               pop   estimate    difference
# b0    11.9230769 11.9179656  0.0051112956
# b1     2.3076923  2.3128462 -0.0051539037
# b2    -1.3461538 -1.3470146  0.0008607601
# SD(e)  4.0191848  4.0154969  0.0036878149
# R2     0.3538462  0.3542002 -0.0003540325

#-------------------------3.3d------------------------------------------#
# Generating for Any Number of Predictors

# Create two functions that combine the predictor data generation with the outcome data generation for methods 1 and 2. Both functions should implement the general matrix expressions given so that they can incorporate any number of predictors. These functions should map onto the following inputs:  function(n, p_x, p_y)

# function to create the covariance matrix
genSigma <- function(params){
    # create a matrix of zeros with variances on diagonal 
    diagVariances = matrix(diag(params$sigmas), nrow = length(params$sigmas))
    # create a matrix with correlations on off diag, 1 on diagonal
    # start with matrix of correct dims with NAs
    corMatrix = matrix(NA, nrow = length(params$sigmas), ncol = length(params$sigmas))
    # put in the diagonal correlations which is always just 1s
    diag(corMatrix) <- rep(1, length(params$sigmas))
    # add in the correlations bw variables on the lower triangle
    corMatrix[lower.tri(corMatrix)] <- params$rhos
    #and the upper triangle
    corMatrix[upper.tri(corMatrix)] <- t(corMatrix)[upper.tri(t(corMatrix))]
    
    # diagVar * corMatrix * diagVar (eq. 20a -- not labelled directly) = covariance matrix  
    covarMatrix = diagVariances %*% corMatrix %*% diagVariances    
    return(covarMatrix)
}

# create function that takes the x_p and y_p inputs and translates them into a list for the next functions like genSigma
transFun <- function(x_p, y_p){
    list(
        # the y params are first in keeping with the partitioning on pg. 15
        mus = c(y_p$mus, x_p$mus),
        sigmas = c(y_p$sigmas, x_p$sigmas),
        rhos = c(y_p$rhos, x_p$rhos)
    )
}

# Method 1, generate from total variance explained
# Test the Method 1 function using five predictors with n = 100000 that are all correlated 0.15 with variances from 1 to 5.Set all regression coefficients to 1,the R2 =0.5, and μy = 25. Set the seed to 6972.

set.seed(6972)
n <- 100000
betas <- c(1,1,1,1,1)
R2 <- .5

x_p <- list(
    mus = c(0,0,0,0,0),
    sigmas = c(sqrt(1), sqrt(2), sqrt(3), sqrt(4), sqrt(5)),
    rhos = c(.15, .15, .15, .15, .15)
)

# *for method 1, there is no sigma or correlations given for y but including here to make it general to be used in method 2
y_p <- list(
    mus = c(25),
    sigmas = NULL,
    rhos = NULL
)

method1Fun <- function(n, x_p, y_p, betas, R2){
    #get params into workable form
    params <- transFun(x_p, y_p)
    # generate the covariance matrix, Sigma
    covarMatrix <- genSigma(params)
    # generate the Xs
    X <- rmvnorm(n, x_p$mus, covarMatrix)
    # add intercept
    X <- cbind(matrix(1, NROW(X)), X)
    # determine the variance of Y given betas and covarMatrix
    varY <- t(betas) %*% covarMatrix %*% betas # eq 22
    # determine the error variance 
    varE <- varY * (1/R2 - 1) 
    # generate y
    # expectation of Y
    Ey <- params$mus[1] - t(params$mus[-1]) %*% betas
    # generate y
    Y <- Ey[1] + X[,-1] %*% betas + rnorm(n, 0, sqrt(sigmaE))
    
    dat <- cbind(Y, X)
    return(dat)
}

test1 <- method1Fun(n, x_p, y_p, betas, R2)
method1Mod <- myLM(test1[,1], test1[,2:7])
#           Estimate          SE    t value p value
# b0    25.0182701 0.012690524 1971.41356       0
# b1     0.9843795 0.013048351   75.44091       0
# b2     1.0166090 0.009238173  110.04438       0
# b3     0.9901229 0.007552943  131.09100       0
# b4     0.9876760 0.006538371  151.05843       0
# b5     1.0042164 0.005881917  170.72944       0
# SD(e)  4.0129881          NA         NA      NA
# R2     0.5893671          NA         NA      NA

# get varE parameter to compare to estimate
params <- transFun(x_p, y_p)
# generate the covariance matrix, Sigma
covarMatrix <- genSigma(params)
# generate the Xs
X <- rmvnorm(n, x_p$mus, covarMatrix)
# add intercept
X <- cbind(matrix(1, NROW(X)), X)
# determine the variance of Y given betas and covarMatrix
varY <- t(betas) %*% covarMatrix %*% betas # eq 22
# determine the error variance 
varE <- varY * (1/R2 - 1) 
# generate y
# expectation of Y
Ey <- params$mus[1] - t(params$mus[-1]) %*% betas

pop <- c(Ey, betas, sqrt(varE), R2)
estimate <- method1Mod[,1]
difference <- pop - estimate

difTab <- cbind(pop, estimate, difference)
#           pop   estimate   difference
# b0    25.000000 25.0182701 -0.018270136
# b1     1.000000  0.9843795  0.015620494
# b2     1.000000  1.0166090 -0.016608976
# b3     1.000000  0.9901229  0.009877059
# b4     1.000000  0.9876760  0.012324021
# b5     1.000000  1.0042164 -0.004216436
# SD(e)  4.825922  4.0129881  0.812933981
# R2     0.500000  0.5893671 -0.089367060

# Method 2, generate from marginal correlations
# Test the Method 1 function using five predictors with n = 100000 that are all cor- related 0.15 with variances from 1 to 5. The correlations with Y ought to be −.15, −.50, .15, 0.30, and 0.20. The mean and variance of Y are μy = 10, σy = 4. Set the seed to 1237
set.seed(1237)
n <- 100000

x_p <- list(
    mus = c(0,0,0,0,0),
    sigmas = c(sqrt(1), sqrt(2), sqrt(3), sqrt(4), sqrt(5)),
    rhos = rep(.15, 10)
)

y_p <- list(
    mus = c(10),
    sigmas = c(4),
    rhos = c(-.15, -.5, .15, .30, .20)
)

method2Fun <- function(n, x_p, y_p, betas){
    #get params into workable form
    params <- transFun(x_p, y_p)
    #gereate covrmtrix for x
    covarMatrix <- genSigma(x_p)
    # generate X
    X <- rmvnorm(n, x_p$mus, covarMatrix)
    #intercept
    X <- cbind(matrix(1, NROW(X)), X)
    # generate the covariance matrix, Sigma
    covarMatrix <- genSigma(params)
    # solve for the betas
    betas <- solve(covarMatrix[-1, -1]) %*% covarMatrix[1, -1]
    # determine the error variance 
    varE <- covarMatrix[1,1] - t(betas) %*%  covarMatrix[1, -1]
    # generate y
    # expectation of Y
    Ey <- params$mus[1] - t(params$mus[-1]) %*% betas
    # generate y
    Y <- Ey[1] + X[,-1] %*% betas + rnorm(n, 0, sqrt(sigmaE))
    
    dat <- cbind(Y, X)
    return(dat)
}

test2 <- method2Fun(100, x_p, y_p, betas = betas)
method2Mod <- myLM(test2[,1], test2[,2:7])
method2Mod

#       Estimate         SE   t value      p value
# b0    10.2086353 0.12571876 81.202162 7.344336e-89
# b1    -0.4859086 0.12179292 -3.989630 1.309803e-04
# b2    -0.7473232 0.09543021 -7.831097 7.208772e-12
# b3     0.1828802 0.07215174  2.534661 1.290787e-02
# b4     0.3673510 0.06622996  5.546599 2.667136e-07
# b5     0.1796883 0.05852594  3.070234 2.796316e-03
# SD(e)  1.2342154         NA        NA           NA
# R2     0.5314705         NA        NA           NA

# calculate R2
# function to create the covariance matrix
params <- transFun(x_p, y_p)
# empty matrix
corMatrix = matrix(NA, nrow = length(params$sigmas), ncol = length(params$sigmas))
    # put in the diagonal correlations which is always just 1s
diag(corMatrix) <- rep(1, length(params$sigmas))
    # add in the correlations bw variables on the lower triangle
corMatrix[lower.tri(corMatrix)] <- params$rhos
    #and the upper triangle
corMatrix[upper.tri(corMatrix)] <- t(corMatrix)[upper.tri(t(corMatrix))]
# get R2 to compare to estimate eq 27
R2 <- corMatrix[-1, 1] %*% solve(corMatrix[-1, -1]) %*% corMatrix[1, -1]

Ey <- params$mus[1] - t(params$mus[-1]) %*% betas
# generate y
Y <- Ey[1] + X[,-1] %*% betas + rnorm(n, 0, sqrt(sigmaE))

#gereate covrmtrix for x
covarMatrix <- genSigma(params)
# solve for the betas
betas <- solve(covarMatrix[-1, -1]) %*% covarMatrix[1, -1]
# determine the error variance 
varE <- covarMatrix[1,1] - t(betas) %*%  covarMatrix[1, -1]

pop <- c(Ey, betas, sqrt(varE), R2)
estimate <- method2Mod[,1]
difference <- pop - estimate

difTab <- cbind(pop, estimate, difference)

#               pop   estimate  difference
# b0    10.0000000 10.4172707 -0.41727069
# b1    -0.7058824 -0.9718173  0.26593491
# b2    -1.6637807 -1.4946464 -0.16913421
# b3     0.4075414  0.3657604  0.04178099
# b4     0.7058824  0.7347021 -0.02881971
# b5     0.4209069  0.3593766  0.06153028
# SD(e)  2.8284271  2.4684309  0.35999626
# R2     0.5000000  0.5314705 -0.03147045
