
# LOAD FUNCTION

#################################################################################
# Simulation data of ALfons 2013, p.235: for p > n
#################################################################################

simDataAlfonsHigh <- function(n = 100, d = 1000, cont = "No", frac=0.1, intercept = T){
  # n: sample size
  # d: number of predictors
  # cont: tipe of contamination
    # "No": clean data
    # "Vert": only for the response
    # "Lev": create "Vert" and also bad leverage points
  # frac: fraction of contamination
  
  # load/install required packages
  suppressWarnings(suppressMessages(library("simFrame")))
  suppressWarnings(suppressMessages(library("mvtnorm")))
  
  # correlated normal predictors
  mean <- rep(0, d-intercept)
  sigma <- matrix(0, nrow = d-intercept, ncol = d-intercept)
  for (i in 1:d-intercept) {
    for (j in 1:d-intercept) {
      sigma[i,j] <- 0.5^abs(i-j)
    }
  }
  # library("corrplot")
  # corrplot(sigma, method = "circle")
  dc <- DataControl(size = n, distribution = rmvnorm,
                    dots = list(mean = mean, sigma = sigma))
  X = as.matrix(generate(dc))
  if (intercept == T){
    X = cbind(rep(1, N), X)
  }
  # coefficients  
  beta = rep(0, d)
  contSet = c(1,7,2,4,11)
  beta[contSet] <- c(1.5, 1.5, 0.5, 1, 1)
  # normal errors (they did not control for SNR=31 here)
  sigma_e_u = 0.5
  err = rnorm(n, mean=0, sd=sigma_e_u)
  # response
  y <-  X %*% beta + err
  #SNR = var(X %*% beta)/var(err) 
  #print(SNR)
  
# generate test set
Xtest = as.matrix(generate(dc))
if (intercept == 1){
  Xtest = cbind(rep(1, n), Xtest)
}
e_test = rnorm(n = n, mean = 0, sd = sigma_e_u)
Ytest = Xtest %*% beta + e_test


  # response (asymmetric) contamination
  if (cont == "Vert" || cont == "Lev"){
    m = round(n*frac)
    subs = 1:m # sample(1:n, size = m)
    y[1:m] = y[1:m] + 20
  } 
  
  # X (asymmetric) contamination -- ONLY on k relevant features, perhaps different than Alfons (2013)
  if (cont == "Lev"){
    dcCont = DataControl(size = m, distribution = rmvnorm,
                         dots = list(mean = rep(20, length(contSet)-intercept), sigma = diag(rep(1, length(contSet)-intercept))))
    dtCont = generate(dcCont)
    dtCont = as.matrix(dtCont)
    if (intercept == 1){
		dtCont = cbind(rep(1, dim(dtCont)[1]), dtCont)
	}
    X[subs, contSet] = dtCont
  }
  
  # put the k relevant features in the first k columns of X
  X = cbind(X[,contSet], X[,-contSet])
  Xtest = cbind(Xtest[,contSet], Xtest[,-contSet])
  beta = c(beta[contSet], beta[-contSet])
  
  result <- list(y = y, X = X, beta = beta, Ytest=Ytest, Xtest=Xtest)
  return(result)
  
}

###########################################################################

intercept = 0
dat = simDataAlfonsHigh(n = N, d = d, cont = "Lev", frac=R, intercept = intercept)
Y = dat$y
X = dat$X
beta_true = dat$beta
Ytest = dat$Ytest
Xtest = dat$Xtest