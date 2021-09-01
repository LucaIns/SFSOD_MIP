##########################################
# Data Generation parameters
########################################## 
  
# ----------------------------------------------------
# Parameters 
# ----------------------------------------------------
# num. of obs.
n = N 	
# fraction of contaminated units
frac_cont = R
# including intercept (T/F)
intercept = 1
# total predictors (including intercept)
p = d
# relevant 'non-zero' predictors (possibly including intercept)
s =  true_s
# irrelevant 'zero' predictors
irr_p = p-s
# True coefficient vector
# draw the first s 'non-zero' coefficients from an uniform(a, b)
# a = 2; b = 3
# beta_true = runif(s, min = a, max = b)
beta_fixed = 2
beta_true = rep(beta_fixed, s)
# attach zeros parameter
beta_true = c(beta_true, rep(0, irr_p))
# Uncontaminated data  -----------------------------------------------
# centrality for X
mu_X_u = rep(0, p - intercept)
# dispersion for X
sigma_X_u = 1
rho_X_u = diag(1, p - intercept) * sigma_X_u^2
# errors' dispersions
# sigma_e_u = 1
# Contamination mechanisms  --------------------------------------------
# Signal to Noise Ratio
if (lowSNR == 0){
	SNR = 5
} else if (lowSNR == 1){
	SNR = 3
}
# contamination mechanisms considered
mech_nam = "relevant"  #"row-wise"
# mean shifting for leverages
mu_X_out = 10
# # dispersion for leverages
#sigma_X_out = 2
# mean shifting for vertical outliers
mu_y_out = -10
# dispersion for vertical outliers
#sigma_y_out = 2
# flag for contaminating also the response
y_cont_too = 1
# ----------------------------------------------------
# Data generation 
# ----------------------------------------------------
if (multicoll == 0){
	X_u = SimDesign::rmvnorm(n, mean = mu_X_u, sigma = rho_X_u)
} else if (multicoll == 1){
	suppressWarnings(suppressMessages(library("simFrame")))
	# correlated normal predictors
	mean <- rep(0, p-intercept)
	sigma <- matrix(0, nrow = p-intercept, ncol = p-intercept)
	for (i in 1:(p-intercept)) {
		for (j in 1:(p-intercept)) {
		  sigma[i,j] <- 0.3^abs(i-j)
		}
	}
	dc <- simFrame::DataControl(size = n, distribution = SimDesign::rmvnorm,
	                dots = list(mean = mean, sigma = sigma))
	X_u = as.matrix(simFrame::generate(dc))
}
if (intercept == 1){
	X_u = cbind(rep(1, n), X_u)
}
sigma_e_u = sqrt(var(X_u %*% beta_true) / SNR)
e_u = rnorm(n = n, mean = 0, sd = sigma_e_u)
y_u = X_u %*% beta_true + e_u

# check = lm(y_u~X_u-1)
# summary(check)

# generate test set
if (multicoll == 0){
	Xtest = SimDesign::rmvnorm(n, mean = mu_X_u, sigma = rho_X_u)
} else if (multicoll == 1){
	# correlated normal predictors
	mean <- rep(0, p-intercept)
	sigma <- matrix(0, nrow = p-intercept, ncol = p-intercept)
	for (i in 1:(p-intercept)) {
	for (j in 1:(p-intercept)) {
	  sigma[i,j] <- 0.5^abs(i-j)
	}
	}
	dc <- simFrame::DataControl(size = n, distribution = rmvnorm,
	                dots = list(mean = mean, sigma = sigma))
	Xtest = as.matrix(simFrame::generate(dc))
}
if (intercept == 1){
	Xtest = cbind(rep(1, n), Xtest)
}
e_test = rnorm(n = n, mean = 0, sd = sigma_e_u)
Ytest = Xtest %*% beta_true + e_test
# ----------------------------------------------------
# Data contamination 
# ----------------------------------------------------
n_out_i = round(n * frac_cont)
if (n_out_i>0) {
	if (mech_nam=="relevant") {
	  k = s 
	} else if (mech_nam=="row-wise"){
	  k = p
	}
	X_c = X_u
	X_c[1:n_out_i, (1+intercept):k] = X_u[1:n_out_i, (1+intercept):k] + mu_X_out
	y_c = y_u
	if (y_cont_too == 1){
	  y_c[1:n_out_i] = y_u[1:n_out_i] + mu_y_out
	}
} else{
	X_c = X_u
	y_c = y_u
}
# (Sanity check)
# plot(X_u[,1+intercept], y_u)
# plot(X_c[,1+intercept], y_c)
#X = cbind(X_c, diag(n))
X = X_c
Y = y_c