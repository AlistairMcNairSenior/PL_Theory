
# Functions to calculate the expected strength of PL based on mean P and NP, the SD of P and NP, and their correlation, rho_PNP. Assumes a bivariate normal distribution

get_L<-function(mu_P, mu_NP, s_P, s_NP, rho_PNP, return_params=F){
	
	# Mean and SD for total intake
	mu_T<-mu_P + mu_NP
	s_T<-sqrt(s_P^2 + s_NP^2 + 2*rho_PNP*s_P*s_NP)
	
	# And the correlation between total intake and protein intake
	rho_PT<-(s_P^2 + rho_PNP*s_P*s_NP) / (s_P * s_T)
	
	# The mean of the proportion of energy from protein (pP) found by integration of the pdf for pP, which is as a function below
	f<-function(w){
		w*as.numeric(f_W(w=w, mu_X=mu_P, mu_Y=mu_T, s_X=s_P, s_Y=s_T, rho=rho_PT))
	}
	mu_pP<-integrate(f, -Inf, Inf)$value
	
	# The sd for the log of pP is found by probability integral transformation of the pdf for pP.
	# First find the expected log of pP
	f<-function(z){
		z * as.numeric(f_W(w=exp(z), mu_X=mu_P, mu_Y=mu_T, s_X=s_P, s_Y=s_T, rho=rho_PT)) * exp(z)
	}
	mu_lnpP<-integrate(f, -10, 10)$value
	
	# Then the sd is
	f<-function(z){
		(z - mu_lnpP)^2 * as.numeric(f_W(w=exp(z), mu_X=mu_P, mu_Y=mu_T, s_X=s_P, s_Y=s_T, rho=rho_PT)) * exp(z)
	}
	s_lnpP<-sqrt(integrate(f, -10, 10)$value)
	
	# Find L and return it
	L<-log(1 + (mu_P - mu_pP*mu_T)/(mu_pP*mu_T)) / s_lnpP^2
	if(return_params == F){
		return(L)
	}else{
		output<-data.frame(L, mu_pP, mu_lnpP, s_lnpP)
		return(output)	
	}
}

# Function for the pdf of the ratio of two non-central correlated, random normal distributions.
# Based on eq. 13 of Pham-Gia et al. 2006. Communications in Statistics - Theory and Methods, 35: 1569-1591
# Uses the fAsianOptions package
# w is the value of the ratio, mu_X and mu_Y are the means of the distributions, s_X and s_Y their SDs, and rho the correlation

f_W<-function(w, mu_X, mu_Y, s_X, s_Y, rho){
	
	require(fAsianOptions) # For function kummerM - Kummer's classical confluent hypergeometric function of first kind - F1 below
	
	K<-1 / (2*pi*s_X*s_Y*sqrt(1-rho^2)) * exp(-((s_Y^2*mu_X^2 - 2*rho*s_X*s_Y*mu_X*mu_Y + mu_Y^2*s_X^2) / (2*(1-rho^2)*s_X^2*s_Y^2)))
	
	theta<-(-s_Y^2*mu_X*w + rho*s_X*s_Y*(mu_Y*w+mu_X)-mu_Y*s_X^2)^2 / (2*s_X^2*s_Y^2*(1-rho^2) * (s_Y^2*w^2 - 2*rho*s_X*s_Y*w+s_X^2))
	
	F1<-kummerM(theta, 1, 1/2)
	
	p<-K * ((2 * (1-rho^2)*s_X^2*s_Y^2) / (s_Y^2*w^2 - 2*rho*s_X*s_Y*w + s_X^2)) * F1
	
	return(p)
	
}

# Function to get the means and variances for P, NP and Total Energy assuming PL, mean and variance in prop energy from p, strength L, a given P parameter and an assumed residual variance

# mu_W = mean proportion energy from protein
# s_W = sd proportion energy from protein
# L = strength of leverage - or the beta in the log-log regression 
# P = P parameter in Pp^L - or the exp(alpha) in the log-log regression
# e = the sd of the residual variance 

get_intakes<-function(mu_W, s_W, L, P, e){
	
	# The PDF for W is
	alpha<-((1 - mu_W) / s_W^2 - (1 / mu_W)) * mu_W^2
	beta<-alpha * (1 / mu_W - 1)
	f_W<-function(w, alpha, beta){
		dbeta(w, alpha, beta)
	}
	
	# The pdf for X = log(W) is
	f_X<-function(x, alpha, beta){
		f_W(exp(x), alpha, beta) * exp(x)
	}
	
	# And the mean and sd of X are
	# And the mean will be
	mu_X<-integrate(function(x, alpha, beta) x * f_X(x, alpha, beta), lower=-Inf, upper=0, alpha=alpha, beta=beta)$value
	s_X<-sqrt(integrate(function(x, mu_X, alpha, beta) (x - mu_X)^2 * f_X(x, alpha, beta), lower=-Inf, upper=0, mu_X=mu_X, alpha=alpha, beta=beta)$value)

	# Which means that yi = a + L*xi + ei  = , as is in PL we have; Y = G + R, 
	# And G will have as transformation the ditribution: f_X(g/L, alpha, beta) * abs(1/L)
	f_G<-function(g, L, alpha, beta){
		f_X(g/L, alpha, beta) * abs(1/L)
	}
	# R will have norm(log(P), e)
	f_R<-function(r) dnorm(r, log(P), e)
	f_Y<-function(y, L, alpha, beta) integrate(function(g,y) f_R(y-g)*f_G(g, L=L, alpha=alpha, beta=beta), 0, Inf, y)$value
	f_Y<-Vectorize(f_Y)

	# With the mean and sd for Y as
	mu_Y<-integrate(function(y) y * f_Y(y, L=L, alpha=alpha, beta=beta), -Inf, Inf)$value
	s_Y<-sqrt(integrate(function(y) (y - mu_Y)^2 * f_Y(y, L=L, alpha=alpha, beta=beta), -Inf, Inf)$value)

	# Now backtransform Y to get Z
	# Which will have density: f_Y(log(z)) * (1/z)
	f_Z<-function(z, L, alpha, beta){
		f_Y(log(z), L, alpha, beta) * (1/z)
	}

	# And mean and SD of
	mu_Z<-integrate(function(z) z * f_Z(z, L, alpha, beta), 0, 10000000)$value
	s_Z<-sqrt(integrate(function(z) (z - mu_Z)^2 * f_Z(z, L, alpha, beta), 0, 10000000)$value)

	# We are 1/3rd of the way there! On to logU - log absolute protein intake
	
	# Assuming X and Y are approx. normal this means that mu_lnU and sd_lnU can be approxiamted as
	mu_lnU<-mu_X + mu_Y
	s_lnU<-sqrt(s_X^2 + s_Y^2 + 2 * L * s_X^2)
	
	# Which backtransform to
	mu_U<-exp(mu_lnU + 0.5*s_lnU^2)
	s_U<-sqrt((exp(s_lnU^2) - 1) * mu_U^2)
	
	# 2/3rds of the way there - that was quick! Now on to V - intake of non-protein
	
	# The mean will be 
	mu_V<-mu_Z - mu_U
	
	# The sd for V is more of a pain
	# Start with the covariance between proportion non-protein energy and total energy, and transform to log scale - saying S = log(1 - W)
	covSY<-log(1 + (-(mu_U - mu_Z*mu_W) / (mu_Z * (1 - mu_W))))
	
	# The pdf for S is:
	f_S<-function(s, alpha, beta){
		dbeta(1-exp(s), shape1=alpha, shape2=beta) * exp(s)
	}
	
	# The mean and sd of S is:
	mu_S<-integrate(function(s, alpha, beta) s * f_S(s, alpha, beta), lower=-Inf, upper=0, alpha=alpha, beta=beta)$value
	s_S<-sqrt(integrate(function(s, mu_S, alpha, beta) (s - mu_S)^2 * f_S(s, alpha, beta), lower=-Inf, upper=0, mu_S=mu_S, alpha=alpha, beta=beta)$value)

	# The approximate the sd in log(V) and backtransform
	s_lnV<-sqrt(s_Y^2 + s_S^2 + 2*covSY)
	s_V<-sqrt((exp(s_lnV^2) - 1) * mu_V^2)
	
	# Phew - parcel up and return
	mus<-c(mu_U, mu_V, mu_Z)
	sds<-c(s_U, s_V, s_Z)
	
	return(list(mus, sds))
	
}


