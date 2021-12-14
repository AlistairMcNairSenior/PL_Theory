
# Clean up R
rm(list=ls())

# Load packages and the header
source("Header.R")
library(MASS)
library(anySim)
library(ggplot2)
library(gridExtra)
library(Cairo)
library(lemon)

####################################################
#### Parameterise P as a % of 8700 kJ/day ##########
####################################################

# A reasonable range of values to test over

# Mean protein and non protein intakes
mu_P<-seq(800, 4800, 100)

# I parameterise from the dispersion index
D_P<-seq(50, 200, 10)
D_NP<-mean(D_P)

# Correlation between P and NP
rho_PNP<-c(0, 0.25, 0.5)

# Get the combinations
conditions<-expand.grid(mu_P, D_P, D_NP, rho_PNP)
names(conditions)<-c("mu_P", "D_P", "D_NP", "rho_PNP")

# Add in mu_NP as 8700 - mu_P
conditions$mu_NP<-8700 - mu_P

# Calculate the %P - useful for later
conditions$Percent_P<-conditions$mu_P / (conditions$mu_P + conditions$mu_NP) * 100

# Add in the SDs by calculating from D
conditions$s_P<-sqrt(conditions$D_P * conditions$mu_P)
conditions$s_NP<-sqrt(conditions$D_NP * conditions$mu_NP)

# Go through each condition and get the expected and simualted value of L
conditions$L<-NA
conditions$L_sim<-NA
conditions$L_sim2<-NA
for(i in 1:nrow(conditions)){
	
	# Pull out the mus and sds
	mus<-c(conditions$mu_P[i], conditions$mu_NP[i])
	sds<-c(conditions$s_P[i], conditions$s_NP[i])
	
	# Find the expected value from the analytical solution
	conditions$L[i]<-get_L(mu_P=mus[1], mu_NP=mus[2], s_P=sds[1], s_NP=sds[2], rho_PNP=conditions$rho_PNP[i])
	
	# Now get a simulated value of L under a bivariate normal distribution
	# Create the vcov matrix
	covar<-conditions$rho_PNP[i]*sds[1]*sds[2]
	vcov<-array(covar, c(2,2))
	diag(vcov)<-sds^2
	
	# Simulate
	data<-as.data.frame(mvrnorm(20000, mu=mus, Sigma=vcov, empirical=TRUE))
	names(data)<-c("P", "NP")
	data$Tot<-data$P + data$NP
	data$pP<-data$P / data$Tot
	
	# Fit the log-log model to get the simulated L
	conditions$L_sim[i]<-lm(log(Tot) ~ log(pP), data=data)$coef[2]
	
	# Now repeat with a log-normal distribution 
	
	# Convert the raw means and sds to the approapriate log values
	ln_mus<-log(mus^2 / sqrt(mus^2 + sds^2))
	ln_sds<-sqrt(log(1 + (sds^2 / mus^2)))
	
	# If we have correlated lnNorm use the EstCorrRVs functions, but for non-correlated, just use mvrnorm and take the exponent
	if(conditions$rho_PNP[i] > 0){
		
		# Package likes this listed up
		DistrParams<-list()
		DistrParams[[1]]<-list(meanlog=ln_mus[1], sdlog=ln_sds[1])
		DistrParams[[2]]<-list(meanlog=ln_mus[2], sdlog=ln_sds[2])
	
		# The package takes a correlation matrix
		CorrelMat<-matrix(c(1,conditions$rho_PNP[i],
					        conditions$rho_PNP[i],1), ncol=2, nrow=2)

		# Get the converted parameters
		paramsRVs<-EstCorrRVs(R=CorrelMat, dist=c("qlnorm", "qlnorm"), params=DistrParams, NatafIntMethod="GH")
		
		# Simulate			    
		data<-as.data.frame(SimCorrRVs(n=20000, paramsRVs=paramsRVs))
		
	}else{
		
		# VCOV with 0 cov
		vcov<-array(0, c(2,2))
		diag(vcov)<-ln_sds^2
		
		# Simulate and take the exponent
		data<-as.data.frame(mvrnorm(20000, mu=ln_mus, Sigma=vcov, empirical=TRUE))
		data<-exp(data)
	}
	
	# Calculate total and pP
	names(data)<-c("P", "NP")
	data$Tot<-data$P + data$NP
	data$pP<-data$P / data$Tot
	
	# Fit the log-log model to get the simulated L
	conditions$L_sim2[i]<-lm(log(Tot) ~ log(pP), data=data)$coef[2]
	
}

# Check the range of L
range(conditions$L)

# Differences between simulated and theoretical values
conditions$diff1<-conditions$L_sim - conditions$L
conditions$diff2<-conditions$L_sim2 - conditions$L

conditions$IDR<-conditions$D_P / conditions$D_NP

# Plot out the different rho's
plots_list<-list()
plots_list2<-list()
plots_list3<-list()
plots_list4<-list()
plots_list5<-list()
for(i in 1:length(rho_PNP)){
	
	plots_list[[i]]<-ggplot(data=conditions[which(conditions$rho_PNP == rho_PNP[i]),], aes(x=Percent_P, y=IDR, fill=L)) + 
	geom_tile() + theme_bw() + 
	labs(x=expression(frac(italic("\U03BC")[italic(U)], (italic("\U03BC")[italic(U)] + italic("\U03BC")[italic(V)])) %*% 100), y=expression(frac(ID[italic(U)], ID[italic(V)])), fill=expression(italic(L)), subtitle=substitute(paste(italic("\U03C1")[italic(UV)], " = ", x), list(x=rho_PNP[i]))) + 
	scale_fill_gradient2(limits=c(-1, 1)) + 
	geom_hline(yintercept=1, size=1) + geom_vline(xintercept=50, size=1) +
	theme(text=element_text(size=13)) +
	theme(legend.title=element_text(size=20), legend.text=element_text(size=20), legend.key.size=unit(1, "cm"))
	
	plots_list2[[i]]<-ggplot(data=conditions[which(conditions$rho_PNP == rho_PNP[i]),], aes(x=L, y=L_sim)) + geom_point(col="firebrick", size=0.5) + 
	geom_abline(intercept=0, slope=1) + 
	geom_vline(xintercept=0) + geom_hline(yintercept=0) +
	theme_bw() + 
	labs(x="Theoretical L (Eqn. 23)", y=expression("Simulated"~~italic(L)~~"(Normal Distribution)"), subtitle=substitute(paste(italic("\U03C1")[italic(UV)], " = ", x), list(x=rho_PNP[i]))) + 
	theme(text=element_text(size=15))
	
	plots_list3[[i]]<-ggplot(data=conditions[which(conditions$rho_PNP == rho_PNP[i]),], aes(x=L, y=L_sim2)) + geom_point(col="firebrick", size=0.5) + 
	geom_abline(intercept=0, slope=1) + 
	geom_vline(xintercept=0) + geom_hline(yintercept=0) +
	theme_bw() + 
	labs(x="Theoretical L (Eqn. 23)", y=expression("Simulated"~~italic(L)~~"(Log-Normal Distribution)"), subtitle=substitute(paste(italic("\U03C1")[italic(UV)], " = ", x), list(x=rho_PNP[i]))) + 
	theme(text=element_text(size=15))


}

# Pull out a legend
legend<-g_legend(plots_list[[1]])

cairo_pdf("Fig_2.pdf", height=7, width=7)

grid.arrange(plots_list[[1]]+theme(legend.position="hidden"), plots_list[[2]]+theme(legend.position="hidden"), plots_list[[3]]+theme(legend.position="hidden"), legend, layout_matrix=(rbind(c(1,2),c(3,4))))

dev.off()

cairo_pdf("Fig_3.pdf", height=9, width=13.5)

grid.arrange(plots_list2[[1]], plots_list2[[2]], plots_list2[[3]], plots_list3[[1]], plots_list3[[2]], plots_list3[[3]], layout_matrix=(rbind(c(1,2,3),c(4,5,6))))

dev.off()
