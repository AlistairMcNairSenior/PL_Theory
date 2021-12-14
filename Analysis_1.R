
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
mu_P<-seq(1000, 4800, 100)

# Add in the SDs - assume constant
s_P<-500
s_NP<-s_P

# Correlation between P and NP
rho_PNP<-c(0, 0.25, 0.5)

# Get the combinations
conditions<-expand.grid(mu_P, s_P, s_NP, rho_PNP)
names(conditions)<-c("mu_P", "s_P", "s_NP", "rho_PNP")

# Add in mu_NP as 8700 - mu_P
conditions$mu_NP<-8700 - mu_P

# Calculate the %P - useful for later
conditions$Percent_P<-conditions$mu_P / (conditions$mu_P + conditions$mu_NP) * 100

# Go through each condition and get the expected and simualted value of L
conditions$L<-NA
for(i in 1:nrow(conditions)){
	
	# Pull out the mus and sds
	mus<-c(conditions$mu_P[i], conditions$mu_NP[i])
	sds<-c(conditions$s_P[i], conditions$s_NP[i])
	
	# Find the expected value from the analytical solution
	conditions$L[i]<-get_L(mu_P=mus[1], mu_NP=mus[2], s_P=sds[1], s_NP=sds[2], rho_PNP=conditions$rho_PNP[i])

}

# Check the range of L
range(conditions$L)

# Plot out the different rho's
plots_list<-list()
	
plots_list[[1]]<-ggplot(data=conditions, aes(x=Percent_P, y=L, color=as.factor(rho_PNP))) + 
	geom_hline(yintercept=0, size=0.5) + geom_vline(xintercept=50, size=0.5) +
	geom_path(size=1.5) + theme_bw() + 
	labs(x=expression(frac(italic("\U03BC")[italic(U)], (italic("\U03BC")[italic(U)] + italic("\U03BC")[italic(V)])) %*% 100), y=expression(italic(L)), color=expression(italic("\U03C1")[italic(UV)])) + 
	theme(text=element_text(size=15)) +
	theme(legend.title=element_text(size=20), legend.text=element_text(size=20), legend.key.size=unit(0.25, "cm"), legend.position=c(0.1, 0.15))



# A reasonable range of values to test over

# now paramterise from the index of dispersion D
D_P<-100
D_NP<-D_P

# Get the combinations
conditions<-expand.grid(mu_P, D_P, D_NP, rho_PNP)
names(conditions)<-c("mu_P", "D_P", "D_NP", "rho_PNP")

# Add in mu_NP as 8700 - mu_P
conditions$mu_NP<-8700 - mu_P

# Add in the SDs from the CVs
conditions$s_P<-sqrt(conditions$D_P * conditions$mu_P)
conditions$s_NP<-sqrt(conditions$D_NP * conditions$mu_NP)

# Calculate the %P - useful for later
conditions$Percent_P<-conditions$mu_P / (conditions$mu_P + conditions$mu_NP) * 100

# Go through each condition and get the expected and simualted value of L
conditions$L<-NA
for(i in 1:nrow(conditions)){
	
	# Pull out the mus and sds
	mus<-c(conditions$mu_P[i], conditions$mu_NP[i])
	sds<-c(conditions$s_P[i], conditions$s_NP[i])
	
	# Find the expected value from the analytical solution
	conditions$L[i]<-get_L(mu_P=mus[1], mu_NP=mus[2], s_P=sds[1], s_NP=sds[2], rho_PNP=conditions$rho_PNP[i])

}

# Check the range of L
range(conditions$L)
	
plots_list[[2]]<-ggplot(data=conditions, aes(x=Percent_P, y=L, color=as.factor(rho_PNP))) + 
	geom_hline(yintercept=0, size=0.5) + geom_vline(xintercept=50, size=0.5) +
	geom_path(size=1.5) + theme_bw() + 
	labs(x=expression(frac(italic("\U03BC")[italic(U)], (italic("\U03BC")[italic(U)] + italic("\U03BC")[italic(V)])) %*% 100), y=expression(italic(L)), color=expression(italic("\U03C1")[italic(UV)])) + 
	theme(text=element_text(size=15)) +
	theme(legend.title=element_text(size=20), legend.text=element_text(size=20), legend.key.size=unit(0.25, "cm"), legend.position="none")

cairo_pdf("Fig 1.pdf", height=7, width=14)

grid.arrange(plots_list[[1]]+labs(title="A."), plots_list[[2]]+labs(title="B"), layout_matrix=array(c(1,2), c(1,2)))

dev.off()

