
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

# The mean value of W - proportion energy from protein, and its SD
mu_W<-c(0.15, 0.2, 0.3)
s_W<-seq(0.025, 0.07, 0.005)

# The value of L - strength of leverage
L<-seq(-0.55, -0.1, 0.025)

# And some residual variance
e<-sqrt(log(8700) * 0.02)

# Get the combinations
conditions<-expand.grid(mu_W, s_W, L, e)
names(conditions)<-c("mu_W", "s_W", "L", "e")

# The value of P in Pp^L - if mean intake at mu_Z wants to be 8700
conditions$P<-8700 / conditions$mu_W^conditions$L

# Go through each condition and get the expected and simualted value of L
conditions$mu_U<-NA
conditions$mu_V<-NA
conditions$mu_Z<-NA
conditions$s_U<-NA
conditions$s_V<-NA
conditions$s_Z<-NA
conditions$Sim_mu_U<-NA
conditions$Sim_mu_V<-NA
conditions$Sim_mu_Z<-NA
conditions$Sim_s_U<-NA
conditions$Sim_s_V<-NA
conditions$Sim_s_Z<-NA
conditions$corUV

# Run through the conditions
for(i in 1:nrow(conditions)){
	
	# Find the expected value from the analytical solution
	intakes<-get_intakes(mu_W=conditions$mu_W[i], s_W=conditions$s_W[i], L=conditions$L[i], P=conditions$P[i], e=conditions$e[i])
	conditions$mu_U[i]<-intakes[[1]][1]
	conditions$mu_V[i]<-intakes[[1]][2]
	conditions$mu_Z[i]<-intakes[[1]][3]
	conditions$s_U[i]<-intakes[[2]][1]
	conditions$s_V[i]<-intakes[[2]][2]
	conditions$s_Z[i]<-intakes[[2]][3]

	# Now get simulated results
	
	# Simulate W
	alpha<-((1 - conditions$mu_W[i]) / conditions$s_W[i]^2 - (1 / conditions$mu_W[i])) * conditions$mu_W[i]^2
	beta<-alpha * (1 / conditions$mu_W[i] - 1)
	W<-rbeta(100000, alpha, beta)
	
	# Simualte the residual 
	R<-rnorm(100000, 0, conditions$e[i])	
	
	# Calculate Y - log energy intake based on PL
	Y<-log(conditions$P[i]) + conditions$L[i]*log(W) + R
	
	# Backtransform
	Z<-exp(Y)
	
	# Calculate U and V - absolute intakes of P and NP
	U<-Z * W
	V<-Z * (1-W)
	
	# Now get the mean and SD of those
	conditions$Sim_mu_U[i]<-mean(U)
	conditions$Sim_mu_V[i]<-mean(V)
	conditions$Sim_mu_Z[i]<-mean(Z)
	conditions$Sim_s_U[i]<-sd(U)
	conditions$Sim_s_V[i]<-sd(V)
	conditions$Sim_s_Z[i]<-sd(Z)
	conditions$corUV[i]<-cor(U,V)
	
}

# Calculate the ratio of ID for P and NP - simulated and calculated
conditions$IDR<-(conditions$s_U^2 / conditions$mu_U) / (conditions$s_V^2 / conditions$mu_V)
conditions$Sim_IDR<-(conditions$Sim_s_U^2 / conditions$Sim_mu_U) / (conditions$Sim_s_V^2 / conditions$Sim_mu_V)

conditions$CVR<-(conditions$s_U / conditions$mu_U) / (conditions$s_V / conditions$mu_V)

# plot(conditions$L, conditions$IDR)
# plot(conditions$IDR, conditions$Sim_IDR)
# abline(a=0, b=1)

# Plot out the different rho's
plots_list<-list()
plots_list2<-list()
plots_list3<-list()
for(i in 1:length(mu_W)){
	
	plots_list[[i]]<-ggplot(data=conditions[which(conditions$mu_W == mu_W[i]),], aes(x=L, y=s_W, fill=IDR)) + 
	geom_tile() + theme_bw() + 
	labs(x=expression(italic(L)), y=expression(italic("\U03C3")[italic(W)]), fill=expression(frac(ID[italic(U)], ID[italic(V)])), subtitle=substitute(paste(italic("\U03BC")[italic(W)], " = ", x), list(x=mu_W[i]))) + 
	scale_fill_gradient2(limits=c(0.05,0.75), midpoint=0.4, low="blue", mid="green", high="red") + 
	theme(text=element_text(size=13)) +
	theme(legend.title=element_text(size=20), legend.text=element_text(size=20), legend.key.size=unit(1, "cm"))
	
	plots_list2[[i]]<-ggplot(data=conditions[which(conditions$mu_W == mu_W[i]),], aes(x=IDR, y=Sim_IDR)) + geom_point(col="firebrick", size=0.5) + 
	geom_abline(intercept=0, slope=1) + 
	theme_bw() + 
	labs(x="Model 2 IDR", y=expression("Simulated IDR"), subtitle=substitute(paste(italic("\U03BC")[italic(W)], " = ", x), list(x=mu_W[i]))) + 
	theme(text=element_text(size=15))

	plots_list3[[i]]<-ggplot(data=conditions[which(conditions$mu_W == mu_W[i]),], aes(x=L, y=s_W, fill=Sim_IDR)) + 
	geom_tile() + theme_bw() + 
	labs(x=expression(italic(L)), y=expression(italic("\U03C3")[italic(W)]), fill=expression(frac(ID[italic(U)], ID[italic(V)])), subtitle=substitute(paste(italic("\U03BC")[italic(W)], " = ", x), list(x=mu_W[i]))) + 
	scale_fill_gradient2(limits=c(0.05,0.75), midpoint=0.4, low="blue", mid="green", high="red") + 
	theme(text=element_text(size=13)) +
	theme(legend.title=element_text(size=20), legend.text=element_text(size=20), legend.key.size=unit(1, "cm"))
}

# Pull out a legend
legend<-g_legend(plots_list[[1]])

cairo_pdf("Fig_4.pdf", height=7, width=7)

grid.arrange(plots_list[[1]]+theme(legend.position="hidden"), plots_list[[2]]+theme(legend.position="hidden"), plots_list[[3]]+theme(legend.position="hidden"), legend, layout_matrix=(rbind(c(1,2),c(3,4))))

dev.off()

cairo_pdf("Fig_5.pdf", height=5, width=15)

grid.arrange(plots_list2[[1]], plots_list2[[2]], plots_list2[[3]], layout_matrix=array(c(1,2,3), c(1,3)))

dev.off()


cairo_pdf("Fig_6.pdf", height=7, width=7)

grid.arrange(plots_list3[[1]]+theme(legend.position="hidden"), plots_list3[[2]]+theme(legend.position="hidden"), plots_list3[[3]]+theme(legend.position="hidden"), legend, layout_matrix=(rbind(c(1,2),c(3,4))))

dev.off()


