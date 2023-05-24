source("scripts/Parameters_MSFR_NR.R")

# Function to compute the number of predators
NR_sim<- function (hr) {
  overlap = 0.18 # fox home range overlap
  Nb_p= 2 #number of predators in H_o
  H_o = (hr*(1-overlap))
  y = (Nb_p/H_o)
  return(y)
}

#Plot Functional response of prey 3 to prey 1 and prey 2 densities
N3=(0:7)

# FR
plot(MSFRp_g(700,0,N3)~N3,ylim=c(0,2),type="l",lwd=4,col="black",bty="n",lty=3,ylab="Number of nests predated per fox per day")
lines(MSFRp_g(0,0,N3)~N3,lwd=4,col="black")
lines(MSFRp_g(0,255,N3)~N3,lwd=4,col="firebrick")
lines(MSFRp_g(700,255,N3)~N3,lwd=4,col="firebrick",lty=3)
dev.copy2pdf(file="FR_sandpipers_densityN2_N3.pdf")
dev.off()

# FR x NR
plot(MSFRp_g(700,0,N3)*NR_sim(18.2)~N3,ylim=c(0,0.4),type="l",lwd=4,col="black",bty="n",lty=3,ylab="Total number of sandpiper nests predated per km2 per day")
lines(MSFRp_g(0,0,N3)*NR_sim(18.2)~N3,lwd=4,col="black")
lines(MSFRp_g(0,255,N3)*NR_sim(10.8)~N3,lwd=4,col="firebrick")
lines(MSFRp_g(700,255,N3)*NR_sim(10.8)~N3,lwd=4,col="firebrick",lty=3)

dev.copy2pdf(file="FRNR_sandpipers_densityN2_N3.pdf")
dev.off()

#N1
N3=(3.1)
N1=(1:700)

plot(MSFRp_g(N1,0,N3)*NR_sim(18.2)~N1,ylim=c(0,0.2),type="l",lwd=4,col="black",bty="n",lty=3,ylab="Total number of sandpiper nests predated per km2 per day")
lines(MSFRp_g(N1,255,N3)*NR_sim(10.8)~N1,lwd=4,col="firebrick")

dev.copy2pdf(file="FRNR_sandpipers_densityN2_N3.pdf")
dev.off()

plot(MSFRp_g(N1,0,N3)~N1,ylim=c(0,1),type="l",lwd=4,col="black",bty="n",lty=3,ylab="Number of nests predated per fox per day")
lines(MSFRp_g(N1,255,N3)~N1,lwd=4,col="firebrick")
dev.copy2pdf(file="FRNR_sandpipers_densityN2_N3.pdf")
dev.off()

#N2
N2=(50:300)
N3=(3.1)

plot(MSFRp_g(0,N2,N3)~N2,ylim=c(0,1),type="l",lwd=4,col="black",bty="n",lty=3,ylab="Number of nests predated per fox per day")
lines(MSFRp_g(700,N2,N3)~N2,lwd=4,col="firebrick")
dev.copy2pdf(file="FRN3_N2.pdf")
dev.off()
