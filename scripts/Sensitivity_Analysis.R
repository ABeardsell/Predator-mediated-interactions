library(RColorBrewer)
library(deSolve)
library(pse)
library(MCMCglmm) #truncated normal distribution
library(MASS)

#--------------------------------------------------------------------
# MSFR parameters
#--------------------------------------------------------------------
## GENERAL Parameter values and distribution
hr_dist  <- rtnorm(n = 1000, mean = 14.5, sd=12,lower=3,upper=48) #based on ARGOS data n = 113
s_dist <- rtnorm(n = 1000, mean = 82, sd=5,lower=64,upper=98) # low and upp 95%CI Predator speed (km/day) - **** mean, lower=min, upper=max
phi_dist <- rtnorm(n = 1000, mean = 0.5,sd=0.06,lower = 0.40, upper = 0.60) # low and upp 95%CI proportion of time spent active

## Parameter values and distribution - LEMMINGS
Tp1_dist <- rtnorm(n = 1000, mean = 0.024/24, sd = 0.031/24, lower = 0.0221/24, upper = 0.0259/24)
Te1_dist <- rtnorm(n = 1000, mean = 33/86400, sd = 26/86400, lower = 26/86400, upper = 40/86400)
e1_dist <- rtnorm(n = 1000, mean = 0.48, lower = 0.429, upper = 0.531) #lower mean - SE upper mean + SE n=97
To1_dist <- rtnorm(n = 1000, mean = 42/86400, sd = 26/86400, lower = 28/86400, upper =46/86400)
o1_dist <- rtnorm(n = 1000, mean = 0.32, lower = 0.273, upper = 0.367)#lower mean - SE upper mean + SE n=97
Tde1_dist <- rtnorm(n = 1000, mean = 337/86400,lower = 169/86400, upper = 505/86400) #+- 50% de la valeur
de1_dist <- rtnorm(n = 1000, mean = 0.20, lower = 0.159, upper = 0.24)#lower mean - SE upper mean + SE n=97
#337/86400 #sec to days
f41_dist <- rtnorm(n = 1000, mean = 0.51, lower = 0.48, upper = 0.54)

## Parameter values and distribution - GEESE
Tp2_dist <- rtnorm(n = 1000, mean = 0.022/24,sd=0.014/24,lower = 0.008/24, upper = 0.036/24) #Chasing time (day/nest) - mean +/- SD (lower and upper limit)
Tm2_dist <- rtnorm(n = 1000, mean = 0.138/24,sd=0.07/24,lower = 0.068/24, upper = 0.208/24) #Manipulation time (day/nest) - mean +/- SD (lower and upper limit)
f42_ua_dist <- rtnorm(n = 1000, mean = 0.93, sd=0.022, lower = 0.908, upper = 0.952) # p_vide (nest unattended)
f42_a_dist <- rtnorm(n = 1000, mean = 0.098, sd=0.011,lower = 0.087, upper = 0.109)## p_adult (nest attended) - Beta distribution donne quasi le meme resultat
p_c_ua_dist <- rtnorm(n = 1000, mean = 0.69,sd=0.115,lower = 0.5744, upper = 0.8056) # complete predation probability on unattended nest
p_c_a_dist <- rtnorm(n = 1000, mean = 0.47,sd=0.128,lower = 0.3412, upper = 0.5988)# complete predation probability on attended nest
## Parameter values and distribution - SANDPIPER
Tm3_dist <- rtnorm(n = 1000, mean = 0.069/24, sd = 0.0525/24, lower = 0.0017/24, upper = 0.15/24) # manipulation time


qdata_hr <- function(p, data) quantile(x=hr_dist, probs=p)
qdata_phi <- function(p, data) quantile(x=phi_dist, probs=p)
qdata_s <- function(p, data) quantile(x=s_dist, probs=p)
## Quantiles functions - LEMMINGS
qdata_Tp1 <- function(p, data) quantile(x=Tp1_dist, probs=p)
qdata_Te1 <- function(p, data) quantile(x=Te1_dist, probs=p)
qdata_e1 <- function(p, data) quantile(x=e1_dist, probs=p)
qdata_To1 <- function(p, data) quantile(x=To1_dist, probs=p)
qdata_o1 <- function(p, data) quantile(x=o1_dist, probs=p)
qdata_Tde1 <- function(p, data) quantile(x=Tde1_dist, probs=p)
qdata_de1 <- function(p, data) quantile(x=de1_dist, probs=p)
qdata_f41 <- function(p, data) quantile(x=f41_dist, probs=p)
## Quantiles functions - GEESE
qdata_Tp2 <- function(p, data) quantile(x=Tp2_dist, probs=p)
qdata_Tm2 <- function(p, data) quantile(x=Tm2_dist, probs=p)
qdata_f42_ua <- function(p, data) quantile(x=f42_ua_dist, probs=p)
qdata_f42_a <- function(p, data) quantile(x=f42_a_dist, probs=p)
qdata_p_c_ua <- function(p, data) quantile(x=p_c_ua_dist, probs=p)
qdata_p_c_a <- function(p, data) quantile(x=p_c_a_dist, probs=p)
## Quantiles functions - LEMMINGS
qdata_Tm3 <- function(p, data) quantile(x=Tm3_dist, probs=p)

area = 70 #spatial scale - outside or inside the goose colony
overlap = 0.18 # 2019: 7 HR in 52 km2 donc 7.4 km2 sans chevauchement - HR average 9.1 km2 - overlap=1-0.82
d1=0.0075 #km
d2_a=0.033 #km
d2_ua=0.11 #km
# Sandpiper
d3=0.085

# Parameters
factors_global <- c("hr","phi","s","f21_f31","Tp1","Te1","e1","To1","o1","Tde1","de1","f41","f22_ua","f32_a","Tp2",
"Tm2","f42_ua","f42_a","p_c_ua","p_c_a","w","f23","Tm3")
# Distribution for each parameters
q_global <- c("qdata_hr","qdata_phi","qdata_s","qunif","qdata_Tp1","qdata_Te1","qdata_e1","qdata_To1","qdata_o1","qdata_Tde1","qdata_de1",
"qdata_f41","qunif","qunif","qdata_Tp2","qdata_Tm2","qdata_f42_ua",
"qdata_f42_a","qdata_p_c_ua","qdata_p_c_a","qunif","qunif","qdata_Tm3")
q.arg_global <- list(list(data=hr_dist),list(data=phi_dist),list(data=s_dist),list(min=0.075,max=0.225),list(data=Tp1_dist),
list(data=Te1_dist),list(data=e1_dist),list(data=To1_dist),list(data=o1_dist),list(data=Tde1_dist),list(data=de1_dist)
,list(data=f41_dist),list(min=0.185,max=0.555),list(min=0.025,max=0.075),list(data=Tp2_dist),
list(data=Tm2_dist),list(data=f42_ua_dist),list(data=f42_a_dist),list(data=p_c_ua_dist),list(data=p_c_a_dist)
,list(min=0.965,max=0.993),list(min=0.015,max=0.059),list(data=Tm3_dist))

##########################################################################################################
mod_MSFR_global <- function (hr,phi,s,f21_f31,Tp1,Te1,e1,To1,o1,Tde1,de1,f41,f22_ua,f32_a,Tp2,Tm2,f42_ua,f42_a,p_c_ua,p_c_a,w,f23,Tm3){
    NS <- array();
    for(i in N2) {
      N1=350
      N3 = 3.5
      Nb_renard = (area/(hr*(1-overlap))) * 2
      alpha_1=s*(2*d1)*f21_f31*f41
      alpha_3=s*(2*d3)*f23
      alpha_2a_complete = s*f32_a*(2*d2_a)*f42_a*p_c_a #Capture rate of a nest by a predator - COMPLETE
      alpha_2a_partial = (s*f32_a*(2*d2_a)*f42_a*(1-p_c_a))/3.7 #Capture rate of a nest by a predator - PARTIAL
      alpha_2a = alpha_2a_complete + alpha_2a_partial # Predation totale
      alpha_2ua_complete = s*f22_ua*(2*d2_ua)*p_c_ua*f42_ua #Capture rate of a nest by a predator - COMPLETE
      alpha_2ua_partial = (s*f22_ua*(2*d2_ua)*f42_ua*(1-p_c_ua))/3.7 #Capture rate of a nest by a predator - PARTIAL
      alpha_2ua = alpha_2ua_complete + alpha_2ua_partial # Predation totale

      h_1=  (Tp1/f41) + ((Te1 * e1) + (To1 * o1) + (Tde1*de1))
      h_2ua= Tp2/(f42_ua * p_c_ua) + Tm2
      h_2a= Tp2/(f42_a * p_c_a) + Tm2
      h_3= Tm3

      AR3 <- (alpha_3 * N3 * phi)/(1+(alpha_1 * h_1 * N1) + (alpha_3 * h_3* N3)
          + (alpha_2ua * h_2ua* ((1-w)*i)) + (alpha_2a * h_2a* (w*i)))
      PR <- AR3 * Nb_renard
      #NS[i] <- 1 - ((PR * 24)/(3.5*52))
      NS[i] <- ifelse((1 - ((PR * 24)/(3.5*70)))>0, 1 - ((PR * 24)/(3.5*70)),0)
    }
      return(NS)
  }
modelRun_MSFR_global <- function (my.data) {
    return(mapply(mod_MSFR_global,my.data[,1], my.data[,2], my.data[,3], my.data[,4],
      my.data[,5],my.data[,6],my.data[,7],my.data[,8],my.data[,9],my.data[,10],my.data[,11],my.data[,12],my.data[,13],my.data[,14],my.data[,15],
      my.data[,16],my.data[,17],my.data[,18],my.data[,19],my.data[,20],my.data[,21],my.data[,22],my.data[,23]))
}


N2=255
res.names <- paste("N2",N2) ##give colnames (which will be used in the plots below)

## Run the model
myLHS <- LHS(modelRun_MSFR_global,factors_global,1000, q_global, q.arg_global, res.names, nboot=40) # Latin Hypercube sampling (LHS)
dfdata <- as.data.frame(get.data(myLHS)) # values of parameters used in simulation
df1 <- as.data.frame(get.results(myLHS)) #Results in data frame
res3  <- df1[colSums(!is.na(df1)) > 0] # Removing the NA in data frame
head(res3)

#Scatterplot
pdf(file="Scatterplot_Geese_255.pdf",width=10,height=10)
plotscatter(myLHS,index.res=c(255),index.data=c(1,2,3,13,14,15,16,21,22),pch=19,cex.axis=1.15,cex.lab=1.15)
dev.off()

#Partial rank coefficients
pdf(file="PRCC_Geese_400.pdf")
plotprcc(myLHS, index.res=c(255),ylab="Partial correlation coefficient")
dev.off()

# plot influential parameters
df = cbind(res3,dfdata)

#Home range size
plot(df$hr,df$V255,pch=19,cex.axis=1.15,ylab="",xlab="Home range size (hr; km2)")
abline(lm(V255~hr,data=df),lwd=4)

#Phi
plot(df$phi,df$V255,pch=19,cex.axis=1.15,ylab="",xlab="Phi")
abline(lm(V255~phi,data=df),lwd=4)

# W
plot(df$w,df$V255,pch=19,cex.axis=1.15,ylab="",xlab="w")
abline(lm(V255~w,data=df),lwd=4)

# s
plot(df$s,df$V255,pch=19,cex.axis=1.15,ylab="",xlab="s")
abline(lm(V255~s,data=df))

# detection probability
plot(df$f22_ua,df$V255,pch=19,cex.axis=1.15,ylab="",xlab="f22_ua")
abline(lm(V255~f22_ua,data=df),lwd=4)

# attack probability
plot(df$f32_a,df$V255,pch=19,cex.axis=1.15,ylab="",xlab="f32_ua")
abline(lm(V255~f32_a,data=df),lwd=4)

# Pursue time prey 2
plot(df$Tp2,df$V255,pch=19,cex.axis=1.15,ylab="",xlab="Tp2")
abline(lm(V255~Tp2,data=df),lwd=4)

# Manipulation time prey 2
plot(df$Tm2,df$V255,pch=19,cex.axis=1.15,ylab="",xlab="Tm2")
abline(lm(V255~Tm2,data=df),lwd=4)

# Detection probability prey 3
plot(df$f23,df$V255,pch=19,cex.axis=1.15,ylab="",xlab="f23")
abline(lm(V255~f23,data=df),lwd=4)
dev.off()
