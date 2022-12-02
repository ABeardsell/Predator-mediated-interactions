library(RColorBrewer)
library(deSolve)
library(ggplot2)

# Load parameters and dd fcts
getwd()
source("scripts/Parameters_MSFR.R")

# Function to compute the number of predators
NR_sim<- function (hr) {
  area = 50 #spatial scale - outside or inside the goose colony
  overlap = 0.18 # fox home range overlap
  y = (area/(hr*(1-overlap))) * 2
  return(y)
}

# -------------------------------------------------------
days <- c(1:24) #Compute average predation rate on the bird incubation period (24 days)

# fct to compute the decreasing of goose nest density in 24 days
dec_gdensity_in <- function(t,state,parameters)
 {with(as.list(c(t,state,parameters)),{  # unpack the parameters
          dAR2 <- MSFRg_g(N1,N2,N3) * NR_sim(hr)
          dN2 <- - N2 + ((n_max_2 - AR2)/plot_size)
    list(c(dN2,dAR2),
        c(N1=N1,N3=N3,n_max_2=n_max_2))
    })}

# fct to compute the average predation rate after 24 days
Ave_pred_ing <- function(t,state,parameters)
 {with(as.list(c(t,state,parameters)),{  # unpack the parameters
    N2 = N2_den(t)
    dAR3 <- MSFRp_g(N1,N2,N3) * NR_sim(hr)
    d3 <- - N3 + ((n_max_3 - AR3)/plot_size)
    list(c(d3,dAR3),
        c(N1=N1,N3=N3,n_max_3=n_max_3,N2=N2))
    })}

#------------------------------------------------------------------------------------------------------
# ---- Goose nests and lemming densities fct of Home range size - Effects on predation rate -----------
#------------------------------------------------------------------------------------------------------
N2_range = c(255)
N1_range_1 = c(42,281,14,384,504,9,2,648,365,253,19,1.6,137) #empirical lemmings densities from 2007 to 2019
#N1_range_2 = c(204,204,204,204,204,204,204,204,204,204,204,204,204)
#N1_range_3 = c(1,322,14,500,504,9,2,648,10,600,19,1.6,20)

hr_range = seq(3.75,25.1,by=0.5) #Range in the presence of a goose colony
hr_range = seq(6.5,48.5,by=0.5) #Range in the absence of a goose colony
hr_range = seq(3.75,48.5,by=4.5) #Complete range
hr_range=c(10.8,18.2)
d_grid = expand.grid(N1_vec = N1_range_1, N3_vec = 3.1, N2_vec= N2_range,hr_vec=hr_range)

Out <- list()
  for (k in seq_along(d_grid$N3_vec)){
              state_2d= c(N2=d_grid$N2_vec[k],AR2=MSFRg_g(d_grid$N1_vec[k],d_grid$N2_vec[k],d_grid$N3_vec[k])
              *NR_sim(d_grid$hr_vec[k]))
              output_2d <- data.frame(ode(y=state_2d,times=days,func=dec_gdensity_in,parms=c(N1=d_grid$N1_vec[k],
              N3=d_grid$N3_vec[k],hr=d_grid$hr_vec[k], plot_size=50,n_max_2= d_grid$N2_vec[k] * 50)))
              N2_density=output_2d$N2
              # Forcing fct for goose nest density
              N2_den <- approxfun(days,N2_density,rule=2)

              state<-c(N3=d_grid$N3_vec[k],AR3 = MSFRp_g(d_grid$N1_vec[k],d_grid$N2_vec[k],d_grid$N3_vec[k])
              *NR_sim(d_grid$hr_vec[k]))
              days <- c(1:24)
              output <- data.frame(ode(y=state,times=days,func=Ave_pred_ing,parms=c(N1=d_grid$N1_vec[k],
                N3=d_grid$N3_vec[k], N2=d_grid$N2_vec[k],hr=d_grid$hr_vec[k], plot_size=50,n_max_3=d_grid$N3_vec[k] * 50)))
              pred_out <- output[24,]$AR3/output[24,]$n_max_3 # predation rate/fox after 24 days
              Out[[k]] <- c(pred_out,d_grid$N1_vec[k],d_grid$N2_vec[k],d_grid$N3_vec[k],d_grid$hr_vec[k])
          }

df_ing <- do.call(rbind,Out)
df_ing <- as.data.frame(df_ing)
colnames(df_ing) <- c("PR","N1","N2","N3","HR")
dff = df_ing
dff

################################################################################
# Mean predation rate - Need to aggregate lemming density
group_mean <- aggregate(dff$PR, list(dff$HR), mean)
ave_PR = as.data.frame(group_mean)
colnames(ave_PR) = c("HR","PR")

#Write csv of the entire range of home range and the associated nesting success for the matrix model
#write.csv(ave_PR,"data/ave_PR.csv", row.names = FALSE)

##################################################################################
#Plot Average Nesting success in presence and absence of goose
ave_PR$NS= 1- ave_PR$PR

#sce_noG=ave_PR
sce_G=ave_PR
sce_noG
sce_G

sce_G$pred_density=NR_sim(sce_G$HR)/50
sce_noG$pred_density=NR_sim(sce_noG$HR)/50

plot(sce_noG$pred_density,sce_noG$NS,ylim=c(0,1),xlim=c(0.05,0.4),bty="n",lwd=3,type="l",col="blue")
lines(sce_G$pred_density,sce_G$NS,col="red",lwd=3)
abline(v=NR_sim(10.81)/50)#mean in
abline(v=NR_sim(18.2)/50)#mean out

plot(sce_noG$HR,sce_noG$NS,ylim=c(0,1),xlim=c(0,50),bty="n",lwd=3,type="l",col="blue")
lines(sce_G$HR,sce_G$NS,col="red",lwd=3)
abline(v=10.81)#mean in
abline(v=18.2)#mean out

dev.copy2pdf(file="NS_noG_G.pdf")
dev.off()

###########################################################################
# PLOT LEMMING TIME SERIES
###########################################################################
N1_range_1 = c(42,281,14,384,504,9,2,648,365,253,19,1.6,137) #2007 à 2019
year_l = c(2007:2019)

plot(year_l,N1_range_1,type="b",pch=19,xaxt="n",ylim=c(0,700),ylab="",xlab="",bty="n",cex.axis=1.1,cex.lab=1.2,lwd=3)
lines(year_l,N1_range_2,type="b",col="lightgrey",lwd=3,pch=19)
lines(year_l,N1_range_3,type="b",col="firebrick",lwd=3)

axis(1,at=seq(from=2007,to=2019,by=1),font.lab=2,cex.axis=1.1,cex.lab=1.2)
dev.copy2pdf(file="Lemmings_timeseries.pdf",width=8, height=6)
dev.off()

########JENSEN INEQUALITY###############################
ave_PR$NS= 1- ave_PR$PR
sce1=ave_PR

plot(sce1$HR,sce1$NS,ylim=c(0,1),xlim=c(0,50),bty="n",lwd=3,type="l",col="blue")
lines(sce2$HR,sce2$NS,col="lightgrey",lwd=3)
lines(sce3$HR,sce3$NS,col="firebrick",lwd=3)
dev.copy2pdf(file="NS_JensenInequality.pdf")
dev.off()
