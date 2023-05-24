library(RColorBrewer)
library(deSolve)
library(ggplot2)

# Load parameters and dd fcts
getwd()
source("scripts/Parameters_MSFR_NR.R")

# -------------------------------------------------------
days <- c(1:24) #Compute average predation rate on the bird incubation period (24 days)

# fct to compute the decreasing of goose nest density in 24 days
Ave_pred_N2<- function(t,state,parameters)
     {with(as.list(c(t,state,parameters)),{  # unpack the parameters
              dAR2 <- MSFRg_g(N1,N2,N3) * NR_sim(hr)
              dN2 <- - dAR2
        list(c(dN2,dAR2),
            c(N1=N1,N3=N3))
        })}

#state_2d= c(N2=255,AR2=MSFRg_g(204,255,3.1)*NR_sim(10.8))
#output_2d <- data.frame(ode(y=state_2d,times=days,func=Ave_pred_N2,parms=c(N1=204,N3=3.1,hr=10.8)))
#output_2d

# fct to compute the average predation rate after 24 days
Ave_pred_N3 <- function(t,state,parameters)
 {with(as.list(c(t,state,parameters)),{  # unpack the parameters
    #N2 = N2_den(t)
    dAR3 <- MSFRp_g(N1,N2,N3) * NR_sim(hr)
    dN3 <- - dAR3
    list(c(dN3,dAR3),
        c(N1=N1,N2=N2))
    })}
#state_2d= c(N3=3.1,AR3=MSFRp_g(204,255,3.1)*NR_sim(10.8))
#output_2d <- data.frame(ode(y=state_2d,times=days,func=Ave_pred_N3,parms=c(N1=204,N2=255,hr=10.8)))

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
              output_2d <- data.frame(ode(y=state_2d,times=days,func=Ave_pred_N2,parms=c(N1=d_grid$N1_vec[k],
              N3=d_grid$N3_vec[k],hr=d_grid$hr_vec[k])))
              N2_density=output_2d$N2
              # Forcing fct for goose nest density
              N2_den <- approxfun(days,N2_density,rule=2)

              state<-c(N3=d_grid$N3_vec[k],AR3 = MSFRp_g(d_grid$N1_vec[k],d_grid$N2_vec[k],d_grid$N3_vec[k])
              *NR_sim(d_grid$hr_vec[k]))
              days <- c(1:24)
              output <- data.frame(ode(y=state,times=days,func=Ave_pred_N3,parms=c(N1=d_grid$N1_vec[k], N2=d_grid$N2_vec[k],hr=d_grid$hr_vec[k])))
              n_max_3=d_grid$N3_vec[k]
              pred_out <- output[24,]$AR3/n_max_3 # predation rate/fox after 24 days
              Out[[k]] <- c(pred_out,d_grid$N1_vec[k],d_grid$N2_vec[k],d_grid$N3_vec[k],d_grid$hr_vec[k])
          }

df_ing <- do.call(rbind,Out)
df_ing <- as.data.frame(df_ing)

colnames(df_ing) <- c("PR","N1","N2","N3","HR")
dff = df_ing

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
sce_noG = ave_PR

sce_G$pred_density=NR_sim(sce_G$HR)
sce_noG$pred_density=NR_sim(sce_noG$HR)
write.csv(sce_G,"data/ave_NS_G.csv", row.names = FALSE)
#write.csv(sce_noG,"data/ave_NS_noG.csv", row.names = FALSE)

plot(sce_G$HR,sce_G$pred_density,ylim=c(0,0.7))

plot(sce_noG$pred_density,sce_noG$NS,ylim=c(0,1),xlim=c(0.05,0.4),bty="n",lwd=3,type="l",col="blue")
lines(sce_G$pred_density,sce_G$NS,col="red",lwd=3)
abline(v=NR_sim(10.81)/50)#mean in
abline(v=NR_sim(18.2)/50)#mean out

plot(sce_noG$HR,sce_noG$NS,ylim=c(0,0.8),xlim=c(0,50),bty="n",lwd=3,type="l",col="black",lty=3,ylab="Average sandpiper nesting success", xlab="Fox home range size (km2)",cex.lab=1.4)
lines(sce_G$HR,sce_G$NS,col="black",lwd=3)
points(10.81,0.444,bg="#d7191cff",pch=22,lwd=3,cex=1.4) #combined effect
points(18.2,0.624,bg="#fee090ff",pch=21,lwd=3,cex=1.4) # FR only
points(10.81,0.363,bg="#f46d43ff",pch=21,lwd=3,cex=1.4) #NR only
points(18.2,0.555, bg="#68a0d9ff",pch=21,lwd=3,cex=1.4) #Absence
dev.copy2pdf(file="NS_Hr.pdf")
dev.off()

abline(v=10.81)#mean in
abline(v=18.2)#mean out

dev.copy2pdf(file="NS_noG_G.pdf")
dev.off()

#### VIOLON PLOT Fig. S3########
library(ggplot2)
df <- read.csv("data/Home_Range_Foxes.csv",sep=",")
df$year = as.factor(df$year)
df = df[which(df$goose_category!="outside with overlap"),]#Remove HR overlapping with goose colony
df$goose_category = as.factor(df$goose_category)

df$area = df$area * 0.7 #correction ARGOS - See Christin et al. 2015 PlOs One

in_colony = df[which(df$goose_category=="inside colony",),]
out_colony= df[which(df$goose_category=="outside colony",),]

hr_range = in_colony$area #Range in the presence of a goose colony
hr_range = out_colony$area #Range in the absence of a goose colony

hist(in_colony$area,breaks=10)
hist(out_colony$area,breaks=10)


N2=c(0)#or 255
d_grid = expand.grid(N1_vec = 204, N3_vec = 3.1, N2_vec= N2,hr_vec=hr_range)

Out <- list()
  for (k in seq_along(d_grid$N3_vec)){
              AR3 = MSFRp_g(d_grid$N1_vec[k],d_grid$N2_vec[k],d_grid$N3_vec[k])*NR_sim(d_grid$hr_vec[k])
              Out[[k]] <- c(AR3,d_grid$N1_vec[k],d_grid$N2_vec[k],d_grid$N3_vec[k],d_grid$hr_vec[k])
          }

df <- do.call(rbind,Out)
df <- as.data.frame(df)
colnames(df) <- c("PR","N1","N2","N3","HR")
max(df$PR)
df_255=df
df_0=df

df = rbind(df_0,df_255)
df$N2 = as.factor(as.character(df$N2))

p <- ggplot(df, aes(x=N2, y=PR)) + geom_violin(trim=TRUE)
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
coord_cartesian(ylim = c(0,0.25))
dev.copy2pdf(file="Violon_Fig.S3.pdf")
dev.off()
