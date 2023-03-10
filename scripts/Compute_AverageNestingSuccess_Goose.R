#Code to compute average nesting success of goose
library(RColorBrewer)
library(deSolve)
library(ggplot2)

# Load parameters and dd fcts
source("scripts/Parameters_MSFR.R")

# general variable
days <- c(1:28) #Compute average predation rate on the bird incubation period (24 days) + laying 4 days

# Function to compute the number of predators
NR_sim<- function (hr) {
  overlap = 0.18 # fox home range overlap
  Nb_p= 2 #number of predators in H_o
  H_o = (hr*(1-overlap))
  y = Nb_p/H_o
  return(y)
}

# -------------------------------------------------------
# -------------------------------------------------------

# fct to compute the decreasing of goose nest density in 24 days
dec_gdensity <- function(t,state,parameters)
 {with(as.list(c(t,state,parameters)),{  # unpack the parameters
          dAR2 <- MSFRg_g(N1,N2,N3) * NR_sim(hr)
          dN2 <- - dAR2
    list(c(dN2,dAR2),
        c(N1=N1,N3=N3))
    })}

#state_2d= c(N2=255,AR2=MSFRg_g(100,255,3.1)*NR_sim(10))
#output_2d <- data.frame(ode(y=state_2d,times=days,func=dec_gdensity_in,parms=c(N1=100,N3=3.1,hr=10, plot_size=50,n_max_2= 255 * 50)))

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
N2_range = c(255)
N1_range_1 = c(42,281,14,384,504,9,2,648,365,253,19,1.6,137) #2007 à 2019
hr_range = 10.8 #average home range size in the goose colony

d_grid = expand.grid(N1_vec = N1_range_1, N3_vec = 3.1, N2_vec= N2_range,hr_vec=hr_range)

Out <- list()
  for (k in seq_along(d_grid$N3_vec)){
              state_2d= c(N2=d_grid$N2_vec[k],AR2=MSFRg_g(d_grid$N1_vec[k],d_grid$N2_vec[k],d_grid$N3_vec[k])
              *NR_sim(d_grid$hr_vec[k]))
              output_2d <- data.frame(ode(y=state_2d,times=days,func=dec_gdensity,parms=c(N1=d_grid$N1_vec[k],
              N3=d_grid$N3_vec[k],hr=d_grid$hr_vec[k])))
              n_max_2=d_grid$N2_vec[k]
              pred_out <- output_2d[28,]$AR2/n_max_2 # predation rate/fox after 24 days
              Out[[k]] <- c(pred_out,d_grid$N1_vec[k],d_grid$N2_vec[k],d_grid$N3_vec[k],d_grid$hr_vec[k])
          }
df_ing <- do.call(rbind,Out)
df_ing <- as.data.frame(df_ing)
colnames(df_ing) <- c("PR","N1","N2","N3","HR")
dff = df_ing

# Mean predation rate - Need to aggregate lemming density
group_mean <- aggregate(dff$PR, list(dff$HR), mean)
ave_PR = as.data.frame(group_mean)
colnames(ave_PR) = c("HR","PR")
ave_PR$NS= 1- ave_PR$PR

ave_PR
##################################################################################
