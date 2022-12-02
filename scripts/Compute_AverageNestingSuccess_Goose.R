#Code to compute average nesting success of goose
library(RColorBrewer)
library(deSolve)
library(ggplot2)

# Load parameters and dd fcts
source("scripts/Parameters_MSFR.R")

# general variable
days <- c(1:28) #Compute average predation rate on the bird incubation period (24 days) + laying 4 days

NR_sim<- function (hr) {
  area = 50 #spatial scale - outside or inside the goose colony
  overlap = 0.18 # 2019: 7 HR in 52 km2 donc 7.4 km2 sans chevauchement - HR average 9.1 km2 - overlap=1-0.82
  y = (area/(hr*(1-overlap))) * 2
  return(y)
}

# -------------------------------------------------------
#  INSIDE - MSFR - WITH GEESE
# -------------------------------------------------------

# fct to compute the decreasing of goose nest density in 24 days
dec_gdensity_in <- function(t,state,parameters)
 {with(as.list(c(t,state,parameters)),{  # unpack the parameters
          dAR2 <- MSFRg_g(N1,N2,N3) * NR_sim(hr)
          dN2 <- - N2 + ((n_max_2 - AR2)/plot_size)
    list(c(dN2,dAR2),
        c(N1=N1,N3=N3,n_max_2=n_max_2))
    })}

#state_2d= c(N2=255,AR2=MSFRg_g(100,255,3.1)*NR_sim(10))
#output_2d <- data.frame(ode(y=state_2d,times=days,func=dec_gdensity_in,parms=c(N1=100,N3=3.1,hr=10, plot_size=50,n_max_2= 255 * 50)))

#------------------------------------------------------------------------------------------------------
# ---- Goose nests and lemming densities fct of Home range size - Effects on predation rate -----------
#------------------------------------------------------------------------------------------------------
N2_range = c(255)
N1_range_1 = c(42,281,14,384,504,9,2,648,365,253,19,1.6,137) #2007 à 2019
hr_range = 10.8 #average home range size in the goose colony

d_grid = expand.grid(N1_vec = N1_range_1, N3_vec = 3.1, N2_vec= N2_range,hr_vec=hr_range)

Out <- list()
  for (k in seq_along(d_grid$N3_vec)){
              state_2d= c(N2=d_grid$N2_vec[k],AR2=MSFRg_g(d_grid$N1_vec[k],d_grid$N2_vec[k],d_grid$N3_vec[k])
              *NR_sim(d_grid$hr_vec[k]))
              output_2d <- data.frame(ode(y=state_2d,times=days,func=dec_gdensity_in,parms=c(N1=d_grid$N1_vec[k],
              N3=d_grid$N3_vec[k],hr=d_grid$hr_vec[k], plot_size=50,n_max_2= d_grid$N2_vec[k] * 50)))
              pred_out <- output_2d[28,]$AR2/output_2d[28,]$n_max_2 # predation rate/fox after 24 days
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
