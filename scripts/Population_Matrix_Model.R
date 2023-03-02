## Set up workspace ----
library(popbio)
library(ggplot2)
options(scipen=10, stringsAsFactors=FALSE)

## Source functions ----
# To calculate mean clutch size:
fn_clutch <- function(p4, p3e, p2e, p1e){
	clutch <- mean(c(
		rep(4,100*p4),
		rep(3,100*(1-p4)*p3e),
		rep(2,100*(1-p4)*p2e),
		rep(1,100*(1-p4)*p1e)))
	clutch
}

################################################################################
## Input values of the population model:

# Probability of returning to breeding area:
		pbreed1 <- 0.67		#  Juv
		pbreed2 <- 0.925	# Yearlings
		pbreed3 <- 1					# 2+
# Rates for first nest attempt:
		p4e1 <- 0.90	# probability that a clutch has 4 eggs
		csurv1 <- 0.71	# probability of a chick surviving from hatch to fledge
# Rates for second nest attempt:
		prnestm <- 0.73 # probability of renesting  if the first nest fails
		p4e2 <- 0.78	# probability that a clutch has 4 eggs
		csurv2 <- 0.23	# probability of a chick surviving from hatch to fledge
# For clutches with <4 eggs, proportion with 1, 2, or 3, respectively:
			p1e <- 0.014
			p2e <- 0.143
			p3e <- 0.843
# Proportion of eggs that hatch:
			phatch <- 0.94
# Sex ratio:
			sexr <- 0.5
# Survival rates:
		sJuv <- 0.44	# juvenile survival

################################################################################
# deterministic simulation
#Code for producing BOXPLOT in figure 3B2:
#############################################################
NS_range=c(0.459, 0.383,0.566,0.633) #Average nesting success for the different combinaisons (numerical response only, functional response only, absence of geese and combined effects)
pnestm_range <- seq(0.8,0.95,by=0.05)# Nesting probability range
as_range = seq(0.722,0.798, by=0.01) # adult survival range

grid = expand.grid(NS = NS_range, as= as_range,pnestm=pnestm_range)

Out <- list()
  for (k in seq_along(grid$NS)){

    fn1 <- grid$pnestm[k]*(grid$NS[k])*
      (fn_clutch(p4e1, p3e, p2e, p1e)*phatch)*
      csurv1*sexr
    fn2 <- grid$pnestm[k]*(1 - (grid$NS[k]))*prnestm*(grid$NS[k])*
      (fn_clutch(p4e2, p3e, p2e, p1e)*
      phatch)*csurv2*sexr

      # Multiply by the probability of returning at each age to calculate the
      # number of male fledglings produced per male in the whole population:
			F1 <- sJuv * pbreed1*(fn1 + fn2)
			F2 <- grid$as[k] * pbreed2*(fn1 + fn2)
			F3 <- grid$as[k] * pbreed3*(fn1 + fn2)

			# Set up the transition matrix:
				mat <- matrix( c(
					F1, F2, F3,
					sJuv, 0, 0,
					0, grid$as[k], grid$as[k]), nrow=3, ncol=3, byrow=TRUE)

        Out[[k]] <- lambda(mat)
      }

sims <- do.call(rbind,Out)
sims <- as.data.frame(cbind(sims, grid$NS,grid$as))
colnames(sims) <- c("Lambda","NS","AS")

sims$fact = cut(sims$NS,4, labels=c(3,4,1,2))
sims$fact=as.numeric(as.character(sims$fact))


group_median <- aggregate(sims$Lambda, list(sims$fact), median)

par(mfrow = c(1,2))
boxplot(sims$Lambda~sims$fact,col=c("#68a0d990","#fee09094","#f46d439d","#d7191ca9"),xaxt="n")
abline(h=1)#lambda

dev.copy2pdf(file="BoxPlot_Lambda.pdf")
dev.off()

##############################################################
#Code for producing PLOT for figure S6
#############################################################
ave_PR <- read.csv("data/ave_PR.csv")
ave_PR$NS=1-(ave_PR$PR) #Nesting success for different values of home range

as_range = seq(0.5,0.95,by=0.05) # adult survival range
NS_range = ave_PR$NS
pnestm_range = 0.8 #probability of nesting in a given year

grid = expand.grid(NS = NS_range, as= as_range,pnestm=pnestm_range)
Out <- list()
  for (k in seq_along(grid$NS)){
    fn1 <- grid$pnestm[k]*(grid$NS[k])*
      (fn_clutch(p4e1, p3e, p2e, p1e)*phatch)*
      csurv1*sexr

    fn2 <- grid$pnestm[k]*(1 - (grid$NS[k]))*prnestm*(grid$NS[k])*
      (fn_clutch(p4e2, p3e, p2e, p1e)*
      phatch)*csurv2*sexr
      # Multiply by the probability of returning at each age to calculate the
      # number of male fledglings produced per male in the whole population:
			F1 <- sJuv * pbreed1*(fn1 + fn2)
			F2 <- grid$as[k] * pbreed2*(fn1 + fn2)
			F3 <- grid$as[k] * pbreed3*(fn1 + fn2)

			# Set up the transition matrix:
				mat <- matrix( c(
					F1, F2, F3,
					sJuv, 0, 0,
					0, grid$as[k], grid$as[k]), nrow=3, ncol=3, byrow=TRUE)
        Out[[k]] <- lambda(mat)
      }

sims <- do.call(rbind,Out)
sims <- as.data.frame(cbind(sims, grid$NS,grid$as))
final = cbind(ave_PR$HR,sims)
colnames(final) = c("HR","L","NS","AS")

#PLOT For Lambda<1
ggplot(final, aes(HR,AS, fill=L)) + geom_tile() + scale_fill_distiller(limits = c(0.55,1),palette = "OrRd",direction=-1)+
theme_bw() +
geom_hline(yintercept=0.76)+
geom_vline(xintercept=10.8)+
coord_cartesian(xlim = c(2.5,47.5)) +
scale_x_continuous(breaks = seq(2.5,47.5, by = 7.5)) +
  theme(axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
		axis.text.y= element_text(face="bold", colour="black",
	                                         size=12),
		axis.text.x= element_text(face="bold", colour="black",
	                                         size=12))
dev.copy2pdf(file="HR_AS_L_decreaseG255.pdf")
dev.off()

#PLOT For Lambda>1
ggplot(final, aes(HR,AS, fill=L)) + geom_tile() + scale_fill_distiller(limits = c(1,1.4),palette = "Blues",direction=1)+
theme_bw() +
coord_cartesian(xlim = c(2.5,47.5)) +
scale_x_continuous(breaks = seq(2.5,47.5, by = 7.5)) +
  theme(axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
		axis.text.y= element_text(face="bold", colour="black",
	                                         size=12),
		axis.text.x= element_text(face="bold", colour="black",
	                                         size=12))

    dev.copy2pdf(file="HR_AS_L_increaseG350.pdf")
    dev.off()

# Addind predator density in x axis
NR_sim<- function (hr) {
  area = 50 #spatial scale - outside or inside the goose colony
  overlap = 0.18 # 2019: 7 HR in 52 km2 donc 7.4 km2 sans chevauchement - HR average 9.1 km2 - overlap=1-0.82
  y = (area/(hr*(1-overlap))) * 2
  return(y)
}
final$pred_density=NR_sim(final$HR)/50
