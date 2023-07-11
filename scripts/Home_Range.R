###########################################################################
#Estimated summer home range size of arctic foxes using telemetry data (Argos) of 113 foxes from 2008 to 2016 on Bylot Island.
###########################################################################
df <- read.csv("data/Home_Range_Foxes.csv",sep=",")
df$year = as.factor(df$year)
df = df[which(df$goose_category!="outside with overlap"),]#Remove HR overlapping with goose colony
df$goose_category = as.factor(df$goose_category)

df$area = df$area * 0.7 #correction ARGOS - See Christin et al. 2015 PlOs One
in_colony = df[which(df$goose_category=="inside colony",),]
out_colony= df[which(df$goose_category=="outside colony",),]

summary(in_colony$area)
summary(out_colony$area)
length(out_colony$area)
length(in_colony$area)

##------------------------------------------------------------------------
# Plot DV and predator density Fig3A
hr = c(3.75:48)

# Function to compute the number of predators
NR_sim<- function (hr) {
  overlap = 0.18 # fox home range overlap
  Nb_p= 2 #number of predators in H_o
  H_o = (hr*(1-overlap))
  y = (Nb_p/H_o)
  return(y)
}
pred_density=NR_sim(hr)

plot(hr,pred_density,type="l",bty="n",xlim=c(0,50),ylim=c(0,0.7),lwd=4)
abline(v=18.2)
abline(v=10.8)
dev.copy2pdf(file="hr_preddensity.pdf")
dev.off()

hist(out_colony$area,xlim=c(0,50),ylim=c(0,25),breaks = seq(min(in_colony$area),max(out_colony$area), length.out = 10),col="#68a0d951")
hist(in_colony$area,breaks = seq(min(in_colony$area),max(out_colony$area), length.out = 10),add=T,col="#ff970041")
dev.copy2pdf(file="hist_both.pdf")
dev.off()
