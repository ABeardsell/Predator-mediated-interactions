#Figure S2 Empirical time series of lemming density on Bylot Island from 2007 to 2019 measured by live-trapping (see methods in Fauteux et al. 2018). The density of brown and collared lemmings was summed.
N1_range = c(42,281,14,384,504,9,2,648,365,253,19,1.6,137) #2007 à 2019
year_l = c(2007:2019)

plot(year_l,N1_range,type="b",pch=19,xaxt="n",ylim=c(0,700),ylab="",xlab="",bty="n",cex.axis=1.1,cex.lab=1.2,lwd=3)
axis(1,at=seq(from=2007,to=2019,by=1),font.lab=2,cex.axis=1.1,cex.lab=1.2)
dev.copy2pdf(file="Lemmings_timeseries.pdf",width=8, height=6)
dev.off()
