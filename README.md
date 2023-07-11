# Predator-mediated-interactions
Code for the article entitled: Predator-mediated interactions through changes in predator home range size can lead to local prey exclusion

#Parameters_MSFR_NR.R
In this file, you will find the parameter values and the main equations of the predation model.

#Compute_AverageNestingSuccess_Goose. R
To estimate the nesting success of snow geese.

#Compute_AverageNestingSuccess_Sandpipers
 We estimated annual nesting success of sandpipers (prey 3) for the whole range of home range sizes observed in presence (N_2 = 255 nests per km^2) and absence (N_2 = 0 nest per km^2) of the goose colony using a set of differential equations. These equations allowed us to calculate the total number of nests predated per km$^2$ over the sandpiper nesting period (i.e., the average duration between the laying date and hatching date) while considering that the density of nests decreases each day. This script produces figure 3B and S3(C2).

#Sensitivity_Analysis.R
We quantified the relative influence of model parameter values on the estimation of sandpiper annual nesting success by using the Latin hypercube sampling technique (an efficient implementation of the Monte Carlo methods; Marino et al. 2008). Each parameter was represented by a probability distribution (uniform or normal truncated) based on the distribution of empirical data. For some parameters, the biological information was limited, so we assigned a uniform distribution allowing for a large range bounded by minimum and maximum values. Latin hypercube sampling was then applied to each distribution (N = 1,000 iterations). For simplicity, the sensitivity analysis was conducted on the predation model including the presence of the goose colony, without density-dependence in parameters $s$ and $\phi_{active}$, and with fixed prey densities ($N_1$ = 204 individuals/km$^2$ , $N_2$ = 255 nests/km$^2$, $N_3$ = 3.1 nests/km $^2$).

This script produces Fig.S4 and S5.

#Home_Range.R
We estimated summer home range size of arctic foxes using telemetry data (Argos) of 113 foxes from 2008 to 2016 on Bylot Island. We estimated the area of the 95% home range contour for each individual-year between May-October using the autocorrelation-informed home range estimation workflow described in Fleming et al. 2015, and implemented in ctmm R package (Calabrese et al. 2016, Dulude et al. in prep.). This script produces Figure 3A.

#Population_Matrix_Model.R
This code contains the sandpiper population matrix model. It generates Figures 3C and S7.
