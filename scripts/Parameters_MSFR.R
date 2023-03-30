#--------------------------------------------------------------------
# MSFR parameters
#--------------------------------------------------------------------
s=93 #Daily distance traveled - km/day
phi=0.5#Average proportion of time the predator spend active in a day

# Lemmings
d1=0.0075 #km
f21_f31=0.15 #detection and attack probabilities
f41=0.51 # success probability of an attack
Tp1 = 0.024/24 #Pursue time - day
Te1= 0.00038 #consumption time - day
e1= 0.48 # consumption probability
To1= 0.000486 # hoarding time- day
o1= 0.32 # hoarding probability
Tde1= 0.00388# delivery time - day
de1= 0.20 #delivery probability

# fcts densité-dépendance Lemmings - From Beardsell et al. Ecology 2022.
s_L_concav <- function(x) {
64 + (93-64)* (exp(-x/150))
}

phi_L_concav <- function(x) {
0.40 + (0.53-0.40)* (exp(-x/180))
}

#N1= c(0:700)
#plot(N1,s_L_concav(N1),type="l",lwd=4,bty="n",ylim=c(60,100))
#plot(N1,phi_L_concav(N1),type="l",lwd=4,bty="n",ylim=c(0.35,0.55))

# Sandpiper
d3=0.085 #detection range - km
f23=0.029 # detection probability
Tm3=0.069/24 #consumption time - day

# Geese
w= 1-0.021 #0.021 # nest attendance prob.
d2_a=0.033 #Attack distance - km
d2_ua=0.11 #Detection distance - km
f22_ua=0.366 # detection probability
f42_a=0.098 # success probability
f42_ua=0.934 # success probability
f32_a=0.05#0.11 # attack probability (nest attended)
p_c_ua=0.69# Complete predation probability of unattended nest
p_c_a=0.47#Complete predation probability of attended nest
Tp2=0.02/24 #Pursue time - day/nest
Tm2 = 0.14/24 #Manipulation time - day/nest
# -------------------------------------------------------
# Model inside the goose colony
# -------------------------------------------------------
MSFRp_g <- function(N1,N2,N3) {
  phi = phi_L_concav(N1)
  s = s_L_concav(N1)

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
      + (alpha_2ua * h_2ua* ((1-w)*N2)) + (alpha_2a * h_2a* (w*N2)))
  return(AR3)
}

MSFRg_g <- function(N1,N2,N3) {
  phi = phi_L_concav(N1)
  s = s_L_concav(N1)

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

  AR2ga <- (alpha_2a * (w*N2) * phi)/(1+(alpha_1 * h_1 * N1) + (alpha_3 * h_3* N3)
        + (alpha_2ua * h_2ua* ((1-w)*N2)) + (alpha_2a * h_2a* (w*N2)))
  AR2gua <- (alpha_2ua * ((1-w)*N2) * phi)/(1+(alpha_1 * h_1 * N1) + (alpha_3 * h_3* N3)
        + (alpha_2ua * h_2ua* ((1-w)*N2)) + (alpha_2a * h_2a* (w*N2)))
  AR2=AR2ga+AR2gua
  return(AR2)
  #list("ARl"=ARl,"ARp"=ARp,"ARs"=ARs,"ARg"=ARg,"ARga"=ARga,"ARgua"=ARgua,"L"=x,"P"=y,"S"=z,"G"=a)
}
