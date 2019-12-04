# Code for the Speciation Algorithm by Bob Jackson
require(ggplot2)
require(reshape2)

# Initial Species Values
z1 = 100# WT homozygotes in niche 1
z2 = 1# heterozygotes in niche 1
z3 = 1# re-arranged homozygotes in niche 1
z4 = 0# WT homozygotes in niche 1
z5 = 0# heterozygotes in niche 2
z6 = 0# re-arranged homozygotes in niche 1
z0 = z1+z2+z3 
z7 = z4+z5+z6 

# Population Maximum Values
n1max = 10000 # maximum allowed in niche 1
n2max = 10000 # maximum allowed in niche 2
capacity = n1max+n2max # max. population size across both niches

# Initial time and time-step
t = 0
dt = 0.1

# Doubling times
DT1 = 2 # doubling time of z1 WT
DT4 = 4 # doubling time of z4 WT
DT = c(DT1,2*DT1,DT1,DT4,2*DT4,DT4)
fac = exp(dt*log(2)/DT)

# Survival probability when leaving Niche 1
sf = c(0.1,0.1,0.1)


# general population death rate - set as the time-step value
LF<-c(dt,dt,dt,dt,dt,dt)/10

# data-frame that will keep track of each time-step
rt<-as.data.frame(array(NA,c(1,9)))
colnames(rt)<-c("Time","Z1","Z2","Z3","Z4","Z5","Z6","T1","T2")
rt[1,]<-c(t,z1,z2,z3,z4,z5,z6,z0,z7)
# run the algorithm
while (c(z0+z7<=capacity) & c(z2>=0.5)){
  
  # production
  z1 = z1*fac[1]
  z2 = z2*fac[2]
  z3 = z3*fac[3]
  z4 = z4*fac[4]
  z5 = z5*fac[5]
  z6 = z6*fac[6]
  
  # death
  z1 = z1-z1*LF[1]
  z2 = z2-z2*LF[2]
  z3 = z3-z3*LF[3]
  z4 = z4-z4*LF[4]
  z5 = z5-z5*LF[5]
  z6 = z6-z6*LF[6]
  
  # calculate totals
  z0 = z1+z2+z3
  z7 = z4+z5+z6
  
  # calculate gametes
  ga1 = 2*z1+z2
  ga2 = z2+2*z3
  ga4 = 2*z4+z5
  ga5 = z5+2*z6
  
  # calculate gamete totals
  ga0 = ga1+ga2
  ga7 = ga4+ga5
  
  # calculate fraction of gametes
  gfrac1 = ga1/ga0
  gfrac2 = ga2/ga0
  # if statement needed to prevent division by 0 and creating NaNs in R
  if(ga7>0){
    gfrac3 = ga4/ga7
    gfrac4 = ga5/ga7
  }else{
    gfrac3 = 0
    gfrac4 = 0
  }
  
  #fertilisation
  z1 = (ga1*gfrac1)/2
  z2 = (ga1*gfrac2 +ga2*gfrac1)/2
  z3 = (ga2*gfrac2)/2
  z4 = (ga4*gfrac3)/2
  z5 = (ga4*gfrac4 + ga5*gfrac3)/2
  z6 = (ga5*gfrac4)/2
  
  # re-calculate totals
  z0 = z1+z2+z3
  z7 = z4+z5+z6
  
  # calculate fractions
  f1 = z1/z0
  f2 = z2/z0
  f3 = z3/z0
  # if statement needed to prevent division by 0 and creating NaNs in R
  if(z7>0){
    f4 = z4/z7
    f5 = z5/z7
    f6 = z6/z7
  }else{
    f4 = 0
    f5 = 0
    f6 = 0
  }
  
  
  # migration from niche 1 to niche 2 if we have reached capacity in niche 1
  
  if (z0>n1max){
    # leaving niche 1
    z1 = z1-(z0-n1max)*f1
    z2 = z2-(z0-n1max)*f2
    z3 = z3-(z0-n1max)*f3
    # entering niche 2 based on surviving the journey 
    z4 = z4+(z0-n1max)*sf[1]*f1
    z5 = z5+(z0-n1max)*sf[2]*f2
    z6 = z6+(z0-n1max)*sf[3]*f3
    
  }
  # keep track of total numbers
  z0 = z1+z2+z3
  z7 = z4+z5+z6 
  
  # increment time and record population
  t = t+dt
  rt<-rbind(rt,c(t,z1,z2,z3,z4,z5,z6,z0,z7))
  
}

rt<-melt(rt,id=c("Time"))
ggplot(rt,aes(Time,value,col=variable,group=variable))+geom_line()+
  theme_bw(base_size=16)+xlab("Time (Days)")+ylab("Population (no. species)")+
  scale_y_log10(lim=c(0.5,capacity))
