#Natural mortality
M.fun=function(AGE,Amax,age.mat,LinF,kk,awt,bwt,Lo)
{
  #STEP 1. calculate M from different methods (see Kenchington 2013)
  
  #1.1. Age-independent
  #Jensen (1996)
  m.Jensen.2=1.65/age.mat
  m.Jensen.2=rep(m.Jensen.2,length(AGE))
  
  #Pauly (1980)  
  #m.Pauly=10^(-0.0066-0.279*log10(LinF)+0.6543*log10(kk)+0.4634*log10(Aver.T)) #not used, updated by Then et al (2015)
  m.Pauly=4.118*(kk^0.73)*LinF^(-0.33)
  m.Pauly=rep(m.Pauly,length(AGE))
  
  #Hoenig (1983), combined teleost and cetaceans    
  #m.Hoenig=exp(1.44-0.982*log(Amax))      #not used, updated by Then et al (2015)
  #m.Hoenig=rep(m.Hoenig,length(AGE))
  
  #Then et al (2015)
  m.Then.1=4.899*Amax^(-0.916)
  m.Then.1=rep(m.Then.1,length(AGE))
  
  #1.2. Age-dependent
  #Peterson and Wroblewski 1984 (dry weight in grams, length in cm)
  Dry.w=0.2   # Cortes (2002)
  TL=Lo+(LinF-Lo)*(1-exp(-kk*AGE))
  wet.weight=1000*awt*TL^bwt
  m.PetWro=1.92*(wet.weight*Dry.w)^-0.25
  m.PetWro[m.PetWro>1]=NA
  
  #Lorenzen 1996 (weight in grams)
  m.Lorenzen=3*wet.weight^-0.288
  m.Lorenzen[m.Lorenzen>1]=NA
  
  #Gislason et al (2010) (length in cm)
  m.Gislason=1.73*(TL^-1.61)*(LinF^1.44)*kk
  m.Gislason[m.Gislason>1]=NA
  
  
  #STEP 2. get mean at age
  nat.mort=data.frame(m.Jensen.2,m.Pauly,m.Then.1,
                      m.PetWro,m.Lorenzen,m.Gislason)  
  #for dogfish, due to their small size, weight-based M estimators highly overestimate M
  if(mean(m.PetWro/m.Gislason,na.rm=T)>5)  nat.mort=data.frame(m.Jensen.2,m.Pauly,m.Then.1,
                                                               m.Gislason)  
  
  #STEP 3. Calculate mean
  if(M.averaging=='mean') MoRt=rowMeans(nat.mort,na.rm=T)
  if(M.averaging=='min') MoRt=apply(nat.mort, 1, function(x) min(x,na.rm=T))
  if(M.averaging=='weighted.mean')
  {
    W8ts=rep(1,ncol(nat.mort))
    MoRt=apply(nat.mort, 1, function(x) weighted.mean(x,W8ts))
  }
  
  if(MoRt[1]<MoRt[2]) MoRt[1]=1.2*MoRt[2]  #for analysed species, M[1] is 1.2 times M[2] on average
  
  return(MoRt)
}