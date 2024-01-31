#Calculate M from different methods (Kenchington 2013, Dureuil et al 2021)

#notes: Age-based methods more robust (Maunder et al 2023)
#       PetWro, Lorenzen & Gislason yield too high M and unrealistic h (below 0.2) for several species if k & Amax are high

M.fun=function(AGE,Amax,age.mat,LinF,kk,awt,bwt,Lo,
               weight.method=c(Dureuil=1),weight.method.at.age=c(Dureuil.at.age=1))  
{
  #---1. Calculate M by method
  TL=Lo+(LinF-Lo)*(1-exp(-kk*AGE))
  
    #1.1. Age-independent
  #Jensen (1996)
  m.Jensen.2=1.65/age.mat
  m.Jensen.2=rep(m.Jensen.2,length(AGE))
  
  #Pauly (1980)  
  #m.Pauly=10^(-0.0066-0.279*log10(LinF)+0.6543*log10(kk)+0.4634*log10(Aver.T)) 
  m.Pauly=4.118*(kk^0.73)*LinF^(-0.33)
  m.Pauly=rep(m.Pauly,length(AGE))
  
  #Hoenig (1983), combined teleost and cetaceans    
  m.Hoenig=exp(1.44-0.982*log(Amax))      
  m.Hoenig=rep(m.Hoenig,length(AGE))
  
  #Then et al (2015)
  m.Then.1=4.899*Amax^(-0.916)
  m.Then.1=rep(m.Then.1,length(AGE))
  
  #Dureuil et al 2021 (adults)     
  P=0.0178  #for elasmos
  tmax=(1/kk)*log((LinF-Lo)/((1-0.99)*LinF))  #Cortes & Taylor 2023
  AAA=Amax
  if(reset.max.Age) AAA=max(Amax,tmax)
  Mr1=exp(1.551-1.066*log(AAA))
  Mr2=-log(P)/AAA
  Mr=mean(Mr1,Mr2)
  Dureuil=rep(Mr,length(AGE))

  #Hamel & Cope 2022 
  Hamel.Cope.tmax=5.40/AAA
  Hamel.Cope.tmax=rep(Hamel.Cope.tmax,length(AGE))
  Hamel.Cope.k=1.55*kk
  Hamel.Cope.k=rep(Hamel.Cope.k,length(AGE))
  
  
    #1.2. Age-dependent
  
  #Dureuil et al 2021 (juveniles)
  x=(2/(log(P)))+1
  ta=((2/(log(P))+1)*AAA)
  E=(AAA-(x*AAA))/2
  Lta=Lo+(LinF-Lo)*(1-exp(-kk*ta))
  Dureuil.at.age=Mr*(Lta/TL)  
  
  #Peterson and Wroblewski 1984 (dry weight in grams, length in cm)
  Dry.w=0.2   # Cortes (2002)
  wet.weight=1000*awt*TL^bwt
  m.PetWro=1.92*(wet.weight*Dry.w)^-0.25
  m.PetWro[m.PetWro>1]=NA
  
  #Lorenzen 1996 (weight in grams)
  m.Lorenzen=3*wet.weight^-0.288
  m.Lorenzen[m.Lorenzen>1]=NA
  
  #Gislason et al (2010) (length in cm)
  m.Gislason=1.73*(TL^-1.61)*(LinF^1.44)*kk  
  m.Gislason[m.Gislason>1]=NA
  
  
  #---2. Combine all estimates
  All.methods=data.frame(m.Pauly,m.Hoenig,m.Then.1,m.Jensen.2,Dureuil,Hamel.Cope.tmax,Hamel.Cope.k,
                         Dureuil.at.age,m.PetWro,m.Lorenzen,m.Gislason)
  nat.mort=All.methods%>%dplyr::select(names(weight.method))
  nat.mort.at.age=All.methods%>%dplyr::select(names(weight.method.at.age))

  if(M.averaging=='mean')
  {
    MoRt=apply(nat.mort, 1, function(x) weighted.mean(x,w=weight.method,na.rm=T))
    MoRt.at.age=apply(nat.mort.at.age, 1, function(x) weighted.mean(x,w=weight.method,na.rm=T))
  }
  if(M.averaging=='min')
  {
    MoRt=apply(nat.mort, 1, function(x) min(x,na.rm=T))
    MoRt.at.age=apply(nat.mort.at.age, 1, function(x) min(x,na.rm=T))
  }
  if(MoRt.at.age[1]<MoRt.at.age[2]) MoRt.at.age[1]=1.2*MoRt.at.age[2]  #for analysed species, M[1] is 1.2 times M[2] on average
  
  
  #---3. Output estimates
  return(list(MoRt=MoRt,MoRt.at.age=MoRt.at.age,nat.mort=All.methods,Amax=Amax,tmax=tmax))
}