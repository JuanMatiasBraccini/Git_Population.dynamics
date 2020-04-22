library(popbio)     #for solving matrices
library(EnvStats)
if("package:VGAM" %in% search()) detach("package:VGAM", unload=TRUE)
library(triangle)
fun.Leslie=function(N.sims,k,Linf,Aver.T,A,first.age,RangeMat,Rangefec,sexratio,Reprod_cycle,bwt,awt,Lo)
{
  fn.draw.samples=function()
  {
    #Max Age
    if(length(A)==1) Max.A=A
    if(length(A)>1)if(A[1]==A[2]) Max.A=A[1]
    if(length(A)>1)if(A[1]<A[2]) Max.A=ceiling(rtriangle(1,a=A[1],b=A[2],c=ceiling((A[1]+A[2])/2)))
    
    #Age vector
    age=first.age:Max.A
    
    #fecundity at age
    if(Rangefec[1]==Rangefec[2]) Meanfec.sim=rep(Rangefec[1],length(age))
    if(Rangefec[1]<Rangefec[2]) Meanfec.sim=rep(ceiling(rtriangle(1,a=Rangefec[1],b=Rangefec[2],
                                                            c=ceiling((Rangefec[1]+Rangefec[2])/2))),length(age)) 
    
    #Age at 50% maturity
    if(RangeMat[1]==RangeMat[2]) age.mat.sim=ceiling(RangeMat[1])
    if(RangeMat[1]<RangeMat[2]) age.mat.sim=ceiling(runif(1,RangeMat[1],RangeMat[2]))  
    
    #Reproductive cycle
    if(length(Reprod_cycle)==1) Rep_cycle.sim=Reprod_cycle else
    {
      if(Reprod_cycle[1]==Reprod_cycle[2]) Rep_cycle.sim=round(Reprod_cycle[1])
      if(Reprod_cycle[1]<Reprod_cycle[2]) Rep_cycle.sim=round(runif(1,Reprod_cycle[1],Reprod_cycle[2]))    
    }
    
    return(list(Max.A=Max.A,age.mat=age.mat.sim,Meanfec=Meanfec.sim,Rep_cycle=Rep_cycle.sim))    
  }

  M.fun=function(Amax,age.mat)
  {
    #STEP 1. calculate M from different methods (see Kenchington 2013)
    
      #1.1. Age-independent
        #Jensen (1996)
    m.Jensen.2=1.65/age.mat
    m.Jensen.2=rep(m.Jensen.2,length(age))
    
        #Pauly (1980)  
    m.Pauly=10^(-0.0066-0.279*log10(Linf)+0.6543*log10(k)+0.4634*log10(Aver.T))
    m.Pauly=rep(m.Pauly,length(age))
    
        #Hoenig (1983), combined teleost and cetaceans    
    m.Hoenig=exp(1.44-0.982*log(Amax))      
    m.Hoenig=rep(m.Hoenig,length(age))
    
        #Then et al (2015)
    m.Then.1=4.899*Amax^(-0.916)
    m.Then.1=rep(m.Then.1,length(age))
    
      #1.2. Age-dependent
        #Peterson and Wroblewski 1984 (dry weight in grams, length in cm)
    Dry.w=0.2   # Cortes (2002)
    TL=Lo+(Linf-Lo)*(1-exp(-k*age))
    wet.weight=1000*awt*TL^bwt
    m.PetWro=1.92*(wet.weight*Dry.w)^-0.25
    m.PetWro[m.PetWro>1]=NA
    
        #Lorenzen 1996 (weight in grams)
    m.Lorenzen=3*wet.weight^-0.288
    m.Lorenzen[m.Lorenzen>1]=NA
    
        #Gislason et al (2010) (weight in grams, length in cm)
    m.Gislason=1.73*(TL^-1.61)*(Linf^1.44)*k
    m.Gislason[m.Gislason>1]=NA
    
    
    #STEP 2. get mean at age
    nat.mort=data.frame(m.Jensen.2,m.Pauly,m.Hoenig,m.Then.1,
                        m.PetWro,m.Lorenzen,m.Gislason)  
    #for dogfish, due to their small size, weight-based M estimators highly overestimate M
    if(mean(m.PetWro/m.Gislason,na.rm=T)>5)  nat.mort=data.frame(m.Jensen.2,m.Pauly,m.Hoenig,m.Then.1,
                                                                   m.Gislason)  
    
    #STEP 3. Calculate mean
    MoRt=rowMeans(nat.mort,na.rm=T)
    #MoRt=apply(nat.mort, 1, function(x) weighted.mean(x, c(1,1,1.5,1.5,1,1,1)))
    
    if(MoRt[1]<MoRt[2]) MoRt[1]=1.2*MoRt[2]  #for analysed species, M[1] is 1.2 times M[2] on average
    
    return(MoRt)
  }
  
  Leslie=function(M,age.mat,Meanfec,CyclE)
  {  
    #survivorship
    S=exp(-M)         
    
    #proportion surviving
    lx=rep(NA,length(age))
    lx[1]=1.0
    for (i in 2:(length(age)))lx[i]=lx[i-1]*S[i]
    
    #reproductive schedules   
    MF=c(rep(0,(age.mat-1)),Meanfec[age.mat:length(Meanfec)])
    mx=MF*sexratio/CyclE
    
    #probability of surviving (for birth-pulse, post-breeding census)
    px=vector(length=length(lx))
    for(i in 2:length(lx)) px[i-1]=(lx[i])/(lx[i-1])
    
    #fertility  (for birth-pulse, post-breeding census)
    bx=mx*px
    
    #projection matrix
    PX=px
    PX=PX[-length(PX)]
    BX=bx
    n=length(BX)
    Data=matrix(0,nrow=n,ncol=n)
    diag(Data)[-nrow(Data)]=PX
    Data=rbind(matrix(0,nrow=1,ncol=n),Data)
    Data=Data[-(n+1),]
    Data[1,]=BX
    rownames(Data)=colnames(Data)=(first.age+1):n
    
    #solve projection matrix
    LAMBDA=lambda(Data)
    r=log(LAMBDA)  
    t2=log(2)/r 
    
    return(list(r=r,t2=t2,lambda=LAMBDA,M=M))  
  }
  
  Store=vector('list',length(N.sims))
  for(i in 1:N.sims)
  {
    a=fn.draw.samples()
    A.sim=a$Max.A
    age=first.age:A.sim
    Age.mat.sim=a$age.mat
    Meanfec.sim=a$Meanfec
    Reprod_cycle.sim=a$Rep_cycle
    M.sim=M.fun(Amax=A.sim,age.mat=Age.mat.sim)
    Store[[i]]=Leslie(M=M.sim,age.mat=Age.mat.sim,Meanfec=Meanfec.sim,CyclE=Reprod_cycle.sim)
    rr=Store[[i]]$r
    
    #avoid negative r
    if(rr<=0)repeat 
    {
      a=fn.draw.samples()
      A.sim=a$Max.A
      age=first.age:A.sim
      Age.mat.sim=a$age.mat
      Meanfec.sim=a$Meanfec
      Reprod_cycle.sim=a$Rep_cycle
      M.sim=M.fun(Amax=A.sim,age.mat=Age.mat.sim)
      Store[[i]]=Leslie(M=M.sim,age.mat=Age.mat.sim,Meanfec=Meanfec.sim,CyclE=Reprod_cycle.sim)
      rr=Store[[i]]$r
      if(rr>0)break
    }
    
  }
  
  r.prior=do.call("c", lapply(Store, "[[", 1))
  M.all=do.call("list", lapply(Store, "[[", 4))
  return(list(r.prior=r.prior,M=M.all))
}