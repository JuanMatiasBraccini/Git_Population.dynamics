#Notes:
#   see Kendal et al 2019 for common issues with Leslie Matrices
#   Leslie matrices project the population from one nominal census date to the next.
#   Census occurs just before breeding (“prebreeding census”) e.g. Cortes 2002
#              or just after breeding (“postbreeding census”) e.g. Aires da Silva 2007
#   birth-pulse populations= short breeding season (i.e. the season in which individuals
#                                 are born or hatched)
#   Age starts counting from birth, so that a reproductively
#       mature individual breeds on or about their birthday: e.g., if the age at
#       first reproduction is 5 years, then an individual has its first offspring on
#       its fifth birthday



library(popbio)     #for solving matrices
library(EnvStats)
if("package:VGAM" %in% search()) detach("package:VGAM", unload=TRUE)
library(triangle)
library(MASS)   #for sampling from multivariate distribution

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/Natural.mortality.R"))

fun.Leslie=function(N.sims,k,Linf,k.sd,Linf.sd,k.Linf.cor,A,first.age,RangeMat,Rangefec,
                    sexratio,Reprod_cycle,bwt,awt,Lo,Resamp)
{
  #univariate distributions
  fn.draw.samples=function()
  {
    #Max Age
    if(length(A)==1) Max.A=A
    if(length(A)>1)if(A[1]==A[2]) Max.A=A[1]
    if(length(A)>1)if(A[1]<A[2]) Max.A=floor(rtriangle(1,a=A[1],b=A[2],c=A[1]))
    
    #Age vector
    age=first.age:Max.A
    
    #Age at 50% maturity
    if(RangeMat[1]==RangeMat[2]) age.mat.sim=floor(RangeMat[1])
    if(RangeMat[1]<RangeMat[2]) age.mat.sim=floor(runif(1,RangeMat[1],RangeMat[2]))  
    
    #fecundity at age
      #single value
    if(Rangefec[1]==Rangefec[2]) Meanfec.sim=rep(Rangefec[1],length(age))
      #triangular distribution
    if(!linear.fec=="YES")
    {
      if(Rangefec[1]<Rangefec[2]) Meanfec.sim=rep(ceiling(rtriangle(1,a=Rangefec[1],b=Rangefec[2],
                                                                    c=ceiling((Rangefec[1]+Rangefec[2])/2))),length(age))
    }
      #linear increase
    if(linear.fec=="YES")
    {
      Meanfec.sim=c(rep(0,age.mat.sim),
                    ceiling(seq(Rangefec[1],Rangefec[2],length.out=1+Max.A-age.mat.sim)))
    }

    
    #Reproductive cycle
    if(length(Reprod_cycle)==1) Rep_cycle.sim=Reprod_cycle else
    {
      if(Reprod_cycle[1]==Reprod_cycle[2]) Rep_cycle.sim=round(Reprod_cycle[1])
      if(Reprod_cycle[1]<Reprod_cycle[2]) Rep_cycle.sim=sample(Reprod_cycle,1,replace=T)    
    }
    
    return(list(Max.A=Max.A,age.mat=age.mat.sim,Meanfec=Meanfec.sim,Rep_cycle=Rep_cycle.sim))    
  }
  
  #multivariate distributions
  fun.multivar=function(mu,stddev,corMat,N)
  {
    covMat <- stddev %*% t(stddev) * corMat
    return(mvrnorm(n = N, mu = mu, Sigma = covMat, empirical = TRUE))
  }
  
  #Leslie
  Leslie=function(M,age.mat,Meanfec,CyclE)
  {  
    #survivorship
    S=exp(-M)         
    
    #proportion surviving
    lx=rep(NA,length(age))
    lx[1]=1.0
    for (i in 2:(length(age)))lx[i]=lx[i-1]*S[i]
    
    #reproductive schedules 
      #knife edge
    MF=c(rep(0,(age.mat-1)),Meanfec[age.mat:length(Meanfec)])
    MF[age.mat]=MF[age.mat]*.5   #because age.mat is actually 50% maturity
      #ogive
    #MF=plogis(age,age.mat,1)
    
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
    if(first.age==0) from=first.age+1
    if(first.age==1) from=first.age
    rownames(Data)=colnames(Data)=from:n
    
    #solve projection matrix
    LAMBDA=lambda(Data)
    r=log(LAMBDA)  
    t2=log(2)/r # doubling time
    G=sum(age*lx*mx*exp(-(r*age))) # generation time Cortes 2002
    
    return(list(r=r,t2=t2,lambda=LAMBDA,M=M,G=G))  
  }
  
  #Monte Carlo simulations
  Growth.sim=fun.multivar(mu <- c(Linf,k),
                          stddev <- c(Linf.sd,k.sd),
                          corMat <- matrix(c(1, k.Linf.cor,k.Linf.cor, 1),ncol = 2),
                          N=N.sims)
  Store=vector('list',length(N.sims))
  Input.pars=Store
  nat.mort.sim=Store
  for(i in 1:N.sims)
  {
    a=fn.draw.samples()
    A.sim=a$Max.A
    age=first.age:A.sim
    Age.mat.sim=a$age.mat
    Meanfec.sim=a$Meanfec
    Reprod_cycle.sim=a$Rep_cycle
    Linf.sim=Growth.sim[i,1]
    k.sim=Growth.sim[i,2]
    M.sim=M.fun(AGE=age,Amax=A.sim,age.mat=Age.mat.sim,LinF=Linf.sim,
                kk=k.sim,awt=awt,bwt=bwt,Lo=Lo)   
    if(GET.all.Ms) nat.mort.sim[[i]]=M.sim$nat.mort
    if(what.M=='age.invariant') M.sim=M.sim$MoRt
    if(what.M=='at.age') M.sim=M.sim$MoRt.at.age
    Store[[i]]=Leslie(M=M.sim,age.mat=Age.mat.sim,Meanfec=Meanfec.sim,CyclE=Reprod_cycle.sim)
    rr=Store[[i]]$r
    
      #avoid negative r
    if(Resamp=="YES")
    {
      if(rr<=0)repeat 
      {
        a=fn.draw.samples()
        A.sim=a$Max.A
        age=first.age:A.sim
        Age.mat.sim=a$age.mat
        Meanfec.sim=a$Meanfec
        Reprod_cycle.sim=a$Rep_cycle
        M.sim=M.fun(AGE=age,Amax=A.sim,age.mat=Age.mat.sim,LinF=Linf.sim,kk=k.sim,awt=awt,bwt=bwt,Lo=Lo)
        nat.mort.sim[[i]]=M.sim$nat.mort
        if(what.M=='age.invariant') M.sim=M.sim$MoRt
        if(what.M=='at.age') M.sim=M.sim$MoRt.at.age
        Store[[i]]=Leslie(M=M.sim,age.mat=Age.mat.sim,Meanfec=Meanfec.sim,CyclE=Reprod_cycle.sim)
        rr=Store[[i]]$r
        if(rr>0)break
      }
    }
    
    Input.pars[[i]]=list(Max.age=A.sim,
                         Age.mat=Age.mat.sim,
                         Meanfec=mean(Meanfec.sim),
                         Reprod_cycle=Reprod_cycle.sim,
                         Linf=Linf.sim,
                         k=k.sim)
  }
  
  r.prior=do.call("c", lapply(Store, "[[", 1))
  M.all=do.call("list", lapply(Store, "[[", 4))
  G.all=do.call("c", lapply(Store, "[[", 5))
  
  if(GET.all.Ms) 
  {
   for(e in 1:length(nat.mort.sim))
    {
     if(first.age==0) to=nrow(nat.mort.sim[[e]])-1
     if(first.age==1) to=nrow(nat.mort.sim[[e]])
     nat.mort.sim[[e]]$Age=first.age:to  
    }
    nat.mort.sim=plyr::compact(nat.mort.sim)
    nat.mort.sim=do.call(rbind,nat.mort.sim)%>%
      tidyr::gather(key = "Method", value = "M", -Age)%>%
      dplyr::mutate(Age=factor(Age,levels=min(Age):max(Age)))%>%
      dplyr::group_by(Age,Method)%>%
      dplyr::summarise(M.mean=mean(M,na.rm=T),
                M.sd=sd(M,na.rm=T))%>%
      dplyr::ungroup()%>%
      dplyr::mutate(M.CV=M.sd/M.mean)
    
  }
  if(GET.all.Ms)
  {
    return(list(r.prior=r.prior,M=M.all,Input.pars=Input.pars,G=G.all,nat.mort.sim=nat.mort.sim))
  }else
  {
    return(list(r.prior=r.prior,M=M.all,Input.pars=Input.pars,G=G.all))
  }
}
fun.rprior.dist=function(Nsims,K,LINF,K.sd,LINF.sd,k.Linf.cor,Amax,MAT,FecunditY,
                         Cycle,BWT,AWT,LO)
{
  Fecu=unlist(FecunditY)
  Rprior=fun.Leslie(N.sims=Nsims,k=K,Linf=LINF,k.sd=K.sd,Linf.sd=LINF.sd,k.Linf.cor=k.Linf.cor,
                    A=Amax,first.age=First.Age,RangeMat=MAT,Rangefec=Fecu,
                    sexratio=0.5,Reprod_cycle=Cycle,
                    bwt=BWT,awt=AWT,Lo=LO,
                    Resamp=RESAMP)  
  
  #get mean and sd from gamma and normal distribution
  normal.pars=suppressWarnings(fitdistr(Rprior$r.prior, "normal"))
  gamma.pars=suppressWarnings(fitdistr(Rprior$r.prior, "gamma"))  
  shape=gamma.pars$estimate[1]        
  rate=gamma.pars$estimate[2]      
  return(list(shape=shape,rate=rate,
              mean=normal.pars$estimate[1],sd=normal.pars$estimate[2],
              M=Rprior$M,
              G=Rprior$G,
              Input.pars=Rprior$Input.pars,
              nat.mort.sim=Rprior$nat.mort.sim))
  
  #get mean and sd from lognormal distribution
  #LogN.pars=fitdistr(Rprior, "lognormal")  
  #log_mean.r=LogN.pars$estimate[1]    #already in log space     
  #log_sd.r=LogN.pars$estimate[2]      #already in log space     
  #return(list(log_mean.r=log_mean.r,log_sd.r=log_sd.r))
}