if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/Natural.mortality.R"))


fun.steepness=function(Nsims,K,LINF,Linf.sd,k.sd,first.age,sel.age,F.mult,Amax,MAT,
                       FecunditY,Cycle,sexratio,spawn.time,AWT,BWT,LO,Resamp,simsout=1e3)
{
  #samples from univariate distribution
  fn.draw.samples=function(A=Amax,RangeMat=MAT,Rangefec=Fecu,Reprod_cycle=Cycle)
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
      if(Reprod_cycle[1]<Reprod_cycle[2]) Rep_cycle.sim=round(runif(1,Reprod_cycle[1],Reprod_cycle[2]))    
    }
    
    return(list(Max.A=Max.A,age.mat=age.mat.sim,Meanfec=Meanfec.sim,Rep_cycle=Rep_cycle.sim))    
  }
  
  #samples from multivariate distribution
  fun.multivar=function(mu,stddev,corMat,N)
  {
    covMat <- stddev %*% t(stddev) * corMat
    return(mvrnorm(n = N, mu = mu, Sigma = covMat, empirical = TRUE))
  }
  
  #Steepness calculation
  Stipns=function(max.age,M,age.mat,Meanfec,CyclE,Sel)
  {  
    #survivorship
    surv=exp(-M)
    
    #fecundity  
    fecundity=Meanfec*sexratio/CyclE
    
    #maturity
    #knife edge
    maturity=ifelse(age>=age.mat,1,0)   
    maturity[which(age==age.mat)]=0.5   #age.mat is actually 50% maturity
    #ogive
    #maturity=plogis(age,age.mat,1)      
    
    # maximum age is plus group
    phi.o=0.0
    cum.survive=1.0
    z=0.0
    for (i in 2:max.age)
    {
      z=M[i] + F.mult*Sel[i]
      z.ts=(M[i]+F.mult*Sel[i])*spawn.time
      phi.o=phi.o+cum.survive*fecundity[i]*maturity[i]*exp(-z.ts)
      cum.survive=cum.survive*exp(-z )
    }
    #plus group
    if(first.age==0) to=max.age+1
    if(first.age==1) to=max.age
    z= M[to] + F.mult*Sel[to]
    z.ts=(M[to]+F.mult*Sel[to])*spawn.time
    phi.o=phi.o + fecundity[to]*maturity[to]*cum.survive*exp(-z.ts)/( 1- exp(-z ) )
    
    #maximum lifetime reproductive rate at low density
    alpha=phi.o*surv[1]
    
    #steepness
    h=alpha/(4+alpha)
    
    #spawning potential ratio at maximum excess recruitment (MER) (Beverton-Holt relationship)
    SPR.mer=1/alpha^0.5
    
    #optimal depletionlevel (i.e.depletion at MER, the proportional reduction from unexploited level)
    Dep.MER=((alpha^0.5)-1)/(alpha-1) 
    
    return(list(steepness=h,alpha=alpha))  
  }
  
  Fecu=unlist(FecunditY)
  
  
  #Monte Carlo simulations
  Growth.sim=fun.multivar(mu <- c(LINF,K),
                          stddev <- c(Linf.sd,k.sd),
                          corMat <- matrix(c(1, k.Linf.cor,k.Linf.cor, 1),ncol = 2),
                          N=Nsims)
  Store=rep(NA,Nsims)
  Store.alfa=Store
  M.all=dummies=vector('list',Nsims)
  for(i in 1:Nsims)
  {
    a=fn.draw.samples()
    A.sim=a$Max.A
    age=first.age:A.sim
    Age.mat.sim=a$age.mat
    Meanfec.sim=a$Meanfec
    Reprod_cycle.sim=a$Rep_cycle
    Linf.sim=Growth.sim[i,1]
    k.sim=Growth.sim[i,2]
    M.sim=M.fun(AGE=age,Amax=A.sim,age.mat=Age.mat.sim,LinF=Linf.sim,kk=k.sim,awt=AWT,bwt=BWT,Lo=LO)
    
    if(what.M=='age.invariant') M.sim=M.sim$MoRt
    if(what.M=='at.age') M.sim=M.sim$MoRt.at.age
    
    Sel=sel.age
    if(length(Sel)<length(age)) Sel=c(Sel,rep(Sel[length(Sel)],length(age)-length(Sel)))
    if(length(Sel)>length(age)) Sel=Sel[1:length(age)]
      
    hh=Stipns(max.age=A.sim,M=M.sim,age.mat=Age.mat.sim,
              Meanfec=Meanfec.sim,CyclE=Reprod_cycle.sim,
              Sel)$steepness
    Alfa=Stipns(max.age=A.sim,M=M.sim,age.mat=Age.mat.sim,
              Meanfec=Meanfec.sim,CyclE=Reprod_cycle.sim,
              Sel)$alpha
    
    #avoid non-sense h
    if(Resamp=="YES")
    {
      if(hh<=0.20)repeat 
      {
        a=fn.draw.samples()
        A.sim=a$Max.A
        age=first.age:A.sim
        Age.mat.sim=a$age.mat
        Meanfec.sim=a$Meanfec
        Reprod_cycle.sim=a$Rep_cycle
        M.sim=M.fun(AGE=age,Amax=A.sim,age.mat=Age.mat.sim,LinF=Linf.sim,kk=k.sim,awt=AWT,bwt=BWT,Lo=LO)
        
        if(what.M=='age.invariant') M.sim=M.sim$MoRt
        if(what.M=='at.age') M.sim=M.sim$MoRt.at.age
        
        Sel=sel.age
        if(length(Sel)<length(age)) Sel=c(Sel,rep(Sel[length(Sel)],length(age)-length(Sel)))
        if(length(Sel)>length(age)) Sel=Sel[1:length(age)]
        
        hh=Stipns(max.age=A.sim,M=M.sim,age.mat=Age.mat.sim,
                  Meanfec=Meanfec.sim,CyclE=Reprod_cycle.sim,
                  Sel)$steepness
        Alfa=Stipns(max.age=A.sim,M=M.sim,age.mat=Age.mat.sim,
                    Meanfec=Meanfec.sim,CyclE=Reprod_cycle.sim,
                    Sel)$alpha
        if(hh>0.2)break
      } 
    }

    Store[i]=hh
    Store.alfa[i]=Alfa
    M.all[[i]]=M.sim
    dummies[[i]]=data.frame(Max.age=A.sim,Fecundity=mean(Meanfec.sim),Rep.cycle=Reprod_cycle.sim,
                            Age.mat=Age.mat.sim,Linf=Linf.sim,k=k.sim,M=mean(M.sim),h=hh)
  }
  if(!Resamp=="YES")  Store=Store[Store>0.2]
 
  dummies=do.call(rbind,dummies)
  dummies=subset(dummies,h>0.2)
  dummies=dummies[sample(1:nrow(dummies),simsout,replace = T),]
    
  #get mean and sd from lognormal distribution
  # normal.pars=suppressWarnings(fitdistr(Store, "normal"))
  # gamma.pars=suppressWarnings(fitdistr(Store, "gamma"))  
  # shape=gamma.pars$estimate[1]        
  # rate=gamma.pars$estimate[2]  
  
  return(list(mean=mean(dummies$h),
              sd=sd(dummies$h),
              M=M.all,
              Alpha=Store.alfa,
              Runs=dummies))
  # return(list(shape=shape,rate=rate,
  #             mean=normal.pars$estimate[1],
  #             sd=normal.pars$estimate[2],
  #             M=M.all,
  #             Alpha=Store.alfa,
  #             Runs=dummies))
}

Alpha.Brooks=function(max.age,M,age.mat,Meanfec,CyclE,sexratio=0.5,spawn.time=0,first.age=First.Age)
{  
  age=first.age:max.age
  
  #survivorship
  M=M[1:length(age)]
  surv=exp(-M)
  
  #fecundity  
  fecundity=rep(Meanfec*sexratio/CyclE,length(age))
  
  #maturity
  #knife edge
  maturity=ifelse(age>=age.mat,1,0)   
  maturity[which(age==age.mat)]=0.5   #age.mat is actually 50% maturity
  #ogive
  #maturity=plogis(age,age.mat,1)      
  
  # maximum age is plus group
  phi.o=0.0
  cum.survive=1.0
  z=0.0
  for (i in 2:(max.age)  )
  {
    z=M[i] 
    z.ts=M[i]*spawn.time
    phi.o=phi.o+cum.survive*fecundity[i]*maturity[i]*exp(-z.ts)
    cum.survive=cum.survive*exp(-z )
  }
  #plus group
  if(first.age==0) to=max.age+1
  if(first.age==1) to=max.age
  z= M[to] 
  z.ts=M[to]*spawn.time
  phi.o=phi.o + fecundity[to]*maturity[to]*cum.survive*exp(-z.ts)/( 1- exp(-z ) )
  
  #maximum lifetime reproductive rate at low density
  alpha=phi.o*surv[1]
  
  return(alpha)  
}
