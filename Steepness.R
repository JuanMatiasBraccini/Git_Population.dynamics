source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/Natural.mortality.R")

fun.steepness=function(Nsims,K,LINF,Linf.sd,k.sd,first.age,sel.age,F.mult,Amax,MAT,
                       FecunditY,Cycle,sexratio,spawn.time,AWT,BWT,LO,Resamp)
{
  Fecu=unlist(FecunditY)
  
  #samples from univariate distribution
  fn.draw.samples=function(A=Amax,RangeMat=MAT,Rangefec=Fecu,Reprod_cycle=Cycle)
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
    maturity=ifelse(age>=age.mat,1,0)   #knife edge
    maturity[which(age==age.mat)]=0.5   #age.mat is actually 50% maturity
    #maturity=plogis(age,age.mat,1)      #ogive
    
    # maximum age is plus group
    phi.o=0.0
    cum.survive=1.0
    z=0.0
    for (i in 2:(max.age)  )
    {
      z=M[i] + F.mult*Sel[i]
      z.ts=(M[i]+F.mult*Sel[i])*spawn.time
      phi.o=phi.o+cum.survive*fecundity[i]*maturity[i]*exp(-z.ts)
      cum.survive=cum.survive*exp(-z )
    }
    #plus group  
    z= M[max.age+1] + F.mult*Sel[max.age+1]
    z.ts=(M[max.age+1]+F.mult*Sel[max.age+1])*spawn.time
    phi.o=phi.o + fecundity[max.age+1]*maturity[max.age+1]*cum.survive*exp(-z.ts)/( 1- exp(-z ) )
    
    #maximum lifetime reproductive rate at low density
    alpha=phi.o*surv[1]
    
    #steepness
    h=alpha/(4+alpha)
    
    #spawning potential ratio at maximum excess recruitment (MER) (Beverton-Holt relationship)
    SPR.mer=1/alpha^0.5
    
    #optimal depletionlevel (i.e.depletion at MER, the proportional reduction from unexploited level)
    Dep.MER=((alpha^0.5)-1)/(alpha-1) 
    
    return(steepness=h)  
  }
  
  #Monte Carlo simulations
  Growth.sim=fun.multivar(mu <- c(LINF,K),
                          stddev <- c(Linf.sd,k.sd),
                          corMat <- matrix(c(1, k.Linf.cor,k.Linf.cor, 1),ncol = 2),
                          N=Nsims)
  Store=rep(NA,Nsims)
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
    
    Sel=sel.age
    if(length(Sel)<length(age))
    {
      Sel=c(Sel,rep(Sel[length(Sel)],length(age)-length(Sel)))
    }
    hh=Stipns(max.age=A.sim,M=M.sim,age.mat=Age.mat.sim,
              Meanfec=Meanfec.sim,CyclE=Reprod_cycle.sim,
              Sel)
    #avoid non-sense h
    if(Resamp=="YES")
    {
      if(hh<0.20)repeat 
      {
        a=fn.draw.samples()
        A.sim=a$Max.A
        age=first.age:A.sim
        Age.mat.sim=a$age.mat
        Meanfec.sim=a$Meanfec
        Reprod_cycle.sim=a$Rep_cycle
        M.sim=M.fun(AGE=age,Amax=A.sim,age.mat=Age.mat.sim,LinF=Linf.sim,kk=k.sim,awt=AWT,bwt=BWT,Lo=LO)
        Sel=sel.age
        if(length(Sel)<length(age))
        {
          Sel=c(Sel,rep(Sel[length(Sel)],length(age)-length(Sel)))
        }
        hh=Stipns(max.age=A.sim,M=M.sim,age.mat=Age.mat.sim,
                  Meanfec=Meanfec.sim,CyclE=Reprod_cycle.sim,
                  Sel)
        if(hh>0.2)break
      } 
    }

    Store[i]=hh
  }
  if(!Resamp=="YES") Store=Store[Store>=0.2]
  
  #get mean and sd from lognormal distribution
  normal.pars=suppressWarnings(fitdistr(Store, "normal"))
  gamma.pars=suppressWarnings(fitdistr(Store, "gamma"))  
  shape=gamma.pars$estimate[1]        
  rate=gamma.pars$estimate[2]  
  
  return(list(shape=shape,rate=rate,
              mean=normal.pars$estimate[1],
              sd=normal.pars$estimate[2]))
}