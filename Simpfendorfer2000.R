###########################################################################
###   WHISKERY SHARK STOCK ASSESSMENT MODEL (SIMPFENDORFER ET AL 2000)  ###
###########################################################################

#notes: . This script performs the whiskery shark stock assessment using the model developed by
#           Simpfendorfer et al 2000.
#       . Age and sex structured model
#       . Model conditioned on catch so TC = Catch
#       . Change in fishing patterns in the mid 1980s: from targeting whiskery to targeting dusky sharks due to
#           decrease in whiskery shark cpue (modelled through 2 different catchabililties). Hence, introduction of 
#           management plan (limited entry, effort controls and gear restrictions).
#       . Management objective (set in 1995):  Bt=0.4Bo by 2010-11     
#       . Regulation by time-gear units
#       . Fishing year starts in June

Simpfendorfer=function(Pth,q.periods)
{
  #--- DATA AND PARAMETERS SECTION ---#
    setPath(Pth)
  
   #inputs
  Catch=(INpts$CATCH.F+INpts$CATCH.M)*1000  #in kg
  CPUE=INpts$CPUE
  Year=yr.start:yr.end  
  age=INpts$age.F
  age.13=match(13,age)    
  n.years=length(Year)
  n.age=length(age)  
  M=INpts$M[1]  
  names(M)=NULL
  sex.ratio=INpts$pup.sx.ratio
  age.mat=age[which(INpts$FEC>0)[1]]  
  mid.Wt.male=ifelse(age<=13,INpts$mid.TwT.M,0)
  mid.Wt.fem=INpts$mid.TwT.F  
  sel.male=ifelse(age<=13,INpts$Selectivity.M,0)
  sel.fem=INpts$Selectivity.F  
  breed.freq=INpts$Breed.freq.min
  fecundity=INpts$FEC*breed.freq
  if(q.periods=="YES") Yr.q.change=INpts$Yr.q.change
  
  
  #--- PROCEDURE SECTION ---#
  
  #1. Define objects to fill in
  N.male=matrix(nrow=n.years,ncol=n.age)
  colnames(N.male)=age
  B.male=N.fem=B.fem=C.male=C.fem=C.male.wt=C.fem.wt=N.male
  C=REC=EGGS=f=B=B.mat=B.exp=vector(length=n.years)
  No.male=No.fem=Per.recruit=N.1975.male=N.1975.fem=vector(length=n.age) 
  
  
  
  #2. Population dynamics
  stock.ass=function(pars)                                    
  {
    PARS=pars
    R.star=PARS[match("R.star",names(PARS))]*1e+6
    z=PARS[match("z",names(PARS))]
    if(q.periods=="YES")
    {
      q.first=PARS[match("q.first",names(PARS))]*1e-6
      q.sec=PARS[match("q.sec",names(PARS))]*1e-6
    }
    if(q.periods=="NO") q.first=PARS[match("q.first",names(PARS))]*1e-6
    
    Fo=PARS[match("Fo",names(PARS))]
    
    
    # 1. Virgin conditions  
    #numbers
    No.male[1]=R.star*sex.ratio     
    for(v in 2:n.age) No.male[v]=No.male[v-1]*exp(-M)
    No.male[age.13]=No.male[age.13]/(1-exp(-M))   #plus group      
    #No.male[age.13]=No.male[age.13]/M   #incorrect parametrisation used in spreadsheet model
    No.male[(age.13+1):n.age]=0
    
    No.fem[1]=R.star*sex.ratio     
    for(v in 2:n.age) No.fem[v]=No.fem[v-1]*exp(-M)
    No.fem[n.age]=No.fem[n.age]/(1-exp(-M))
    #No.fem[n.age]=No.fem[n.age]/M    #incorrect parametrisation used in spreadsheet model 
    
    #eggs
    egg.o=No.fem*fecundity #at age
    Egg.o=sum(egg.o)  #all 
    
    
    # 2. Initial state (this accounts for impact of fishing from 1940 to 1974)
    
    #penalty to avoid Fo<0    
    minF=0.001 #minimum accepted exploitation rate
    pen=0
    pen=sum((Fo<minF)*100*(Fo-minF)^2)
    
    #n per recruit
    Per.recruit[1]=1
    for(p in 2:n.age) Per.recruit[p]=Per.recruit[p-1]*exp(-(M+Fo))
    Per.recruit[n.age]=Per.recruit[n.age]/(1-exp(-M-Fo)) #plus group      
    #Per.recruit[n.age]=Per.recruit[n.age]/(M+Fo)  #incorrect parametrisation used in spreadsheet model     
    
    #eggs per recruit
    egg.per.rec.exploit=Per.recruit*fecundity*sex.ratio
    Egg.per.rec=sum(egg.per.rec.exploit)  #X0 in Simpfendorfer's notation
    
    
    #penalty for z<Z.max (Simpfendorfer 2000)
    z.max=Egg.o/(4*R.star+Egg.o) #maximimum possible z
    pen=pen+sum((z>z.max)*10*(z-z.max)^2)
    
    #penalty for z<Z.min (Simpfendorfer 2000)
    z.min=0.205
    pen=pen+sum((z<z.min)*10*(z-z.min)^2)
    
    
    #initial recruitment
    #note: R.star sets the scale of initial recruitment  
    a=Egg.o/R.star*(1-(z-0.2)/(0.8*z))
    b=(z-0.2)/(0.8*z*R.star)
    Rec.pre.fishing=(Egg.per.rec-a)/(b*Egg.per.rec)   #Rec.pre.fishing = R0 in Simpfendorfer's notation
    
    #numbers 1975 
    N.1975.male[1]=Rec.pre.fishing*sex.ratio
    for(n in 2:n.age) N.1975.male[n]=N.1975.male[n-1]*exp(-(M+Fo))
    N.1975.male[age.13]=N.1975.male[age.13]/(1-exp(-M-Fo)) #plus group      
    #N.1975.male[age.13]=N.1975.male[age.13]/(M+Fo) #incorrect parametrisation used in spreadsheet model      
    N.1975.male[(age.13+1):n.age]=0  
    
    N.1975.fem[1]=Rec.pre.fishing*sex.ratio
    for(n in 2:n.age) N.1975.fem[n]=N.1975.fem[n-1]*exp(-(M+Fo))
    N.1975.fem[n.age]=N.1975.fem[n.age]/(1-exp(-M-Fo))
    #N.1975.fem[n.age]=N.1975.fem[n.age]/(M+Fo)  #incorrect parametrisation used in spreadsheet model 
    
    # store quantities
    Bo=sum(mid.Wt.male*No.male+mid.Wt.fem*No.fem)
    Bo.exp=sum(mid.Wt.male*No.male*sel.male+mid.Wt.fem*No.fem*sel.fem)
    B.1975=sum(mid.Wt.male*N.1975.male+mid.Wt.fem*N.1975.fem)
    Bo.mature=sum(ifelse(age>=age.mat,mid.Wt.fem*No.fem,0))
    
    
    # 3. Post 1974  
    for (i in 1:n.years)
    {
      #fill in numbers
      if(Year[i]==1975) 
      {
        N.fem[i,]=N.1975.fem
        N.male[i,]=N.1975.male
      }else
      {
        for(j in 2:(n.age-1)) N.fem[i,j]=(N.fem[i-1,j-1]-C.fem[i-1,j-1])*exp(-M)      
        for(j in n.age) N.fem[i,j]=(N.fem[i-1,j-1]-C.fem[i-1,j-1]+N.fem[i-1,j]-C.fem[i-1,j])*exp(-M)            
        
        for(j in 2:(n.age-(1+2))) N.male[i,j]=(N.male[i-1,j-1]-C.male[i-1,j-1])*exp(-M)      
        for(j in n.age-2) N.male[i,j]=(N.male[i-1,j-1]-C.male[i-1,j-1]+N.male[i-1,j]-C.male[i-1,j])*exp(-M)
        N.male[i,(age.13+1):n.age]=0    #set males age 14 and 15 to 0
      }
      
      #eggs
      EGGS[i]=sum(N.fem[i,]*fecundity,na.rm=T)
      
      #recruitment
      REC[i]=EGGS[i]/(a+b*EGGS[i])
      N.male[i,1]=(1-sex.ratio)*REC[i]
      N.fem[i,1]=(1-sex.ratio)*REC[i]
      
      #total biomass    
      B[i]=sum(N.fem[i,]*mid.Wt.fem)+sum(N.male[i,]*mid.Wt.male,na.rm=T)
      
      #mature female biomass
      B.mat[i]=sum(ifelse(age>=age.mat,N.fem[i,]*mid.Wt.fem,0))
      
      #exploitable biomass
      B.exp[i]=sum(N.fem[i,]*mid.Wt.fem*sel.fem)+sum(N.male[i,]*mid.Wt.male*sel.male,na.rm=T)
      
      #fishing mortality
      f[i]=ifelse(Catch[i]>=B.exp[i],1,Catch[i]/B.exp[i])   #set f to 1 if catch>biomass
      
      #predicted catch (numbers)
      C.fem[i,]=N.fem[i,]*f[i]*sel.fem
      C.male[i,]=N.male[i,]*f[i]*sel.male
      
      #predicted catch (weight)
      C.fem.wt[i,]=C.fem[i,]*mid.Wt.fem
      C.male.wt[i,]=C.male[i,]*mid.Wt.male
    }
    
    # total catch (weight)
    TC=rowSums(C.fem.wt)+rowSums(C.male.wt,na.rm=T)  
    
    
    # Observation model
        #use different catchabilities according to period
    if(q.periods=="YES") cpue.pred=ifelse(Year<=Yr.q.change,q.first*B.exp,q.sec*B.exp) 
    if(q.periods=="NO") cpue.pred=q.first*B.exp
    
    #penalty to avoid cpue.pred<0    
    min.cpue=0.001 #minimum cpue
    pen=pen+sum((cpue.pred<min.cpue)*.1*(cpue.pred-min.cpue)^2)
    
    
    # Objective function
    residual=log(CPUE+0.000001)-log(cpue.pred+0.000001)
    
    if(q.periods=="YES")
    {
      range.1=match(Year[1]:Yr.q.change,Year)
      range.2=match((Yr.q.change+1):Year[length(Year)],Year) 
      var.1=var(residual[range.1])
      var.2=var(residual[range.2])
      correct.res=ifelse(Year<=Yr.q.change,residual/(var.1)^0.5,residual/(var.2)^0.5)
    }
    if(q.periods=="NO") 
    {
      var.1=var(residual)
      correct.res=residual/(var.1)^0.5
    }
    
    
    #  ssq=sum((correct.res)^2)     #original model
    ssq=sum((residual)^2)         #current, assuming single variance and uncorrected res
    d.f.=n.years-length(pars)-1
    resid.vars=ssq/d.f.
    
    # Penalised neg log Likelihoods
    Log.Like=-(-0.5*(n.years)*(log(2*pi)+log(resid.vars)+1)-100*pen)  #neg log like 
 
    if(is.na(Log.Like)) Log.Like=1e3    #Penalty. Setting very high if getting NA values (e.g.negative biomass)
    
    return=list(Bo=Bo,Bo.exp=Bo.exp,Bo.mature=Bo.mature,B=B,B.MAT=B.mat,B.EXP=B.exp,
                cpue=cpue.pred,TC=TC,LL=Log.Like,residual=residual,std.res=correct.res)
  }
  
    
  #--- MAIN SECTION ---#
  
  #Param estimation through MLE
  metodos=c("Nelder-Mead","BFGS","CG", "L-BFGS-B", "SANN")
  fn_obj=function(pars)stock.ass(pars)$LL  #objfun to MAXIMIZE 
  MLE.1=optim(pars,fn_obj,control=list(fnscale=1,trace=T,maxit=5000),hessian=T) #default is Nelder-Mead algorithm
  MLE=optim(MLE.1$par,method=metodos[2],fn_obj,control=list(fnscale=1,trace=T,maxit=5000),hessian=T)
  
  #Error estimation 
  v_ob=solve(MLE$hessian)  #variance covariance matrix 
  std_ob=sqrt(diag(v_ob))
  CV=100*std_ob/MLE$par
  R_ob=v_ob/(std_ob%o%std_ob)  	#correlation
  
  #Calculate dynamics using MLE
  est.dynamics=stock.ass(MLE$par)

  
  
  #--- REPORT SECTION ---#  
  
  #Goodness of fit (e.g. page 120 O'Neill et al 2005)
  tiff(file="goodness.of.fit.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
  
  par(mfcol=c(2,2),las=1,mai=c(0.7,0.7,0.01,0.01), omi=c(0.10,0.1,0.10,.05),mgp=c(2,0.7,0))
  
  #predicted and observed cpue
  plot(Year,CPUE,col="red",ylab="cpue",xlab="year",ylim=c(0,max(c(CPUE,est.dynamics$cpue))),pch=19)
  lines(Year,est.dynamics$cpue,col="blue",lwd=2)
  #points(Year,est.dynamics$cpue,col="blue",pch=19,cex=0.75)
  legend("topright",c("observed", "predicted"),lty=c(0,1),bty="n",col=c("red","blue"),cex=1,pch=19)
  #title('Observed vs Predicted cpue',line = -0.8,cex.main =1)
  
  #standardised residuals
  #histogram of standardised residuals
  hist(est.dynamics$std.res,breaks=25, xlim=c(-3,3),col="gray", ylab="Frequency", 
       xlab="Standardised residuals",main="")
  title("Standardised residuals",line = -0.8,cex.main =1.05)
  box()
  
  #predicted cpue vs standardised residuals
  plot(est.dynamics$cpue,est.dynamics$std.res, xlab="Predicted cpue", ylab="Standardised residuals",
       ylim=c(-3,3))
  abline(0,0)
  title("Fitted values",line = -0.8,cex.main =1.05)
  
  #quantile-quantile plot of standardised residuals
  qqnorm(est.dynamics$std.res,main="",xlab="Theoretical quantiles",ylab="Sample quantiles")
  qqline(est.dynamics$std.res, col = 2)
  title("Normality plot",line = -0.8,cex.main =1.05)
  
  dev.off()
  
  return(est.dynamics)
}