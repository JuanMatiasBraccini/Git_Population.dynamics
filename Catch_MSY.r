##################### Catch-based MSY estimation #############################
fn.fig=function(NAME,Width,Height)
{
  if(Do.tiff=="YES") tiff(file=paste(NAME,".tiff",sep=""),width=Width,height=Height,units="px",res=300,compression="lzw")
  if(Do.jpeg=="YES") jpeg(file=paste(NAME,".jpeg",sep=""),width=Width,height=Height,units="px",res=300)
}
fn.frmt=function(x) format(x,digits=3)

#note: implementation of Martell & Froese 2012.
#       use total catches (i.e. all sources of F)
library(future.apply)
plan(multiprocess)   #for parallel processing

Catch_MSY=function(ct,yr,r.prior,user,k.max,startbio,finalbio,res,n,sigR,ct.future,yr.future)
{
  set.seed(999)  ## for same random sequence
  
  nyr  <- length(yr)    ## number of years in the time series
  
  #-r prior
  #user defined pars
  if(user=="Yes") 
  {
    start_r=r.prior
  }else
    #defaults
  {
    start_r  <- if(res == "Very low")            
    {c(0.015, 0.1)}else if
    (res == "Low")
    {c(0.05,0.5)}else if
    (res == "High")
    {c(0.6,1.5)}else
    {c(0.2,1)} 
  }
  
  #-range of possible K values
  start_k=c(max(ct),k.max*max(ct))
  
  #- boundaries of initial and final depletion ranges if not user-defined (as fraction of k) 
  if(is.na(startbio[1]))
  {
    #biomass range at start of time series
    #note: this assumes C trend is Abundance trend proxy, it ignores targeting changes, etc
    B.i.def.fn=function(DAT)if(DAT[1]/max(DAT)<0.5)B.i=c(.5,.9)else{B.i=c(.3,.6)}
    startbio=B.i.def.fn(ct)  
  }
  if(is.na(finalbio[1]))
  {
    #biomass range after last catch 
    B.f.def.fn=function(DAT)if(DAT[length(DAT)]/max(DAT)>0.5)B.f=c(.3,.7)else{B.f=c(.01,.4)}
    finalbio=B.f.def.fn(ct)
  }
  
  
  #---PROCEDURE SECTION-----
  
  ## FUNCTIONS
  #Surplus production model for calculating biomass
  .schaefer  <- function(theta)
  {
    with(as.list(theta), {  ## for all combinations of ri & ki
      bt=vector()
      ell = 0  ## initialize ell
      for (j in startbt)
      {
        if(ell == 0) 
        {
          bt[1]=j*k*exp(rnorm(1,0, sigR))  ## set biomass in first year
          for(i in 1:nyr) ## for all years in the time series
          {
            xt=rnorm(1,0, sigR)
            bt[i+1]=(bt[i]+r*bt[i]*(1-bt[i]/k)-ct[i])*exp(xt) ## calculate biomass as function of previous year's biomass plus net production minus catch
          }
          
          #Bernoulli likelihood, assign 0 or 1 to each combination of r and k
          ell = 0
          if(bt[nyr+1]/k>=lam1 && bt[nyr+1]/k <=lam2 && min(bt) > 0 && max(bt) <=k && bt[which(yr==interyr)]/k>=interbio[1] && bt[which(yr==interyr)]/k<=interbio[2]) 
            ell = 1
        }	
      }
      return(list(ell=ell,bt=bt,bt.rel=bt/k,KTCH=ct))
    })
  }
  
  #This function conducts the stock reduction analysis for N trials
  #args:
  #  theta - a list object containing:
  #		r (lower and upper bounds for r)
  #		k (lower and upper bounds for k)
  #		lambda (limits for current depletion)
  sraMSY	<-function(theta, N)
  {
    with(as.list(theta), 
         {
           ## get N values of r, assign to ri
           #lognormal r  
           if(user=="Yes")
           {
             ri = rgamma(N, r[1], r[2])  #gamma
             #ri = rlnorm(N, r[1], r[2])  #lognormal
           } else ri = exp(runif(N, log(r[1]), log(r[2])))  #uniform r  
           
           ## get N values between k[1] and k[2], assing to ki
           ki = exp(runif(N, log(k[1]), log(k[2])))  
           
           ## assign ri, ki, and final biomass range to itheta
           itheta=cbind(r=ri,k=ki, lam1=lambda[1],lam2=lambda[2], sigR=sigR) 
           
           ## call Schaefer function with parameters in itheta
           M = future_apply(itheta,1,.schaefer) #parallel processing
           i=1:N
           ## prototype objective function
           get.ell=function(i) M[[i]]$ell
           ell = sapply(i, get.ell) 
           
           get.bt=function(i) M[[i]]$bt
           bt = sapply(i, get.bt) 
           
           get.bt.rel=function(i) M[[i]]$bt.rel
           bt.rel = sapply(i, get.bt.rel) 
           
           
           get.ct=function(i) M[[i]]$KTCH
           KTCH = sapply(i, get.ct) 
           
           U=KTCH/bt[1:nrow(KTCH),]
           
           return(list(r=ri,k=ki, ell=ell,bt=bt,bt.rel=bt.rel,U=U))	
         })
  }
  
  #Main sections
  interyr 	<- yr[2]   ## interim year within time series for which biomass estimate is available; set to yr[2] if no estimates are available
  interbio 	<- c(0, 1) ## biomass range for interim year, as fraction of k; set to 0 and 1 if not available
  startbt     <- seq(startbio[1], startbio[2], by = 0.05) ## apply range of start biomass in steps of 0.05	
  parbound <- list(r = start_r, k = start_k, lambda = finalbio, sigR)
  
  
  #print out some data and assumption (checks)
  cat("Last year =",max(yr),", last catch (tons)=",ct[nyr],";  ","Resilience =",res,";  ",
      "Process error =", sigR,";  ","Assumed initial biomass (B/k) =", startbio[1],"-", startbio[2], " k",";  ",
      "Assumed intermediate biomass (B/k) in", interyr, " =", interbio[1],"-",interbio[2]," k",";  ",
      "Assumed final biomass (B/k) =", parbound$lambda[1],"-",parbound$lambda[2]," k",";  ",
      "Initial bounds for k (tons)=", parbound$k[1], "-", parbound$k[2],"\n")
  if(user=="Yes") cat("r prior (logmean and logSD) =", parbound$r,"\n")	
  if(!user=="Yes") cat("Initial bounds for r =", parbound$r[1], "-", parbound$r[2],"\n")  
  
  
  ## Execute stock reduction analysis
  R1 = sraMSY(parbound, n)
  
  ## Get statistics on r, k, MSY and determine new bounds for r and k
  r1 	<- R1$r[R1$ell==1]
  k1 	<- R1$k[R1$ell==1]
  msy1  <- r1*k1/4
  mean_msy1 <- exp(mean(log(msy1))) 
  
  if(user=="Yes") max_k1 <- max(k1[r1*k1/4<mean_msy1])  else   
  {
    max_k1a  <- min(k1[r1<1.1*parbound$r[1]]) ## smallest k1 near initial lower bound of r
    max_k1b  <- max(k1[r1*k1/4<mean_msy1]) ## largest k1 that gives mean MSY
    max_k1 <- if(max_k1a < max_k1b) {max_k1a} else {max_k1b}
  }
  
  #Repeat stock reduction analysis analysis with improved initial pars
  if(length(r1)<10)    cat("Too few (", length(r1), ") possible r-k combinations, check input parameters","\n")
  if(length(r1)>=10)
  {
    ## set new upper bound of r to 1.2 max r1
    if(!user=="Yes") parbound$r[2] <- 1.2*max(r1)
    
    ## set new lower bound for k to 0.9 min k1 and upper bound to max_k1 
    parbound$k 	  <- c(0.9 * min(k1), max_k1)
    
    cat("First MSY =", mean_msy1,";  ","First r =", fn.frmt(exp(mean(log(r1)))),";  ",
        "New range for k (tons)=", parbound$k[1], "-", parbound$k[2],"\n")
    if(!user=="Yes") cat("New upper bound for r =", fn.frmt(parbound$r[2]),"\n")	
    
    ## Repeat analysis with new r-k bounds
    R1 = sraMSY(parbound, n)
    
    ## Get statistics on r, k and msy
    r = R1$r[R1$ell==1]
    k = R1$k[R1$ell==1]
    msy = r * k / 4
    mean_ln_msy = mean(log(msy))
    bt=R1$bt[,R1$ell==1]
    bt.rel=R1$bt.rel[,R1$ell==1]
    
    U=R1$U[,R1$ell==1]
    U.neg <- apply(U, 2, function(col) any(col < 0)) 
    if(length(which(U.neg))>0)U=U[,-which(U.neg)]
    U.more.1 <- apply(U, 2, function(col) any(col > 1))
    if(length(which(U.more.1))>0)U=U[,-which(U.more.1)]
    Fish.mort=apply(U, 2, function(x) -log(1-x)) 
    Mean.MSY=exp(mean(log(msy)))  #geometric mean
    
    
    ## Plot statistics
    Mean.MSY_UP=exp(mean_ln_msy + 1.96 * sd(log(msy)))  
    Mean.MSY_LOW=exp(mean_ln_msy - 1.96 * sd(log(msy)))
    
    #Catch and MSY
    fn.fig("CatchMSY_Plots",2000,2400) 
    par(mai=c(1,1.1,.1,.1),las=1,mgp=c(3,.5,0))
    plot(yr, ct, type="l", ylim = c(0, max(ct)), xlab = "Financial year", ylab = "Total catch (tonnes)", main = "",lwd=2,
         cex.lab=2,cex.axis=1.5)
    all.yrs=c(yr[1]-2,yr,yr[length(yr)]+2)
    polygon(c(all.yrs,rev(all.yrs)),  c(rep(Mean.MSY_LOW,length(all.yrs)),rep(Mean.MSY_UP,length(all.yrs))),
            col=rgb(.1,.1,.1,alpha=.2),border="transparent") 
    abline(h=Mean.MSY,col="orange", lwd=2.5)
    legend("bottom",c("MSY (±1.96 SE)"),bty='n',col=c("orange"),lty=1,lwd=2.5,cex=1.5)
    dev.off()
    
    #Multiplot
    fn.fig("CatchMSY_Multi",2400,2400) 
    par(mfcol=c(2,3),mar=c(4,4,1,1),cex.lab=1.5,mgp=c(2.25,.7,0))    
    plot(yr, ct, type="l", ylim = c(0, max(ct)), xlab = "Financial year", ylab = "Total catch (tonnes)",lwd=2)
    abline(h=Mean.MSY,col="black", lwd=2)
    abline(h=Mean.MSY_LOW,col="red")
    abline(h=Mean.MSY_UP,col="red")
    
    hist(r, freq=F, xlim=c(0, 1.2 * max(r)), main = "",col="grey70")
    abline(v=exp(mean(log(r))),col="black",lwd=2)
    abline(v=exp(mean(log(r))-1.96*sd(log(r))),col="red")
    abline(v=exp(mean(log(r))+1.96*sd(log(r))),col="red")
    box()
    
    plot(r1, k1, xlab="r", ylab="k (tonnes)")
    
    hist(k, freq=F, xlim=c(0, 1.2 * max(k)), xlab="k (tonnes)", main = "",col="grey70")
    abline(v=exp(mean(log(k))),col="black", lwd=2)	
    abline(v=exp(mean(log(k))-1.96*sd(log(k))),col="red")
    abline(v=exp(mean(log(k))+1.96*sd(log(k))),col="red")
    box()
    
    plot(log(r), log(k),xlab="ln(r)",ylab="ln(k)")
    abline(v=mean(log(r)))
    abline(h=mean(log(k)))
    abline(mean(log(msy))+log(4),-1, col="black",lwd=2)
    abline(mean(log(msy))-1.96*sd(log(msy))+log(4),-1, col="red")
    abline(mean(log(msy))+1.96*sd(log(msy))+log(4),-1, col="red")
    
    hist(msy, freq=F, xlim=c(0, 1.2 * max(msy)), xlab="MSY (tonnes)",main = "",col="grey70")
    abline(v=exp(mean(log(msy))),col="black", lwd=2)
    abline(v=exp(mean_ln_msy - 1.96 * sd(log(msy))),col="red")
    abline(v=exp(mean_ln_msy + 1.96 * sd(log(msy))),col="red")
    box()
    
    dev.off()
    
    
    #Forward projections, starting at current year  
    nyr.future  <- length(yr.future)    ## number of years in the time series
    ## assign selected r and k combinations
    theta.sel=cbind(r=r,k=k,bt.current=bt[nrow(bt),],sigR=sigR) 
    schaefer.future  <- function(theta)
    {
      with(as.list(theta), {  ## for all combinations of ri & ki
        bt=vector()
        for(i in 1:nyr.future) ## for all years in the time series
        {
          xt=rnorm(1,0, sigR)
          if(i==1) bt[i]=(bt.current+r*bt.current*(1-bt.current/k)-ct[nyr])*exp(xt)else
            bt[i]=(bt[i-1]+r*bt[i-1]*(1-bt[i-1]/k)-ct.future[i-1])*exp(xt) ## calculate biomass as function of previous year's biomass plus net production minus catch
        }
        
        return(list(bt=bt,bt.rel=bt/k))
      })
    }
    Projections = apply(theta.sel,1,schaefer.future)   
    bt.future=matrix(nrow=length(Projections[[1]]$bt),ncol=nrow(theta.sel)) 
    bt.rel.future=bt.future
    for(qq in 1:nrow(theta.sel))
    {
      bt.future[,qq]=Projections[[qq]]$bt
      bt.rel.future[,qq]=Projections[[qq]]$bt.rel
    }
    return(list(r=r,k=k,msy=msy,bt=bt[1:nyr,],bt.rel=bt.rel[1:nyr,],Fish.mort=Fish.mort,
                bt.future=bt.future,bt.rel.future=bt.rel.future,
                "Possible combinations r-k"=length(r),
                "geom. mean r"=fn.frmt(exp(mean(log(r)))),  
                "r +/- 1.96 SD"=paste(fn.frmt(exp(mean(log(r))-1.96*sd(log(r)))),"-",fn.frmt(exp(mean(log(r))+1.96*sd(log(r))))),
                "geom. mean k (tons)"=fn.frmt(exp(mean(log(k)))),
                "k +/- 1.96 SD (tons)"=paste(fn.frmt(exp(mean(log(k))-1.96*sd(log(k)))),"-",fn.frmt(exp(mean(log(k))+1.96*sd(log(k))))),
                "geom. mean MSY (tons)"=fn.frmt(Mean.MSY),
                "MSY +/- 1.96 SD (tons)"=paste(fn.frmt(exp(mean(log(msy))-1.96*sd(log(msy)))),"-",fn.frmt(exp(mean(log(msy))+1.96*sd(log(msy)))))
    ))
  }
}




