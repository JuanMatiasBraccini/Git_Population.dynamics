
library(MASS)

#Prior for r (biomass dynamics) using Leslie matrix
source.hnld="C:/Matias/Analyses/SOURCE_SCRIPTS/Population dynamics/"
fn.source=function(script)source(paste(source.hnld,script,sep=""))
fn.source("Leslie.matrix.R")    
fn.source("fn.fig.R")
fn.source("Catch_MSY.R")


#Functions
fun.rprior.dist=function(Nsims,K,LINF,Temp,Amax,MAT,FecunditY,Cycle)
{
  Fecu=unlist(FecunditY)
  Rprior=fun.Leslie(N.sims=Nsims,k=K,Linf=LINF,Aver.T=Temp,
                    A=Amax,first.age=0,RangeMat=MAT,Rangefec=Fecu,
                    sexratio=0.5,Reprod_cycle=Breed.cycle,Hoenig.only="NO")  
  
  #get mean and sd from lognormal distribution
  gamma.pars=suppressWarnings(fitdistr(Rprior, "gamma"))  
  shape=gamma.pars$estimate[1]        
  rate=gamma.pars$estimate[2]      
  return(list(shape=shape,rate=rate))
  
  #get mean and sd from lognormal distribution
  #LogN.pars=fitdistr(Rprior, "lognormal")  
  #log_mean.r=LogN.pars$estimate[1]    #already in log space     
  #log_sd.r=LogN.pars$estimate[2]      #already in log space     
  #return(list(log_mean.r=log_mean.r,log_sd.r=log_sd.r))
}

density.fun2=function(what,B.ref,CEX)
{
  #Prob above ref point
  f=ecdf(what)
  Prob.below=f(B.ref)
  SEQ=seq(0,1,0.001)
  f.range=f(SEQ)
  id=which.min(abs(SEQ - B.ref))
  X=SEQ[1:id]
  Y=f.range[1:id]
  plot(SEQ,f.range,main='',ylab="",xlab="",type='l',lwd=2,cex.axis=1)
  polygon(c(X,rev(X)),c(Y,rep(0,length(Y))),col=CL)
  Prob.above=round(1-Prob.below,3)
  Prob.below=round(Prob.below,3)
  if(Prob.below<=0.2 & Prob.below>0.01)SRT=55
  if(Prob.below<=0.01) SRT=70 else SRT=0
  
  abline(v=B.ref,lty=2,lwd=2,col="grey50")
  if(Prob.below>=0.01|Prob.below==0)text(B.ref*.99,mean(f.range),fn.exp.less(B.ref,"=",Prob.below),pos=2,cex=CEX,col="brown4",srt=SRT)else
    text(B.ref*.99,mean(f.range),fn.exp.less(B.ref,"<","0.01"),pos=2,cex=CEX,col="brown4",srt=SRT)
  
  if(Prob.above>=0.01|Prob.above==0)text(B.ref*1.05,mean(f.range),fn.exp.more(B.ref,"=",Prob.above),pos=4,cex=CEX,col="brown4",srt=0)else
    text(B.ref*1.05,mean(f.range),fn.exp.more(B.ref,"<","0.01"),pos=4,cex=CEX,col="brown4",srt=0)
}
fn.exp.less=function(x,symb,y) bquote("P"["<"][.(x)]~.(paste0(symb,y)))
fn.exp.more=function(x,symb,y) bquote("P"[">"][.(x)]~.(paste0(symb,y)))
Probs.ref.point=function(what)
{
  f=ecdf(what)
  
  P.below.target=f(B.target)
  P.below.threshold=f(B.threshold)
  P.below.limit=f(B.limit)
  
  P.above.target=1-P.below.target
  P.above.threshold=1-P.below.threshold
  P.above.limit=1-P.below.limit
  
  P.between.thre.tar=P.below.target-P.below.threshold
  P.between.lim.thre=P.below.threshold-P.below.limit
  
  return(data.frame(P.above.target=P.above.target,
                    P.between.thre.tar=P.between.thre.tar,
                    P.between.lim.thre=P.between.lim.thre,
                    P.below.limit=P.below.limit))
}


setPath=function(hndl)setwd(paste(hndl,AssessYr,"/",sep="")) #set paths

#Get r prior
r.prior.dist=fun.rprior.dist(Nsims=10000,K=Growth.F$k,LINF=Growth.F$FL_inf,Temp=TEMP,Amax=Max.age.F,
           MAT=unlist(Age.50.mat),FecunditY=Fecundity,Cycle=Breed.cycle)

for(sc in 1:length(ktch_msy_scen)) if(!is.na(ktch_msy_scen[[sc]]$r.prior)) ktch_msy_scen[[sc]]$r.prior=unlist(r.prior.dist)

#Run catch_msy function                   
setPath(hndl=paste("C:/Matias/Analyses/Population dynamics/",Spec," shark/",sep=''))

if(!file.exists(file.path(getwd(), "/Catch_MSY"))) dir.create(file.path(getwd(), "/Catch_MSY"))   
setwd(file.path(getwd(), "/Catch_MSY"))
Path.ktch_msy=getwd()
Ktch_MSY=ktch_msy_scen
for(sc in 1:length(ktch_msy_scen))
{
  Folder=names(Ktch_MSY)[sc]
  if(!file.exists(paste(Path.ktch_msy,Folder,sep="/"))) dir.create(paste(Path.ktch_msy,Folder,sep="/"))   
  setwd(paste(Path.ktch_msy,Folder,sep="/"))
  Yrs=ct$finyear             
  Tot.Ktch=ct$Total.ktch
  Ktch_MSY[[sc]]=Catch_MSY(ct=Tot.Ktch,
                           yr=Yrs,
                           r.prior=ktch_msy_scen[[sc]]$r.prior,
                           user=ktch_msy_scen[[sc]]$user,
                           k.max=ktch_msy_scen[[sc]]$k.max,
                           startbio=ktch_msy_scen[[sc]]$startbio,
                           finalbio=ktch_msy_scen[[sc]]$finalbio,
                           res=ktch_msy_scen[[sc]]$res,
                           n=ktch_msy_scen[[sc]]$niter,
                           sigR=ktch_msy_scen[[sc]]$sigR,
                           ct.future=ct.future,           
                           yr.future=yr.future
                           )
  
  #Export outputs
  Table1_ktch_MSY=with(Ktch_MSY[[sc]],data.frame(`geom. mean r`,`r +/- 1.96 SD`,`geom. mean k (tons)`,`k +/- 1.96 SD (tons)`,
                                                 `geom. mean MSY (tons)`,`MSY +/- 1.96 SD (tons)`))
  write.csv(Table1_ktch_MSY,"Table1_ktch_MSY.csv",row.names=F)
  write.csv(cbind(Yrs=c(Yrs,yr.future),ct=c(Tot.Ktch,ct.future)),"ct.future.csv",row.names=F)
}


#Outputs
setwd(Path.ktch_msy)
smart.par=function(n.plots,MAR,OMA,MGP) return(par(mfrow=n2mfrow(n.plots),mar=MAR,oma=OMA,las=1,mgp=MGP))
CL='forestgreen'
CL.mean='transparent'

#Output table of scenarios
Tabl.scen.Ktch.MSY=vector('list',length(ktch_msy_scen))
for(i in 1:length(ktch_msy_scen))
{
  dummy=ktch_msy_scen[[i]]
  for(a in 1:length(ktch_msy_scen[[i]])) if(length(dummy[[a]])>1) dummy[[a]]=paste(dummy[[a]],collapse=";")
  Tabl.scen.Ktch.MSY[[i]]=unlist(dummy)
}
Tabl.scen.Ktch.MSY=do.call(rbind,Tabl.scen.Ktch.MSY)
row.names(Tabl.scen.Ktch.MSY)=names(ktch_msy_scen)
write.csv(Tabl.scen.Ktch.MSY,"Scenarios.csv")


#Plot r prior dist
fn.fig("Prior_r", 2000, 2000)
par(las=1,mai=c(1,1.15,.1,.15),mgp=c(3.5,.75,0))
plot(density(rgamma(10000, shape = r.prior.dist$shape, rate = r.prior.dist$rate)),
     lwd=3,main="",xlab=expression(paste(plain("Intrinsic rate of increase (year") ^ plain("-1"),")",sep="")),
     cex.lab=2,cex.axis=1.25,col=1)
# plot(density(rlnorm(10000, meanlog = ktch_msy_scen$`Base case`$r.prior[1], sdlog = ktch_msy_scen$`Base case`$r.prior[2])),
#       lwd=3,main="",xlab=expression(paste(plain("Intrinsic rate of increase (years") ^ plain("-1"),")",sep="")),
#       cex.lab=2,cex.axis=1.25,col=1)
dev.off()


#Total biomass trend
Yrs=c(Yrs,yr.future)  
indx.ftur=(length(Yrs)-years.futures+1):length(Yrs)

Low.percentile=function(Nper,DAT) apply(DAT, 1, function(x) quantile(x, (0+Nper)/100))   #get percentiles
High.percentile=function(Nper,DAT) apply(DAT, 1, function(x) quantile(x, (100-Nper)/100))
fn.cons.po=function(low,up) c(low, tail(up, 1), rev(up), low[1])  #construct polygon

#colfunc <- colorRampPalette(c("grey90","grey50"))
colfunc <- colorRampPalette(c("aliceblue","lightblue3"))
COLS=colfunc(3)
colfunc.f <- colorRampPalette(c("white","burlywood3"))
COLS.f=colfunc.f(3)
fn.plot.percentile=function(DAT,YR,ADD.prob)
{
  #50% of data
  Nper=(100-50)/2
  LOW.50=Low.percentile(Nper,DAT)
  UP.50=High.percentile(Nper,DAT)
  
  #75% of data
  Nper=(100-75)/2
  LOW.75=Low.percentile(Nper,DAT)
  UP.75=High.percentile(Nper,DAT)
  
  #100% of data
  Nper=(100-100)/2
  LOW.100=Low.percentile(Nper,DAT)
  UP.100=High.percentile(Nper,DAT)
  
  #construct polygons
  Year.Vec <-  fn.cons.po(YR[-indx.ftur],YR[-indx.ftur])
  Biom.Vec.50 <- fn.cons.po(LOW.50[-indx.ftur],UP.50[-indx.ftur]) 
  Biom.Vec.75 <- fn.cons.po(LOW.75[-indx.ftur],UP.75[-indx.ftur]) 
  Biom.Vec.100 <-fn.cons.po(LOW.100[-indx.ftur],UP.100[-indx.ftur]) 
  
  id.futr=c((indx.ftur[1]-1),indx.ftur)
  Year.Vec.f <-  fn.cons.po(YR[id.futr],YR[id.futr])
  Biom.Vec.50.f <- fn.cons.po(LOW.50[id.futr],UP.50[id.futr]) 
  Biom.Vec.75.f <- fn.cons.po(LOW.75[id.futr],UP.75[id.futr]) 
  Biom.Vec.100.f <-fn.cons.po(LOW.100[id.futr],UP.100[id.futr]) 
  
  
  #plot
  plot(YR,UP.100,ylim=c(0,max(UP.100)),type="l",ylab="",xlab="",xaxt='n',col='transparent',cex.axis=1.25)
  
  polygon(Year.Vec, Biom.Vec.100, col = COLS[3], border = "grey20")
  polygon(Year.Vec, Biom.Vec.75, col = COLS[2], border = "grey20")
  polygon(Year.Vec, Biom.Vec.50, col = COLS[1], border = "grey20")
  
  polygon(Year.Vec.f, Biom.Vec.100.f, col = COLS.f[3], border = "grey20")
  polygon(Year.Vec.f, Biom.Vec.75.f, col = COLS.f[2], border = "grey20")
  polygon(Year.Vec.f, Biom.Vec.50.f, col = COLS.f[1], border = "grey20")
  
  
  #add probs
  if(ADD.prob=="YES")
  {
    add.probs(id.yr=match(Current,YR),YR,DAT,UP.100,LOW.100)
    add.probs(id.yr=length(YR),YR,DAT,UP.100,LOW.100)
    
    abline(h=B.target,lwd=2,col='grey30',lty=2)
    text(YR[4],B.target,"Target",pos=3,cex=1.1)
    abline(h=B.threshold,lwd=2,col='grey30',lty=2)
    text(YR[4],B.threshold,"Threshold",pos=3,cex=1.1)
    abline(h=B.limit,lwd=2,col='grey30',lty=2)
    text(YR[4],B.limit,"Limit",pos=3,cex=1.1)
  }
  axis(1,at=YR,labels=F,tck=-0.01)
  axis(1,at=seq(YR[1],YR[length(YR)],5),labels=seq(YR[1],YR[length(YR)],5),tck=-0.02,cex.axis=1.25)
}
add.probs=function(id.yr,YR,DAT,UP.100,LOW.100)
{
  f=ecdf(DAT[id.yr,])
  P.below.target=f(B.target)
  P.below.threshold=f(B.threshold)
  P.below.limit=f(B.limit)
  P.above.target=1-P.below.target
  P.above.threshold=1-P.below.threshold
  P.above.limit=1-P.below.limit
  P.between.thre.tar=P.below.target-P.below.threshold
  P.between.lim.thre=P.below.threshold-P.below.limit
  if(P.above.target>0)
  {
    segments(YR[id.yr],B.target,YR[id.yr],UP.100[id.yr],col="forestgreen",lwd=8,lend="butt")  
    text(YR[id.yr],B.target*1.5,paste(round(100*P.above.target,1),"%",sep=""),
         col="black",cex=1.1,srt=45,pos=2,font=2)
  }
  if(P.between.thre.tar>0)
  {
    segments(YR[id.yr],B.target,YR[id.yr],B.threshold,col="yellow",lwd=8,lend="butt")  
    text(YR[id.yr],mean(c(B.target,B.threshold))*1.1,paste(round(100*P.between.thre.tar,1),"%",sep=""),
         col="black",cex=1.1,srt=45,pos=2,font=2)
  }
  if(P.between.lim.thre>0)
  {
    segments(YR[id.yr],B.threshold,YR[id.yr],B.limit,col="orange",lwd=8,lend="butt")
    text(YR[id.yr],mean(c(B.threshold,B.limit))*1.1,paste(round(100*P.between.lim.thre,1),"%",sep=""),
         col="black",cex=1.1,srt=45,font=2,pos=2)
  }
  if(P.below.limit>0)
  {
    segments(YR[id.yr],B.limit,YR[id.yr],LOW.100[id.yr],col="red",lwd=8,lend="butt")
    text(YR[id.yr],B.limit*0.8,paste(round(100*P.below.limit,1),"%",sep=""),
         col="black",cex=1.1,srt=45,pos=2,font=2)
  }
}


  #All scenarios
fn.fig("Biomass",2000,2400)
smart.par(n.plots=length(ktch_msy_scen),MAR=c(3,5,1,1),OMA=rep(.5,4),MGP=c(1,.6,0))
for(sc in 1:length(ktch_msy_scen))
{
  Ktch_MSY_Rel.bio=rbind(Ktch_MSY[[sc]]$bt,Ktch_MSY[[sc]]$bt.future)   

  #get geometric mean
  Ktch_MSY_rel_bt_mean=Ktch_MSY_rel_bt_lowSE=Ktch_MSY_rel_bt_upSE=nrow(Ktch_MSY_Rel.bio)
  for(nr in 1:nrow(Ktch_MSY_Rel.bio))
  {
    Ktch_MSY_rel_bt_mean[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])))
    Ktch_MSY_rel_bt_upSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) + 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
    Ktch_MSY_rel_bt_lowSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) - 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
  }
  YMAX=max( Ktch_MSY_rel_bt_upSE,na.rm=T)
  plot(Yrs,Ktch_MSY_rel_bt_mean,cex=1.5,pch=19,col=CL.mean,lwd=3,ylim=c(0,YMAX),xaxt='n',xlab="",ylab="",cex.axis=1.25)
  segments(Yrs,Ktch_MSY_rel_bt_lowSE,Yrs,Ktch_MSY_rel_bt_upSE,col=CL)
  points(Yrs[indx.ftur],Ktch_MSY_rel_bt_mean[indx.ftur],pch=19,col=CL.mean,cex=1.65)  #highlight future projections
  segments(Yrs[indx.ftur],Ktch_MSY_rel_bt_lowSE[indx.ftur],Yrs[indx.ftur],Ktch_MSY_rel_bt_upSE[indx.ftur],col="brown4")
  axis(1,Yrs,labels=F,tck=-0.015)
  axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
  legend("bottomleft",names(Ktch_MSY)[sc],bty='n',cex=1.5)
}
mtext("Total biomass (tonnes)",2,cex=1.75,las=3,line=-1.5,outer=T)
mtext("Financial year",1,line=-1,cex=1.75,outer=T)
dev.off()


  #Base case only
fn.fig("Biomass_Base case",2000,2400)
par(mar=c(3,5,1,1),oma=rep(.5,4),mgp=c(1,.6,0),las=1)
sc=match("Base case",names(ktch_msy_scen))
Ktch_MSY_Rel.bio=rbind(Ktch_MSY[[sc]]$bt,Ktch_MSY[[sc]]$bt.future)   

#get geometric mean
Ktch_MSY_rel_bt_mean=Ktch_MSY_rel_bt_lowSE=Ktch_MSY_rel_bt_upSE=nrow(Ktch_MSY_Rel.bio)
for(nr in 1:nrow(Ktch_MSY_Rel.bio))
{
  Ktch_MSY_rel_bt_mean[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])))
  Ktch_MSY_rel_bt_upSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) + 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
  Ktch_MSY_rel_bt_lowSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) - 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
}
YMAX=max( Ktch_MSY_rel_bt_upSE,na.rm=T)
plot(Yrs,Ktch_MSY_rel_bt_mean,cex=1.5,pch=19,col=CL.mean,lwd=3,ylim=c(0,YMAX),xaxt='n',xlab="",ylab="",cex.axis=1.25)
segments(Yrs,Ktch_MSY_rel_bt_lowSE,Yrs,Ktch_MSY_rel_bt_upSE,col=CL)
points(Yrs[indx.ftur],Ktch_MSY_rel_bt_mean[indx.ftur],pch=19,col=CL.mean,cex=1.65)  #highlight future projections
segments(Yrs[indx.ftur],Ktch_MSY_rel_bt_lowSE[indx.ftur],Yrs[indx.ftur],Ktch_MSY_rel_bt_upSE[indx.ftur],col="brown4")
axis(1,Yrs,labels=F,tck=-0.015)
axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
mtext("Total biomass (tonnes)",2,cex=1.75,las=3,line=-1,outer=T)
mtext("Financial year",1,line=-1,cex=1.75,outer=T)
dev.off()



#Relative biomass trend

  #All scenarios
fn.fig("Biomass_relative",2000,2400)
smart.par(n.plots=length(ktch_msy_scen),MAR=c(3,4,1,1),OMA=rep(.5,4),MGP=c(1,.6,0))
for(sc in 1:length(ktch_msy_scen))
{
  Ktch_MSY_Rel.bio=rbind(Ktch_MSY[[sc]]$bt.rel,Ktch_MSY[[sc]]$bt.rel.future)   

  #Percentile   
  fn.plot.percentile(DAT=Ktch_MSY_Rel.bio,YR=Yrs,ADD.prob="YES")
  
  #Geometric mean
  # Ktch_MSY_rel_bt_mean=Ktch_MSY_rel_bt_lowSE=Ktch_MSY_rel_bt_upSE=nrow(Ktch_MSY_Rel.bio)
  # for(nr in 1:nrow(Ktch_MSY_Rel.bio))
  # {
  #   zz=Ktch_MSY_Rel.bio[nr,]
  #   zz[zz<0]=1e-6
  #   Ktch_MSY_rel_bt_mean[nr]=exp(mean(log(zz)))
  #   Ktch_MSY_rel_bt_upSE[nr]=exp(mean(log(zz)) + 1.96 * sd(log(zz)))
  #   Ktch_MSY_rel_bt_lowSE[nr]=exp(mean(log(zz)) - 1.96 * sd(log(zz)))
  # }
  # plot(Yrs,Ktch_MSY_rel_bt_mean,cex=1.5,pch=19,col=CL.mean,lwd=3,ylim=c(0,1),xaxt='n',xlab="",ylab="",cex.axis=1.25)
  # segments(Yrs,Ktch_MSY_rel_bt_lowSE,Yrs,Ktch_MSY_rel_bt_upSE,col=CL)
  # points(Yrs[indx.ftur],Ktch_MSY_rel_bt_mean[indx.ftur],pch=19,col=CL.mean,cex=1.65)  #highlight future projections
  # segments(Yrs[indx.ftur],Ktch_MSY_rel_bt_lowSE[indx.ftur],Yrs[indx.ftur],Ktch_MSY_rel_bt_upSE[indx.ftur],col="brown4")
  # abline(h=B.target,lwd=2,col='grey30',lty=2)
  # text(Yrs[3],B.target,"Target",pos=3,cex=1.1)
  # abline(h=B.threshold,lwd=2,col='grey30',lty=2)
  # text(Yrs[5],B.threshold,"Threshold",pos=3,cex=1.1)
  # abline(h=B.limit,lwd=2,col='grey30',lty=2)
  # text(Yrs[3],B.limit,"Limit",pos=3,cex=1.1)
  # axis(1,Yrs,labels=F,tck=-0.015)
  # axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
  legend("bottomleft",names(Ktch_MSY)[sc],bty='n',cex=1.5)
}
mtext("Relative biomass",2,cex=1.75,las=3,line=-1.5,outer=T)
mtext("Financial year",1,line=-1,cex=1.75,outer=T)
dev.off()


  #Base case only
fn.fig("Biomass_relative_Base case",2000,2400)
par(mar=c(3,3.5,1,1),oma=rep(.5,4),mgp=c(1,.6,0),las=1)
sc=match("Base case",names(ktch_msy_scen))
Ktch_MSY_Rel.bio=rbind(Ktch_MSY[[sc]]$bt.rel,Ktch_MSY[[sc]]$bt.rel.future)   

  #Percentile   
fn.plot.percentile(DAT=Ktch_MSY_Rel.bio,YR=Yrs,ADD.prob="YES")

#   #Geometric mean
Ktch_MSY_rel_bt_mean=Ktch_MSY_rel_bt_lowSE=Ktch_MSY_rel_bt_upSE=nrow(Ktch_MSY_Rel.bio)
for(nr in 1:nrow(Ktch_MSY_Rel.bio))
{
  Ktch_MSY_rel_bt_mean[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])))
  Ktch_MSY_rel_bt_upSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) + 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
  Ktch_MSY_rel_bt_lowSE[nr]=exp(mean(log(Ktch_MSY_Rel.bio[nr,])) - 1.96 * sd(log(Ktch_MSY_Rel.bio[nr,])))
}
# plot(Yrs,Ktch_MSY_rel_bt_mean,cex=1.5,pch=19,col=CL.mean,lwd=3,ylim=c(0,1),xaxt='n',xlab="",ylab="",cex.axis=1.25)
# segments(Yrs,Ktch_MSY_rel_bt_lowSE,Yrs,Ktch_MSY_rel_bt_upSE,col=CL)
# points(Yrs[indx.ftur],Ktch_MSY_rel_bt_mean[indx.ftur],pch=19,col=CL.mean,cex=1.65)  #highlight future projections
# segments(Yrs[indx.ftur],Ktch_MSY_rel_bt_lowSE[indx.ftur],Yrs[indx.ftur],Ktch_MSY_rel_bt_upSE[indx.ftur],col="brown4")
# 
# abline(h=B.target,lwd=2,col='grey30',lty=2)
# text(Yrs[3],B.target,"Target",pos=3,cex=1.25)
# abline(h=B.threshold,lwd=2,col='grey30',lty=2)
# text(Yrs[5],B.threshold,"Threshold",pos=3,cex=1.25)
# abline(h=B.limit,lwd=2,col='grey30',lty=2)
# text(Yrs[3],B.limit,"Limit",pos=3,cex=1.25)
# axis(1,Yrs,labels=F,tck=-0.015)
# axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
legend("bottomleft",c("50%","75%","100%"),fill=COLS,bty='n',title="Percentile",cex=1.25)
mtext("Relative biomass",2,cex=1.75,las=3,line=-1,outer=T)
mtext("Financial year",1,line=-1,cex=1.75,outer=T)
dev.off()

write.csv(cbind(Yrs=Yrs,geom.mean=Ktch_MSY_rel_bt_mean,
         Low95CI=Ktch_MSY_rel_bt_lowSE,Up95CI=Ktch_MSY_rel_bt_upSE),"rel.biom.base.case.csv",row.names=F)


#Current depletion of total biomass

  #Base case only
sc=match("Base case",names(ktch_msy_scen))
Ktch_MSY_Rel.bio=rbind(Ktch_MSY[[sc]]$bt.rel,Ktch_MSY[[sc]]$bt.rel.future)
Ktch_MSY_current_yr=Ktch_MSY_Rel.bio[match(Current,Yrs),]

fn.fig("Posterior_current_year_depletion_Base case",1200,2400)
par(mfcol=c(3,1),mar=c(3.5,3.5,.1,1),oma=rep(.5,4),mgp=c(1,.6,0),las=1)
density.fun2(what=Ktch_MSY_current_yr,B.ref=B.target,CEX=1.5) 
legend('bottomright',"Target",bty='n',cex=2)
density.fun2(what=Ktch_MSY_current_yr,B.ref=B.threshold,CEX=1.5) 
legend('bottomright',"Threshold",bty='n',cex=2)
density.fun2(what=Ktch_MSY_current_yr,B.ref=B.limit,CEX=1.5) 
legend('bottomright',"Limit",bty='n',cex=2)
mtext("Probability",2,cex=1.75,las=3,line=-1.75,outer=T)
mtext(paste(Current,"relative biomass"),1,line=-1,cex=1.75,outer=T)
dev.off()

#export probabilities for Weight of Evidence
write.csv(Probs.ref.point(what=Ktch_MSY_current_yr),"Consequence_likelihood_total.csv",row.names=F)

#export relative biomass runs
write.csv(Ktch_MSY_Rel.bio,"Base case/Ktch_MSY_Rel.bio.csv",row.names=F)


#Future depletion of total biomass   

#Base case only
Ktch_MSY_future=Ktch_MSY_Rel.bio[length(Yrs),]

fn.fig("Posterior_future_year_depletion_Base case",1200,2400)
par(mfcol=c(3,1),mar=c(3.5,3.5,.1,1),oma=rep(.5,4),mgp=c(1,.6,0),las=1)
density.fun2(what=Ktch_MSY_future,B.ref=B.target,CEX=1.5) 
legend('bottomright',"Target",bty='n',cex=2)
density.fun2(what=Ktch_MSY_future,B.ref=B.threshold,CEX=1.5) 
legend('bottomright',"Threshold",bty='n',cex=2)
density.fun2(what=Ktch_MSY_future,B.ref=B.limit,CEX=1.5) 
legend('bottomright',"Limit",bty='n',cex=2)
mtext("Probability",2,cex=1.75,las=3,line=-1.75,outer=T)
mtext(paste(Yrs[length(Yrs)],"relative biomass"),1,line=-1,cex=1.75,outer=T)
dev.off()

#export probabilities for Weight of Evidence
write.csv(Probs.ref.point(what=Ktch_MSY_future),"Consequence_likelihood_total_future.csv",row.names=F)



#Fishing mortality trend
Yrs=Yrs[-match(yr.future,Yrs)]              

  #All scenarios
fn.fig("Fishing_mortality",2000,2400)
smart.par(n.plots=length(ktch_msy_scen),MAR=c(3,5,1,1),OMA=rep(.5,4),MGP=c(1,.6,0))
for(sc in 1:length(ktch_msy_scen))
{
  Fish.mort=Ktch_MSY[[sc]]$Fish.mort
  
  #Percentiles
  fn.plot.percentile(DAT=Fish.mort,YR=Yrs,ADD.prob="NO")

  # #get geometric mean
  # Fish.mort_mean=Fish.mort_lowSE=Fish.mort_upSE=nrow(Fish.mort)
  # for(nr in 1:nrow(Fish.mort))
  # {
  #   Fish.mort_mean[nr]=exp(mean(log(Fish.mort[nr,])))
  #   Fish.mort_upSE[nr]=exp(mean(log(Fish.mort[nr,])) + 1.96 * sd(log(Fish.mort[nr,])))
  #   Fish.mort_lowSE[nr]=exp(mean(log(Fish.mort[nr,])) - 1.96 * sd(log(Fish.mort[nr,])))
  # }
  # plot(Yrs,Fish.mort_mean,cex=1.5,pch=19,col=CL.mean,lwd=3,ylim=c(0,.7),xaxt='n',xlab="",ylab="",cex.axis=1.25)
  # segments(Yrs,Fish.mort_lowSE,Yrs,Fish.mort_upSE,col=CL)
  # points(Yrs[indx.ftur],Fish.mort_mean[indx.ftur],pch=19,col=CL.mean,cex=1.65)  #highlight future projections
  # segments(Yrs[indx.ftur],Fish.mort_lowSE[indx.ftur],Yrs[indx.ftur],Fish.mort_upSE[indx.ftur],col="brown4")
  # axis(1,Yrs,labels=F,tck=-0.015)
  # axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
  legend("topleft",names(Ktch_MSY)[sc],bty='n',cex=1.5)
}
mtext(expression(paste(plain("Fishing mortality (year") ^ plain("-1"),")",sep="")),
      2,cex=1.75,las=3,line=-2.5,outer=T)   
mtext("Financial year",1,line=-1,cex=1.75,outer=T)
dev.off()


  #Base case only
fn.fig("Fishing_mortality_Base case",2000,2400)
par(mar=c(3,5,1,1),oma=rep(.5,4),mgp=c(1,.6,0),las=1)
sc=match("Base case",names(ktch_msy_scen))
Fish.mort=Ktch_MSY[[sc]]$Fish.mort

  #Percentiles
fn.plot.percentile(DAT=Fish.mort,YR=Yrs,ADD.prob="NO")

#   #Geometric mean
# Fish.mort_mean=Fish.mort_lowSE=Fish.mort_upSE=nrow(Fish.mort)
# for(nr in 1:nrow(Fish.mort))
# {
#   Fish.mort_mean[nr]=exp(mean(log(Fish.mort[nr,])))
#   Fish.mort_upSE[nr]=exp(mean(log(Fish.mort[nr,])) + 1.96 * sd(log(Fish.mort[nr,])))
#   Fish.mort_lowSE[nr]=exp(mean(log(Fish.mort[nr,])) - 1.96 * sd(log(Fish.mort[nr,])))
# }
# plot(Yrs,Fish.mort_mean,cex=1.5,pch=19,col=CL.mean,lwd=3,ylim=c(0,max(Fish.mort_upSE)),xaxt='n',xlab="",ylab="",cex.axis=1.25)
# segments(Yrs,Fish.mort_lowSE,Yrs,Fish.mort_upSE,col=CL)
# points(Yrs[indx.ftur],Fish.mort_mean[indx.ftur],pch=19,col=CL.mean,cex=1.65)  #highlight future projections
# segments(Yrs[indx.ftur],Fish.mort_lowSE[indx.ftur],Yrs[indx.ftur],Fish.mort_upSE[indx.ftur],col="brown4")
# axis(1,Yrs,labels=F,tck=-0.015)
# axis(1,Yrs[seq(1,length(Yrs),5)],labels=Yrs[seq(1,length(Yrs),5)],tck=-0.030,cex.axis=1.25)
mtext(expression(paste(plain("Fishing mortality (year") ^ plain("-1"),")",sep="")),
      2,cex=1.75,las=3,line=-2,outer=T)
mtext("Financial year",1,line=-1,cex=1.75,outer=T)
legend("topleft",c("50%","75%","100%"),fill=COLS,bty='n',title="Percentile",cex=1.25)
dev.off()

