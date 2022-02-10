#Zhou et al 2018 Optimized Catch-Only Assessment Method
#source: https://academic.oup.com/icesjms/article/75/3/964/4772849#supplementary-data

# ----------------Optimized catch-only method function----------------
##  Predictors for BRT model
# First load the following function catchParam.
# This function derives 56 predictors based on catch history.
# It can be used for batch processing of multiple stocks.
# Data format must be 3 columns: stock, yr, catch.

require(segmented)  
require(gbm)
require(dismo)
require(fGarch) # for skewed normal dist

catchParam = function(catchData) {
  sid = unique(as.character(catchData$stock))
  n.stock = length(sid)
  para = matrix(NA, n.stock, 57)
  n.ab = 7  # number of years at the beginning and ending time series
  a=b=a0=b0=matrix(0, n.stock, n.ab)
  for(i in 1:n.stock) {
    stk = sid[i]
    dat = subset(catchData, stock==stk)
    yr = dat$yr-min(dat$yr) +1
    midYr = mean(yr)
    C = dat$catch/max(dat$catch)   # regressions are based on scaled catch
    C05max = sum(C[-length(C)]>0.5)
    Cmean = mean(C)
    nyr = length(yr)
    nyrToCmax = min(yr[C==max(C)]-yr[1]+1)
    nyrToCmaxR = nyrToCmax/nyr
    nyrAfterCmax = nyr-nyrToCmax
    yr.cent = yr-midYr           # regressions are based on centred year
    line0 = lm(C~yr.cent)
    yr1.cent = yr[1:nyrToCmax]-mean(yr[1:nyrToCmax])
    line1 = lm(C[1:nyrToCmax]~yr1.cent)
    yr2.cent = yr[nyrToCmax:nyr]-mean(yr[nyrToCmax:nyr])
    line2 = lm(C[nyrToCmax:nyr]~yr2.cent)
    aa0 = summary(line0)$coeff[1]; bb0 = summary(line0)$coeff[2]    #all yr
    aa1 = summary(line1)$coeff[1]; bb1 = summary(line1)$coeff[2]    #before Cmax
    aa2 = summary(line2)$coeff[1]; bb2 = summary(line2)$coeff[2]    #after Cmax
    for (j in 1:n.ab) { # periodical regressions
      yrLast.cent = yr[(nyr-j):nyr]-mean(yr[(nyr-j):nyr])
      l.last = lm(C[(nyr-j):nyr]~yrLast.cent)                   # last several years
      a[i,j] = summary(l.last)$coeff[1];
      b[i,j]=summary(l.last)$coeff[2]
      yrBegin.cent = yr[1:(j+1)]-mean(yr[1:(j+1)])
      l.begin = lm(C[1:(j+1)] ~ yrBegin.cent)                  # beginning several years
      a0[i,j] = summary(l.begin)$coeff[1]
      b0[i,j] = summary(l.begin)$coeff[2] }
    # segmented regression: use yr and breakpoint is between 0-1
    f = tryCatch(segmented(line0, seg.Z=~yr, psi=list(yr=median(yr))) , error=function(err) {})
    if(is.null(f)) {
      a.spline = NA; b1.spline = NA; b2.spline = NA; breakPoint = NA
    } else {
      a.spline = summary(f)$coef[1] 
      b1.spline = summary(f)$coef[2]
      slp= slope(f) 
      b2.spline = slp$yr[2]
      breakPoint = (round(f$psi.history[[5]],0)-yr[1] +1)/nyr
    }
    para[i,] = c(aa0, aa1, aa2, bb0, bb1, bb2, a[i,], b[i,], a0[i,], b0[i,], a.spline, b1.spline,b2.spline, breakPoint, C[1:n.ab], C[(nyr-n.ab): nyr], Cmean, nyrToCmaxR, nyr, C05max)
  }   # end params
  colnames(para) = c('a.AllYr', 'a.BfMax', 'a.AfMax', 'b.AllYr', 'b.BfMax', 'b.AfMax', paste('a.LsY',1:n.ab, sep=""), paste('b.LsY',1:n.ab, sep=""), paste("a.BgY", 1:n.ab, sep = ""), paste('b.BgY', 1:n.ab, sep=""), 'a.seg', 'b1.seg', 'b2.seg','breakPoint',  paste('c.BgY',1:n.ab, sep=""), paste('c.LsY', n.ab:0, sep=""), 'Cmean','nyrToCmaxR','nyr', 'C05max')
  para = data.frame(stock=sid, para)
  return(para) 
} 

#### biomass dynamics model: one parameter r, using optimize
BDM = function(K, r, S, b, C) {  
  nyr = length(C)
  B = vector()   
  B[1] = K*b  
  for (i in 1:nyr) {
    B[i+1] = max(min(B[i]+r*B[i]*(1-B[i]/K)-C[i], K), 0) 
  }
  if (all(B[-nyr]>C) & all(B<=K)) abs(B[nyr]/K-S)  else max(K)*10^4 
}

### derive S = B/K prior bwt [0,1]
Sdistrib = function(n, s_mean) {
  nv = 0 ; n.redo = 0
  while(nv < n) {
    n.redo = n.redo+1       
    if(s_mean<=0.5) {
      si1 = rsnorm(n*n.redo, mean=max(s_mean,0)-0.072, sd=0.189, xi=0.763)
      si = si1[si1>0 & si1<1] 
    } else if(s_mean>0.5) {
      si1 = rsnorm(n*n.redo, mean=max(s_mean,0)+0.179, sd=0.223, xi=0.904)
      si = si1[si1>0 & si1<1] }
    if(length(si)>n) si = sample(si,n); 
    nv = length(si) }
  return (si) 
}

### draw biomass trajectories
drawBt = function(cdat, oc0){
  B = Bmed = vector()
  n.sim = 100
  oc1 = oc0[oc0$obj<0.01 & oc0$k>1.01*min(oc0$k) & oc0$k<0.99*max(oc0$k), c(2:5)] #  eliminate bordering effect
  oc2 = oc1[oc1$k>quantile(oc1$k, 0.25) & oc1$k<quantile(oc1$k, 0.75),]
  smp = sample(1:nrow(oc2), n.sim, replace=T)
  if(min(oc2$k)>1000) { 
    k = oc2[smp,1]/1000
    C = cdat$catch/1000
  } else { 
    k = oc2[smp,1]
    C = cdat$catch }  
  r = oc2[smp,2]
  kmed = median(k); rmed = median(r)
  stock = cdat$stock[1]     
  yr=cdat$yr; nyr = length(yr)
  
  plot(rep(0, nyr)~yr, type='n',  ylim=c(0, max(k)*1.2),xlab='', ylab='', las=1, yaxs='i')
  legend("topright", legend=stock, bty='n', cex=0.8)
  Bmed[1] = kmed
  for (j in 1:n.sim){
    B[1] = k[j]
    for(t in 1:(nyr-1)) {
      B[t+1] = (B[t])+r[j]*B[t]*(1-B[t]/k[j])-C[t]
      Bmed[t+1] = (Bmed[t])+rmed*Bmed[t]*(1-Bmed[t]/kmed)-C[t]
    }
    lines(yr, B, col=rgb(runif(1,0,j)/n.sim,(n.sim-runif(1,0,j))/n.sim, 1/(n.sim+100), alpha=0.6) )
  }
  lines(yr, Bmed, lwd=2)    
  lines(C~yr, type='h', lwd=2 )
  
  mtext('Year', side=1, line=2, outer=T)
  if(min(oc2$k)>1000) mtext(expression('Biomass (' %*% '1000 ton)'), side=2, line=2, outer=T) 
  else  mtext(expression('Biomass (ton)'), side=2, line=2, outer=T) 
} 



# ----------------Optimized Catch-Only Method----------------
require(fGarch) # for skewed normal dist
require(plyr)

## Import catch and M data
## catch data must have three columns: stock, yr, catch
## natural mortality data has two coulumns: stock, M   
## Example: SESS stocks in Australia
Dir =("E:/")
catchData =read.table(paste0(Dir, 'catchData.txt'), head=T)
Mdata =read.table(paste0(Dir, 'sessM.txt'), head=T)  

## if M is unkonwn, derive M from other life history parameters, e.g.: 
M = vector()
M[1] = 4.899*Tmax^-0.916             # Tmax = max age
M[2] = 4.11*k^0.73*Linf^-0.33        # k and Linf: Bertalanffy growth parameters
M[3] = 1.82*k
M[4] = 1.65/Tmat                     # Tmax = maturation age
M = mean(M)

## derive prior for stock saturation S = B/K
## load saved BRT models: two alternative models. You can use 1 or both models.
## model 1: BRTmodelP8   (using 8 predictors )
## model 2: BRTmodelP38  (using 38 predictors)
load(file = paste0(Dir, "BRTmodelP8.RData"))   # model 1
load(file = paste0(Dir, "BRTmodelP38.RData"))  # model 2

## derive predictors from catch history
sessPar = catchParam(catchData)

## centering the first 37 predictors   
sessParCent = scale(sessPar[,2:38], center=BRTmodelP8$parMean, scale=F)

## construct predition data  
stockName = unique(as.character(catchData$stock))
nstk = length(stockName)
predDat = data.frame(stock=stockName, sessParCent, sessPar[,39:57])

## estimating saturation S 
s8 = predict.gbm(BRTmodelP8$model, predDat, n.trees=BRTmodelP8$model$gbm.call$best.trees, type="response")
s38 = predict.gbm(BRTmodelP38$model, predDat, n.trees=BRTmodelP38$model$gbm.call$best.trees, type="response")  

## Bias correction
## for model 1
slr8.a = BRTmodelP8$slr[[1]]; slr8.b = BRTmodelP8$slr[[2]]; 
sBC8 = (s8-slr8.a)/slr8.b
sBC8[sBC8<=0] = 0.01 

## for model 2
slr38.a = BRTmodelP38$slr[[1]]; slr38.b = BRTmodelP38$slr[[2]]; 
sBC38 = (s38-slr38.a)/slr38.b
sBC38[sBC38<=0] = 0.01

## Output S 
Sdata = data.frame(stock=stockName, S8=sBC8, S38=sBC38, S = (sBC8+sBC38)/2)


# ----------------OCOM implementation----------------
nsim = 10000 
summ = array(NA, dim=c(5,4,nstk)) 

for (i in 1:nstk) {  # i=1
  stk = stockName[i]
  dat = subset(catchData, stock==stk)
  C = dat$catch ; 
  yr = dat$yr
  ## derive r prior with lognorm dist 
  M = Mdata$M[Mdata$stock==stk] 
  r_median = 2*0.866*M;              # for teleosts
  r_sig2 = (2*0.866)^2*(0.0012+0.23)
  #   r_median = 2*0.41*M;              # for chondrichthyans
  #      r_sig2 = (2*0.41)^2*(0.0012+0.2   
  r_sig = sqrt(r_sig2)          
  r_mu = log(r_median)    
  ri = rlnorm(nsim, r_mu, r_sig)
  
  ## derive S prior btw 0 and 1
  s_mean = Sdata$S[Sdata$stock==stk] # column 2 = BRTmodelP8
  si = Sdistrib(nsim, s_mean)
  
  rs = cbind(r=ri, s=si)
  k.low = max(C); k.up=max(C)*200
  opt = apply(rs, 1, function(x) { optimize(BDM, c(k.low, k.up), r=x[["r"]], S=x[["s"]], b=1, C=C) }) 
  
  ki = sapply(opt, '[[', 1)
  msy = ki*ri/4
  obji = sapply(opt,'[[',2)
  kr = data.frame(k=ki, r=ri, msy, s=si, obj=obji)
  write.csv(kr, paste0(Dir, "test/krms", i, ".csv", sep=""))
  
  ## summary 
  kr2 = kr
  kr2[kr2$k<1.01*k.low | kr2$k>0.99*k.up | kr2$obj>0.01,] = NA # eliminate bordering effect
  kr2 = na.omit(kr2)
  # plot(log(kr2$k), log(kr2$r), xlab="K", ylab='r')  
  # abline(h=mean(log(kr$r)), v=mean(log(kr$k)), lty=2, col=2)
  summ[,,i] = apply(kr2, 2, function(x) quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)))[,1:4]
}

## summary result
summ2 = adply(summ, 3) 
colnames(summ2)=c('stockID', 'k', 'r', 'MSY', 'S')
sumOut = data.frame(percent=rep(c(5, 25, 50, 75, 95),nstk), summ2) 
write.csv(sumOut, paste0(Dir, "result.csv"))

#####################################
## distribution 
par(mfrow=c(2,2), mai=c(.8,.8,0.1,0.1))
hist(kr2$k, xlab='K', main='', las=1); abline(v=median(kr$k), lwd=2); 
hist(kr2$r, xlab='r', main='', las=1); abline(v=median(kr$r), lwd=2);
hist(kr2$msy, xlab='MSY', main='', las=1); abline(v=median(kr$msy), lwd=2);
hist(kr2$s, xlab='S', main='', las=1); abline(v=median(kr$s), lwd=2);

#####################################
## draw 100 biomass trajectories

windows()
par(mfrow=c(5,3), mai=c(.2,.25,0.1,0.1), omi=c(0.5,0.5,0,0))

for (i in 1:nstk) { # i=1
  cdat = subset(catchData, stock==stockName[i])
  oc0 = read.csv(paste0(Dir,'test\\krms', i, '.csv', sep=''))
  drawBt(cdat, oc0)
}

