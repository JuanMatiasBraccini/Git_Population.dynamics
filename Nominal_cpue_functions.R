
##=================================================================================================
## Ainslie Denham 
## March 2018
##
## This code contains functions to 
## (1) Calculate season-specific qualification level (QL)
## (2) Calculate and plot means (and 95% CIs) of data set via 4 calculations
##
##     Method 1: Ratio Estimator
##     CPUE =  average catch  =  total catch / n  = total catch
##             --------------    ----------------   ------------ 
##             average effort    total effort / n   total effort
##     with CI calculated using SE =sqrt(Var) where
##     Var = 1/n * (  My^2*Sx^2/(Mx^4) + Sy^2/(Mx^2) - 2*My*r*Sx*Sy/(Mx^3)  )
##     see various survey statistics texts.
##
##     Method 2: Arithmetic Mean with CI calculated as mean +/- 1.96*sigma/sqrt(n)
##
##     Method 3: Lognormal Mean with CI calculated via Cox method (Zhou & Gao 1997)
##
##     Method 4: Delta Lognormal Mean with CI calculated via modified Cox method (Fletcher 2008)
##
##=================================================================================================

library(dplyr)
library(fields)
if("plyr" %in% (.packages())) detach("package:plyr", unload=TRUE)  #remove plyr if loaded as it interfer with dplyr



##=================================================================================================
## Function to calculate season-specific qualification level to explain x% (default 90%) of annual catch 
##=================================================================================================

CalcQL = function(dat, catch.column="catch", prop.column="prop", season.column="season", prop.catch=0.9)
{
    eval(parse(text=paste("dat$prop = dat[, which(names(dat)=='", prop.column, "')]", sep='')))
    eval(parse(text=paste("dat$catch = dat[, which(names(dat)=='", catch.column, "')]", sep='')))
    dat = dat[order(dat$season, 1-dat$prop), ]
    seasons = sort(unique(dat$season))
    outdat = data.frame(season = seasons, ql=NA)
    for (i in seasons)
    {
      sdat = subset(dat, season==i & !is.na(prop))
      sdat$cumsum = cumsum(sdat$catch)
      dx = min(which(cumsum(sdat$catch) >= prop.catch * sum(sdat$catch)))
      outdat$ql[outdat$season==i] = sdat[dx, ]$prop
    }
    return(outdat)
  }


##=================================================================================================
## Function to calculate mean cpue (+95% CI) using ratio estimator, lognormal and delta-lognormal mean
##=================================================================================================

CalcMeanCPUE = function(cpuedata, catch.column="catch", effort.column="effort",
                        plot.title, cpue.units, draw.plot=TRUE, show.legend=FALSE,PaR,showLNMean)
  {
    eval(parse(text=paste("cpuedata$catch = cpuedata[, which(names(cpuedata)=='", catch.column, "')]", sep='')))
    eval(parse(text=paste("cpuedata$effort = cpuedata[, which(names(cpuedata)=='", effort.column, "')]", sep='')))
    cpuedata$cpue = cpuedata$catch / cpuedata$effort
    
    ## Nominal CPUE
    out1 = cpuedata %>%
      group_by(season) %>%
      summarise(My = mean(catch),
                Mx = mean(effort),
                Sy = sd(catch),
                Sx = sd(effort),
                r = cor(catch, effort),
                n = length(catch)) %>%
      as.data.frame
    out1$r[is.na(out1$r)] = 0
    out1 = out1 %>%
      mutate(mean=My/Mx,
             se =  sqrt(1/n*(My^2*Sx^2/(Mx^4) + Sy^2/(Mx^2) - 2*My*r*Sx*Sy/(Mx^3))),
             lowCL = mean - 1.96*se,
             uppCL = mean + 1.96*se) %>%
      as.data.frame
    out1$method = "Nominal"
    #head(out1,2)
    
    ## Arithmetic Mean CPUE
    out2 = cpuedata %>%
      group_by(season) %>%
      summarise(mean = mean(cpue),
                n = length(cpue),
                sd = sd(cpue)) %>%
      mutate(lowCL = mean - 1.96*sd/sqrt(n),
             uppCL = mean + 1.96*sd/sqrt(n)) %>%
      as.data.frame
    #head(out2,2)
    out2$method = "Mean"
    
    ## Lognormal CPUE
    out3 = cpuedata %>%
      filter(cpue>0) %>%
      group_by(season) %>%
      summarise(ymean = mean(log(cpue)),
                n = length(cpue),
                ysigma = sd(log(cpue))) %>%
      mutate(mean=exp(ymean + ysigma^2/2) ,
             ySE = sqrt(ysigma^2/n+ysigma^4/(2*(n-1))),
             lowCL = exp(ymean + ysigma^2/2 - 1.96*ySE),
             uppCL = exp(ymean + ysigma^2/2 + 1.96*ySE)) %>%
      as.data.frame
   # head(out3,2)
    out3$method = "LnMean"
    
    ## Delta - Lognormal CPUE
    out4 = cpuedata %>%
      group_by(season) %>%
      summarise(n = length(cpue),
                m = length(cpue[cpue>0]),
                mean.lognz = mean(log(cpue[cpue>0])),
                sd.lognz = sd(log(cpue[cpue>0]))) %>%
      mutate(p.nz = m/n,
             theta = log(p.nz)+mean.lognz+sd.lognz^2/2,
             c = (1-p.nz)^(n-1),
             d = 1+(n-1)*p.nz,
             vartheta = ((d-c)*(1-c*d)-m*(1-c)^2)/(m*(1-c*d)^2)+
               sd.lognz^2/m+sd.lognz^4/(2*(m+1)),
             mean = exp(theta),
             lowCL = exp(theta - 1.96*sqrt(vartheta)),
             uppCL = exp(theta + 1.96*sqrt(vartheta))) %>%
      as.data.frame
    #head(out4,2)
    out4$method = "DLnMean"
    
    ## merge results into single data frame
    cpue.results = merge(out1, out2, all=TRUE)
    cpue.results = merge(cpue.results, out3, all=TRUE)
    cpue.results = merge(cpue.results, out4, all=TRUE)
    cpue.results = subset(cpue.results, select = c(season, method, mean, lowCL, uppCL))
    cpue.results = cpue.results[order(cpue.results$method, cpue.results$season), ]
    
    ymax = max(c(out1$uppCL, out2$uppCL, out3$uppCL, out4$uppCL),na.rm=T)
    if (ymax>20) ymax = ceiling(ymax/10)*10
    if (ymax>200) ymax = ceiling(ymax/100)*100
    
    if (draw.plot)
    {
      tc=0.1
      if(PaR==1)par(mfrow=c(1,1), mar=c(3,3,1,1),oma=c(1,1,1,1),mgp=c(1.75,.5,0))

      ## Plot mean and CI using ratio estimator
      plot(out1$season-1.5*tc, out1$mean, "o", pch=16, lty=2, col="black", 
           xlim=range(out1$season), ylim=c(0,ymax), las=1, xlab="Season", ylab=paste("CPUE (", cpue.units, ")", sep=""), 
           main=plot.title)
      arrows(x0=out1$season-1.5*tc, y0=out1$lowCL, 
             x1=out1$season-1.5*tc, y1=out1$uppCL, 
             code=3, angle=90, length=0.05, col="black")
      
      ## And assuming normal distribution...
      points(out2$season-0.5*tc, out2$mean, "o", pch=16, lty=2, col="grey")
      arrows(x0=out2$season-0.5*tc, y0=out2$lowCL, 
             x1=out2$season-0.5*tc, y1=out2$uppCL, 
             code=3, angle=90, length=0.05, col="grey")
      
      ## And assuming lognormal distribution...
      if(showLNMean=="YES")
      {
        points(out3$season+0.5*tc, out3$mean, "o", pch=16, lty=2, col="blue")
        arrows(x0=out3$season+0.5*tc, y0=out3$lowCL, 
               x1=out3$season+0.5*tc, y1=out3$uppCL, 
               code=3, angle=90, length=0.05, col="blue")
      }
      
      ## And assuming delta-lognormal distribution...
      points(out4$season+1.5*tc, out4$mean, "o", pch=16, lty=2, col="cyan")
      arrows(x0=out4$season+1.5*tc, y0=out4$lowCL, 
                           x1=out4$season+1.5*tc, y1=out4$uppCL, 
                           code=3, angle=90, length=0.05, col="cyan")
      if (show.legend)
      {
        if(showLNMean=="YES")
        {
          legend("topright", legend=c("Ratio", "Mean", "LnMean", "DLnMean"), 
               pch=16, lty=1, col=c("black", "grey", "blue", "cyan"), bty="n")
        }else
        {
          legend("topright", legend=c("Ratio", "Mean","DLnMean"), 
                        pch=16, lty=1, col=c("black", "grey", "cyan"), bty="n")
        }
      }
    }
  
    ## return cpue dataframe
    rownames(cpue.results) = NULL
    return(cpue.results)
}


##=================================================================================================
## END OF CODE
##=================================================================================================
