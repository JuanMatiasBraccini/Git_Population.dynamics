library(AquaticLifeHistory)
library(BayesGrowth)
Geraghty.dusky=read.csv(handl_OneDrive('Data/Age and growth/Dusky FINAL AGE Final (9)_Geraghty et al 2013.csv'))
Geraghty.sandbar=read.csv(handl_OneDrive('Data/Age and growth/Sandbar FINAL AGE Final (9)_Geraghty et al 2013.csv'))
plot.Bayes.growth=function(gc,DD)
{
  p=gc%>%
    ggplot(aes(Age, LAA))+
    geom_point(data = DD, aes(Age, Length), alpha = .3)+
    geom_lineribbon(aes( ymin = .lower, ymax = .upper, fill = factor(.width)), size = .8) +
    labs(y = "Fork length (cm)", x = "Age (yrs)")+
    scale_fill_brewer(palette="BuPu", direction=-1,name = "Credibility interval")+
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(expand = c(0,0), breaks = seq(0,max(DD$Age),1))+
    theme_bw()+
    theme(text = element_text(size = 14),
          legend.position = c(0.8,0.325),
          legend.background = element_rect(colour = "transparent"))
  return(p)
}
add.Frequentist.plot=function(p,Estimates)
{
  Estimates=Estimates%>%
    mutate(LAA=AVG)
  p <- p+
    geom_line(data=Estimates, aes(Age, AVG, col = Model),size = 1) + 
    geom_line(data=Estimates, aes(Age, low, col = Model),size = 1,alpha=0.4,linetype = "dotted") +
    geom_line(data=Estimates, aes(Age, upp, col = Model),size = 1,alpha=0.4,linetype = "dotted")
#    geom_ribbon(data=Estimates,aes(ymin = low, ymax = upp, fill = Model), alpha = 0.4)  
  return(p)

}
fit.growth.curve=function(SP,dat,LH)
{
  for(i in 1:length(SP))
  {
    print(paste("refitting growth for ----",names(SP)[i]))
    ii=match(names(SP)[i],names(dat))
    dd=dat[[ii]]$age_length%>%
      rename(Length=FL)%>%
      mutate(Sex=ifelse(Sex=='Male','M',ifelse(Sex=='Female','F',NA)))%>%
      dplyr::select(Age,Length,Sex)
    birth_size=LH[[ii]]$Lzero
    birth_size_se=1
    max_size_se=5
    wd=paste(handl_OneDrive("Analyses/Population dynamics/1."),capitalize(names(SP)[i]),"/",AssessYr,sep='')
    if(names(SP)[i]%in%c("dusky shark","sandbar shark"))
    {
      if(names(SP)[i]=="dusky shark") add.geraghty=Geraghty.dusky
      if(names(SP)[i]=="sandbar shark") add.geraghty=Sandbar.dusky
      add.geraghty=add.geraghty%>%
                    filter(!Readability..amended.==1)%>%
                    dplyr::select(AgeAgree,FL,Sex)%>%
                    rename(Age=AgeAgree,
                           Length=FL)
      
      rbind(add.geraghty%>%mutate(Data.set='Geraghty'),dd%>%mutate(Data.set='WA'))%>%
        mutate(Age=as.numeric(Age))%>%
        ggplot(aes(Age,Length,color=Data.set))+
        geom_point()+
        facet_wrap(~Sex,ncol = 1)
      ggsave(paste(wd,"/1_Inputs/Visualise data/refit_growth_compare_WA_Geraghty.tiff",sep=''), 
             width = 8,height = 8, dpi = 300, compression = "lzw")

      ns=table(dd$Sex)
      n.fem=ns[match('F',names(ns))]
      n.mal=ns[match('M',names(ns))]
      
      dd.gera.f=add.geraghty%>%filter(Sex=='F')
      dd.gera.f=dd.gera.f[sample(1:nrow(dd.gera.f),n.fem,replace = T),]
      dd.gera.m=add.geraghty%>%filter(Sex=='M')
      dd.gera.m=dd.gera.m[sample(1:nrow(dd.gera.m),n.mal,replace = T),]
      
      dd=rbind(dd,dd.gera.f,dd.gera.m)%>%mutate(Age=as.numeric(Age))

    }
    #Females
    dd1=dd%>%filter(Sex=='F')
    max_size=max(dd1$Length)
    #frequentist
    Fem.growth=Estimate_Growth(dd1,Birth.Len=birth_size,plots = FALSE)
    #bayesian
    Fem.growth.B <- Estimate_MCMC_Growth(data = dd1%>%dplyr::select(Length,Age), 
                                         Model = "VB" ,iter = 5e3,Linf = max_size*1.1,Linf.se = max_size_se,
                                         L0 = birth_size,sigma.max = 100,L0.se = birth_size_se,k.max = 2)
    #plots and estimates
    p=plot.Bayes.growth(gc= Calculate_MCMC_growth_curve(Fem.growth.B, Model = "VB",max.age = max(dd1$Age), probs = c(.5,.95)),DD=dd1)
    add.Frequentist.plot(p=p,Estimates=Fem.growth$Estimates)
    ggsave(paste(wd,"/1_Inputs/Visualise data/refit_growth_female.tiff",sep=''), 
           width = 8,height = 8, dpi = 300, compression = "lzw")
    
    out=rbind(Fem.growth$VonB%>%data.frame%>%mutate(Model='VonB'),
              Fem.growth$Logistic%>%data.frame %>%mutate(Model='Logistic'),
              Fem.growth$Gompertz%>%data.frame%>%mutate(Model='Gompertz'))
    fit.summary=summary(Fem.growth.B)
    out=rbind(out,
              fit.summary$summary%>%
                data.frame%>%
                dplyr::select(mean,sd)%>%
                rename( Parameter=mean,SE=sd)%>%
                mutate(Model='Bayesian.VonB'))
    write.csv(out,paste(wd,"/1_Inputs/Visualise data/refit_growth_female.csv",sep=''))
    
    #Males
    dd1=dd%>%filter(Sex=='M')
    max_size=max(dd1$Length)
    #frequentist
    Male.growth=Estimate_Growth(dd1,Birth.Len=birth_size,plots = FALSE)
    #bayesian
    Male.growth.B <- Estimate_MCMC_Growth(data = dd1%>%dplyr::select(Length,Age), 
                                          Model = "VB" ,iter = 5e3,Linf = max_size*1.1,Linf.se = max_size_se,
                                          L0 = birth_size,sigma.max = 100,L0.se = birth_size_se,k.max = 2)
    #plots
    p=plot.Bayes.growth(gc=Calculate_MCMC_growth_curve(Male.growth.B, Model = "VB",max.age = max(dd1$Age), probs = c(.5,.95)),DD=dd1)
    add.Frequentist.plot(p=p,Estimates=Male.growth$Estimates)
    ggsave(paste(wd,"/1_Inputs/Visualise data/refit_growth_male.tiff",sep=''), 
           width = 8,height = 8, dpi = 300, compression = "lzw")
    
    out=rbind(Male.growth$VonB%>%data.frame%>%mutate(Model='VonB'),
              Male.growth$Logistic%>%data.frame %>%mutate(Model='Logistic'),
              Male.growth$Gompertz%>%data.frame%>%mutate(Model='Gompertz'))
    fit.summary=summary(Male.growth.B)
    out=rbind(out,
              fit.summary$summary%>%
                data.frame%>%
                dplyr::select(mean,sd)%>%
                rename( Parameter=mean,SE=sd)%>%
                mutate(Model='Bayesian.VonB'))
    write.csv(out,paste(wd,"/1_Inputs/Visualise data/refit_growth_male.csv",sep=''))
  }
}
fit.growth.curve(SP=Indicator.species,dat=Species.data,LH=List.sp)



