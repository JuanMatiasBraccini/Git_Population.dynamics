if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Smart_par.R"))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Population dynamics/fn.fig.R"))
Do.jpeg="YES"
Do.tiff="NO"
fun.cons.like.mat=function(TAB,CX,Species)
{
  TAB$Risk=TAB$Consequence*TAB$Likelihood
  TAB$Clr=with(TAB,ifelse(Risk>0 & Risk<=2,"deepskyblue3",
        ifelse(Risk>2 & Risk<=4,"darkolivegreen3",
        ifelse(Risk>4 & Risk<=8,"yellow",
        ifelse(Risk>8 & Risk<=12 & Consequence <4,"orange",
        ifelse(Risk>=12 & Consequence ==4,"red",'transparent'))))))
  
  plot(TAB$Likelihood,TAB$Consequence,xaxt='n',yaxt='n',ylab="",xlab="",cex=CX,pch=19,
       col=TAB$Clr,xlim=c(0.25,4.25),ylim=c(0.25,4.25))
  axis(1,1:4,F)
  axis(2,1:4,F)
  text(TAB$Likelihood,TAB$Consequence,TAB$Risk,cex=2,col="black")
  legend("topright",Species,cex=1.5,bty='n')
}