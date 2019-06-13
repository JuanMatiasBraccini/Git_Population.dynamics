#Function for calculating Pearson's residuals
Pearson.res=function(OBS,PRED,Ns,BINS,YRS,YLAB,scaler,CL.pos,CL.neg,SeqMin,SeqMax)
{
  #remove years with no observations
  ID=rowSums(OBS)
  names(ID)=1:length(ID)
  ID=subset(ID,ID>0)
  ID=as.numeric(names(ID))
  
  ID.c=colSums(OBS)
  names(ID.c)=1:length(ID.c)
  ID.c=subset(ID.c,ID.c>0)
  ID.c=as.numeric(names(ID.c))
  
  OBS=OBS[ID,ID.c]
  PRED=PRED[ID,ID.c]
  Ns=Ns[ID]
  YRS=YRS[ID]
  BINS=BINS[ID.c]
  OBS[OBS<1e-4]=0
  OBS[OBS<=0]=NA
  Pear.res=(OBS - PRED)/sqrt(PRED*(1-PRED)/Ns) 
  Col.bg.mat=ifelse(Pear.res>0,CL.pos,CL.neg)
  
  par(mfcol=c(1,1),mar=c(4.5, 5,3, 1), xpd=TRUE,mgp = c(2.8,0.5, 0))
  plot(1,1,'n',xlim=range(YRS),ylim=range(BINS),ylab=YLAB,xaxt='n',yaxt='n',xlab='Financial year',bty='n',cex.lab=2)  
  
  for(n in 1:nrow(OBS)) points(rep(YRS[n],length(BINS)),BINS,pch=21,cex=scaler*abs(Pear.res[n,]),
                               col=1, bg=Col.bg.mat[n,],lwd=1.5)
  LEG=seq(SeqMin,SeqMax,2)
  legend('top',paste(c(-rev(LEG),LEG)),pch=21,col=1,pt.bg=c(rep(CL.neg,length(LEG)),rep(CL.pos,length(LEG))),
         bty='n',horiz=T,pt.cex=scaler*c(rev(LEG),LEG),cex=1.75,inset=c(0,-0.11))
  
  axis(1,at=YRS,las=3)
  axis(2,at=BINS,las=1,cex.axis=.85)
  box()  
}