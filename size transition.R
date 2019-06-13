#Size transition matrix
#Guide: Haddon 2001 page 222
#Use Faben's version of vonB to estiamte average length increment for a given size-class.
# then use the statisical distribution chosen (gamma, normal, lognormal, etc) to describe
# how individuals would be distributed arond the average increment

size.trans.mat=function(Linf,K,SD,bin.low,bin.up)
{
  Size.bin=data.frame(low=bin.low,up=bin.up)
  Size.bin$mid=rowMeans(Size.bin)
  Size.bin$Pred.annual.inc=(Linf-Size.bin$mid)*(1-exp(-K))  #Faben's model
  Size.bin$Pred.annual.inc=with(Size.bin,ifelse(Pred.annual.inc<0,0,Pred.annual.inc))
  
  STM=matrix(ncol=length(Size.bin$mid),nrow=length(Size.bin$mid))
  colnames(STM)=rownames(STM)=Size.bin$mid
  
  for(j in 1:ncol(STM))
  {
    Size=as.numeric(colnames(STM)[j])
    for(k in 1:length(STM[,j]))
    {
      Size.next=as.numeric(rownames(STM)[k])
      dummy=0
      if(Size<Linf) if(Size.next>=Size) dummy=pnorm(Size.next-Size+1,Size.bin$Pred.annual.inc[j],SD)
      if(Size>=Linf & Size.next>=Size)  dummy=1
      STM[k,j]=dummy
    }
  }
  
  #normalise
  STM.n=STM
  nn=nrow(STM)
  STM.n[1,]=STM[1,]/STM[nn,]
  for(k in 2:nn) STM.n[k,]=(STM[k,]-STM[k-1,])/STM[nn,]
  
  return(STM.n)
}


#Sadovy et al 2007. Stock assessment approach for the Napoleon fish, Chilinus undulatus, in Indonesia.
Sadovy.size.trans.mat=function(Linf,K,SD,bin.low,bin.up,Truncate)
{
  Size.bin=data.frame(low=bin.low,up=bin.up)
  Mids=rowMeans(Size.bin)
  
  STM=matrix(ncol=length(Mids),nrow=length(Mids))
  colnames(STM)=rownames(STM)=Mids
  
  prob.density=function(l.start) exp(-(Mids-(Linf*(1-exp(-K))+(l.start*exp(-K))))^2/(2*SD^2))
  for(qq in 1:ncol(STM))
  {
    dummy=prob.density(Mids[qq])
    if(Truncate=="YES") dummy[1:(qq-1)]=0
    STM[,qq]=dummy
  }
    
   #normalise
  STM.n=STM
  for(qq in 1:ncol(STM)) STM.n[,qq]=STM.n[,qq]/sum(STM.n[,qq])
   
  return(STM.n)
}


