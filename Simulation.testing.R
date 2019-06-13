#Funciton for Simulation testing

Simul.Eval=function(MODEL,Spec,PARS,n.sims,Jitr)
{
  #create simulation folder
  setPath(Scenarios[match(MODEL,Scenarios$Model),]$Model)
  A=getwd()
  if(!file.exists(file.path(A, "Simulation.Testing"))) dir.create(file.path(A, "Simulation.Testing"))  
  setwd(file.path(getwd(), "Simulation.Testing"))
  B=getwd()
  
  PINPARS=Pin.pars[[match(MODEL,names(Pin.pars))]]
  
  #copy files
  file.copy(paste(A,paste("/",Spec,".tpl",sep=""),sep=""), paste(Spec,".tpl",sep=""),overwrite =T)
  file.copy(paste(A,paste("/",Spec,".pin",sep=""),sep=""), paste(Spec,".pin",sep=""),overwrite =T)
  
  #create .dat
  #note: change phases so params are not estimated
  n=Inputs[[match(MODEL,names(Inputs))]]
  n$Phases[1]=-n$Phases[1]    #don't estimate
  n$Phases[2:length(n$Phases)]=-abs(n$Phases[2:length(n$Phases)])
  FILE=paste(Spec,".dat",sep="")
  nzones=n$nzone
  ModDims=unlist(c(yr.start,yr.end,nzones))
  Hdr="#Basic model dimensions (yr.start, yr.end, nzones)"
  write(Hdr,file = FILE)
  write(ModDims,file = FILE,sep = "\t",append=T)
  for(k in (length(ModDims)+1):length(n))
  {
    nn=n[[k]]
    if(is.data.frame(nn)|is.matrix(nn))
    {
      Hdr=paste("#",paste(c(names(n)[k],"(",names(nn),")"),collapse=' '))
      write(Hdr,file = FILE,append=T)      
      write.table(nn,file = FILE,row.names=F,col.names=F,append=T)
    }else
    {
      Hdr=paste("#",names(n)[k],sep='')
      write(Hdr,file = FILE,append=T)
      write(n[[k]],file = FILE,sep = "\t",append=T)
    }
  }
  
  #run .tpl
  args=paste(paste("./",Spec, " -ind ", paste(Spec,".dat",sep=""), " -est",sep=""), sep="")
  clean_admb(Spec)
  compile_admb(Spec,verbose=T)
  system(args)  
  
  #Extract CPUE with no param estimation
  Pred.CpuE=reptoRlist(paste(Spec,".rep",sep="")) 
  Pred.CpuE=Pred.CpuE$CPUE_out
  #Pred.CpuE=Pred.CpuE$Est_CPUE
  
  
  #run sim eval
  Tab=matrix(nrow=n.sims,ncol=length(PARS))
  colnames(Tab)=PARS
  for(nnnn in 1:n.sims)
  {
    #create pin file with jitter
    if(file.exists(paste(Spec,".pin",sep="")))file.remove(paste(Spec,".pin",sep=""))
    These.pars=sapply(PINPARS,fn.ji)
    par.nms=names(These.pars)
    FILE=paste(Spec,".pin",sep="")
    write("# Input parameters",file = FILE)
    for(k in 1:length(These.pars))
    {
      Hdr=paste("#",par.nms[k])
      write(Hdr,file = FILE,append=T)
      write(These.pars[k],file = FILE,sep = "\t",append=T)
    }
    
    #create sim data 
    #note: update CPUE
    n=Inputs[[match(MODEL,names(Inputs))]]
    FILE=paste(Spec,".dat",sep="")
    nzones=n$nzone
    
    if(is.numeric(Pred.CpuE)) cpue.plus.error=jitter(Pred.CpuE,Jitr)
    if(is.matrix(Pred.CpuE))cpue.plus.error=apply(Pred.CpuE, 2,function(x) jitter(x,Jitr))
    n$CPUE=abs(cpue.plus.error)
    ModDims=unlist(c(yr.start,yr.end,nzones))
    Hdr="#Basic model dimensions (yr.start, yr.end, nzones)"
    write(Hdr,file = FILE)
    write(ModDims,file = FILE,sep = "\t",append=T)
    for(k in (length(ModDims)+1):length(n))
    {
      nn=n[[k]]
      if(is.data.frame(nn)|is.matrix(nn))
      {
        Hdr=paste("#",paste(c(names(n)[k],"(",names(nn),")"),collapse=' '))
        write(Hdr,file = FILE,append=T)      
        write.table(nn,file = FILE,row.names=F,col.names=F,append=T)
      }else
      {
        Hdr=paste("#",names(n)[k],sep='')
        write(Hdr,file = FILE,append=T)
        write(n[[k]],file = FILE,sep = "\t",append=T)
      }
    }
    
    #run .tpl
    args=paste(paste("./",Spec, " -ind ", paste(Spec,".dat",sep=""), " -est",sep=""), sep="")
    clean_admb(Spec)
    compile_admb(Spec,verbose=T)
    system(args)  
    
    if(file.exists(paste(Spec,".par",sep="")))
    {      
      a=scan(paste(Spec,".par",sep=""), what="", sep="\n")
      Dat=as.numeric(a[which(!substr(a,1,1)=="#")])
      Tab[nnnn,]=c(Dat)
    }
    
  }
  return(Tab)
}
# Simul.Eval=function(MODEL,Spec,PARS,n.sims,Jitr)
# {
#   #create simulation folder
#   setPath(Scenarios[match(MODEL,Scenarios$Model),]$Model)
#   A=getwd()
#   
#   
#   if(!file.exists(file.path(A, "Simulation.Testing"))) dir.create(file.path(A, "Simulation.Testing"))  
#   setwd(file.path(getwd(), "Simulation.Testing"))
#   B=getwd()
#   
#   PHASES=Par.phases[[match(MODEL,names(Par.phases))]]
#   PINPARS=Pin.pars[[match(MODEL,names(Pin.pars))]]
#   
#   #copy files
#   file.copy(paste(A,paste("/",Spec,".tpl",sep=""),sep=""), paste(Spec,".tpl",sep=""),overwrite =T)
#   file.copy(paste(A,paste("/",Spec,".pin",sep=""),sep=""), paste(Spec,".pin",sep=""),overwrite =T)
#     
#   #extract predicted cpue
#     #change phases
#   fn.scan=function(What,Chng,Chng.to)
#   {
#   setwd(A)
#   y <- scan(paste(Spec,What,sep=""), what="", sep="\n")  
#   setwd(B)
#   if(file.exists(paste(Spec,What,sep="")))file.remove(paste(Spec,What,sep=""))
#   StRings=which(substr(y,1,1)=="#")
#   Dat=which(!substr(y,1,1)=="#")
#   D=vector('list',length(StRings))
#   for(k in 1:(length(D)-1))
#   {
#     if(diff(c(StRings[k],StRings[k+1]))==2) test=y[StRings[k]+1]
#     if(diff(c(StRings[k],StRings[k+1]))>2) test=y[seq(StRings[k]+1,StRings[k+1]-1)]
#     D[[k]]=as.numeric(unlist(strsplit(test, "[[:space:]]+")))
#   }
#   D[[length(D)]]=y[length(y)]
#   
#   D[[which(y[StRings]==Chng)]]=Chng.to  
#   
#   FILE=paste(Spec,".dat",sep="")
#   for(k in 1:length(D)) 
#   {
#     write(y[StRings][k],file = FILE,append=T)
#     write(D[[k]],file = FILE,sep = "\t",append=T)
#   }
#   
# }
#   fn.scan(What=".dat",Chng="#Phases",Chng.to=-PHASES)
# 
#   #run .tpl
#   args=paste(paste("./",Spec, " -ind ", paste(Spec,".dat",sep=""), " -est",sep=""), sep="")
#   clean_admb(Spec)
#   compile_admb(Spec,verbose=T)
#   system(args)  
#   Pred.CpuE=scan(paste(Spec,".rep",sep=""), what="", sep="\n")
#   Pred.CpuE=Pred.CpuE[match("Est_CPUE",Pred.CpuE)+1]
#   Pred.CpuE=as.numeric(unlist(strsplit(Pred.CpuE, "[[:space:]]+")))
#   Pred.CpuE=Pred.CpuE[-is.na(Pred.CpuE)]
# 
#   
#   #run sim eval
#   Tab=matrix(nrow=n.sims,ncol=length(PARS))
#   colnames(Tab)=PARS
#   for(n in 1:n.sims)
#   {
#     #create pin file with jitter
#     if(file.exists(paste(Spec,".pin",sep="")))file.remove(paste(Spec,".pin",sep=""))
#     These.pars=sapply(PINPARS,fn.ji)
#     par.nms=names(These.pars)
#     FILE=paste(Spec,".pin",sep="")
#     write("# Input parameters",file = FILE)
#     for(k in 1:length(These.pars))
#     {
#       Hdr=paste("#",par.nms[k])
#       write(Hdr,file = FILE,append=T)
#       write(These.pars[k],file = FILE,sep = "\t",append=T)
#     }
#     
#     
#     #create sim data 
#     cpue.plus.error=sapply(Pred.CpuE,function(x) jitter(x,Jitr))
#     fn.scan(What=".dat",Chng="#CPUE",Chng.to=cpue.plus.error)
#     
#     #run .tpl
#     args=paste(paste("./",Spec, " -ind ", paste(Spec,".dat",sep=""), " -est",sep=""), sep="")
#     clean_admb(Spec)
#     compile_admb(Spec,verbose=T)
#     system(args)  
#   
#     if(file.exists(paste(Spec,".par",sep="")))
#     {      
#       a=scan(paste(Spec,".par",sep=""), what="", sep="\n")
#       a=a[-1]      
#       Dat=as.numeric(a[which(!substr(a,1,1)=="#")])
#       Tab[n,]=c(Dat)
#     }
#   }
#   
#   return(Tab)
# }