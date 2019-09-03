#Script for creating table of input paramters and source
library(htmlTable)
library(xtable)
#source all input parameters
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/Organise input parameters.R")

specify_decimal <- function(x, k) format(round(x, k), nsmall=k)
is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
zeros.after.point=function(x)
{
  n=1:10
  n[(!floor(x*10^n)==0)][1]
}
                    
fn.par.ref.tbl=function(SP,add.growth.cv,add.Mrt.age)
{
  ParS=fn.input.pars(SP,add.growth.cv,add.Mrt.age)
  
  #remove empty objects
  cond <- sapply(ParS[[1]], function(x) length(x) > 0)
  ParS[[1]]=ParS[[1]][cond]
  ParS[[2]]=ParS[[2]][cond]
  
  a=ParS
  Nms=names(a$pars)
  N=length(Nms)
  Tabl=vector('list',N)
  names(Tabl)=Nms
  for(n in 1:N)
  {    
    Par=a[[1]][n]
    
    if(!Nms[n]%in%c("Reporting","Prop.males.in.ktch"))
    {
      
      Par=unlist(Par)
      Nm=names(Par)
      nn=length(Par)
      Par.s=paste(unique(unlist(a[[2]][n])), collapse = ",")
      RR=rep(NA,nn)
      dumy=data.frame(Parameter=RR,Value=RR,Comment=NA,Source=NA)
      for(r in 1:nn)
      {
        if(length(Nm)==1)
        {
          dumy$Parameter[r]=Nms[n]
          dumy$Comment[r]=  substr(Nm, nchar(Nms[n])+2, 100) 
        }
        if(length(Nm)>1) dumy$Parameter[r]=Nm[r]
        
        if(is.wholenumber(Par[r])) dumy$Value[r]=Par[r]
        if(!is.wholenumber(Par[r]))
        {
          zeros=zeros.after.point(Par[r])
          if(zeros<=3) dumy$Value[r]=specify_decimal(Par[r],3)
          if(zeros>3) dumy$Value[r]=specify_decimal(Par[r],zeros+2)
         
        }
          
        if(r==1)dumy$Source[r]=Par.s[r]
      }
      
    }
    
    if(Nms[n]=="Reporting")
    {
      Par=Par[[1]]
      Par1=c(Par$Zn1,Par$Zn2,Par$WC)
      names(Par1)=paste(as.character(Par$FinYear),rep(c("Zn1","Zn2","WC"),each=length(Par$FinYear)))

     Nm=names(Par1)
     nn=length(Par1)
     Par.s=paste(unique(unlist(a[[2]][n])), collapse = ",")
     RR=rep(NA,nn)
     dumy=data.frame(Parameter=RR,Value=RR,Comment=NA,Source=NA)
     for(r in 1:nn)
     {
       dumy$Parameter[r]=Nms[n]
       dumy$Comment[r]=Nm[r]
       dumy$Value[r]=specify_decimal(Par1[r],3)
       if(r==1)dumy$Source[r]=Par.s[r]
     }

    }
    if(Nms[n]=="Prop.males.in.ktch")
    {
      Par1=Par[[1]]
      Nm=names(Par1)
      nn=length(Par1)
      Par.s=paste(unique(unlist(a[[2]][n])), collapse = ",")
      RR=rep(NA,nn)
      dumy=data.frame(Parameter=RR,Value=RR,Comment=NA,Source=NA)
      for(r in 1:nn)
      {
        dumy$Parameter[r]=Nms[n]
        dumy$Comment[r]=Nm[r]
        dumy$Value[r]=specify_decimal(Par1[r],3)
        if(r==1)dumy$Source[r]=Par.s[r]
      }
      
    }
      
    Tabl[[n]]=dumy
  }
  Tabl=do.call(rbind,Tabl)
  Tabl[is.na(Tabl)]=""
  return(Tabl)
}