#----- R Monte Carlo simulation evaluations of X.tpl model-----#
#notes:   This code runs the X.tpl script to simulate data and fit model for testing  model precision and bias
#         The X.tpl must have a simulation() function
Simul.test=function(MODEL,n.sims,PARS,DATA)
{
  #create input pars
  x.log=paste(MODEL,"Simpar.log",sep="")
  if(file.exists(x.log)) file.remove(x.log)  
  write(paste(PARS,"\t",sep="",collapse=" "), file=x.log)
  
  #create random seeds
  seeds=seq(999, by=2, length=n.sims)
  args=paste(paste("./",MODEL, " -ind ", DATA," -est -sim ",sep=""), seeds, sep="")
  
  #run model with random seeds
  for(sim in args) system(sim)

  #Read in and plot X.log
  br=read.table(x.log, header=T)
  boxplot(br, ylim=c(-2, 2))
}