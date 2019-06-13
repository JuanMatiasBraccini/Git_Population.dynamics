#Script for running MCMC in ADMB from R
#note: use -mcgrope n to add fatter tail to inverse hessian
fn.run.MCMC=function(MODEL,DATA,niter,thinning)
{
  args=paste("./",MODEL, " -ind ", DATA," -mcmc ",niter, " -mcsave ", thinning," -mcscale",sep="")
  system(args)
  args1=paste("./",MODEL, " -ind ", DATA," -mceval",sep="")
  system(args1)
}