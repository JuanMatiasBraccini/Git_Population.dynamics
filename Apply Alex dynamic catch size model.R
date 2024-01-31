#--------------- Data-poor dynamic catch + size model --------------
#Developed by Alex Hesp, 1/7/2021

# The framework is that of a single-area, single sex, size-structured integrated
# model. Selectivity is assumed to be known. Estimated parameters
# include R0 (Initial unfished recruitment). To allow for some uncertainty in 
# stock productivity, natural mortality and steepness (i.e. of the Beverton-Holt 
# stock recruitment relationship) are 'estimated', constrained using normal,
# penalty prior function. To allow for further uncertainty associated with recruitment,
# a parameter describing the mean deviation in annual recruitment (i.e. from
# those estimated from the stock recruitment relationship) is also 'estimated'
# using a normal penalty prior function, assuming a mean of zero and specified value
# for the standard deviation of the natural logarithms of recruitment.

library(MASS)
library(tidyverse)
library("Rcpp")
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')
sourceCpp(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/DynamicCatchLenModfunctions_Matias.cpp"))


apply.Alex.catch.length=function(Init_F,InitRec,NatMort,NatMort_sd,Steepness,Steepness_sd,lnRecDev,
                                 Lo,Linf,vbK,CVLenAtAge,SDGrowthRec,
                                 MaxLen,LenInc,MaxAge, MatL50,MatL95,PropFemAtBirth,
                                 wtlen_a,wtlen_b,SelAtLength,SelL50, SelL95,
                                 Len_SimYr,n_SimLen,ObsLenComp,Katch,
                                 nsims,UnfishRec,lnSigmaR) 
{
  ObsAnnCatch=Katch$Tonnes
  
  params = c(log(InitRec),log(NatMort),log(Steepness),lnRecDev)
  NatMort_mean=NatMort
  Steepness_mean=Steepness
  
  # ages
  Ages = 1:MaxAge
  nAges = length(Ages)
  
  # length bins
  lbnd = seq(0,MaxLen - LenInc, LenInc)
  ubnd = lbnd + LenInc
  midpt = lbnd + (LenInc/2)
  nLenCl = length(midpt)
  
  n_SimYrs = length(unique(Len_SimYr))
  nLen_SimYrs = length(Len_SimYr)
  nYears = length(ObsAnnCatch)
  
  
  #objects for rcpp
  # for rccp
  ObsLenCompVec1=ObsLenCompVec2=ObsLenCompVec3=rep(0,nLenCl)
  ObsLenCompVec1 = as.numeric(unlist(ObsLenComp))
  if(exists('ObsLenComp2'))ObsLenCompVec2 = ObsLenComp2
  if(exists('ObsLenComp3'))ObsLenCompVec3 = ObsLenComp3
  
  
  #Fit model
  res=GetBiologyAndFisheryParams_cpp(Ages, nAges, nLenCl, midpt, ubnd, lbnd, 
                                     Lo, Linf, vbK, CVLenAtAge, wtlen_a, wtlen_b, 
                                     MatL50, MatL95, SelL50, SelL95)
  EstLenAtAge=res$EstLenAtAge
  RecLenDist=res$RecLenDist
  LTM=res$LTM
  WtatLen=res$WtatLen
  MatAtLen=res$MatAtLen
  ObjFunc <- function(params)
  {
    res=ObjectiveFunc_cpp(params, nYears, nLen_SimYrs, Len_SimYr, ObsLenCompVec1, ObsLenCompVec2, ObsLenCompVec3,
                          ObsAnnCatch, MaxAge, nLenCl, midpt, LTM, WtatLen, MatAtLen, RecLenDist, Init_F,
                          PropFemAtBirth, SelAtLength, NatMort_mean, NatMort_sd, Steepness_mean, Steepness_sd, lnSigmaR)
    return(res)
  }
  nlmb <- nlminb(params, ObjFunc, gradient = NULL, hessian = TRUE)
  
  #Get variance-covariance matrix from fitted model and distribution of estim pars
  hess.out = optimHess(nlmb$par, ObjFunc)
  vcov.Params = solve(hess.out)
  sims = data.frame(mvrnorm(n = nsims, nlmb$par, vcov.Params))
  names(sims) = c("ln_R0", "ln_M","ln_h","lnRecDev")
  Table.estimates=data.frame(Parameter=c("ln_R0", "ln_M","ln_h","lnRecDev"),
                             Median=nlmb$par,
                             SE = sqrt(diag(vcov.Params)))
  
  # observed vs expected length composition
  res=AssessmentModel_cpp(nlmb$par, nYears, ObsAnnCatch, MaxAge, nLenCl, midpt, 
                          LTM, WtatLen, MatAtLen, RecLenDist, Init_F, PropFemAtBirth, SelAtLength, lnSigmaR)
  Predicted.LenComp=Observed.LenComp=ObsLenComp
  for (i in 1:length(Len_SimYr))
  {
    Observed.LenComp[i,]=ObsLenComp[i,]/sum(ObsLenComp[i,])
    Predicted.LenComp[i,] = res$CatchN[Len_SimYr[i],] / sum(res$CatchN[Len_SimYr[i],])
  }
  
    
  #Get biomass and fishing mortality predictions
  BiomEst = data.frame(matrix(nrow=nsims,ncol=nYears))
  colnames(BiomEst) = 1:nYears
  FMortEst = BiomEst
  RelBiomEst = BiomEst
  
  for (i in 1:nsims)
  {
    param=as.vector(unlist(sims[i,1:4]))
    res=AssessmentModel_cpp(param, nYears, ObsAnnCatch, MaxAge, nLenCl, midpt, 
                            LTM, WtatLen, MatAtLen, RecLenDist, Init_F, PropFemAtBirth, SelAtLength, lnSigmaR)
    
    BiomEst[i,] = res$FemSpBiom
    FMortEst[i,] = res$Ann_FMort
    RelBiomEst[i,] = res$FemSpBiom / (res$Unfish_FemSpBiomPerRec * exp(param[1]))
  }
  
  rel.biom=data.frame(year=Katch$finyear,
                      median=apply(RelBiomEst[,], MARGIN=2, function(x) quantile(x, 0.5)),
                      upper.95=apply(RelBiomEst[,], MARGIN=2, function(x) quantile(x, 0.975)),
                      lower.95=apply(RelBiomEst[,], MARGIN=2, function(x) quantile(x, 0.025)))
  f.series=data.frame(year=Katch$finyear,
                      median=apply(FMortEst[,], MARGIN=2, function(x) quantile(x, 0.5)),
                      upper.95=apply(FMortEst[,], MARGIN=2, function(x) quantile(x, 0.975)),
                      lower.95=apply(FMortEst[,], MARGIN=2, function(x) quantile(x, 0.025)))
  
  return(list(rel.biom=rel.biom,f.series=f.series,Table.estimates=Table.estimates,
              Year=Katch$finyear,Observed.catch=ObsAnnCatch,Predicted.catch=res$Catch_Biom,
              Observed.LenComp=Observed.LenComp, Predicted.LenComp=Predicted.LenComp,
              midpt=midpt))
  
}
