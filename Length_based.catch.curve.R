CalcMeanSizeAtAge <- function(Lo,Linf, vbK) {
  EstLenAtAge=(Lo+(Linf-Lo)*(1-exp(-vbK*Ages)))
  return(EstLenAtAge)
}

CalcSizeDistOfRecruits <- function(GrowthCurveResults, CVLenAtAge) {
  # Mean length and SD of 0+ recruits, at 1 year of age
  MeanLenAtAge = GrowthCurveResults
  MeanLenRec <- MeanLenAtAge[1]
  SDAgeOneRecruits = MeanLenRec * CVLenAtAge
  
  RecLenDist = rep(0,nLenCl)
  RecLenDist = pnorm(ubnd, mean=MeanLenRec, sd=SDAgeOneRecruits, lower.tail = T) -
    pnorm(lbnd, mean=MeanLenRec, sd=SDAgeOneRecruits,lower.tail = T)
  
  RecLenDist = RecLenDist / sum(RecLenDist)
  
  return(RecLenDist)
}

CalcLTM <- function(Linf, vbK, CVLenAtAge, midpt) {
  # set up data frame for length-transition matrix
  LTM <- data.frame(matrix(ncol = nLenCl, nrow = nLenCl)) # to, from
  colnames(LTM) <- seq(1,nLenCl) 
  
  for (ii in 1:nLenCl) {  # starting length class
    
    LenClMidPt <- midpt[ii]
    MeanEndingLength = LenClMidPt + (Linf - LenClMidPt) * 
      (1 - exp(-vbK)) # expected 
    # cat("ii",ii,"midpt[ii",midpt[ii],"MeanEndingLength",MeanEndingLength,'\n')
    
    tempsum = 0
    for (i in 1:nLenCl) { # ending length class
      
      temp_SD <- MeanEndingLength * CVLenAtAge
      # temp_SD = SDLenAtAge
      
      if (i == 1) {
        # x+1 is upper bound of current 1 mm length class
        LTM[i,ii] = pnorm(ubnd[i], mean=MeanEndingLength,
                          sd=temp_SD,lower.tail = T)
      } # if
      else {
        # upper bound - lower bound
        tempUpperbnd = pnorm(ubnd[i], mean=MeanEndingLength,
                             sd=temp_SD,lower.tail = T)
        tempLowerbnd = pnorm(lbnd[i], mean=MeanEndingLength,
                             sd=temp_SD,lower.tail = T)
        LTM[i,ii] = tempUpperbnd - tempLowerbnd
        
      } # else
      
      tempsum = tempsum + LTM[i,ii]
      
    } #i
    
    # ensure each row of the LTM sums to 1
    for (i in 1:nLenCl) { # ending length class
      LTM[i,ii] = LTM[i,ii] / tempsum
    } #i
  } #ii
  
  return(LTM)
}

CatchCurveModel <- function(params) {
  
  FishMort = exp(params[1])
  # cat("CatchCurveModel: FishMort",FishMort,'\n')
  
  # per recruit numbers surviving after natural mortality
  Fish_FemNPerRec <- rep(0,nLenCl)
  Fish_FemNPerRecAtAge <- data.frame(matrix(nrow = MaxAge, ncol = nLenCl))
  colnames(Fish_FemNPerRecAtAge) <- midpt
  Fish_FemNPerRecAtAge[1:length(Ages),1:nLenCl] = 0
  Fish_FemSurvPerRecAtAge = Fish_FemNPerRecAtAge
  Catch = Fish_FemNPerRecAtAge
  CatchAtLen = rep(0,nLenCl)
  
  # get length distribution of recruits
  Fish_FemNPerRecAtAge[1,] = RecLenDist
  Fish_FemNPerRec = Fish_FemNPerRec + Fish_FemNPerRecAtAge[1,]
  
  # fishing mortality, calcuated based on selectivity associated with research fishing
  FAtLen = SelAtLength * FishMort
  
  # total mortality experienced by the stock, due to commercial fishing
  ZAtLen = NatMort + FAtLen
  
  # calculate catch at length for timestep (in numbers)
  Catch[1,] = Fish_FemNPerRecAtAge[1,] * (FAtLen / ZAtLen) * (1 - exp(-ZAtLen))  
  CatchAtLen = CatchAtLen + Catch[1,]
  
  # apply mortality to calculate survival
  for (t in 2:MaxAge) {
    if (t < MaxAge) {
      Fish_FemSurvPerRecAtAge[t,] = Fish_FemNPerRecAtAge[t-1,] * exp(-ZAtLen)
    } else {
      Fish_FemSurvPerRecAtAge[t,] = Fish_FemNPerRecAtAge[t-1,] * exp(-ZAtLen) / 
        (1 - exp(-ZAtLen))
    }
    
    # apply growth - looping - slower
    # for (ii in 1:nLenCl) {  # starting length class
    #   for (i in 1:nLenCl) { # ending length class
    #     Fish_FemNPerRecAtAge[t,i] = Fish_FemNPerRecAtAge[t,i] + (Fish_FemSurvPerRecAtAge[t,ii] * LTM[i,ii])
    #   }
    # }
    
    # apply growth - matrix multiplication - faster
    M1=t(t(Fish_FemSurvPerRecAtAge[t,]))
    M2=t(LTM[,])
    Fish_FemNPerRecAtAge[t,] = as.numeric(M1 %*% M2)
    Fish_FemNPerRec = Fish_FemNPerRec + Fish_FemNPerRecAtAge[t,]
    
    # calculate catch at length for timestep (in numbers)
    Catch[t,] = Fish_FemNPerRecAtAge[t,] * (FAtLen / ZAtLen) * (1 - exp(-ZAtLen))  
    
    CatchAtLen = CatchAtLen + Catch[t,]
    
  } # t
  
  if (!is.numeric(sum(Catch))) {
    cat("Problem - OperatingModel. sum(Catch) not numeric")
  }
  
  Res=list(Catch=Catch,
           CatchAtLen=CatchAtLen,
           Fish_FemNPerRec = Fish_FemNPerRec)
  
  return(Res)
}

ObjectiveFunc=function(params) {
  
  # to speed up, have commented these out - but must all be known in global variables
  # MeanSizeAtAge = CalcMeanSizeAtAge(Lo,Linf,vbK)
  # RecLenDist = CalcSizeDistOfRecruits(MeanSizeAtAge, CVLenAtAge)
  # SelAtLengthResults = CalcGillnetSelectivity(Theta1, Theta2, MSHES, midpt)
  # SelAtLength = SelAtLengthResults$Sel.combined
  # LTM = CalcLTM(Linf, vbK, CVLenAtAge, midpt)
  
  CatchCurveResults = CatchCurveModel(params)
  
  CatchAtLen = CatchCurveResults$CatchAtLen
  ExpCatchAtLen = CatchAtLen / sum(CatchAtLen)
  NLL = CalcMarginalNLL(ObsCatchFreqAtLen, ExpCatchAtLen)
  cat("exp(params)",exp(params),"NLL",NLL,"\n")
  return(NLL)
  
}

CalcMarginalNLL <- function(ObsCatchFreqAtLen, ExpCatchAtLen) {
  
  # calculate marginal length composition likelihood
  ExpPropInLenClass = ExpCatchAtLen / sum(ExpCatchAtLen)
  MarLengthNLL = - sum(ObsCatchFreqAtLen * log(ExpPropInLenClass + 1E-4))
  
  return(MarLengthNLL)
}

GetExpCatchAtLen <- function(params) {
  
  MeanSizeAtAge = CalcMeanSizeAtAge(Lo,Linf,vbK)
  RecLenDist = CalcSizeDistOfRecruits(MeanSizeAtAge, CVLenAtAge)
  #SelAtLengthResults = CalcGillnetSelectivity(Theta1, Theta2, MSHES, midpt)
  #SelAtLength = SelAtLengthResults$Sel.combined
  LTM = CalcLTM(Linf, vbK, CVLenAtAge, midpt)
  
  CatchCurveResults = CatchCurveModel(params)
  # CatchCurveResults = CatchCurveModel(params, Linf, vbK, tzero, theta1, theta2, nMeshes, 
  #                                     CVLenAtAge, MaxAge, nLenCl, midpt, MeanSizeAtAge, 
  #                                     RecLenDist, SelAtLength, LTM)
  
  CatchAtLen = CatchCurveResults$CatchAtLen
  ExpCatchAtLen = CatchAtLen / sum(CatchAtLen)
  return(ExpCatchAtLen)
  
}
################################



# Length-based catch curve analysis
# Alex Hesp, 14/9/2020
# For Matias Braccini - Exploring alternative assessment methods for data poor species

############
caca=FALSE
if(caca)
{
  # calculate mean size at age ------------------------------------------------------------------------
  # 2 par
  CalcMeanSizeAtAge <- function(Lo,Linf, vbK) {
    EstLenAtAge=(Lo+(Linf-Lo)*(1-exp(-vbK*Ages)))
    return(EstLenAtAge)
  }
  
  
  
  
  # 3 par
  # CalcMeanSizeAtAge <- function(Linf, vbK, tzero) {
  #   EstLenAtAge = Linf * (1 - exp (-vbK * (Ages - tzero)))
  #   return(EstLenAtAge)
  # }
  
  # calculate gillnet selectivity Kirkwood and Walker (1986) model--------------------------------------------
  # CalcGillnetSelectivity <- function(theta1, theta2, nMeshes, nLenCl, midpt) {
  #   
  #   # calculate alpha and beta parameters for selectivity schedule
  #   alpha_beta = rep(NA,nMeshes)
  #   alpha_beta = MeshSize_mm * theta1
  #   beta = -0.5 * (alpha_beta - sqrt(alpha_beta * alpha_beta + 4 * theta2))
  #   alpha = alpha_beta / beta 
  #   
  #   SelAtLengthForMesh = data.frame(matrix(nrow=nMeshes, ncol=nLenCl))
  #   colnames(SelAtLengthForMesh) = midpt
  #   
  #   SumSelAtLengthForMesh = rep(0,nLenCl)
  #   for (m in 1:nMeshes) {
  #     SelAtLengthForMesh[m,] = ((midpt / alpha_beta[m]) ^ alpha[m]) * 
  #       exp(alpha[m] - midpt / beta[m])
  #     
  #     SumSelAtLengthForMesh = SumSelAtLengthForMesh + SelAtLengthForMesh[m,]
  #   }
  #   SelAtLength = SumSelAtLengthForMesh / max(SumSelAtLengthForMesh)
  #   
  #   SelectResults = list(SelAtLengthForMesh = SelAtLengthForMesh,
  #                        SelAtLength = SelAtLength)
  #   
  #   return(SelectResults)
  #   
  # }
  
  # calculate size distribution of recruits--------------------------------------------
  CalcSizeDistOfRecruits <- function(GrowthCurveResults, CVLenAtAge) {
    
    # Mean length and SD of 0+ recruits, at 1 year of age
    MeanLenAtAge = GrowthCurveResults
    MeanLenRec <- MeanLenAtAge[1]
    SDAgeOneRecruits = MeanLenRec * CVLenAtAge
    
    RecLenDist = rep(0,nLenCl)
    RecLenDist = pnorm(ubnd, mean=MeanLenRec, sd=SDAgeOneRecruits, lower.tail = T) -
      pnorm(lbnd, mean=MeanLenRec, sd=SDAgeOneRecruits,lower.tail = T)
    
    RecLenDist = RecLenDist / sum(RecLenDist)
    
    return(RecLenDist)
    
  }
  
  # Length transition matrix--------------------------------------------
  CalcLTM <- function(Linf, vbK, CVLenAtAge, midpt) {
    
    
    # set up data frame for length-transition matrix
    LTM <- data.frame(matrix(ncol = nLenCl, nrow = nLenCl)) # to, from
    colnames(LTM) <- seq(1,nLenCl) 
    
    ii = 5
    for (ii in 1:nLenCl) {  # starting length class
      
      LenClMidPt <- midpt[ii]
      MeanEndingLength = LenClMidPt + (Linf - LenClMidPt) * 
        (1 - exp(-vbK)) # expected 
      # cat("ii",ii,"midpt[ii",midpt[ii],"MeanEndingLength",MeanEndingLength,'\n')
      
      tempsum = 0
      for (i in 1:nLenCl) { # ending length class
        
        temp_SD <- MeanEndingLength * CVLenAtAge
        # temp_SD = SDLenAtAge
        
        if (i == 1) {
          # x+1 is upper bound of current 1 mm length class
          LTM[i,ii] = pnorm(ubnd[i], mean=MeanEndingLength,
                            sd=temp_SD,lower.tail = T)
        } # if
        else {
          # upper bound - lower bound
          tempUpperbnd = pnorm(ubnd[i], mean=MeanEndingLength,
                               sd=temp_SD,lower.tail = T)
          tempLowerbnd = pnorm(lbnd[i], mean=MeanEndingLength,
                               sd=temp_SD,lower.tail = T)
          LTM[i,ii] = tempUpperbnd - tempLowerbnd
          
        } # else
        
        tempsum = tempsum + LTM[i,ii]
        
      } #i
      
      # ensure each row of the LTM sums to 1
      for (i in 1:nLenCl) { # ending length class
        LTM[i,ii] = LTM[i,ii] / tempsum
      } #i
    } #ii
    
    return(LTM)
    
  }
  
  # Catch Curve Model--------------------------------------------
  CatchCurveModel <- function(params, Linf, vbK, tzero, theta1, theta2, nMeshes, 
                              CVLenAtAge, MaxAge, nLenCl, midpt, MeanSizeAtAge, 
                              RecLenDist, SelAtLength, LTM) {
    
    FishMort = exp(params[1])
    # cat("CatchCurveModel: FishMort",FishMort,'\n')
    
    # per recruit numbers surviving after natural mortality
    Fish_FemNPerRec <- rep(0,nLenCl)
    Fish_FemNPerRecAtAge <- data.frame(matrix(nrow = MaxAge, ncol = nLenCl))
    colnames(Fish_FemNPerRecAtAge) <- midpt
    Fish_FemNPerRecAtAge[1:length(Ages),1:nLenCl] = 0
    Fish_FemSurvPerRecAtAge = Fish_FemNPerRecAtAge
    Catch = Fish_FemNPerRecAtAge
    CatchAtLen = rep(0,nLenCl)
    
    # get length distribution of recruits
    # RecLenDist = CalcSizeDistOfRecruits(MeanSizeAtAge, CVLenAtAge)
    Fish_FemNPerRecAtAge[1,] = RecLenDist
    Fish_FemNPerRec = Fish_FemNPerRec + Fish_FemNPerRecAtAge[1,]
    
    # fishing mortality, calcuated based on selectivity associated with research fishing
    FAtLen = SelAtLength * FishMort
    
    # total mortality experienced by the stock, due to commercial fishing
    ZAtLen = NatMort + FAtLen
    
    # calculate catch at length for timestep (in numbers)
    Catch[1,] = Fish_FemNPerRecAtAge[1,] * (FAtLen / ZAtLen) * (1 - exp(-ZAtLen))  
    CatchAtLen = CatchAtLen + Catch[1,]
    
    # apply mortality to calculate survival
    for (t in 2:MaxAge) {
      if (t < MaxAge) {
        Fish_FemSurvPerRecAtAge[t,] = Fish_FemNPerRecAtAge[t-1,] * exp(-ZAtLen)
      } else {
        Fish_FemSurvPerRecAtAge[t,] = Fish_FemNPerRecAtAge[t-1,] * exp(-ZAtLen) / 
          (1 - exp(-ZAtLen))
      }
      
      # apply growth - looping - slower
      # for (ii in 1:nLenCl) {  # starting length class
      #   for (i in 1:nLenCl) { # ending length class
      #     Fish_FemNPerRecAtAge[t,i] = Fish_FemNPerRecAtAge[t,i] + (Fish_FemSurvPerRecAtAge[t,ii] * LTM[i,ii])
      #   }
      # }
      
      # apply growth - matrix multiplication - faster
      M1=t(t(Fish_FemSurvPerRecAtAge[t,]))
      M2=t(LTM[,])
      Fish_FemNPerRecAtAge[t,] = as.numeric(M1 %*% M2)
      Fish_FemNPerRec = Fish_FemNPerRec + Fish_FemNPerRecAtAge[t,]
      
      # calculate catch at length for timestep (in numbers)
      Catch[t,] = Fish_FemNPerRecAtAge[t,] * (FAtLen / ZAtLen) * (1 - exp(-ZAtLen))  
      
      CatchAtLen = CatchAtLen + Catch[t,]
      
    } # t
    
    if (!is.numeric(sum(Catch))) {
      cat("Problem - OperatingModel. sum(Catch) not numeric")
    }
    
    Res=list(Catch=Catch,
             CatchAtLen=CatchAtLen,
             Fish_FemNPerRec = Fish_FemNPerRec)
    
    return(Res)
    
  }
  
  # multinomial likelihood--------------------------------------------
  CalcMarginalNLL <- function(ObsCatchFreqAtLen, ExpCatchAtLen) {
    
    # calculate marginal length composition likelihood
    ExpPropInLenClass = ExpCatchAtLen / sum(ExpCatchAtLen)
    MarLengthNLL = - sum(ObsCatchFreqAtLen * log(ExpPropInLenClass + 1E-4))
    
    return(MarLengthNLL)
    
  }
  
  # catch curve objective function--------------------------------------------
  ObjectiveFunc <- function(params) {
    # cat("ObjectiveFunc: params",params,'\n')
    
    #MeanSizeAtAge = CalcMeanSizeAtAge(Lo,Linf,vbK)
    #MeanSizeAtAge = CalcMeanSizeAtAge(Linf, vbK, tzero)
    RecLenDist = CalcSizeDistOfRecruits(MeanSizeAtAge, CVLenAtAge)
    #SelAtLengthResults = CalcGillnetSelectivity(theta1, theta2, nMeshes, nLenCl, midpt)
    #SelAtLength = SelAtLengthResults$SelAtLength
    LTM = CalcLTM(Linf, vbK, CVLenAtAge, midpt)
    
    CatchCurveResults = CatchCurveModel(params, Linf, vbK, tzero, theta1, theta2, nMeshes, 
                                        CVLenAtAge, MaxAge, nLenCl, midpt, MeanSizeAtAge, 
                                        RecLenDist, SelAtLength, LTM)
    
    CatchAtLen = CatchCurveResults$CatchAtLen
    ExpCatchAtLen = CatchAtLen / sum(CatchAtLen)
    NLL = CalcMarginalNLL(ObsCatchFreqAtLen, ExpCatchAtLen)
    cat("NLL",NLL,"\n")
    return(NLL)
    
  }
  
  # Get expected catch at length--------------------------------------------
  GetExpCatchAtLen <- function(params) {
    
    #MeanSizeAtAge =CalcMeanSizeAtAge(Lo,Linf,vbK)
    #MeanSizeAtAge = CalcMeanSizeAtAge(Linf, vbK, tzero)
    RecLenDist = CalcSizeDistOfRecruits(MeanSizeAtAge, CVLenAtAge)
    #SelAtLengthResults = CalcGillnetSelectivity(theta1, theta2, nMeshes, nLenCl, midpt)
    #SelAtLength = SelAtLengthResults$SelAtLength
    LTM = CalcLTM(Linf, vbK, CVLenAtAge, midpt)
    
    CatchCurveResults = CatchCurveModel(params, Linf, vbK, tzero, theta1, theta2, nMeshes, 
                                        CVLenAtAge, MaxAge, nLenCl, midpt, MeanSizeAtAge, 
                                        RecLenDist, SelAtLength, LTM)
    
    CatchAtLen = CatchCurveResults$CatchAtLen
    ExpCatchAtLen = CatchAtLen / sum(CatchAtLen)
    return(ExpCatchAtLen)
    
  }
}  
