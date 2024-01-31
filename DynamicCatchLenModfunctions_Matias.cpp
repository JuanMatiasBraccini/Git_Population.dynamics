#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double x_is_gt_bound_cpp(const double x, const double bound, const double slope) {
  // Returns 0 if x << bound and 1 if x >> bound.
  // The parameter, slope, determines the range over
  // which the return value changes. The function
  // used is the logistic function, where the return
  // value is 0.5 when x == bound and 0.95
  // when x95 == bound + log(19)/slope
  double result;
  result = 1.0 / (1.0 + exp(slope * (bound - x)));
  
  return result;
}

// [[Rcpp::export]]
double x_is_lt_bound_cpp(const double x, const double bound, const double slope) {
  // Returns 0 if x >> bound and 1 if x << bound.
  // The parameter, slope, determines the range over
  // which the return value changes. The function
  // used is the logistic function, where the return
  // value is 0.5 when x == bound and 0.05
  // when x05 == bound + log(19)/slope.
  double result;
  result = 1.0 - 1.0 / (1.0 + exp(slope * (bound - x)));
  
  return result;
}

// [[Rcpp::export]]
double posfun_cpp(const double x, const double eps, const double pen) {
  // Assume eps = log(19)/slope
  double slope;
  double bound;
  double result;
  slope = log(19.0) / eps;
  bound = eps;  // Check whether it is sufficient to set this to zero?
  result = x *  x_is_gt_bound_cpp(x, bound, slope);
  
  return result;
}

// [[Rcpp::export]]
double penfun_cpp(const double x, const double eps) {
  // Assume eps = log(19)/slope
  double slope;
  double bound;
  double result;
  slope = log(19.0) / eps;
  bound = eps;  // Check whether it is sufficient to set this to zero?
  result = (x-eps) * (x-eps) * x_is_lt_bound_cpp(x, bound, slope);
  
  return result;
}


// [[Rcpp::export]]
NumericVector CalcMeanSizeAtAge_cpp(const double Lo, const double Linf, 
                                     const double vbK, const NumericVector Ages, 
                                     const int nAges) {
  NumericVector EstLenAtAge(nAges);
  
  EstLenAtAge=(Lo+(Linf-Lo)*(1-exp(-vbK*Ages)));
  // std::cout << "EstLenAtAge " << EstLenAtAge << std::endl;
  
  return EstLenAtAge;
}

// [[Rcpp::export]]
NumericVector CalcSizeDistOfRecruits_cpp(const double CVLenAtAge, const NumericVector EstLenAtAge,
                                         const NumericVector ubnd, const NumericVector lbnd,
                                         const int nLenCl) {
  
  double MeanLenRec;
  double SDAgeOneRecruits;
  NumericVector RecLenDist(nLenCl);

  MeanLenRec = EstLenAtAge(0);
  SDAgeOneRecruits = MeanLenRec * CVLenAtAge;
  
  // calculate the probability of different sizes, given the mean and SD, assuming
  // lengths are normally distributed
  RecLenDist = pnorm(ubnd, MeanLenRec, SDAgeOneRecruits) -
    pnorm(lbnd, MeanLenRec, SDAgeOneRecruits);
  //std::cout << "RecLenDist " << RecLenDist << std::endl;
  
  return (RecLenDist);
}

// [[Rcpp::export]]
NumericMatrix CalcLTM_cpp(const double Linf, const double vbK, const double CVLenAtAge, 
                          const NumericVector midpt, const NumericVector ubnd, const NumericVector lbnd,
                          const int nLenCl) {
  
  NumericMatrix LTM(nLenCl, nLenCl);
  NumericVector MeanEndingLength(nLenCl);
  NumericVector SDAtLen(nLenCl);
  NumericVector temp(nLenCl);
  int i;
  int ii;
  double tempsum;

  MeanEndingLength = midpt + (Linf - midpt) * (1 - exp(-vbK));
  SDAtLen = MeanEndingLength * CVLenAtAge;
    
    for (ii=0; ii<nLenCl; ii++) { // starting length class, upper bound - lower bound
      
      temp = pnorm(ubnd, MeanEndingLength(ii), SDAtLen(ii)) - 
        pnorm(lbnd, MeanEndingLength(ii), SDAtLen(ii));
      
      tempsum = 0;
      for (i=0; i<nLenCl; i++) { // starting length class
        LTM(i,ii) = temp(i);
        tempsum = tempsum + LTM(i,ii);
      }
      // ensure each row of the LTM sums to 1
      for (i=0; i<nLenCl; i++) { // ending length class
        LTM(i,ii) = LTM(i,ii) / tempsum;
      } // i
    }
    
    return LTM;
}

// [[Rcpp::export]]
NumericVector WtatLen_cpp(const double wtlen_a, const double wtlen_b, 
                          const NumericVector midpt, const int nLenCl) {
  
  //weight-length relationship
  NumericVector WtatLen(nLenCl);

  //WtatLen = (wtlen_a * midpt ^ wtlen_b) / 1000;  
  WtatLen = (wtlen_a * pow(midpt, wtlen_b)) / 1000;  
  // std::cout << "WtatLen " << WtatLen << std::endl;
  return WtatLen;
}

// [[Rcpp::export]]
NumericVector MatAtLen_cpp(const double MatL50, const double MatL95, 
                           const NumericVector midpt, const int nLenCl) {
  // maturity at length relationship
  NumericVector MatAtLen(nLenCl);
  
  MatAtLen = 1 / (1 + exp(-log(19) * (midpt - MatL50) / (MatL95 - MatL50)));
  // std::cout << "MatAtLen " << MatAtLen << std::endl;
  return MatAtLen;
}



// [[Rcpp::export]]
List CalcUnFishedFemBiomassPerRecruit_cpp(const int MaxAge, const int nLenCl, const NumericVector midpt, 
                                       const double NatMort, const NumericMatrix LTM, const NumericVector WtatLen, 
                                       const NumericVector MatAtLen, const NumericVector RecLenDist, 
                                       const double PropFemAtBirth) {
  
  // per recruit numbers surviving after natural mortality
  int t;
  int i;
  int ii;
  NumericVector Unfish_FemNPerRec(nLenCl); // numbers at length, across all ages
  NumericMatrix Unfish_FemNPerRecAtAge(MaxAge,nLenCl);
  NumericMatrix Unfish_FemSurvPerRecAtAge(MaxAge,nLenCl);
  double Unfish_FemSpBiomPerRec;
  
  for (t=0; t<MaxAge; t++) {
    if (t==0) {
      for (i=0; i<nLenCl; i++) {
        Unfish_FemNPerRecAtAge(t,i) = RecLenDist(i) * PropFemAtBirth;
      } // for
    } // t = 0
    
    if (t>0) {
      for (i=0; i<nLenCl; i++) {
        if (t < MaxAge-1) { // survival
          Unfish_FemSurvPerRecAtAge(t,i) = Unfish_FemNPerRecAtAge(t-1,i) * exp(-NatMort);
        } // if
        if (t == MaxAge-1) {
          Unfish_FemSurvPerRecAtAge(t,i) = Unfish_FemNPerRecAtAge(t-1,i) * exp(-NatMort) / 
            (1 - exp(-NatMort));
        } // if
      } // for
      
      // apply growth
      for (ii=0; ii<nLenCl; ii++) {  // starting length class
        for (i=0; i<nLenCl; i++) { // ending length class
          Unfish_FemNPerRecAtAge(t,i) = Unfish_FemNPerRecAtAge(t,i) + 
            (Unfish_FemSurvPerRecAtAge(t,ii) * LTM(i,ii));
        } // for
      } // for
    } // t > 0
  } // t
  
  for (t=0; t<MaxAge; t++) {
    for (i=0; i<nLenCl; i++) {
      Unfish_FemNPerRec(i) = Unfish_FemNPerRec(i) + Unfish_FemNPerRecAtAge(t,i);
    }
  }
  
  // unfished biomass per recruit (kg)
  Unfish_FemSpBiomPerRec = sum(Unfish_FemNPerRec * WtatLen * MatAtLen);
  
  return List::create(Named("Unfish_FemSpBiomPerRec") = Unfish_FemSpBiomPerRec, 
                      Named("Unfish_FemNPerRec") = Unfish_FemNPerRec,
                      Named("Unfish_FemNPerRecAtAge") = Unfish_FemNPerRecAtAge);

} 

// [[Rcpp::export]]
List CalcFishedFemBiomassPerRecruit_cpp(const int MaxAge, const int nLenCl, const NumericVector midpt, 
                                       const double NatMort, const NumericMatrix LTM, const NumericVector WtatLen, 
                                       NumericVector const MatAtLen, NumericVector const RecLenDist,
                                       NumericVector const SelAtLength, const double Init_F, const double PropFemAtBirth) {
  
  // per recruit numbers surviving after natural mortality
  int t;
  int i;
  int ii;
  NumericVector Fish_FemNPerRec(nLenCl); // numbers at length, across all ages
  NumericMatrix Fish_FemNPerRecAtAge(MaxAge,nLenCl);
  NumericMatrix Fish_FemSurvPerRecAtAge(MaxAge,nLenCl);
  double Fish_FemSpBiomPerRec;
  double temp_ZAtLen;
  
  
  for (t=0; t<MaxAge; t++) {
    if (t==0) {
      for (i=0; i<nLenCl; i++) {
        Fish_FemNPerRecAtAge(t,i) = RecLenDist(i) * PropFemAtBirth;
      } // for
    } // t = 0
    
    if (t>0) {
      for (i=0; i<nLenCl; i++) {
        temp_ZAtLen = NatMort + (SelAtLength(i) * Init_F);
        if (t < MaxAge-1) { // survival
          Fish_FemSurvPerRecAtAge(t,i) = Fish_FemNPerRecAtAge(t-1,i) * exp(-temp_ZAtLen);
        } // if
        if (t == MaxAge-1) {
          Fish_FemSurvPerRecAtAge(t,i) = Fish_FemNPerRecAtAge(t-1,i) * exp(-temp_ZAtLen) / 
            (1 - exp(-temp_ZAtLen));
        } // if
      } // for
      
      // apply growth
      for (ii=0; ii<nLenCl; ii++) {  // starting length class
        for (i=0; i<nLenCl; i++) { // ending length class
          Fish_FemNPerRecAtAge(t,i) = Fish_FemNPerRecAtAge(t,i) + 
            (Fish_FemSurvPerRecAtAge(t,ii) * LTM(i,ii));
        } // for
      } // for
    } // t > 0
  } // t
  
  for (t=0; t<MaxAge; t++) {
    for (i=0; i<nLenCl; i++) {
      Fish_FemNPerRec(i) = Fish_FemNPerRec(i) + Fish_FemNPerRecAtAge(t,i);
    }
  }
  
  // fished biomass per recruit (kg)
  Fish_FemSpBiomPerRec = sum(Fish_FemNPerRec * WtatLen * MatAtLen);
  
  return List::create(Named("Fish_FemSpBiomPerRec") = Fish_FemSpBiomPerRec, 
                      Named("Fish_FemNPerRec") = Fish_FemNPerRec,
                      Named("Fish_FemNPerRecAtAge") = Fish_FemNPerRecAtAge);

} 

// [[Rcpp::export]]
List CalcStockRecruitParams_cpp(const double Unfish_FemSpBiomPerRec, const double InitRec, 
                            const double Steepness) {
  // calculate B-H stock recruitment parameters
  double SR_alpha;
  double SR_beta;
  
  SR_alpha = (Unfish_FemSpBiomPerRec / InitRec) * ((1 - Steepness) / (4 * Steepness));
  
  SR_beta = (Steepness - 0.2) / (0.8 * Steepness * InitRec);
  
  return List::create(Named("SR_alpha") = SR_alpha, 
                      Named("SR_beta") = SR_beta);
}

// [[Rcpp::export]]
List CalcAnnCatch_cpp(int t, const double tempF, const double NatMort, NumericMatrix N, 
                        const NumericVector SelAtLength, const NumericVector WtatLen, const int nLenCl) {
  // calculate annual catch, given value of F
  
  int i;
  NumericVector temp_FAtLen;
  NumericVector temp_ZAtLen;
  NumericVector tempCatch_N(nLenCl);
  double tempAnnCatch;
  
  
  tempAnnCatch = 0;
  temp_FAtLen = SelAtLength * tempF;
  temp_ZAtLen = NatMort + temp_FAtLen;
  for (i=0; i<nLenCl; i++) {
    tempCatch_N(i) = N(t,i) * (temp_FAtLen(i) / temp_ZAtLen(i)) * (1 - exp(-temp_ZAtLen(i)));
    tempAnnCatch = tempAnnCatch + (tempCatch_N(i) * WtatLen(i));
  }
  
  return List::create(Named("tempCatch_N") = tempCatch_N, 
                      Named("tempAnnCatch") = tempAnnCatch);
}

// [[Rcpp::export]]
List CalcAnnualF_cpp(const int t, NumericVector const SelAtLength, NumericMatrix N, 
                 const NumericVector WtatLen, NumericVector ObsAnnCatch, 
                 const int nLenCl, const double NatMort, double FMort_pen) {
  
 
  // Use Newton's algorithm to calculate annual fishing mortality for the current year
  // by modifying F to match the expected catch to the observed catch
  int r;
  double tempObsCatch;
  double tempF;
  double eps;
  double ExpCatch1;
  double ExpCatch2;
  double Deriv;
  double NewF;
  double tempAnnCatch;
  double AnnCatch;
  NumericVector temp_FAtLen;
  NumericVector temp_ZAtLen;
  NumericVector tempCatch_N;
  List AnnCatchRes;
  
  // get observed catch for the year
  tempObsCatch = ObsAnnCatch(t);
  
  // penalty calculation to ensure fishing mortality does not get too high.
  // The expected catch is calculated for each year, from a very high value of F (0.99 y-1)
  // and if the observed catch for that year is greater than this value, the objective
  // function is penalized, signaling that the initial recruitment is too low. If so,
  // the observed catch for that year is reset to ensure the model is stable.
  tempF=0.99;
  temp_FAtLen = SelAtLength * tempF;
  temp_ZAtLen = NatMort + temp_FAtLen;
  
  AnnCatchRes=CalcAnnCatch_cpp(t, tempF, NatMort, N, SelAtLength, WtatLen, nLenCl);
  tempAnnCatch = as<double>(AnnCatchRes["tempAnnCatch"]);
  AnnCatch = tempAnnCatch;
  
  eps = 0.001;
  FMort_pen = FMort_pen + penfun_cpp(AnnCatch - tempObsCatch, eps);
  tempObsCatch = AnnCatch - posfun_cpp(AnnCatch - tempObsCatch, 0.01, FMort_pen);
  
  //std::cout << "CalcAnnualF_cpp: t " << t << " ObsAnnCatch(t) " << ObsAnnCatch(t) << 
  //  " tempObsCatch " << tempObsCatch << std::endl;
 
  tempF = 0.05;
  for (r=0; r<5; r++) {
    if (r > 0)  {
      tempF = NewF;
    } // if
    
    AnnCatchRes=CalcAnnCatch_cpp(t, tempF, NatMort, N, SelAtLength, WtatLen, nLenCl);
    tempAnnCatch = as<double>(AnnCatchRes["tempAnnCatch"]);
    ExpCatch1 = tempAnnCatch;
    
    tempF = tempF + 0.00000001;
    AnnCatchRes=CalcAnnCatch_cpp(t, tempF, NatMort, N, SelAtLength, WtatLen, nLenCl);
    tempAnnCatch = as<double>(AnnCatchRes["tempAnnCatch"]);
    ExpCatch2 = tempAnnCatch;
    
    Deriv = (ExpCatch2 - ExpCatch1) / 0.00000001;
    NewF = tempF - ((ExpCatch1 - tempObsCatch) / Deriv);
    //std::cout << "CalcAnnualF_cpp: t " << t << " tempObsCatch " << tempObsCatch << 
    //  " ExpCatch1 " << ExpCatch1 << " ExpCatch2 " << ExpCatch2 << " NewF " << NewF << std::endl;
    
  } // r
 
  tempCatch_N = as<NumericVector>(AnnCatchRes["tempCatch_N"]);

  return List::create(Named("NewF") = NewF, 
                      Named("tempAnnCatch") = tempAnnCatch,
                      Named("tempCatch_N") = tempCatch_N,
                      Named("FMort_pen") = FMort_pen);
}

// lnRecDev = params(3);
// Rec_Var = lnSigmaR * lnSigmaR;

// [[Rcpp::export]]
List AssessmentModel_cpp(NumericVector params, const int nYears, const NumericVector ObsAnnCatch, 
                        const int MaxAge, const int nLenCl, const NumericVector midpt, 
                        const NumericMatrix LTM, const NumericVector WtatLen, 
                        const NumericVector MatAtLen, const NumericVector RecLenDist, 
                        const double Init_F, const double PropFemAtBirth, const NumericVector SelAtLength,
                        const double lnSigmaR) {
  
  // reset penalty for estimated parameters, which penalises objective
  // function if parameters are outside specified bounds
  int i;
  int ii;
  int t;
  double param_pen;
  double FMort_pen;
  double NatMort_pen;
  double Steepness_pen;
  double InitRec;
  double NatMort;
  double Steepness;
  double eps;
  double lwbound;
  double uppbound;
  double temp_parm;
  double Unfish_FemSpBiomPerRec;
  double Fish_FemSpBiomPerRec;
  double SR_alpha;
  double SR_beta;
  double NewF;
  double tempAnnCatch;
  double lnRecDev;
  NumericMatrix N(nYears+1,nLenCl);
  NumericMatrix Surv(nYears+1,nLenCl);
  NumericMatrix CatchN(nYears,nLenCl); // numbers by size for all years
  NumericVector tempCatch_N(nLenCl); // numbers by size, for current year
  NumericVector AnnRecruit(nYears);
  NumericVector FemSpBiom(nYears);
  NumericVector VulnPopnBiom(nYears);
  NumericVector TotNum(nYears); // sum across length classes for each timestep
  NumericVector Catch_Biom(nYears); 
  NumericVector Ann_FMort(nYears);
  NumericVector Fish_FemNPerRec(nLenCl); // numbers at length, across all ages
  NumericVector temp_FAtLen(nLenCl);
  NumericVector temp_ZAtLen(nLenCl);
  
  List UnfishBPRres;
  List FishBPRres;
  List SRres;
  List AnnFres;
  
  param_pen = 0;
  NatMort_pen = 0;
  Steepness_pen = 0;
  FMort_pen = 0;
  
  // back log-transform parameters
  InitRec = exp(params(0));
  
  // constrain NatMort to values between 0.01 and 0.99
  eps = 0.001;
  lwbound = 0.01;
  uppbound = 0.99;
  temp_parm = exp(params[1]);
  NatMort_pen = NatMort_pen + penfun_cpp(uppbound - temp_parm, eps);
  temp_parm = uppbound - posfun_cpp(uppbound - temp_parm, eps, NatMort_pen);
  NatMort_pen = NatMort_pen + penfun_cpp(temp_parm - lwbound, eps);
  temp_parm = lwbound + posfun_cpp(temp_parm - lwbound, eps, NatMort_pen);
  NatMort = temp_parm;
  param_pen = param_pen + NatMort_pen;
  
  // constrain steepness to values between 0.2 and 1.0
  lwbound = 0.21;
  uppbound = 0.99;
  temp_parm = exp(params[2]);
  Steepness_pen = Steepness_pen + penfun_cpp(uppbound - temp_parm, eps);
  temp_parm = uppbound - posfun_cpp(uppbound - temp_parm, eps, Steepness_pen);
  Steepness_pen = Steepness_pen + penfun_cpp(temp_parm - lwbound, eps);
  temp_parm = lwbound + posfun_cpp(temp_parm - lwbound, eps, Steepness_pen);
  Steepness = temp_parm;
  param_pen = param_pen + Steepness_pen;
  
  lnRecDev = params(3);
  
  // calc unfished and fished biomass per recruit
  UnfishBPRres=CalcUnFishedFemBiomassPerRecruit_cpp(MaxAge, nLenCl, midpt, NatMort, LTM, 
                                                    WtatLen, MatAtLen, RecLenDist, PropFemAtBirth);
  Unfish_FemSpBiomPerRec = as<double>(UnfishBPRres["Unfish_FemSpBiomPerRec"]);
  
  FishBPRres=CalcFishedFemBiomassPerRecruit_cpp(MaxAge, nLenCl, midpt, NatMort, LTM, 
                                                WtatLen, MatAtLen, RecLenDist, SelAtLength, Init_F, PropFemAtBirth);
  Fish_FemSpBiomPerRec = as<double>(FishBPRres["Fish_FemSpBiomPerRec"]);   
  Fish_FemNPerRec = as<NumericVector>(FishBPRres["Fish_FemNPerRec"]); 
  
  // stock recruitment params
  SRres=CalcStockRecruitParams_cpp(Unfish_FemSpBiomPerRec, InitRec, Steepness);
  SR_alpha = as<double>(SRres["SR_alpha"]);  
  SR_beta = as<double>(SRres["SR_beta"]);  
  
  for (t=0; t<nYears; t++) {
    
    // get numbers and female spawning biomass for first time step
    if (t == 0) {
      AnnRecruit(t) = (Fish_FemSpBiomPerRec - SR_alpha) / (SR_beta * Fish_FemSpBiomPerRec) * 
           exp(lnRecDev - 0.5 * lnSigmaR * lnSigmaR);
      for (i=0; i<nLenCl; i++) {
        N(t,i) = InitRec * Fish_FemNPerRec(i); // 1000s
        FemSpBiom(t) = FemSpBiom(t) + (N(t,i) * MatAtLen(i) * WtatLen(i)); // tonnes
      }
    } 
    
    // calculate recruitment
    if (t > 0) {
      for (i=0; i<nLenCl; i++) {
        FemSpBiom(t) = FemSpBiom(t) + (N(t,i) * MatAtLen(i) * WtatLen(i));
      }
      // annual recruitment
      AnnRecruit(t) = FemSpBiom(t) / (SR_alpha + SR_beta * FemSpBiom(t)) * 
        exp(lnRecDev - 0.5 * lnSigmaR * lnSigmaR); 
    }
    
    // get fishing mortality
    AnnFres=CalcAnnualF_cpp(t, SelAtLength, N, WtatLen, ObsAnnCatch, 
                            nLenCl, NatMort, FMort_pen); // Newton's method
    
    NewF = as<double>(AnnFres["NewF"]);
    tempAnnCatch = as<double>(AnnFres["tempAnnCatch"]);
    tempCatch_N = as<NumericVector>(AnnFres["tempCatch_N"]);
    
    Ann_FMort(t) = NewF; // estimate of F
    temp_FAtLen = SelAtLength * Ann_FMort(t);
    temp_ZAtLen = NatMort + temp_FAtLen;
    
    for (i=0; i<nLenCl; i++) {
      CatchN(t,i) = tempCatch_N(i); // store catches (numbers by year, size)
      
      Surv(t+1,i) = N(t,i) * exp(-temp_ZAtLen(i)); // survival
    }
    
    // apply growth
    for (ii=0; ii<nLenCl; ii++) {  // starting length class
      for (i=0; i<nLenCl; i++) { // ending length class
        N(t+1,i) = N(t+1,i) + (Surv(t+1,ii) * LTM(i,ii));
      } // for
    } // for
    
    // add recruits
    for (i=0; i<nLenCl; i++) {
      N(t+1,i) = N(t+1,i) + (AnnRecruit(t) * RecLenDist(i));
      
      // get annual total population number
      TotNum(t) = TotNum(t) + N(t+1,i);
      
      // calculate popn biomass for timestep
      VulnPopnBiom(t) = VulnPopnBiom(t) + (N(t+1,i) * SelAtLength(i) * WtatLen(i));
    }
    Catch_Biom(t) = tempAnnCatch;
    //std::cout << "AssessmentModel_cpp: t " << t << " Catch_Biom(t) " << Catch_Biom(t) << std::endl;
  } // t
  
  // get estimate catch number and biomass
  FMort_pen = as<double>(AnnFres["FMort_pen"]);
  
  return List::create(Named("Unfish_FemSpBiomPerRec") = Unfish_FemSpBiomPerRec,
                      Named("Fish_FemSpBiomPerRec") = Fish_FemSpBiomPerRec,
                      Named("ObsAnnCatch") = ObsAnnCatch, 
                      Named("Catch_Biom") = Catch_Biom,
                      Named("CatchN") = CatchN,
                      Named("FemSpBiom") = FemSpBiom,
                      Named("AnnRecruit") = AnnRecruit,
                      Named("Ann_FMort") = Ann_FMort,
                      Named("param_pen") = param_pen,
                      Named("FMort_pen") = FMort_pen,
                      Named("NatMort_pen") = NatMort_pen,
                      Named("Steepness_pen") = Steepness_pen);
}

// [[Rcpp::export]]
double ObjectiveFunc_cpp(NumericVector params, const int nYears, const int nLen_SimYrs, 
                         NumericVector Len_SimYr, const NumericVector ObsLenCompVec1, 
                         const NumericVector ObsLenCompVec2, 
                         const NumericVector ObsLenCompVec3, 
                         const NumericVector ObsAnnCatch, const int MaxAge, const int nLenCl, 
                         const NumericVector midpt, const NumericMatrix LTM, const NumericVector WtatLen, 
                         const NumericVector MatAtLen, const NumericVector RecLenDist, 
                         const double Init_F, const double PropFemAtBirth, const NumericVector SelAtLength,
                         const double NatMort_mean, const double NatMort_sd,
                         const double Steepness_mean, const double Steepness_sd, const double lnSigmaR) {
  
  List OpModres;
  int y;
  int t;
  int i;
  double NLL;
  double temp;
  double MPen;
  double hPen;
  double RecPen;
  double NatMort;
  double NatMort_Var;
  double dev;
  double Steepness;
  double Steepness_Var;
  double lnRecDev;
  double Rec_Var;
  double param_pen;
  double FMort_pen;
  double ObjFunc;
  double tempsumCatch;
  NumericMatrix CatchN(nYears,nLenCl); // numbers by size for all years
  NumericVector ExpCatchAtLen(nLenCl);
  NumericVector ExpPropInLenClass(nLenCl);
  
  // do population dynamics
  OpModres=AssessmentModel_cpp(params, nYears, ObsAnnCatch, MaxAge, nLenCl, midpt, 
                              LTM, WtatLen, MatAtLen, RecLenDist, Init_F, PropFemAtBirth, 
                              SelAtLength, lnSigmaR);
  
  CatchN = as<NumericMatrix>(OpModres["CatchN"]);
  param_pen = as<double>(OpModres["param_pen"]);
  FMort_pen = as<double>(OpModres["FMort_pen"]);
  
  //multinomial likelihood associated with length composition(s)

  temp = 0;
  for (y=0; y<nLen_SimYrs; y++) {
    t = Len_SimYr(y) - 1;
    
    tempsumCatch = 0;
    for (i=0; i<nLenCl; i++) {
      tempsumCatch = tempsumCatch + CatchN(t,i);
    }
    
    for (i=0; i<nLenCl; i++) {
      ExpCatchAtLen(i) = CatchN(t,i);
      ExpPropInLenClass(i) = ExpCatchAtLen(i) / tempsumCatch;
      if (y==0) {
        temp = temp + (ObsLenCompVec1(i) * log(ExpPropInLenClass(i) + 1E-4));        
      }
      if (y==1) {
        temp = temp + (ObsLenCompVec2(i) * log(ExpPropInLenClass(i) + 1E-4));        
      }
      if (y==2) {
        temp = temp + (ObsLenCompVec3(i) * log(ExpPropInLenClass(i) + 1E-4));        
      }
    }
  }
  NLL = -temp;
  
  // calculate values associated with penalty priors for NatMort and Steepness
  MPen=0;
  NatMort = exp(params(1));
  NatMort_Var = NatMort_sd * NatMort_sd;
  dev = NatMort - NatMort_mean;
  MPen = 0.5 * log(NatMort_Var) + 0.5 * log(2*3.141592654) + (dev * dev) / (2 * NatMort_Var);

  hPen=0;
  Steepness = exp(params(2));
  Steepness_Var = Steepness_sd * Steepness_sd;
  dev = Steepness - Steepness_mean;
  hPen = 0.5 * log(Steepness_Var) + 0.5 * log(2*3.141592654) + (dev * dev) / (2 * Steepness_Var);

  RecPen=0;
  lnRecDev = params(3);
  Rec_Var = lnSigmaR * lnSigmaR;
  dev = lnRecDev - 0;
  RecPen = 0.5 * log(Rec_Var) + 0.5 * log(2*3.141592654) + (dev * dev) / (2 * Rec_Var);
  
  ObjFunc = NLL + (10.0 * param_pen) + (10.0 * MPen) + (10.0 * hPen) + (10.0 * RecPen) + (0.001 * FMort_pen);
  
  std::cout << "ObjFunc " << ObjFunc << " NLL " << NLL << " param_pen " << 10.0 * param_pen << 
     " MPen " << 10.0 * MPen << " hPen " << 10.0 * hPen << " RecPen " << 10.0 * RecPen << 
       " FMort_pen " << 0.001 * FMort_pen << std::endl;
   std::cout << " InitRec " <<  exp(params(0)) << " NatMort " <<  exp(params(1)) << 
     " Steepness " << exp(params(2)) << std::endl;
    
    // return List::create(Named("ObjFunc") = ObjFunc);
    return ObjFunc;
}


// [[Rcpp::export]]
List GetBiologyAndFisheryParams_cpp(const NumericVector Ages, const int nAges, const int nLenCl, const NumericVector midpt, 
                                    const NumericVector ubnd, const NumericVector lbnd, 
                                    const double Lo, const double Linf, const double vbK, 
                                    const double CVLenAtAge, const double  wtlen_a, const double wtlen_b, 
                                    const double MatL50, const double MatL95,
                                    const double SelL50, const double SelL95) {
  
  // These do not change as model is being fitted, so this function can be called once
  // prior to the model being fitted (rather then within the model, to (slightly) increase speed)
  List Selres;
  NumericVector EstLenAtAge(nAges);
  NumericVector RecLenDist(nLenCl);
  NumericMatrix LTM(nLenCl,nLenCl);
  NumericVector WtatLen(nLenCl);
  NumericVector MatAtLen(nLenCl);
  
  EstLenAtAge=CalcMeanSizeAtAge_cpp(Lo, Linf, vbK, Ages, nAges); // growth curve
  RecLenDist=CalcSizeDistOfRecruits_cpp(CVLenAtAge, EstLenAtAge, ubnd, lbnd, nLenCl); // size distribution of 1+ recruits
  LTM = CalcLTM_cpp(Linf, vbK, CVLenAtAge, midpt, ubnd, lbnd, nLenCl); // length transition matrix
  WtatLen = WtatLen_cpp(wtlen_a, wtlen_b, midpt, nLenCl); // weight-length relationship
  MatAtLen = MatAtLen_cpp(MatL50, MatL95, midpt, nLenCl); // maturity
  return List::create(Named("EstLenAtAge") = EstLenAtAge,
                      Named("RecLenDist") = RecLenDist,
                      Named("LTM") = LTM, 
                      Named("WtatLen") = WtatLen,
                      Named("MatAtLen") = MatAtLen,
                      Named("MatAtLen") = MatAtLen);
}


