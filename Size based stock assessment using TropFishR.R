library(TropFishR)
data("synLFQ7")
#source: https://cran.r-project.org/web/packages/TropFishR/vignettes/tutorial.html

#ELEFAN (ELectronic LEngth Frequency ANalysis) methods to estimate Linf and K from LFQ data

  # set seed value for reproducible results
set.seed(1)

  # adjust bin size
synLFQ7a <- lfqModify(synLFQ7, bin_size = 4)

  # plot raw and restructured LFQ data
lfqbin <- lfqRestructure(synLFQ7a, MA = 5, addl.sqrt = FALSE)
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfqbin, Fname = "catch", date.axis = "modern")
plot(lfqbin, Fname = "rcounts", date.axis = "modern")
par(opar)
#For synLFQ7, a bin size of 4 cm and a moving average of 5 seems appropriate and will be used in the following

# Get a first estimate of Linf using the Powell-Wetherall method (Wetherall, Polovina, and Ralston 1987).
res_PW <- powell_wetherall(param = synLFQ7a,
                           catch_columns = 1:ncol(synLFQ7a$catch),
                           reg_int = c(10,28))
paste("Linf =",round(res_PW$Linf_est), "±", round(res_PW$se_Linf))

# get a quick K value corresponding to the Linf
# ELEFAN with K-Scan
res_KScan <- ELEFAN(synLFQ7a, Linf_fix = res_PW$Linf_est,
                    MA=5, addl.sqrt = TRUE, hide.progressbar = TRUE)

# show results
res_KScan$par; res_KScan$Rn_max


# Response surface analyss
res_RSA <- ELEFAN(synLFQ7a, Linf_range = seq(119,139,1), MA = 5,
                  K_range = seq(0.01,2,0.1), addl.sqrt = TRUE,
                  hide.progressbar = TRUE, contour=5)

# show results
res_RSA$par; res_RSA$Rn_max

# find 3 highest score values
n <- length(res_RSA$score_mat)
best_scores <- sort(res_RSA$score_mat,partial=n-0:2)[n-0:2]
ind <- arrayInd(which(res_RSA$score_mat %in% best_scores),
                dim(res_RSA$score_mat))
Ks <- as.numeric(rownames(res_RSA$score_mat)[ind[,1]])
Linfs <- as.numeric(colnames(res_RSA$score_mat)[ind[,2]])

res_loop <- vector("list", 3)
for(i in 1:3){
  tmp <- ELEFAN(synLFQ7a,
                Linf_range = seq(Linfs[i]-2, Linfs[i]+2, 0.2),
                K_range = seq(Ks[i]-0.1, Ks[i]+0.1, 0.05),
                MA = 5,
                addl.sqrt = TRUE,
                hide.progressbar = TRUE,
                contour=5)
  res_loop[[i]] <- cbind(Rn_max=tmp$Rn_max, t(as.matrix(tmp$par)))
}
results <- do.call(rbind, res_loop)


# run ELEFAN with simulated annealing
res_SA <- ELEFAN_SA(synLFQ7a, SA_time = 60*0.5, SA_temp = 6e5,
                    MA = 5, seasonalised = TRUE, addl.sqrt = FALSE,
                    init_par = list(Linf = 129, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                    low_par = list(Linf = 119, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                    up_par = list(Linf = 139, K = 1, t_anchor = 1, C = 1, ts = 1))
# show results
res_SA$par; res_SA$Rn_max


JK <- vector("list", length(synLFQ7a$dates))
for(i in 1:length(synLFQ7a$dates)){
  loop_data <- list(dates = synLFQ7a$dates[-i],
                    midLengths = synLFQ7a$midLengths,
                    catch = synLFQ7a$catch[,-i])
  tmp <- ELEFAN_SA(loop_data, SA_time = 60*0.5, SA_temp = 6e5,
                   MA = 5, addl.sqrt = TRUE,
                   init_par = list(Linf = 129, K = 0.5, t_anchor = 0.5, C=0.5, ts = 0.5),
                   low_par = list(Linf = 119, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                   up_par = list(Linf = 139, K = 1, t_anchor = 1, C = 1, ts = 1),
                   plot = FALSE)
  JK[[i]] <- unlist(c(tmp$par,list(Rn_max=tmp$Rn_max)))
}
JKres <- do.call(cbind, JK)
# mean
JKmeans <- apply(as.matrix(JKres), MARGIN = 1, FUN = mean)
# confidence intervals
JKconf <- apply(as.matrix(JKres), MARGIN = 1, FUN = function(x) quantile(x, probs=c(0.025,0.975)))
JKconf <- t(JKconf)
colnames(JKconf) <- c("lower","upper")


# run ELEFAN with genetic algorithm
res_GA <- ELEFAN_GA(synLFQ7a, MA = 5, seasonalised = TRUE, maxiter = 50, addl.sqrt = FALSE,
                    low_par = list(Linf = 119, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                    up_par = list(Linf = 139, K = 1, t_anchor = 1, C = 1, ts = 1),
                    monitor = FALSE)
# show results
res_GA$par; res_GA$Rn_max

# plot LFQ and growth curves
plot(lfqbin, Fname = "rcounts",date.axis = "modern", ylim=c(0,130))
lt <- lfqFitCurves(synLFQ7a, par = list(Linf=123, K=0.2, t_anchor=0.25, C=0.3, ts=0),
                   draw = TRUE, col = "grey", lty = 1, lwd=1.5)
# lt <- lfqFitCurves(synLFQ7, par = res_RSA$par,
#                    draw = TRUE, col = "goldenrod1", lty = 1, lwd=1.5)
lt <- lfqFitCurves(synLFQ7a, par = res_SA$par,
                   draw = TRUE, col = "darkblue", lty = 1, lwd=1.5)
lt <- lfqFitCurves(synLFQ7a, par = res_GA$par,
                   draw = TRUE, col = "darkgreen", lty = 1, lwd=1.5)

# assign estimates to the data list
synLFQ7a <- c(synLFQ7a, res_SA$par)
class(synLFQ7a) <- "lfq"


# estimation of M
Ms <- M_empirical(Linf = res_SA$par$Linf, K_l = res_SA$par$K, method = "Then_growth")
synLFQ7a$M <- as.numeric(Ms)
# show results
paste("M =", as.numeric(Ms))
#> [1] "M = 0.286"

# summarise catch matrix into vector and add plus group which is smaller than Linf
synLFQ7b <- lfqModify(synLFQ7a, vectorise_catch = TRUE, plus_group = 118)
# run catch curve
res_cc <- catchCurve(synLFQ7b, reg_int = c(8,26), calc_ogive = TRUE)
# assign estimates to the data list
synLFQ7b$Z <- res_cc$Z
synLFQ7b$FM <- as.numeric(synLFQ7b$Z - synLFQ7b$M)
synLFQ7b$E <- synLFQ7b$FM/synLFQ7b$Z


###Stock size and status ####
synLFQ7c <- synLFQ7b

# assign length-weight parameters to the data list
synLFQ7c$a <- 0.015
synLFQ7c$b <- 3
# run CA
vpa_res <- VPA(param = synLFQ7c, terminalF = synLFQ7c$FM,
               analysis_type = "CA",
               plot=TRUE, catch_corFac = (1+4/12))
# stock size
sum(vpa_res$annualMeanNr, na.rm =TRUE) / 1e3
#> [1] 342.5927
# stock biomass
sum(vpa_res$meanBiomassTon, na.rm = TRUE)
#> [1] 1186801
# assign F per length class to the data list
synLFQ7c$FM <- vpa_res$FM_calc


# Thompson and Bell model with changes in F
TB1 <- predict_mod(synLFQ7c, type = "ThompBell",
                   FM_change = seq(0,1.5,0.05),  stock_size_1 = 1,
                   curr.E = synLFQ7c$E, plot = FALSE, hide.progressbar = TRUE)
# Thompson and Bell model with changes in F and Lc
TB2 <- predict_mod(synLFQ7c, type = "ThompBell",
                   FM_change = seq(0,1.5,0.1), Lc_change = seq(25,50,0.1),
                   stock_size_1 = 1,
                   curr.E = synLFQ7c$E, curr.Lc = res_cc$L50,
                   s_list = list(selecType = "trawl_ogive",
                                 L50 = res_cc$L50, L75 = res_cc$L75),
                   plot = FALSE, hide.progressbar = TRUE)
# plot results
par(mfrow = c(2,1), mar = c(4,5,2,4.5), oma = c(1,0,0,0))
plot(TB1, mark = TRUE)
mtext("(a)", side = 3, at = -1, line = 0.6)
plot(TB2, type = "Isopleth", xaxis1 = "FM", mark = TRUE, contour = 6)
mtext("(b)", side = 3, at = -0.1, line = 0.6)

# Biological reference levels
TB1$df_Es
#>   Fmsy F05     Emsy       E05
#> 1 0.45 0.2 0.611413 0.4115226
# Current yield and biomass levels
TB1$currents
#>   curr.Lc curr.tc    curr.E    curr.F    curr.C   curr.Y curr.V   curr.B
#> 1      NA      NA 0.4449351 0.2292551 0.2749782 1253.466      0 8781.828
