library(ldbounds)
library(tidyverse)
library(xtable)
library("scales")


########################### Global simulation settings
N<-500 #per treatment group, observations come in in pairs
rep<-100000
sigma<-1
alpha<-0.05
eff<-c(0,0.1,0.2,0.3,0.4)


#################################### SEQUENTIAL INFERENCE FRAMEWORKS

######################### EPPO version of Generalized Always Valid Inference (GAVI)
GAVI<-function(i=500,alpha=0.05,phi=500,sigma2=1, est=1){
  rho<-phi/(
    log(log(exp(1)*alpha^(-2)))-2*log(alpha)
  )
  
  ci<-(sqrt(2*sigma2)/i) * sqrt(
    (i+rho)*log(
      (i+rho)/(rho*alpha^2)
    )
  )
  is_sig<-(est-ci)>0
  return(is_sig)
}

######################### Netflix version of Always Valid F-test (mSPRT)
av_analysis_no_cov = function(std_err_tau,std_res,estimate, alpha=0.05, phi=1){
  s = std_res
  z2 = (s/std_err)^2
  radius = std_err*sqrt(log((phi+z2)/(phi*alpha^2))*(phi + z2)/(z2))
  is_sig = (estimate-radius)>0
  return(is_sig)
}

######################## Simple Sequential Testing (SST)
SST <- function(n, sigma, alpha, est){
  return(est > qnorm(1-alpha/2) * sqrt(n) * sigma)
}

library(mvtnorm)
library(stats)

ell <- function(h, k, rho) {
  return (
    1 - (pnorm(h) + pnorm(k))
    + pmvnorm(upper=c(h, k), mean=c(0, 0), sigma=matrix(c(1, rho, rho, 1), nrow=2))
  )
}

compute_convolution <- function(s, sigma, b_1, b_2) {
  # Computes the convolution of a centered normal and a truncated normal at b_2.
  #
  # s - the std of the centered normal
  # sigma - the std of the untruncated "version" of the truncated normal
  # b_1 - the (upper) truncation (cutoff) point
  #
  # See:
  # 1. https://www.erikdrysdale.com/NTS/
  # 2. https://www.kss.or.kr/jounalDown.php?IDX=831
  zeta <- sqrt(sigma^2 + s^2)
  rho <- sigma / sqrt(sigma^2 + s^2)
  h <- b_2 / zeta
  k <- b_1 / sigma
  return (
    1 - (1 - pnorm(h) - ell(h, k, rho)) / pnorm(k)
  )
}

estimate_false_detection_rate_bound <- function(total_events, thresholds, std_signed_value) {
  stopifnot(length(total_events) > 0)
  stopifnot(length(total_events) == length(thresholds))
  
  false_detection_rate <- 2 * (1 - pnorm(thresholds[1] / sqrt(total_events[1]) / std_signed_value))
  
  for (i in 2:length(total_events)) {
    false_detection_rate <- false_detection_rate +
      2 * (pnorm(thresholds[i - 1] / sqrt(total_events[i - 1]) / std_signed_value) *
             (1 - compute_convolution(sqrt(total_events[i] - total_events[i - 1]) * std_signed_value,
                                      sqrt(total_events[i - 1]) * std_signed_value,
                                      thresholds[i - 1],
                                      thresholds[i])))
  }
  
  return(false_detection_rate)
}

compute_pSST_thresholds <- function(expected_cumulative_n, actual_cumulative_n, sigma, alpha){
  thresholds = qnorm(1-alpha/2) * sqrt(expected_cumulative_n) * sigma
  while (estimate_false_detection_rate_bound(expected_cumulative_n, thresholds, sigma) > alpha) {
    thresholds<-thresholds*1.005
  }
  cat("FPR bound =", estimate_false_detection_rate_bound(expected_cumulative_n, thresholds, sigma))
  cat("\nThresholds:", thresholds)
  th<-rep(NA, actual_cumulative_n[length(actual_cumulative_n)])
  th[1:actual_cumulative_n[1]] = thresholds[1]
  for (i in 2:length(actual_cumulative_n)) { 
    th[(actual_cumulative_n[i-1]+1):actual_cumulative_n[i]] = thresholds[i]
  }
  return(th)
}

pSST7 <- compute_pSST_thresholds(round((1:7)*(N/7)), round((1:7)*(N/7)), sqrt(2)*sigma, alpha)
pSST14 <- compute_pSST_thresholds(round((1:14)*(N/14)), round((1:14)*(N/14)), sqrt(2)*sigma, alpha)

######################### Group Sequential Test (GST) 
s14<-seq(1/14,1,1/14)
s28<-seq(1/28,1,1/28)
s42<-seq(1/42,1,1/42)
s56<-seq(1/56,1,1/56)
# For over sampling
s7<-seq(1/7,1,1/7)
s21<-seq(1/21,1,1/21)
# For under sampling
s84<-seq(1/84,1,1/84)
s112<-seq(1/112,1,1/112)


#Accurate sampling
gst1<-ldBounds(t=1,iuse=3, phi=1, alpha=alpha,sides=1)$upper.bounds
gst14<-ldBounds(t=s14,iuse=3, phi=1, alpha=alpha,sides=1)$upper.bounds
gst28<-ldBounds(t=s28,iuse=3, phi=1, alpha=alpha,sides=1)$upper.bounds
gst42<-ldBounds(t=s42,iuse=3, phi=1, alpha=alpha,sides=1)$upper.bounds
gst56<-ldBounds(t=s56,iuse=3, phi=1, alpha=alpha,sides=1)$upper.bounds

gst1phi3<-ldBounds(t=1,iuse=3, phi=3, alpha=alpha,sides=1)$upper.bounds
gst14phi3<-ldBounds(t=s14,iuse=3, phi=3, alpha=alpha,sides=1)$upper.bounds
gst28phi3<-ldBounds(t=s28,iuse=3, phi=3, alpha=alpha,sides=1)$upper.bounds
gst42phi3<-ldBounds(t=s42,iuse=3, phi=3, alpha=alpha,sides=1)$upper.bounds
gst56phi3<-ldBounds(t=s56,iuse=3, phi=3, alpha=alpha,sides=1)$upper.bounds


# Under sampling - Calculate bounds as if there would be twice as many checks as it will be
gst14u<-ldBounds(t=s28,iuse=3, phi=1, alpha=alpha,sides=1)$upper.bounds[1:14]
gst28u<-ldBounds(t=s56,iuse=3, phi=1, alpha=alpha,sides=1)$upper.bounds[1:28]
gst42u<-ldBounds(t=s84,iuse=3, phi=1, alpha=alpha,sides=1)$upper.bounds[1:42]
gst56u<-ldBounds(t=s112,iuse=3, phi=1, alpha=alpha,sides=1)$upper.bounds[1:56]

gst14uphi3<-ldBounds(t=s28,iuse=3, phi=3, alpha=alpha,sides=1)$upper.bounds[1:14]
gst28uphi3<-ldBounds(t=s56,iuse=3, phi=3, alpha=alpha,sides=1)$upper.bounds[1:28]
gst42uphi3<-ldBounds(t=s84,iuse=3, phi=3, alpha=alpha,sides=1)$upper.bounds[1:42]
gst56uphi3<-ldBounds(t=s112,iuse=3, phi=3, alpha=alpha,sides=1)$upper.bounds[1:56]

#Over sampling - Calculate bounds using the over-sample trick suggested in literature
# Helpfunction GST over sample
iter_gstb<-function(actual=14, expected=7){
  counter<-0
  bound<-numeric()
  for(s in (actual-expected+1):actual){
    counter<-counter+1
    bound[counter]<-ldBounds(t=seq(1/s,1,1/s),iuse=3, phi=1, alpha=alpha,sides=1)$upper.bounds[s]
  }
  return(bound)
}

gst14o<-c(ldBounds(t=s7,iuse=3, phi=1, alpha=alpha, sides=1)$upper.bounds, iter_gstb(actual=14, expected=7))
gst28o<-c(ldBounds(t=s14,iuse=3, phi=1, alpha=alpha,sides=1)$upper.bounds, iter_gstb(actual=28, expected=14))
gst42o<-c(ldBounds(t=s21,iuse=3, phi=1, alpha=alpha,sides=1)$upper.bounds, iter_gstb(actual=42, expected=21))
gst56o<-c(ldBounds(t=s28,iuse=3, phi=1, alpha=alpha,sides=1)$upper.bounds, iter_gstb(actual=56, expected=28))

gst14ophi3<-c(ldBounds(t=s7,iuse=3, phi=3, alpha=alpha, sides=1)$upper.bounds, iter_gstb(actual=14, expected=7))
gst28ophi3<-c(ldBounds(t=s14,iuse=3, phi=3, alpha=alpha,sides=1)$upper.bounds, iter_gstb(actual=28, expected=14))
gst42ophi3<-c(ldBounds(t=s21,iuse=3, phi=3, alpha=alpha,sides=1)$upper.bounds, iter_gstb(actual=42, expected=21))
gst56ophi3<-c(ldBounds(t=s28,iuse=3, phi=3, alpha=alpha,sides=1)$upper.bounds, iter_gstb(actual=56, expected=28))

######################### STATSIG sequential test (CAA)
Z_crit_alpha = qnorm(1-alpha)
zcrit_ss<-Z_crit_alpha/((1:N)/N) # Statsig crit value


######################### Help function for finding stopping n
detN<-function(input, nseq,N){
  w<-which(input==TRUE)[1]
  return(ifelse(is.na(w), N, ceiling(nseq[w])))
}

###################################################################################
############################# Simulation 1 - Comparing FPR and Power of all methods

set.seed(8163)
nrows<-81
tib<-matrix(NA, rep*length(eff)*nrows,6)
counter<-1
for(r in 1:rep){
  for(e in 1:length(eff)){
    x<-rnorm(N,1,sigma)
    y<-rnorm(N,1+eff[e],sigma)
    diff<-(cumsum(y)/(1:N)-cumsum(x)/(1:N))
    std_res=sigma #for mSPRT approach (known res var)
    std_err = sqrt(sigma^2/(1:N) + sigma^2/(1:N))
    zdiff<-diff / std_err
    #Assume known variance to not have to deal with t-dist
    
    #Bonferroni
    Bonferroni_stream = qnorm(1-alpha/N)
    Bonferroni_14 = qnorm(1-alpha/14)
    Bonferroni_28 = qnorm(1-alpha/28)
    Bonferroni_42 = qnorm(1-alpha/42)
    Bonferroni_56 = qnorm(1-alpha/56)
    
    ### mSPRT AW Ftest
    sig_mSPRT_phi100<-av_analysis_no_cov(std_err_tau=std_err,std_res=1,estimate=diff, alpha=alpha*2, phi=1/eff[2]^2)
    sig_mSPRT_phi25<-av_analysis_no_cov(std_err_tau=std_err,std_res=1,estimate=diff, alpha=alpha*2, phi=1/eff[3]^2)
    sig_mSPRT_phi11<-av_analysis_no_cov(std_err_tau=std_err,std_res=1,estimate=diff, alpha=alpha*2, phi=1/eff[4]^2)
    
    ### Eppo TUNNC
    GAVI750<-GAVI(1:N,alpha=2*alpha,phi=750,sigma2=sigma^2, est=diff)
    GAVI500<-GAVI(1:N,alpha=2*alpha,phi=500,sigma2=sigma^2, est=diff)
    GAVI250<-GAVI(1:N,alpha=2*alpha,phi=250,sigma2=sigma^2, est=diff)
    
    ### Zalando SST
    sig_SST<-SST(N, sqrt(2)*sigma, alpha, est=diff*(1:N))
    sig_pSST7<-(diff*(1:N) > pSST7)
    sig_pSST14<-(diff*(1:N) > pSST14)
    
    ### Gather data
    tib[counter:(counter+(nrows-1)),]<-rbind(
      cbind("SST",
            c(any(sig_SST),any(sig_SST[N*s14]),any(sig_SST[N*s28]),any(sig_SST[N*s42]),any(sig_SST[N*s56])),
            r,
            c("stream", "14","28","42","56"),
            eff[e],
            c(detN(sig_SST, 1:N,N), detN(sig_SST[N*s14], N*s14,N), detN(sig_SST[N*s28], N*s28,N), detN(sig_SST[N*s42], N*s42,N),detN(sig_SST[N*s56], N*s56,N)
            )),
      cbind("pSST7",
            c(any(sig_pSST7),any(sig_pSST7[N*s14]),any(sig_pSST7[N*s28]),any(sig_pSST7[N*s42]),any(sig_pSST7[N*s56])),
            r,
            c("stream", "14","28","42","56"),
            eff[e],
            c(detN(sig_pSST7, 1:N,N), detN(sig_pSST7[N*s14], N*s14,N), detN(sig_pSST7[N*s28], N*s28,N), detN(sig_pSST7[N*s42], N*s42,N),detN(sig_pSST7[N*s56], N*s56,N)
            )),
      cbind("pSST14",
            c(any(sig_pSST14),any(sig_pSST14[N*s14]),any(sig_pSST14[N*s28]),any(sig_pSST14[N*s42]),any(sig_pSST14[N*s56])),
            r,
            c("stream", "14","28","42","56"),
            eff[e],
            c(detN(sig_pSST14, 1:N,N), detN(sig_pSST14[N*s14], N*s14,N), detN(sig_pSST14[N*s28], N*s28,N), detN(sig_pSST14[N*s42], N*s42,N),detN(sig_pSST14[N*s56], N*s56,N)
            )),
      cbind("GAVI250", 
            c(any(GAVI250),any(GAVI250[N*s14]),any(GAVI250[N*s28]),any(GAVI250[N*s42]),any(GAVI250[N*s56])),
            r,
            c("stream", "14","28","42","56"),
            eff[e],
            c(detN(GAVI250, 1:N,N), detN(GAVI250[N*s14], N*s14,N), detN(GAVI250[N*s28], N*s28,N), detN(GAVI250[N*s42], N*s42,N),detN(GAVI250[N*s56], N*s56,N)
            )),
      cbind("GAVI500", 
            c(any(GAVI500),any(GAVI500[N*s14]),any(GAVI500[N*s28]),any(GAVI500[N*s42]),any(GAVI500[N*s56])),
            r,
            c("stream", "14","28","42","56"),
            eff[e],
            c(detN(GAVI500, 1:N,N), detN(GAVI500[N*s14], N*s14,N), detN(GAVI500[N*s28], N*s28,N), detN(GAVI500[N*s42], N*s42,N),detN(GAVI500[N*s56], N*s56,N)
            )),
      cbind("GAVI750", 
            c(any(GAVI750),any(GAVI750[N*s14]),any(GAVI750[N*s28]),any(GAVI750[N*s42]),any(GAVI750[N*s56])),
            r,
            c("stream", "14","28","42","56"),
            eff[e],
            c(detN(GAVI750, 1:N,N), detN(GAVI750[N*s14], N*s14,N), detN(GAVI750[N*s28], N*s28,N), detN(GAVI750[N*s42], N*s42,N),detN(GAVI750[N*s56], N*s56,N)
            )),
      cbind("mSPRTphi100", 
            c(any(sig_mSPRT_phi100),any(sig_mSPRT_phi100[N*s14]),any(sig_mSPRT_phi100[N*s28]),any(sig_mSPRT_phi100[N*s42]),any(sig_mSPRT_phi100[N*s56])),
            r,
            c("stream", "14","28","42","56"),
            eff[e],
            c(detN(sig_mSPRT_phi100, 1:N,N), detN(sig_mSPRT_phi100[N*s14], N*s14,N), detN(sig_mSPRT_phi100[N*s28], N*s28,N), detN(sig_mSPRT_phi100[N*s42], N*s42,N),detN(sig_mSPRT_phi100[N*s56], N*s56,N)
            )),
      cbind("mSPRTphi25", 
            c(any(sig_mSPRT_phi25),any(sig_mSPRT_phi25[N*s14]),any(sig_mSPRT_phi25[N*s28]),any(sig_mSPRT_phi25[N*s42]),any(sig_mSPRT_phi25[N*s56])),
            r,
            c("stream", "14","28","42","56"),
            eff[e],
            c(detN(sig_mSPRT_phi25, 1:N,N), detN(sig_mSPRT_phi25[N*s14], N*s14,N), detN(sig_mSPRT_phi25[N*s28], N*s28,N), detN(sig_mSPRT_phi25[N*s42], N*s42,N),detN(sig_mSPRT_phi25[N*s56], N*s56,N)
            )),
      cbind("mSPRTphi11", 
            c(any(sig_mSPRT_phi11),any(sig_mSPRT_phi11[N*s14]),any(sig_mSPRT_phi11[N*s28]),any(sig_mSPRT_phi11[N*s42]),any(sig_mSPRT_phi11[N*s56])),
            r,
            c("stream", "14","28","42","56"),
            eff[e],
            c(detN(sig_mSPRT_phi11, 1:N,N), detN(sig_mSPRT_phi11[N*s14], N*s14,N), detN(sig_mSPRT_phi11[N*s28], N*s28,N), detN(sig_mSPRT_phi11[N*s42], N*s42,N),detN(sig_mSPRT_phi11[N*s56], N*s56,N)
            )),
      cbind("Bonferroni", 
            c(any(zdiff>Bonferroni_stream),any(zdiff[N*s14]>Bonferroni_14),any(zdiff[N*s28]>Bonferroni_28),any(zdiff[N*s42]>Bonferroni_42),any(zdiff[N*s56]>Bonferroni_56)),
            r,
            c("stream", "14","28","42","56"),
            eff[e],
            c(detN(zdiff>Bonferroni_stream, 1:N,N), detN(zdiff[N*s14]>Bonferroni_14, N*s14,N), detN(zdiff[N*s28]>Bonferroni_28, N*s28,N), detN(zdiff[N*s42]>Bonferroni_42, N*s42,N),detN(zdiff[N*s56]>Bonferroni_56, N*s56,N)
            )),
      cbind("GST", 
            c(zdiff[N]>gst1, any(zdiff[N*s14]>gst14), any(zdiff[N*s28]>gst28), any(zdiff[N*s42]>gst42), any(zdiff[N*s56]>gst56)),
            r,
            c("1", "14", "28", "42", "56"),
            eff[e],
            c(N, detN(zdiff[N*s14]>gst14, N*s14,N), detN(zdiff[N*s28]>gst28, N*s28,N), detN(zdiff[N*s42]>gst42, N*s42,N), detN(zdiff[N*s56]>gst56, N*s56,N)
            )),
      cbind("GSToversampled", 
            c(any(zdiff[N*s14]>gst14o), any(zdiff[N*s28]>gst28o), any(zdiff[N*s42]>gst42o), any(zdiff[N*s56]>gst56o)),
            r,
            c("14", "28", "42", "56"),
            eff[e],
            c(detN(zdiff[N*s14]>gst14o, N*s14,N), detN(zdiff[N*s28]>gst28o, N*s28,N), detN(zdiff[N*s42]>gst42o, N*s42,N), detN(zdiff[N*s56]>gst56o, N*s56,N)
            )),
      cbind("GSTundersampled", 
            c(any(zdiff[N*s14]>gst14u), any(zdiff[N*s28]>gst28u), any(zdiff[N*s42]>gst42u), any(zdiff[N*s56]>gst56u)),
            r,
            c("14", "28", "42", "56"),
            eff[e],
            c(detN(zdiff[N*s14]>gst14u, N*s14,N), detN(zdiff[N*s28]>gst28u, N*s28,N), detN(zdiff[N*s42]>gst42u, N*s42,N), detN(zdiff[N*s56]>gst56u, N*s56,N)
            )),
      cbind("STATSIG", 
            c(any(zdiff[-1]>zcrit_ss[-1]), any(zdiff[N*s14]>zcrit_ss[N*s14]), any(zdiff[N*s28]>zcrit_ss[N*s28]), any(zdiff[N*s42]>zcrit_ss[N*s42]), any(zdiff[N*s56]>zcrit_ss[N*s56])),
            r,
            c("stream", "14", "28", "42", "56"),
            eff[e],
            c(detN(zdiff>zcrit_ss, 1:N, N), detN(zdiff[N*s14]>zcrit_ss[N*s14], N*s14, N), detN(zdiff[N*s28]>zcrit_ss[N*s28], N*s28, N), detN(zdiff[N*s42]>zcrit_ss[N*s42], N*s42, N), detN(zdiff[N*s56]>zcrit_ss[N*s56], N*s56, N)
            )),
      cbind("GSTphi3", 
            c(zdiff[N]>gst1phi3, any(zdiff[N*s14]>gst14phi3), any(zdiff[N*s28]>gst28phi3), any(zdiff[N*s42]>gst42phi3), any(zdiff[N*s56]>gst56phi3)),
            r,
            c("1", "14", "28", "42", "56"),
            eff[e],
            c(N, detN(zdiff[N*s14]>gst14phi3, N*s14,N), detN(zdiff[N*s28]>gst28phi3, N*s28,N), detN(zdiff[N*s42]>gst42phi3, N*s42,N), detN(zdiff[N*s56]>gst56phi3, N*s56,N)
            )),
      cbind("GSToversampledphi3", 
            c(any(zdiff[N*s14]>gst14ophi3), any(zdiff[N*s28]>gst28ophi3), any(zdiff[N*s42]>gst42ophi3), any(zdiff[N*s56]>gst56ophi3)),
            r,
            c("14", "28", "42", "56"),
            eff[e],
            c(detN(zdiff[N*s14]>gst14ophi3, N*s14,N), detN(zdiff[N*s28]>gst28ophi3, N*s28,N), detN(zdiff[N*s42]>gst42ophi3, N*s42,N), detN(zdiff[N*s56]>gst56ophi3, N*s56,N)
            )),
      cbind("GSTundersampledphi3", 
            c(any(zdiff[N*s14]>gst14uphi3), any(zdiff[N*s28]>gst28uphi3), any(zdiff[N*s42]>gst42uphi3), any(zdiff[N*s56]>gst56uphi3)),
            r,
            c("14", "28", "42", "56"),
            eff[e],
            c(detN(zdiff[N*s14]>gst14uphi3, N*s14,N), detN(zdiff[N*s28]>gst28uphi3, N*s28,N), detN(zdiff[N*s42]>gst42uphi3, N*s42,N), detN(zdiff[N*s56]>gst56uphi3, N*s56,N)
            ))
    )
    counter<-counter+nrows
  }
  print(r)
}



#################### RESULTS 

tib1<-as_tibble(tib)
names(tib1)<-c("type","significant","rep","number_checks","true_effect","detectionN")
tib1$significant<-tib1$significant=="TRUE"
tib1$true_effect<-as.numeric(tib1$true_effect)
tib1$detectionN<-as.numeric(tib1$detectionN)

tib2<-tib1 %>% 
  group_by(type, number_checks, true_effect) %>% 
  summarise(power=mean(significant,na.rm=TRUE), saving=1-mean(detectionN,na.rm=TRUE)/N)

#FPR
fpr_tib<-filter(tib2, true_effect==0, number_checks!=1) %>% select(-saving, -true_effect) %>%  spread(number_checks, power)
print(fpr_tib)

#Power
power_tib_1<-filter(tib2, true_effect==0.1, number_checks!=1,!type%in%c("GSToversampledphi3","GSToversampled","STATSIG","mSPRTtau-2","mSPRTtau-4")) %>% select(-saving, -true_effect) %>%  spread(number_checks, power)
print(power_tib_1)

power_tib_2<-filter(tib2, true_effect==0.2, number_checks!=1,!type%in%c("GSToversampledphi3","GSToversampled","STATSIG","mSPRTtau-2","mSPRTtau-4")) %>% select(-saving, -true_effect) %>%  spread(number_checks, power)
print(power_tib_2)

power_tib_3<-filter(tib2, true_effect==0.3, number_checks!=1,!type%in%c("GSToversampledphi3","GSToversampled","STATSIG","mSPRTtau-2","mSPRTtau-4")) %>% select(-saving, -true_effect) %>%  spread(number_checks, power)
print(power_tib_3)

#Average Savings

n_tib_1<-filter(tib2, true_effect==0.1, number_checks!=1,!type%in%c("GSToversampledphi3","GSToversampled","STATSIG","mSPRTtau-2","mSPRTtau-4")) %>% select(-power, -true_effect) %>%  spread(number_checks, saving)
print(n_tib_1)

n_tib_2<-filter(tib2, true_effect==0.2, number_checks!=1,!type%in%c("GSToversampledphi3","GSToversampled","STATSIG","mSPRTtau-2","mSPRTtau-4")) %>% select(-power, -true_effect) %>%  spread(number_checks, saving)
print(n_tib_2)

n_tib_3<-filter(tib2, true_effect==0.3, number_checks!=1,!type%in%c("GSToversampledphi3","GSToversampled","STATSIG","mSPRTtau-2","mSPRTtau-4")) %>% select(-power, -true_effect) %>%  spread(number_checks, saving)
print(n_tib_3)
