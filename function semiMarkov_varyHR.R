##########################################################################################
##========================================================================================
## START OF SEMIMARKOV_VARYHR FUNCTION
##========================================================================================
##########################################################################################

##========================================================================================
## THIS FUNCTION BUILDS A SEMI-MARKOV MULTI-STATE MODEL AND RETURNS THE ASSOCIATED STATE
## OCCUPANCY PROBABILITIES (PROBABILITIES OVER TIME OF BEING IN EACH STATE)

## IT ALLOWS THE TREATMENT EFFECT IN THE EXTRAPOLATION PERIOD TO VARY FROM THAT IN THE 
## OBSERVED PERIOD. IT DOES THIS BY CHANGING THE HAZARD RATIO IN THE EXTRAPOLATION PERIOD.
## THIS AMENDED HR IS APPLIED TO THE ONE TREATMENT GROUP TO OBTAIN THE HAZARD FOR THE OTHER 
## TREATMENT GROUP

## TREATMENT MUST BE GIVEN AS THE FIRST ARGUMENT FOR EACH TRANSITION FOR covs

## MODIFIBLE ARGUMENTS:
## ntrans       NUMBER OF TRANSITIONS IN THE MODELLING
## ncovs        NUMBER OF COVARIATE PARAMETERS IN THE MODEL (EXCLUDING INTERCEPT) FOR EACH
##              TRANSITION I.E. ONE FOR EACH BINARY AND CONTINUOUS VARIABLE AND K-1 FOR  
##              EACH CATEGORICAL VARIABLE, WHERE K= NO OF CATEGORIES
## covs         VARIABLE NAMES FOR THE COVARIATES
## coveval      VALUE AT WHICH TO EVALUATE EACH COVARIATE
## dist         DISTRIBUTION TO USE FOR EACH TRANSITION. OPTIONS ARE:
##              wei FOR WEIBULL, exp FOR EXPONENTIAL, gom FOR GOMPERTZ,
##              logl FOR LOGLOGISTIC,logn FOR LOGNORMAL AND gam FOR GENERALISED GAMMA.
##              IF FITTING THE SAME MODEL OVER THE OBSERVED AND EXTRAPOLATED PERIOD 
##              DISTRIBUTION WILL SPAN THE WHOLE TIME HORIZON. IF FITTING A DIFFERENT 
##              MODEL OVER THE OBSERVED PERIOD THAN FOR THE EXTRAPOLATION, DISTRIBUTION
##              WILL SPAN THE OBSERVED PERIOD ONLY.
## dist2        ONLY APPLICABLE IF FITTING A DIFFERENT MODEL OVER THE EXTRAPOLATION PERIOD
##              THAN THE OBSERVED PERIOD.DISTRIBUTION TO USE FOR EACH TRANSITION OVER THE
##              EXTRAPOLATION PERIOD. OPTIONS ARE:
##              wei FOR WEIBULL, exp FOR EXPONENTIAL, gom FOR GOMPERTZ,
##              logl FOR LOGLOGISTIC,logn FOR LOGNORMAL AND gam FOR GENERALISED GAMMA.
## timeseq      THE TIME POINTS TO USE FOR PREDICTIONS OVER THE OBSERVED PERIOD OF THE STUDY.
##              THE FIRST ARGUMENT OF seq SHOULD BE THE START TIME, THE SECOND ARGUMENT THE 
##              END TIME AND THE THIRD ARGUMENT THE TIME INCREMENT. 
## timeseq_ext  THE TIME POINTS TO USE FOR PREDICTIONS OVER THE EXTRAPOLATION PERIOD.
##              THE FIRST ARGUMENT OF seq SHOULD BE THE START TIME, THE SECOND ARGUMENT THE 
##              END TIME AND THE THIRD ARGUMENT THE TIME INCREMENT. 
## data         DATASET TO USE FOR MODELLING
## seedno       NUMBER TO USE TO SET THE RANDOM NUMBER GENERATOR SO THAT SIMULATIONS CAN
##              BE REPLICATED EXACTLY. IF NOT REQUIRED SET TO NULL
## M            NUMBER OF SIMULATIONS USED TO CALCULATE PROBABILITIES 
## trans        TRANSITION MATRIX
## predinitial  PREDICT FROM INITIAL STATE? EITHER TRUE OR FALSE
## predfrom     IF NOT PREDICTING FROM INITIAL STATE, NUMBER OF STATE IN WHICH TO START
##              PREDICTION 
## varyHR       ALTER THE HAZARD RATIO SO IT IS DIFFERENT IN THE EXTRAPOLATION PERIOD THAN
##              IN THE OBSERVED PERIOD?. EITHER TRUE OR FALSE.
## HRtreat      ONLY APPLICABLE IF varyHR=TRUE. THE HAZARD RATIOS TO USE IN THE
##              EXTRAPOLATION PERIOD FOR EACH TRANSITION 
## adjust       ONLY APPLICABLE IF varyHR=TRUE. NUMBER OF TIME POINTS TO ADJUST FOR EACH
##              TRANSITION TO FIX PROBLEMS WITH THE GAP IN HAZARDS (DIMINISH TREATMENT EFFECT
##              FOR THE NUMBER OF TIME POINTS GIVEN BEFORE STOPPING THE TREATMENT EFFECT
##              COMPLETELY
##========================================================================================
semiMarkov_varyHR<-function(ntrans=3, ncovs=c(1,1,2), 
                     covs=rbind("covariate1", "covariate1",c("covariate1", "covariate2")),
                     coveval=rbind(0,0,c(0,1)),  
                     dist=cbind("wei", "wei", "wei"),
                     dist2=cbind(NA, NA, NA),
                           timeseq=seq(0,4,1/12),
                     timeseq_ext=c(seq(49/12,118/12,1/12), seq(118/12+1/144, 12, 1/144), 
                                   seq(12+1/600,15, 1/600)),
              data=msmcancer, seedno=12345, M=100,
              trans=tmat,predinitial=TRUE, predfrom=2,
              varyHR=FALSE, HRtreat=c(1,1,1),adjust=c(0,0,1)){
 
  #### set up required lists 
  models<-vector("list", ntrans)
  models_ext<-vector("list", ntrans)
  fmla<-vector("list", ntrans)
  covars<-vector("list", ntrans)
  datasub<-vector("list", ntrans)
  basecaseHR<-rep(NA, ntrans)
  inc<-rep(NA, ntrans)
  lp<-vector("list", ntrans)
  lp_ext<-vector("list", ntrans)
  coeffs<-vector("list", ntrans)
  coeffs_ext<-vector("list", ntrans)
  lp2<-vector("list", ntrans)
  lp3<-vector("list", ntrans)
  coeffs2<-vector("list", ntrans)
  coeffs3<-vector("list", ntrans)
  coeffs4<-vector("list", ntrans)
  x<-vector("list", ntrans)
  cumHaz<-vector("list", ntrans) 
  cumHaz_ext<-vector("list", ntrans) 
  kappa<-vector("list", ntrans) 
  gamma<-vector("list", ntrans) 
  mu<-vector("list", ntrans) 
  sigma<-vector("list", ntrans) 
  z<-vector("list", ntrans) 
  u<-vector("list", ntrans) 
  
  #### create the timepoints
  tt<-timeseq
  tt2<-timeseq_ext
  for (i in 1:ntrans) {
  
  #### coefficients from modelling of each transition
 
  covars[[i]]<-covs[i,1:ncovs[i]]
  fmla[[i]]<-as.formula(paste("Surv(time,status)~ ",paste(covars[[i]],
                          collapse= "+"))) 
  datasub[[i]]<-subset(data,trans==i) 
  x[[i]]<-coveval[i,1:ncovs[i]]
    if (dist[i]=="wei") {
      models[[i]]<-phreg(fmla[[i]],dist="weibull", data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz[[i]]<-exp(-exp(models[[i]][ncovs[i]+2])*models[[i]][ncovs[i]+1]+lp[[i]])*
      tt^exp(models[[i]][ncovs[i]+2])    
    }
  
  if ((is.na(dist2[i])==FALSE & dist2[i] =="wei")|(is.na(dist2[i])==TRUE & dist[i] =="wei")) {
    models_ext[[i]]<-phreg(fmla[[i]],dist="weibull", data=datasub[[i]])$coeff
    coeffs_ext[[i]]<-models_ext[[i]][1:ncovs[i]] 
    
    if (varyHR==FALSE){
      lp_ext[[i]]<-sum(coeffs_ext[[i]]* x[[i]] )}
    if (varyHR==TRUE){
      coeffs_ext[[i]][1]<-log(HRtreat[i])
    lp_ext[[i]]<-sum(coeffs_ext[[i]]* x[[i]] )}
    
    #### cumulative hazards for each transition
    cumHaz_ext[[i]]<-exp(-exp(models_ext[[i]][ncovs[i]+2])*models_ext[[i]][ncovs[i]+1]+lp_ext[[i]])*
      tt2^exp(models_ext[[i]][ncovs[i]+2]) 
    
    if (varyHR==TRUE){
      basecaseHR[i]<-exp(models_ext[[i]][1])
      if (adjust[i]!=0) inc[i]<-(basecaseHR[i]-HRtreat[i])/(adjust[i]+1)
      if (adjust[i]!=0)  coeffs3[[i]][1:adjust[i]]<- log(seq(basecaseHR[i]-inc[i],HRtreat[i], -inc[i])) 
      #coeffs4[[i]]<-coeffs2[[i]]
      if (adjust[i]!=0) { for (j in 1:adjust[i]) {
        coeffs4[[i]][[j]]<-coeffs_ext[[i]]
        coeffs4[[i]][[j]][1]<-coeffs3[[i]][j]
        lp3[[i]][j]<-sum(coeffs4[[i]][[j]]* x[[i]])
      }}
      if (adjust[i]==0) lp3[[i]]<-lp_ext[[i]]
      
      if (adjust[i]!=0) cumHaz_ext[[i]][1:adjust[i]]<-exp(-exp(models_ext[[i]][ncovs[i]+2])*models_ext[[i]][ncovs[i]+1]+lp3[[i]][1:adjust[i]])*
          tt2[1:adjust[i]]^exp(models_ext[[i]][ncovs[i]+2]) 
      }
  }
  
    if (dist[i]=="exp") {
      models[[i]]<-phreg(fmla[[i]],dist="weibull", shape=1, data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz[[i]]<-exp(-models[[i]][ncovs[i]+1]+lp[[i]])*tt    
    } 
  
  if ((is.na(dist2[i])==FALSE & dist2[i] =="exp")|(is.na(dist2[i])==TRUE & dist[i] =="exp")) {
    models_ext[[i]]<-phreg(fmla[[i]],dist="weibull", shape=1,  data=datasub[[i]])$coeff
    coeffs_ext[[i]]<-models_ext[[i]][1:ncovs[i]] 
    
    if (varyHR==FALSE){
      lp_ext[[i]]<-sum(coeffs_ext[[i]]* x[[i]] )}
    if (varyHR==TRUE){
      coeffs_ext[[i]][1]<-log(HRtreat[i])
      lp_ext[[i]]<-sum(coeffs_ext[[i]]* x[[i]] )}
    
    #### cumulative hazards for each transition
    cumHaz_ext[[i]]<-exp(-models_ext[[i]][ncovs[i]+1]+lp_ext[[i]])*tt2    
    
    if (varyHR==TRUE){
      basecaseHR[i]<-exp(models_ext[[i]][1])
      if (adjust[i]!=0) inc[i]<-(basecaseHR[i]-HRtreat[i])/(adjust[i]+1)
      if (adjust[i]!=0)  coeffs3[[i]][1:adjust[i]]<- log(seq(basecaseHR[i]-inc[i],HRtreat[i], -inc[i])) 
      #coeffs4[[i]]<-coeffs2[[i]]
      if (adjust[i]!=0) { for (j in 1:adjust[i]) {
        coeffs4[[i]][[j]]<-coeffs_ext[[i]]
        coeffs4[[i]][[j]][1]<-coeffs3[[i]][j]
        lp3[[i]][j]<-sum(coeffs4[[i]][[j]]* x[[i]])
      }}
      if (adjust[i]==0) lp3[[i]]<-lp_ext[[i]]
      
      if (adjust[i]!=0) cumHaz_ext[[i]][1:adjust[i]]<-exp(-models_ext[[i]][ncovs[i]+1]+lp3[[i]][1:adjust[i]])*tt2[1:adjust[i]]  
      
    }
  }
  

  if (dist[i]=="gom") {
    models[[i]]<-phreg(fmla[[i]],dist="gompertz", param="rate", data=datasub[[i]])$coeff
    coeffs[[i]]<-models[[i]][1:ncovs[i]] 
    lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
    #### cumulative hazards for each transition
    cumHaz[[i]]<-exp(models[[i]][ncovs[i]+2]+lp[[i]])*(1/models[[i]][ncovs[i]+1])*
      (exp(models[[i]][ncovs[i]+1]*tt) -1) 
    
  }
  
  if ((is.na(dist2[i])==FALSE & dist2[i] =="gom")|(is.na(dist2[i])==TRUE & dist[i] =="gom")) {
    models_ext[[i]]<-phreg(fmla[[i]],dist="gompertz", param="rate", data=datasub[[i]])$coeff
    coeffs_ext[[i]]<-models_ext[[i]][1:ncovs[i]] 
    
    if (varyHR==FALSE){
      lp_ext[[i]]<-sum(coeffs_ext[[i]]* x[[i]] )}
    if (varyHR==TRUE){
      coeffs_ext[[i]][1]<-log(HRtreat[i])
      lp_ext[[i]]<-sum(coeffs_ext[[i]]* x[[i]] )}
    
    #### cumulative hazards for each transition
    cumHaz_ext[[i]]<-exp(models_ext[[i]][ncovs[i]+2]+lp_ext[[i]])*(1/models_ext[[i]][ncovs[i]+1])*
      (exp(models_ext[[i]][ncovs[i]+1]*tt2) -1) 
    
    if (varyHR==TRUE){
      basecaseHR[i]<-exp(models_ext[[i]][1])
      if (adjust[i]!=0) inc[i]<-(basecaseHR[i]-HRtreat[i])/(adjust[i]+1)
      if (adjust[i]!=0)  coeffs3[[i]][1:adjust[i]]<- log(seq(basecaseHR[i]-inc[i],HRtreat[i], -inc[i])) 
      #coeffs4[[i]]<-coeffs2[[i]]
      if (adjust[i]!=0) { for (j in 1:adjust[i]) {
        coeffs4[[i]][[j]]<-coeffs_ext[[i]]
        coeffs4[[i]][[j]][1]<-coeffs3[[i]][j]
        lp3[[i]][j]<-sum(coeffs4[[i]][[j]]* x[[i]])
      }}
      if (adjust[i]==0) lp3[[i]]<-lp_ext[[i]]
      
      if (adjust[i]!=0) cumHaz_ext[[i]][1:adjust[i]]<-exp(models_ext[[i]][ncovs[i]+2]+lp3[[i]][1:adjust[i]])*(1/models_ext[[i]][ncovs[i]+1])*
          (exp(models_ext[[i]][ncovs[i]+1]*tt2[1:adjust[i]]) -1) 
      
    }
  }

    if (dist[i]=="logl") {
      models[[i]]<-aftreg(fmla[[i]],dist="loglogistic",  data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz[[i]]<- -log(1/(1+(exp(-(models[[i]][ncovs[i]+1]-lp[[i]]))*tt)^
                              (1/(exp(-models[[i]][ncovs[i]+2])))))
    }
  
  if ((is.na(dist2[i])==FALSE & dist2[i] =="logl")|(is.na(dist2[i])==TRUE & dist[i] =="logl")) {
    models_ext[[i]]<-aftreg(fmla[[i]],dist="loglogistic",  data=datasub[[i]])$coeff
    coeffs_ext[[i]]<-models_ext[[i]][1:ncovs[i]] 
    
    if (varyHR==FALSE){
      lp_ext[[i]]<-sum(coeffs_ext[[i]]* x[[i]] )}
    if (varyHR==TRUE){
      coeffs_ext[[i]][1]<-log(HRtreat[i])
      lp_ext[[i]]<-sum(coeffs_ext[[i]]* x[[i]] )}
    
    #### cumulative hazards for each transition
    cumHaz_ext[[i]]<--log(1/(1+(exp(-(models_ext[[i]][ncovs[i]+1]-lp_ext[[i]]))*tt2)^
                               (1/(exp(-models_ext[[i]][ncovs[i]+2])))))
    
    if (varyHR==TRUE){
      basecaseHR[i]<-exp(models_ext[[i]][1])
      if (adjust[i]!=0) inc[i]<-(basecaseHR[i]-HRtreat[i])/(adjust[i]+1)
      if (adjust[i]!=0)  coeffs3[[i]][1:adjust[i]]<- log(seq(basecaseHR[i]-inc[i],HRtreat[i], -inc[i])) 
      #coeffs4[[i]]<-coeffs2[[i]]
      if (adjust[i]!=0) { for (j in 1:adjust[i]) {
        coeffs4[[i]][[j]]<-coeffs_ext[[i]]
        coeffs4[[i]][[j]][1]<-coeffs3[[i]][j]
        lp3[[i]][j]<-sum(coeffs4[[i]][[j]]* x[[i]])
      }}
      if (adjust[i]==0) lp3[[i]]<-lp_ext[[i]]
      
      if (adjust[i]!=0) cumHaz_ext[[i]][1:adjust[i]]<--log(1/(1+(exp(-(models_ext[[i]][ncovs[i]+1]-lp3[[i]][1:adjust[i]]))*tt2[1:adjust[i]])^
                                                                (1/(exp(-models_ext[[i]][ncovs[i]+2])))))
      
    }
  }
  
    if (dist[i]=="logn") {
      models[[i]]<-aftreg(fmla[[i]],dist="lognormal",  data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz[[i]]<- -log(1-pnorm((log(tt)-(models[[i]][ncovs[i]+1]-lp[[i]]))/
                                   (exp(-models[[i]][ncovs[i]+2]))))  
    }
  
  
  if ((is.na(dist2[i])==FALSE & dist2[i] =="logn")|(is.na(dist2[i])==TRUE & dist[i] =="logn")) {
    models_ext[[i]]<-aftreg(fmla[[i]],dist="lognormal",  data=datasub[[i]])$coeff
    coeffs_ext[[i]]<-models_ext[[i]][1:ncovs[i]] 
    
    if (varyHR==FALSE){
      lp_ext[[i]]<-sum(coeffs_ext[[i]]* x[[i]] )}
    if (varyHR==TRUE){
      coeffs_ext[[i]][1]<-log(HRtreat[i])
      lp_ext[[i]]<-sum(coeffs_ext[[i]]* x[[i]] )}
    
    #### cumulative hazards for each transition
    cumHaz_ext[[i]]<--log(1-pnorm((log(tt2)-(models_ext[[i]][ncovs[i]+1]-lp_ext[[i]]))/
                                    (exp(-models_ext[[i]][ncovs[i]+2]))))  
    
    if (varyHR==TRUE){
      basecaseHR[i]<-exp(models_ext[[i]][1])
      if (adjust[i]!=0) inc[i]<-(basecaseHR[i]-HRtreat[i])/(adjust[i]+1)
      if (adjust[i]!=0)  coeffs3[[i]][1:adjust[i]]<- log(seq(basecaseHR[i]-inc[i],HRtreat[i], -inc[i])) 
      #coeffs4[[i]]<-coeffs2[[i]]
      if (adjust[i]!=0) { for (j in 1:adjust[i]) {
        coeffs4[[i]][[j]]<-coeffs_ext[[i]]
        coeffs4[[i]][[j]][1]<-coeffs3[[i]][j]
        lp3[[i]][j]<-sum(coeffs4[[i]][[j]]* x[[i]])
      }}
      if (adjust[i]==0) lp3[[i]]<-lp_ext[[i]]
      
      if (adjust[i]!=0) cumHaz_ext[[i]][1:adjust[i]]<--log(1-pnorm((log(tt2[1:adjust[i]])-(models_ext[[i]][ncovs[i]+1]-lp3[[i]][1:adjust[i]]))/
                                                                     (exp(-models_ext[[i]][ncovs[i]+2]))))  
      
    }
  }

    if (dist[i]=="gam") {
      models[[i]]<-flexsurvreg(fmla[[i]],dist="gengamma",data=datasub[[i]])$res
      kappa[[i]]<- models[[i]][3,1]
      gamma[[i]]<-(abs(kappa[[i]]))^(-2)
      coeffs[[i]]<-models[[i]][4:(ncovs[i]+3),1]
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      mu[[i]]<- models[[i]][1,1]  +lp[[i]] 
      sigma[[i]]<-  models[[i]][2,1]
      z[[i]]<-rep(0,length(tt))
      z[[i]]<- sign(kappa[[i]])*((log(tt)-mu[[i]])/sigma[[i]])
      u[[i]]<-gamma[[i]]*exp((abs(kappa[[i]]))*z[[i]])
      if(kappa[[i]]>0){
        cumHaz[[i]]<--log(1-pgamma(u[[i]],gamma[[i]]))
      }
      if(kappa[[i]]==0){
        cumHaz[[i]]<--log(1-pnorm(z[[i]]))
      }
      if(kappa[[i]]<0){
        cumHaz[[i]]<--log(pgamma(u[[i]],gamma[[i]]))
      }
    }
  
  
  if ((is.na(dist2[i])==FALSE & dist2[i] =="gam")|(is.na(dist2[i])==TRUE & dist[i] =="gam")) {
    models_ext[[i]]<-flexsurvreg(fmla[[i]],dist="gengamma",data=datasub[[i]])$res
    kappa[[i]]<- models_ext[[i]][3,1]
    gamma[[i]]<-(abs(kappa[[i]]))^(-2)
    coeffs_ext[[i]]<-models_ext[[i]][4:(ncovs[i]+3),1]

    if (varyHR==FALSE){
      lp_ext[[i]]<-sum(coeffs_ext[[i]]* x[[i]] )}
    if (varyHR==TRUE){
      coeffs_ext[[i]][1]<-log(HRtreat[i])
      lp_ext[[i]]<-sum(coeffs_ext[[i]]* x[[i]] )}
    mu[[i]]<- models_ext[[i]][1,1]  +lp_ext[[i]] 
    sigma[[i]]<-  models_ext[[i]][2,1]
    z[[i]]<-rep(0,length(tt2))
    z[[i]]<- sign(kappa[[i]])*((log(tt2)-mu[[i]])/sigma[[i]])
    u[[i]]<-gamma[[i]]*exp((abs(kappa[[i]]))*z[[i]])
    if(kappa[[i]]>0){
      cumHaz_ext[[i]]<--log(1-pgamma(u[[i]],gamma[[i]]))
    }
    if(kappa[[i]]==0){
      cumHaz_ext[[i]]<--log(1-pnorm(z[[i]]))
    }
    if(kappa[[i]]<0){
      cumHaz_ext[[i]]<--log(pgamma(u[[i]],gamma[[i]]))
    }
    
    if (varyHR==TRUE){
      basecaseHR[i]<-exp(models_ext[[i]][4,1])
      if (adjust[i]!=0) inc[i]<-(basecaseHR[i]-HRtreat[i])/(adjust[i]+1)
      if (adjust[i]!=0)  coeffs3[[i]][1:adjust[i]]<- log(seq(basecaseHR[i]-inc[i],HRtreat[i], -inc[i])) 
      #coeffs4[[i]]<-coeffs2[[i]]
      if (adjust[i]!=0) { for (j in 1:adjust[i]) {
        coeffs4[[i]][[j]]<-coeffs_ext[[i]]
        coeffs4[[i]][[j]][1]<-coeffs3[[i]][j]
        lp3[[i]][j]<-sum(coeffs4[[i]][[j]]* x[[i]])
      
      }}
      if (adjust[i]==0) lp3[[i]]<-lp_ext[[i]]
      if (adjust[i]!=0){
      mu[[i]][1:adjust[i]]<- models_ext[[i]][1,1]  +lp3[[i]] [1:adjust[i]]
      sigma[[i]]<-  models_ext[[i]][2,1]
      z[[i]][1:adjust[i]]<-rep(0,length(tt2[1:adjust[i]]))
      z[[i]][1:adjust[i]]<- sign(kappa[[i]])*((log(tt2[1:adjust[i]])-mu[[i]][1:adjust[i]])/sigma[[i]])
      u[[i]][1:adjust[i]]<-gamma[[i]]*exp((abs(kappa[[i]]))*z[[i]])
      }
      if(kappa[[i]]>0 & adjust[i]!=0){
        cumHaz_ext[[i]][1:adjust[i]]<--log(1-pgamma(u[[i]][1:adjust[i]],gamma[[i]]))
      }
      if(kappa[[i]]==0 & adjust[i]!=0){
        cumHaz_ext[[i]][1:adjust[i]]<--log(1-pnorm(z[[i]][1:adjust[i]]))
      }
      if(kappa[[i]]<0 & adjust[i]!=0 ){
        cumHaz_ext[[i]][1:adjust[i]]<--log(pgamma(u[[i]][1:adjust[i]],gamma[[i]]))
      }
      
    }
  }
      cumHaz[[i]]<-c(cumHaz[[i]], cumHaz_ext[[i]])
  } 
    Haz<-unlist(cumHaz)
    if (varyHR==TRUE) tt<- c(timeseq,timeseq_ext)
    tt<-c(timeseq,timeseq_ext)
    newtrans<-rep(1:ntrans,each=length(tt))
    time<-rep(tt,ntrans)
   Haz<-cbind(time=as.vector(time),Haz=as.vector(Haz),trans=as.vector(newtrans))
    Haz<-as.data.frame(Haz)
  ## state occupancy probabilities
    set.seed(seedno)
    if (predinitial==TRUE) {
    stateprobs <- mssample(Haz=Haz,trans=trans,tvec=tt,clock="reset", M=M)
    }
    if (predinitial==FALSE) {
      stateprobs<-mssample(Haz=Haz,trans=trans, tvec=tt,clock="reset", M=M,
                       history=list(state=predfrom,time=0,tstate=NULL))  
    }  
  return(stateprobs) 
    }  
##########################################################################################
##========================================================================================
## END OF SEMIMARKOV_VARYHR FUNCTION
##========================================================================================
##########################################################################################



















