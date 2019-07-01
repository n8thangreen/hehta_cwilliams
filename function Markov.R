##########################################################################################
##========================================================================================
## START OF MARKOV FUNCTION
##========================================================================================
##########################################################################################

##========================================================================================
## THIS FUNCTION BUILDS A MARKOV MULTI-STATE MODEL AND RETURNS THE ASSOCIATED STATE
## OCCUPANCY PROBABILITIES (PROBABILITIES OVER TIME OF BEING IN EACH STATE)

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
## trans        TRANSITION MATRIX
##========================================================================================

Markov<-function(ntrans=3, ncovs=c(1,1,2), 
                      covs=rbind("covariate1", "covariate1",c("covariate1", "covariate2")) ,
                     coveval=rbind(0,0,c(0,1)),  
                 dist=cbind("wei", "wei", "wei"),
                 dist2=cbind(NA, NA, NA),
                 timeseq=seq(0,4,1/12),
                 timeseq_ext=seq(49/12,15,1/12),
                      data=msmcancer, trans=tmat){
  #### set up required lists 
  models<-vector("list", ntrans)
  fmla<-vector("list", ntrans)
  covars<-vector("list", ntrans)
  datasub<-vector("list", ntrans)
  lp<-vector("list", ntrans)
  coeffs<-vector("list", ntrans)
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
  tt2<-timeseq_ext
  for (i in 1:ntrans) {
    
    if (is.na(dist2[i])==TRUE) tt<-c(timeseq,timeseq_ext)
    if (is.na(dist2[i])==FALSE) tt<-timeseq
  
  #### coefficients from modelling of each transition
 
    covars[[i]]<-covs[i,1:ncovs[i]]
    fmla[[i]]<-as.formula(paste("Surv(Tstart,Tstop,status)~ ",paste(covars[[i]], collapse= "+"))) 
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
    if (is.na(dist2[i])==FALSE & dist2[i] =="wei") {
      models[[i]]<-phreg(fmla[[i]],dist="weibull", data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz_ext[[i]]<-exp(-exp(models[[i]][ncovs[i]+2])*models[[i]][ncovs[i]+1]+lp[[i]])*
        tt2^exp(models[[i]][ncovs[i]+2])  
    }
    if (dist[i]=="exp") {
      models[[i]]<-phreg(fmla[[i]],dist="weibull", shape=1, data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz[[i]]<-exp(-models[[i]][ncovs[i]+1]+lp[[i]])*tt    
    } 
    
    if (is.na(dist2[i])==FALSE & dist2[i] =="exp") {
      models[[i]]<-phreg(fmla[[i]],dist="weibull",shape=1, data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz_ext[[i]]<-exp(-models[[i]][ncovs[i]+1]+lp[[i]])*tt2    
    }
    
    if (dist[i]=="gom") {
      models[[i]]<-phreg(fmla[[i]],dist="gompertz", param="rate", data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz[[i]]<-exp(models[[i]][ncovs[i]+2]+lp[[i]])*(1/models[[i]][ncovs[i]+1])*
        (exp(models[[i]][ncovs[i]+1]*tt) -1)   
    }
    
    if (is.na(dist2[i])==FALSE & dist2[i] =="gom") {
      models[[i]]<-phreg(fmla[[i]],dist="gompertz", param="rate", data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz_ext[[i]]<-exp(models[[i]][ncovs[i]+2]+lp[[i]])*(1/models[[i]][ncovs[i]+1])*
        (exp(models[[i]][ncovs[i]+1]*tt2) -1)  
    }
    
    if (dist[i]=="logl") {
      models[[i]]<-aftreg(fmla[[i]],dist="loglogistic",  data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz[[i]]<- -log(1/(1+(exp(-(models[[i]][ncovs[i]+1]-lp[[i]]))*tt)^
                              (1/(exp(-models[[i]][ncovs[i]+2])))))
    }
    
    if (is.na(dist2[i])==FALSE & dist2[i] =="logl") {
      models[[i]]<-aftreg(fmla[[i]],dist="loglogistic",  data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz_ext[[i]]<--log(1/(1+(exp(-(models[[i]][ncovs[i]+1]-lp[[i]]))*tt2)^
                                 (1/(exp(-models[[i]][ncovs[i]+2])))))
    }
    
    
    if (dist[i]=="logn") {
      models[[i]]<-aftreg(fmla[[i]],dist="lognormal",  data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz[[i]]<- -log(1-pnorm((log(tt)-(models[[i]][ncovs[i]+1]-lp[[i]]))/
                                   (exp(-models[[i]][ncovs[i]+2]))))  
    }
    
    if (is.na(dist2[i])==FALSE & dist2[i] =="logn") {
      models[[i]]<-aftreg(fmla[[i]],dist="lognormal",  data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz_ext[[i]]<--log(1-pnorm((log(tt2)-(models[[i]][ncovs[i]+1]-lp[[i]]))/
                                      (exp(-models[[i]][ncovs[i]+2])))) 
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
    
    if (is.na(dist2[i])==FALSE & dist2[i] =="gam") {
      models[[i]]<-flexsurvreg(fmla[[i]],dist="gengamma",data=datasub[[i]])$res
      kappa[[i]]<- models[[i]][3,1]
      gamma[[i]]<-(abs(kappa[[i]]))^(-2)
      coeffs[[i]]<-models[[i]][4:(ncovs[i]+3),1]
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      mu[[i]]<- models[[i]][1,1]  +lp[[i]] 
      sigma[[i]]<-  models[[i]][2,1]
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
    }
    if (is.na(dist2[i])==FALSE) cumHaz[[i]]<-c(cumHaz[[i]], cumHaz_ext[[i]]) 
  } 
    Haz<-unlist(cumHaz)
    tt<-c(timeseq,timeseq_ext)
    newtrans<-rep(1:ntrans,each=length(tt))
    time<-rep(tt,ntrans)
    Haz<-cbind(time=as.vector(time),Haz=as.vector(Haz),trans=as.vector(newtrans))
    Haz<-as.data.frame(Haz)
  ##### state occupancy probabilities
  msf<-list(Haz,trans)
  msf$trans<-trans
  msf$Haz<-Haz
  class(msf)<-"msfit"
  stateprobs <- probtrans(msf, predt = 0,variance=FALSE)
  return(stateprobs)
}
##########################################################################################
##========================================================================================
## END OF MARKOV FUNCTION
##========================================================================================
##########################################################################################




















