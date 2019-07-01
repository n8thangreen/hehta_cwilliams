##########################################################################################
##========================================================================================
## START OF SEMIMARKOV_NOIPD FUNCTION
##========================================================================================
##########################################################################################

##========================================================================================
## THIS FUNCTION RETURNS THE STATE OCCUPANCY PROBABILITIES (PROBABILITIES OVER TIME OF 
## BEING IN EACH STATE) BASED ON THE MODEL COEFFICIENTS OR MORTALITY RATES GIVEN AS INPUTS  

## MODIFIBLE ARGUMENTS:
## ntrans       NUMBER OF TRANSITIONS IN THE MODELLING
## type         SOURCE OF INPUT. EITHER "parametric model" OR "mortality hazard".
##              REQUIRED FOR EACH TRANSITION
## morhaz       APPLICABLE IF TYPE="mortality hazard". 
## ncovs        NUMBER OF COVARIATE PARAMETERS IN THE MODEL (EXCLUDING INTERCEPT) FOR EACH
##              TRANSITION I.E. ONE FOR EACH BINARY AND CONTINUOUS VARIABLE AND K-1 FOR  
##              EACH CATEGORICAL VARIABLE, WHERE K= NO OF CATEGORIES
## covs         VARIABLE NAMES FOR THE COVARIATES
## covcoeff     COEFFICIENTS RELATING TO EACH OF THE COVARIATES FOR EACH OF THE TRANSITIONS
## coveval      VALUE AT WHICH TO EVALUATE EACH COVARIATE
## dist         DISTRIBUTION TO USE FOR EACH TRANSITION. OPTIONS ARE:
##              wei FOR WEIBULL, exp FOR EXPONENTIAL, gom FOR GOMPERTZ,
##              logl FOR LOGLOGISTIC,logn FOR LOGNORMAL AND gam FOR GENERALISED GAMMA 
## logscale     COEFFICIENT FOR LOG(SCALE) FOR EACH OF THE TRANSITIONS - APPLICABLE TO
##              WEIBULL, EXPONENTIAL, LOGLOGISTIC, LOGNORMAL
## logshape     COEFFICIENT FOR LOG(SHAPE) FOR EACH OF THE TRANSITIONS - APPLICABLE TO 
##              WEIBULL, LOGLOGISTIC,LOGNORMAL
## rate         COEFFICIENT FOR RATE - APPLICABLE TO GOMPERTZ
## loglevel     COEFFICIENT FOR LOG(LEVEL) - APPLICABLE TO GOMPERTZ
## mu           COEFFICIENT FOR MU - APPLICABLE TO GENERALISED GAMMA
## sigma        COEFFICIENT FOR SIGMA - APPLICABLE TO GENERALISED GAMMA
## Q            COEFFICIENT FOR Q - APPLICABLE TO GENERALISED GAMMA
## alt1end      THIS SHOULD BE SET TO TRUE IF AN ADJUSTMENT IS REQUIRED TO TIMEPOINTS AT
##              THE END OF SEQUENCE (SEE README FILE). OTHERWISE IT SHOULD BE FALSE
## alt2beg      THIS SHOULD BE SET TO TRUE IF AN ADJUSTMENT IS REQUIRED TO TIMEPOINTS AT
##              THE BEGINNING OF THE SEQUENCE (SEE README FILE).OTHERWISE IT SHOULD BE
##              FALSE
## ttend        APPLICABLE IF alt1end and alt2beg are FALSE. TIMEPOINT AT WHICH TO END 
##              EXTRAPOLATION
## ttinc        APPLICABLE IF alt1end and alt2beg are FALSE. GAP BETWEEN TIMEPOINTS
## tt1end       APPLICABLE IF alt1end = TRUE. THE ENDPOINT OF TIME BEFORE ADJUSTMENT 
##              TO THE LAST PART OF THE TIME SEQUENCE IS REQUIRED
## tt1inc       APPLICABLE IF alt1end = TRUE. THE GAP BETWEEN TIMEPOINTS TO USE
##              BEFORE ADJUSTMENT TO THE LAST PART OF THE TIME SEQUENCE
## tt1end2      APPLICABLE IF alt1end = TRUE. THE TIME HORIZON I.E. THE ENDPOINT OF 
##              THE ADJUSTMENT TO THE LAST PART OF THE TIME SEQUENCE
## tt1inc2      APPLICABLE IF alt1end = TRUE. THE GAP BETWEEN TIMEPOINTS TO USE FOR
##              THE ADJUSTMENT TO THE LAST PART OF THE TIME SEQUENCE (WHICH SHOULD BE 
##              SMALLER THAN tt1inc)
## tt2beg       APPLICABLE IF alt2beg = TRUE. THE TIMEPOINTS TO BE ADDED TO THE
##              BEGINNING OF THE TIME SEQUENCE
## tt2start     APPLICABLE IF alt2beg = TRUE. POINT AT WHICH TO START THE NEXT SECTION OF
##              THE TIME SEQUENCE (AFTER THE ADJUSTMENT OF ADDITIONAL TIMEPOINTS TO THE
##              BEGINNING OF THE SEQUENCE)
## tt2inc       APPLICABLE IF alt2beg = TRUE and alt1end = FALSE. THE GAP BETWEEN 
##              TIMEPOINTS TO USE FOR THE LAST PART OF THE TIME SEQUENCE
## tt2end       APPLICABLE IF alt2beg = TRUE and alt1end = FALSE. THE TIME HORIZON 
##              I.E. THE  FINAL POINT OF THE TIME SEQUENCE
## data         DATASET TO USE FOR MODELLING
## seedno       NUMBER TO USE TO SET THE RANDOM NUMBER GENERATOR SO THAT SIMULATIONS CAN
##              BE REPLICATED EXACTLY. IF NOT REQUIRED SET TO NULL
## M            NUMBER OF SIMULATIONS USED TO CALCULATE PROBABILITIES 
## trans        TRANSITION MATRIX
## predinitial  PREDICT FROM INITIAL STATE? EITHER TRUE OR FALSE
## predfrom     IF NOT PREDICTING FROM INITIAL STATE, NUMBER OF STATE IN WHICH TO START
##              PREDICTION 
##========================================================================================

## NB. THIS FUNCTION USES THE SAME PARAMETERISATIONS AS phreg AND aftreg IN THE eha PACKAGE
## AND flexsurvreg IN THE flexsurv PACKAGE


semiMarkov_noipd<-function(ntrans=3, 
                     type=c("parametric model","parametric model", "mortality hazard"),
                     morhaz=c(0,0, 0.05),
                     ncovs=c(1,1,2), 
                     covs=rbind("covariate1", "covariate1",c("covariate1", "covariate2")),
                     covcoeff=rbind(-0.537, -0.402, c(0.487, -0.687)),     
                     coveval=rbind(0,0,c(0,1)),  
                     dist=cbind("gom", "gom", "gom"),
                     logscale=c(1.289, 3.997, 0.174 ),
                     logshape=c(0.441, -0.236, -0.524 ),
                     rate=c(0.474, 0, -0.487),
                     loglevel=c( -2.187,0,  -2.825),
                     mu=c(1.300, 3.160, 1.370),
                     sigma=c(0.623, 5.780, 0.806),
                     Q=c(1.050, -3.280, 1.160 ),              
                     alt1end=FALSE, alt2beg=FALSE,
                     ttend=15, ttinc=1/12,
                     tt1end=118/12, tt1inc=1/12,tt1end2=15, tt1inc2=0.007,
                     tt2beg=c(0,0.00001,0.005), tt2start=1/12, tt2inc=1/12,tt2end=15,
                     data=msmcancer, seedno=12345, M=100,
                     trans=tmat,predinitial=TRUE, predfrom=2){
  #### create the timepoints
  if(alt1end==FALSE & alt2beg==FALSE){
    tt<-seq(0,ttend,ttinc)
  }
  if(alt1end==TRUE & alt2beg==FALSE){
    tt<-c(seq(0,tt1end,tt1inc), seq(tt1end+tt1inc2, tt1end2,tt1inc2))
  }
  if(alt2beg==TRUE & alt1end==FALSE){
    tt<-c(tt2beg, seq(tt2start, tt2end,tt2inc))
  }
  if(alt2beg==TRUE & alt1end==TRUE){
    tt<-c(tt2beg, seq(tt2start, tt1end,tt1inc), seq(tt1end+tt1inc2, tt1end2,tt1inc2))
  }
  #### set up required lists 
  covcoeffs<-vector("list", ntrans)
  lp<-vector("list", ntrans)
  x<-vector("list", ntrans)
  cumHaz<-vector("list", ntrans) 
  gamma<-vector("list", ntrans) 
  mu2<-vector("list", ntrans) 
  z<-vector("list", ntrans) 
  u<-vector("list", ntrans) 

   for (i in 1:ntrans) {
    #### calculations based on mortality hazards
    
    if (type[i]=="mortality hazard"){
      cumHaz[[i]] <-morhaz[i]*tt
    }
    #### calculations based on coefficients from each transition 
  x[[i]]<-coveval[i,1:ncovs[i]]
  covcoeffs[[i]]<-covcoeff[i,1:ncovs[i]] 
  lp[[i]]<-sum(covcoeffs[[i]]* x[[i]] ) 
  #### cumulative hazards for each transition - distribution specific
  
  if (type[i]== "parametric model" &  dist[i]=="wei") {
    cumHaz[[i]]<-exp(-exp(logshape[i])*logscale[i]+lp[[i]])*tt^exp(logshape[i])
  }
  if (type[i]== "parametric model" & dist[i]=="exp") {
    cumHaz[[i]]<-exp(-logscale[i]+lp[[i]])*tt
  }
  if (type[i]== "parametric model" & dist[i]=="gom") {
    cumHaz[[i]]<-exp(loglevel[i]+lp[[i]])*(1/rate[i])*(exp(rate[i]*tt) -1)   
  }  
  if (type[i]== "parametric model" & dist[i]=="logl") {
    cumHaz[[i]]<- -log(1/(1+(exp(lp[[i]]-logscale[i])*tt)^
                            (1/exp(-logshape[i]))))
  }
  if (type[i]== "parametric model" & dist[i]=="logn") {
    cumHaz[[i]]<- -log(1-pnorm((log(tt)-(logscale[i]-lp[[i]]))/
                                 (exp(-logshape[i]))))  
  }
  if (type[i]== "parametric model" & dist[i]=="gam") {
    gamma[[i]]<-(abs(Q[i]))^(-2)
    mu2[[i]]<- mu[i] +lp[[i]] 
    z[[i]]<-rep(0,length(tt))
    z[[i]]<- sign(Q[i])*((log(tt)-mu2[[i]])/sigma[i])
    u[[i]]<-gamma[[i]]*exp((abs(Q[i]))*z[[i]])
    if(Q[i]>0){
      cumHaz[[i]]<--log(1-pgamma(u[[i]],gamma[[i]]))
    }
    if(Q[i]==0){
      cumHaz[[i]]<--log(1-pnorm(z[[i]]))
    }
    if(Q[i]<0){
      cumHaz[[i]]<--log(pgamma(u[[i]],gamma[[i]]))
    }
  }
  
  }
    Haz<-unlist(cumHaz)
    newtrans<-rep(1:ntrans,each=length(tt))
    time<-rep(tt,ntrans)
    Haz<-cbind(time=as.vector(time),Haz=as.vector(Haz),trans=as.vector(newtrans))
    Haz<-as.data.frame(Haz)
  ##### state occupancy probabilities
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
## END OF SEMIMARKOV_NOIPD FUNCTION
##========================================================================================
##########################################################################################



















