##########################################################################################
##========================================================================================
## START OF MODELPARAM FUNCTION
##========================================================================================
##########################################################################################

##========================================================================================
## THIS FUNCTION DISPLAYS THE MODEL PARAMETERS FOR A TRANSITION

## MODIFIBLE ARGUMENTS:
## Markov       EITHER TRUE OR FALSE, FOR A MARKOV OR SEMI-MARKOV MODEL RESPECTIVELY
## covs         VARIABLE NAMES FOR THE COVARIATES
## transnum     TRANSITION NUMBER OF INTEREST
## dist         DISTRIBUTION TO USE. OPTIONS ARE:
##              wei FOR WEIBULL, exp FOR EXPONENTIAL, gom FOR GOMPERTZ,
##              logl FOR LOGLOGISTIC,logn FOR LOGNORMAL AND gam FOR GENERALISED GAMMA           
## data         DATASET TO USE FOR MODELLING
##========================================================================================
modelparam<-function(Markov=FALSE,covs="treat", transnum=3, dist="wei", data=msmcancer) {
  ##### set up models
    if (Markov==TRUE){
      fmla<-as.formula(paste("Surv(Tstart,Tstop,status)~ ",paste(covs, collapse= "+"))) 
    }
    if (Markov==FALSE){
      fmla<-as.formula(paste("Surv(time,status)~ ",paste(covs, collapse= "+"))) 
    }
    datasub<-subset(data,trans==transnum) 
    
    if (dist=="wei") model<-phreg(fmla,dist="weibull", data=datasub)
    if (dist=="exp") model<-phreg(fmla,dist="weibull", shape=1, data=datasub)
    if (dist=="gom") model<-phreg(fmla,dist="gompertz", param="rate", data=datasub)
    if (dist=="logl") model<-aftreg(fmla,dist="loglogistic", data=datasub)
    if (dist=="logn") model<-aftreg(fmla,dist="lognormal", data=datasub)
    if (dist=="gam") model<-flexsurvreg(fmla,dist="gengamma", data=datasub)
    
    return(model)
}
    
##########################################################################################
##========================================================================================
## END OF MODELPARAM FUNCTION
##========================================================================================
##########################################################################################    
    

