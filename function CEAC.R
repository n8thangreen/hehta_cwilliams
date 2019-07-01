##########################################################################################
##========================================================================================
## START OF CEAC FUNCTION
##========================================================================================
##########################################################################################

##========================================================================================
## THIS FUNCTION PLOTS A COST-EFFECTIVENESS ACCEPTABILITY CURVE TO HELP VISUALISE THE 
## RESULTS OF THE PROBABILISTIC SENSITIVITY ANALYSIS

## MODIFIBLE ARGUMENTS:
## cRatiosim   RANGE OF WILLINGNESS TO PAY CEILING RATIOS TO PLOT ALONG THE HORIZONTAL 
##             AXIS 
## nruns       NUMBER OF DRAWS USED FOR THE PROBABILISTIC SENSITIVITY ANALYSIS (SHOULD BE
##             THE SAME AS USED IN THE PSAprob FUNCTION)
## xlower      LOWER LIMIT FOR THE HORIZONTAL AXIS 
## xupper      UPPER LIMIT FOR THE HORIZONTAL AXIS 
## ylower      LOWER LIMIT FOR THE VERTICAL AXIS
## yupper      UPPER LIMIT FOR THE VERTICAL AXIS
## QALY1       NAME OF VARIABLE CONTAINING THE INCREMENTAL QALY FOR TREATMENT ARM 1
##             (DERIVED USING THE PSAQALY FUNCTION)
## QALY2       NAME OF VARIABLE CONTAINING THE INCREMENTAL QALY FOR TREATMENT ARM 2
##             (DERIVED USING THE PSAQALY FUNCTION)
## cost1       NAME OF VARIABLE CONTAINING THE TOTAL EXPECTED COSTS FOR TREATMENT ARM 1
## cost2       NAME OF VARIABLE CONTAINING THE TOTAL EXPECTED COSTS FOR TREATMENT ARM 2
##========================================================================================
CEAC<-function(cRatiosim=c(seq(0,1000,100), seq(1500,5000,500),
                            seq(6000,30000,1000), seq(35000,100000,5000)),
               nruns=1000,nruns2=1000, xlower=0, xupper=100000, ylower=0, yupper=1,
               QALY1=QALY_RFC_dis, QALY2=QALY_FC_dis,
               cost1=total_cost_RFC_dis, cost2=total_cost_FC_dis, 
               secondcurve=TRUE, QALY1_2=QALY_RFC_dis, QALY2_2=QALY_FC_dis,
               cost1_2=total_cost_RFC_dis, cost2_2=total_cost_FC_dis){
  ### store incremental net monetary benefit for each draw
  incNMBsim<-matrix(nrow=nruns, ncol=length(cRatiosim))
  ### store cost-effective indicator for each draw
  is.cost.effsim<-matrix(rep(0, nruns*length(cRatiosim)),
                         nrow=nruns, ncol=length(cRatiosim))
  ### store mean of indicators
  ### i.e. the probability treatment is cost effective
  meanicesim<-rep(NA,length(cRatiosim))
  
  
  for (i in 1:length(cRatiosim)) {
    incNMBsim[,i]<-(QALY1*cRatiosim[i] - cost1)- 
      (QALY2*cRatiosim[i] - cost2)
  }
  
  for (j in 1:length(cRatiosim)) {
    for (i in 1:nruns) {
      is.cost.effsim[i,j][incNMBsim[i,j]>0] <-1
    }}
  
  for (i in 1:length(cRatiosim)) {
    meanicesim[i]<-mean(is.cost.effsim[,i])  
  }
  
  options(scipen = 3)
  plot(cRatiosim,meanicesim, type="l", xlim=c(xlower,xupper), ylim=c(ylower,yupper),
       xlab="Value of ceiling ratio (£)", 
       ylab="Probability of being cost-effective ")  
 
  if (secondcurve==TRUE){
    
    ### store incremental net monetary benefit for each draw
    incNMBsim2<-matrix(nrow=nruns2, ncol=length(cRatiosim))
    ### store cost-effective indicator for each draw
    is.cost.effsim2<-matrix(rep(0, nruns2*length(cRatiosim)),
                           nrow=nruns2, ncol=length(cRatiosim))
    ### store mean of indicators
    ### i.e. the probability treatment is cost effective
    meanicesim2<-rep(NA,length(cRatiosim))
    
    
    for (i in 1:length(cRatiosim)) {
      incNMBsim2[,i]<-(QALY1_2*cRatiosim[i] - cost1_2)- 
        (QALY2_2*cRatiosim[i] - cost2_2)
    }
    
    for (j in 1:length(cRatiosim)) {
      for (i in 1:nruns2) {
        is.cost.effsim2[i,j][incNMBsim2[i,j]>0] <-1
      }}
    
    for (i in 1:length(cRatiosim)) {
      meanicesim2[i]<-mean(is.cost.effsim2[,i])  
    }
    
    lines(cRatiosim,meanicesim2, type="l",lty=2) 
    
  }
  
}
##########################################################################################
##========================================================================================
## END OF CEAC FUNCTION
##========================================================================================
##########################################################################################