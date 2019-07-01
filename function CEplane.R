##########################################################################################
##========================================================================================
## START OF CEplane FUNCTION
##========================================================================================
##########################################################################################

##========================================================================================
## THIS FUNCTION PLOTS A COST-EFFECTIVENESS PLANE TO HELP VISUALISE THE RESULTS OF THE 
## PROBABILISTIC SENSITIVITY ANALYSIS

## MODIFIBLE ARGUMENTS:
## x           NAME OF VARIABLE CONTAINING INCREMENTAL QALY TO PLOT ALONG THE HORIZONTAL 
##             AXIS (DERIVED USING PSAQALY FUNCTION)
## y           NAME OF VARIABLE CONTAINING INCREMENTAL COST TO PLOT ALONG THE VERTICAL 
##             AXIS
## xlower      LOWER LIMIT FOR THE HORIZONTAL AXIS 
## xupper      UPPER LIMIT FOR THE HORIZONTAL AXIS 
## ylower      LOWER LIMIT FOR THE VERTICAL AXIS
## yupper      UPPER LIMIT FOR THE VERTICAL AXIS
## ICER        WILLINGNESS TO PAY THRESHOLD
## text        TEXT TO ACCOMPANY WILLINGNESS TO PAY THRESHOLD  
##========================================================================================

CEplane<-function(x=incQALY, y=incCost, xlower=-2, xupper=2,
                  ylower=0, yupper=20000, ICER=30000, text="ICER= £30,000"){
  plot(x,y,xlim=c(xlower,xupper), ylim=c(ylower,yupper),
    xlab="Incremental QALY", ylab="Incremental Cost (£)")
  abline(v=0,lty=3)
  abline(h=0,lty=3)
  if(is.na(ICER)==0) abline(0,ICER)
  if (is.na(ICER)==0)text(xupper-0.75, yupper, text)
}
##########################################################################################
##========================================================================================
## END OF CEplane FUNCTION
##========================================================================================
##########################################################################################

