##########################################################################################
##========================================================================================
## START OF MEANLY FUNCTION
##========================================================================================
##########################################################################################

##========================================================================================
## THIS FUNCTION RETURNS THE MEAN LIFE YEARS IN A GIVEN STATE FOR A GIVEN TREATMENT ARM
## BY CALCULATING THE AREA UNDER THE APPROPRIATE CURVE USING THE TRAPEZOIDAL RULE 

## MODIFIBLE ARGUMENTS:
## Markov         EITHER TRUE OR FALSE, FOR BEING BASED ON A MARKOV OR SEMI-MARKOV MODEL
##                RESPECTIVELY
## object         OBJECT RETURNED FROM THE Markov OR semiMarkov FUNCTION AS APPROPRIATE
## state          STATE OF INTEREST
## instate        EITHER TRUE OR FALSE. TRUE GIVES THE MEAN LIFE YEARS SPENT IN THE STATE. 
##                FALSE GIVES THE MEAN LIFE YEARS OUTWITH THE STATE
## discounted     DISCOUNTED? EITHER TRUE OR FALSE.
## dis1yronwards  START DISCOUNTING FROM YEAR 1 INSTEAD OF TIME 0?  EITHER TRUE OR FALSE.
## rate           IF DISCOUNTED, RATE EXPRESSED AS A PROPORTION I.E. BETWEEN 0 AND 1

##========================================================================================
meanLY<-function(Markov=FALSE,object, state=1, instate=TRUE,
                       discounted=TRUE, dis1yronwards=TRUE, rate=0.035){
    if (Markov==FALSE) {
  
      if (instate==TRUE & discounted==TRUE) {
        object[,(ncol(object)+1)]<-object[,state+1]/(1+rate)^object$time
        if (dis1yronwards==TRUE){
        for (i in 1: length(object$time)){
          if(object$time[i]<1) object[i,(ncol(object))]<-object[i,state+1]}}  
        meanLY<- trapz(object$time, object[,(ncol(object))])
        }
      
      if (instate==FALSE & discounted==TRUE) {
        object[,ncol(object)+1]<-(1-object[,state+1])/(1+rate)^object$time
        if (dis1yronwards==TRUE){
        for (i in 1: length(object$time)){
          if(object$time[i]<1) object[i,ncol(object)]<-1- object[i,state+1] }} 
          meanLY<- trapz(object$time, object[,ncol(object)])
        }
        
      if (instate==TRUE & discounted==FALSE) {
         meanLY<-trapz(object$time,object[,state+1])
      }
      if (instate==FALSE & discounted==FALSE) {
         meanLY<-trapz(object$time,1-object[,state+1])
      }
    }
  
    if (Markov==TRUE) {
      
      if (instate==TRUE & discounted==TRUE) {
        object[[1]][,ncol(object[[1]])+1]<-object[[1]][,state+1]/(1+rate)^object[[1]]$time
        if (dis1yronwards==TRUE){
        for (i in 1: length(object[[1]]$time)){
          if(object[[1]]$time[i]<1) object[[1]][i,ncol(object[[1]])]<-object[[1]][i,state+1] }}  
          meanLY<- trapz(object[[1]]$time, object[[1]][,ncol(object[[1]])])
        } 
        
        
      if (instate==FALSE & discounted==TRUE) {
        object[[1]][,ncol(object[[1]])+1]<-(1-object[[1]][,state+1])/(1+rate)^object[[1]]$time
        if (dis1yronwards==TRUE){
        for (i in 1: length(object[[1]]$time)){
          if(object[[1]]$time[i]<1) object[[1]][i,ncol(object[[1]])]<-1-object[[1]][i,state+1] }}  
          meanLY<- trapz(object[[1]]$time, object[[1]][,ncol(object[[1]])])
        }
        
      if (instate==TRUE & discounted==FALSE) {
        meanLY<-trapz(object[[1]]$time,object[[1]][,state+1])
      }
      if (instate==FALSE & discounted==FALSE) {
        meanLY<-trapz(object[[1]]$time,1-object[[1]][,state+1])
      }
    }
  
  return(meanLY) 
}

##########################################################################################
##========================================================================================
## END OF MEANLY FUNCTION
##========================================================================================
##########################################################################################


