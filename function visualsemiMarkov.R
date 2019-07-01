##########################################################################################
##========================================================================================
## START OF VISUALSEMIMARKOV FUNCTION
##========================================================================================
##########################################################################################

##========================================================================================
## THIS FUNCTION ASSESSES THE FIT FROM SEMI-MARKOV MULTI-STATE MODELS BY COMPARING THEIR 
## PREDICTIONS OF BEING IN A GIVEN STATE OVER TIME WITH AN APPROPRIATE OBSERVED PROPORTION
## OF BEING IN THAT STATE

## MODIFIBLE ARGUMENTS:
## nobjects     NUMBER OF DIFFERENT PREDICTIONS TO PLOT
## objects      NAME(S) GIVEN TO OBJECT(S) RESULTING FROM RUNNING THE semiMarkov FUNCTION
## tteach       NUMBER OF TIMEPOINTS USED IN THE MODELLING FROM EACH OF THE objects 
## state        STATE TO USE FOR PREDICTION
## instate      EITHER TRUE OR FALSE. TRUE GIVES THE PREDICTION OF BEING IN THE STATE 
##              OVER TIME. FALSE GIVES THE PREDICTION OF NOT BEING IN THE STATE (USED
##              FOR ABSORBING STATES)
## absorb       NUMBER OF STATES THE ABSORBING STATE IS SPLIT INTO. CAN BE 1 OR 2
##              ASSUMED TO BE THE LAST OR TWO LAST STATES 
## ylab         VERTICAL AXIS LABEL
## xlab         HORIZONTAL AXIS LABEL
## ylim         RANGE OF VALUES FOR VERTICAL AXIS
## xlim         RANGE OF VALUES FOR HORIZONTAL AXIS
## lwd          WIDTH OF LINES IN PLOT
## col          COLOUR OF LINES IN PLOT
## lty          TYPE OF LINES IN PLOT
## obsdata      DATASET TO USE FOR OBSERVED PROPORTION IN STATE I.E. DATASET TO CREATE KM
##              CURVE, COMPETING RISKS CUMULATIVE INCIDENCE CURVE OR OBJECT RETURNED FROM
##              proportionfunc FUNCTION
## observed     EITHER "KM", "1-KM", "cuminc" OR "proportion" FOR KAPLAN-MEIER CURVE,
##              1- KAPLAN-MEIER CURVE, COMPETING RISKS CUMULATIVE INCIDENCE CURVE OR 
##              OBSERVED PROPORTION RESPECTIVELY. COULD ALSO BE FALSE IF OBSERVED CURVE
##              UNAVAILABLE/UNNECESSARY."KM" IS APPROPRIATE FOR AN INITIAL OR ABSORBING
##              STATE AND PLOTS THE PROBABILITY OF STAYING IN AN INITIAL STATE OR OF NOT
##              BEING IN AN ABSORBING STATE. "1-KM" PLOTS THE PROBABILITY OF NOT BEING IN
##              AN INITIAL STATE OR OF BEING IN AN ABSORBING STATE."cuminc" IS APPROPRIATE
##              FOR COMPETING RISKS I.E. WHEN MORE THAN ONE STATE CAN BE ENTERED FROM A        
##              (USUALLY THE INITIAL) STATE. "proportion" IS APPROPRIATE FOR AN INTERMEDIATE
##              STATE E.G. ONE WITH FLOW IN AND OUT OF IT: SEE proportionfunc FUNCTION
## CI           WHETHER TO INCLUDE A 95% CONFIDENCE INTERVAL FOR THE OBSERVED CURVE.
##              EITHER TRUE OR FALSE.       
## KMtime       ONLY APPLICABLE IF observed="KM" OR "1-KM". TIME VARIABLE TO CREATE KM CURVE
## KMstatus     ONLY APPLICABLE IF observed="KM" OR "1-KM". STATUS VARIABLE TO CREATE KM CURVE
## CRtime       ONLY APPLICABLE IF observed="cuminc". TIME VARIABLE TO CREATE COMPETING RISKS 
##              CUMULATIVE INCIDENCE CURVE.
## CRstatus     ONLY APPLICABLE IF observed="cuminc". STATUS VARIABLE TO CREATE COMPETING RISKS
##              CUMULATIVE INCIDENCE CURVE.
## ncr          ONLY APPLICABLE IF observed="cuminc". ORDINAL NUMBER OF COMPETING 
##              RISK OF INTEREST
## legendpos    POSITION ARGUMENT FOR LEGEND
## legendncol   NUMBER OF COLUMNS ARGUMENT FOR LEGEND
## legendcex    RELATIVE SIZE ARGUMENT FOR LEGEND
## legendcurves NAMES OF CURVES FOR LEGEND
## legendlty    LINE TYPE ARGUMENT FOR LEGEND
## legendbty    BOX TYPE ARGUMENT FOR LEGEND
## legendlwd    LINE WIDTH ARGUMENT FOR LEGEND
## legendcol    COLOUR ARGUMENT FOR LEGEND
## main         TITLE OF PLOT. IF NOT REQUIRED SET TO NULL
## cex.main     RELATIVE SIZE OF TITLE OF PLOT
##========================================================================================
visualsemiMarkov<-function(nobjects=6, 
  objects=rbind(weisimRFC, expsimRFC, gomsimRFC, loglsimRFC, lognsimRFC, gamsimRFC),
  tteach=c(181,181, 857, 181, 181, 181),
  state=1, instate=TRUE, absorb=1,
  ylab = "Probability of being progression-free", 
  xlab="Years since start of study", ylim=c(0,1),xlim=c(0,15),lwd=2,
  col=rep("black",6),
  lty=c("51","dotted","21", "longdash","twodash", "16"),
  obsdata=RFCdata,                         
  observed="KM", CI=TRUE,
  KMtime=RFCdata$progdeath_ty, KMstatus=RFCdata$progdeath,  
  CRtime=RFCdata$progdeath_ty, CRstatus=RFCdata$crstatus, ncr=2,
  legendpos="topright",legendncol=1,legendcex=1,
  legendcurves=c("Kaplan-Meier", "95% CI Kaplan-Meier","Exponential", "Lognormal",
               "Loglogistic", "Weibull", "Generalised gamma", "Gompertz"), 
  legendlty=c("solid","solid","dotted","twodash","longdash","51","16", "21" ),
  legendbty=c("n","o","n","n","n","n","n","n") ,legendlwd=c(2,15,2,2,2,2,2,2),
  legendcol=c("black","lightgrey", "black", "black", "black", "black", "black","black" ),
  main="Probability of being progression-free", cex.main=1.25)
  {
  if (instate==TRUE) {
  #### plot the curve corresponding to the first semiMarkov object
  plot(objects$time[(1:tteach[1])],objects[(1:tteach[1]),state+1],type='l',
       lwd=lwd, lty=lty[1], col=col[1], xlab="", ylab="", ylim=ylim, xlim=xlim,
       font=2, font.axis=2, font.lab=2, cex.axis=1.5, cex.lab=1.5, las=1)

  #### add the observed proportion (if applicable)
  if (is.na(observed)==0 & observed=="proportion") {
    if (CI==TRUE){
    polygon(c(obsdata[,1], rev(obsdata[,1])), c(obsdata[,4],
             rev(obsdata[,3])), col = "lightgrey", border = NA)
    lines(objects$time[(1:tteach[1])],objects[(1:tteach[1]),state+1],type='l',
          lwd=lwd, lty=lty[1])}
    lines(obsdata[,1], obsdata[,2], type="l", lwd=3,lty=1)  
  }
  
  #### add the observed curves (if applicable)
  if (is.na(observed)==0 & observed=="KM" & CI==FALSE) {
    lines(survfit(Surv(KMtime,KMstatus)~ 1,data=obsdata),
          conf.int=FALSE,mark.time=F,lwd=lwd)
  }
  if (is.na(observed)==0 & observed=="KM" & CI==TRUE) {
    temp<-survfit(Surv(KMtime, KMstatus)~1,data=obsdata )
    polygon(c(temp$time, rev(temp$time)), c(temp$upper,
        rev(temp$lower)), col = "lightgrey", border = NA)
    lines(survfit(Surv(KMtime,KMstatus)~ 1,data=obsdata),
          conf.int=FALSE,mark.time=F,lwd=lwd)
  }
  
  if (is.na(observed)==0 & observed=="1-KM" & CI==FALSE) {
    lines(survfit(Surv(KMtime,KMstatus)~ 1,fun="event",data=obsdata),
          conf.int=FALSE,mark.time=F,lwd=lwd)
  }
  if (is.na(observed)==0 & observed=="1-KM" & CI==TRUE) {
    temp<-survfit(Surv(KMtime, KMstatus)~1, data=obsdata )
    polygon(c(temp$time, rev(temp$time)), c(1-temp$upper,
      rev(1-temp$lower)), col = "lightgrey", border = NA)
    lines(survfit(Surv(KMtime,KMstatus)~ 1,data=obsdata),
    fun="event", conf.int=FALSE,mark.time=F,lwd=lwd)
  }
  
  if (is.na(observed)==0 & observed=="cuminc" & CI==FALSE) {
    CRcuminc<-cuminc(CRtime,CRstatus) 
    lines(CRcuminc[[ncr]]$time,CRcuminc[[ncr]]$est, type="s", lwd=lwd, lty=1)
  }
  if (is.na(observed)==0 & observed=="cuminc" & CI==TRUE) {
    CRcuminc<-cuminc(CRtime,CRstatus)
    upper<-CRcuminc[[ncr]]$est+1.96*sqrt(CRcuminc[[ncr]]$var)
    lower<-CRcuminc[[ncr]]$est-1.96*sqrt(CRcuminc[[ncr]]$var)
    polygon(c(CRcuminc[[ncr]]$time, rev(CRcuminc[[ncr]]$time)), c(upper,
      rev(lower)), col = "lightgrey", border = NA)
    lines(CRcuminc[[ncr]]$time,CRcuminc[[ncr]]$est, type="s", lwd=lwd, lty=1)
  }
  
  #### add the curve corresponding to the first semiMarkov object again
  lines(objects$time[(1:tteach[1])],objects[(1:tteach[1]),state+1],type='l',
       lwd=lwd, lty=lty[1],col=col[1],xlab="", ylab="", ylim=ylim, xlim=xlim)
  
  #### add the curves corresponding to the other semiMarkov objects
  for (i in 2:nobjects) {
    lines(objects$time[(sum(tteach[1:i-1])+1):sum(tteach[1:i])],
       objects[(sum(tteach[1:i-1])+1):sum(tteach[1:i]),state+1], lty=lty[i], col=col[i],lwd=lwd)
  }
  
  }
  if (instate==FALSE & absorb==1) {
    #### plot the curve corresponding to the first semiMarkov object
    plot(objects$time[(1:tteach[1])],1-(objects[(1:tteach[1]),state+1]),type='l',
          lwd=lwd, lty=lty[1],col=col[1],xlab="", ylab="", ylim=ylim, xlim=xlim,
         font=2, font.axis=2, font.lab=2, cex.axis=1.5, cex.lab=1.5, las=1)
    #### add the observed curve
    if (is.na(observed)==0 & observed=="KM" & CI==FALSE) {
      lines(survfit(Surv(KMtime,KMstatus)~ 1,data=obsdata),
            conf.int=FALSE,mark.time=F,lwd=lwd)
    }
    if (is.na(observed)==0 & observed=="KM" & CI==TRUE) {
      temp<-survfit(Surv(KMtime,KMstatus)~1)
      polygon(c(temp$time, rev(temp$time)), c(temp$upper,
    rev(temp$lower)), col = "lightgrey", border = NA)
      lines(survfit(Surv(KMtime,KMstatus)~ 1),
            conf.int=FALSE,mark.time=F,lwd=lwd)
    }
    if (is.na(observed)==0 & observed=="1-KM" & CI==FALSE) {
      lines(survfit(Surv(KMtime,KMstatus)~ 1,fun="event",data=obsdata),
            conf.int=FALSE,mark.time=F,lwd=lwd)
    }
    if (is.na(observed)==0 & observed=="1-KM" & CI==TRUE) {
      temp<-survfit(Surv(KMtime, KMstatus)~1, data=obsdata )
      polygon(c(temp$time, rev(temp$time)), c(1-temp$upper,
        rev(1-temp$lower)), col = "lightgrey", border = NA)
      lines(survfit(Surv(KMtime,KMstatus)~ 1,data=obsdata),
            fun="event", conf.int=FALSE,mark.time=F,lwd=lwd)
    }
    #### add the curve corresponding to the first semiMarkov object again
    lines(objects$time[(1:tteach[1])],1-(objects[(1:tteach[1]),state+1]),type='l',
         lwd=lwd, lty=lty[1],col=col[1],xlab="", ylab="", ylim=ylim, xlim=xlim)
    
    #### add the curves corresponding to the other semiMarkov objects
    for (i in 2:nobjects) {
      lines(objects$time[(sum(tteach[1:i-1])+1):sum(tteach[1:i])],
        1-(objects[(sum(tteach[1:i-1])+1):sum(tteach[1:i]),state+1]), lty=lty[i], col=col[i],lwd=lwd)
    } 
  }
  if (instate==FALSE & absorb==2) {
    #### plot the curve corresponding to the first semiMarkov object
    plot(objects$time[(1:tteach[1])],1-(objects[(1:tteach[1]),state+1] +
                        objects[(1:tteach[1]),state]),type='l',
         lwd=lwd, lty=lty[1],col=col[1],xlab="", ylab="", ylim=ylim, xlim=xlim,
         font=2, font.axis=2, font.lab=2, cex.axis=1.5, cex.lab=1.5, las=1)
    #### add the curves corresponding to the other semiMarkov objects
    for (i in 2:nobjects) {
      lines(objects$time[(sum(tteach[1:i-1])+1):sum(tteach[1:i])],
        1-(objects[(sum(tteach[1:i-1])+1):sum(tteach[1:i]),state+1] + 
        objects[(sum(tteach[1:i-1])+1):sum(tteach[1:i]),state]  ), lty=lty[i], col=col[i],lwd=lwd)
    }
    #### add the observed curve
    if (is.na(observed)==0 & observed=="KM" & CI==FALSE) {
      lines(survfit(Surv(KMtime,KMstatus)~ 1,data=obsdata),
            conf.int=FALSE,mark.time=F,lwd=lwd)
    }
    if (is.na(observed)==0 & observed=="KM" & CI==TRUE) {
      temp<-survfit(Surv(KMtime, KMstatus)~1,data=obsdata )
      polygon(c(temp$time, rev(temp$time)), c(temp$upper,
        rev(temp$lower)), col = "lightgrey", border = NA)
      lines(survfit(Surv(KMtime, KMstatus)~ 1,data=obsdata),
            conf.int=FALSE,mark.time=F,lwd=lwd)
    }
    if (is.na(observed)==0 & observed=="1-KM" & CI==FALSE) {
      lines(survfit(Surv(KMtime,KMstatus)~ 1,fun="event",data=obsdata),
            conf.int=FALSE,mark.time=F,lwd=lwd)
    }
    if (is.na(observed)==0 & observed=="1-KM" & CI==TRUE) {
      temp<-survfit(Surv(KMtime, KMstatus)~1, data=obsdata )
      polygon(c(temp$time, rev(temp$time)), c(1-temp$upper,
      rev(1-temp$lower)), col = "lightgrey", border = NA)
      lines(survfit(Surv(KMtime,KMstatus)~ 1,data=obsdata),
            fun="event", conf.int=FALSE,mark.time=F,lwd=lwd)
    }
   
  } 
  #### add the legend
  legend(legendpos, legendcurves, lty=legendlty,bty=legendbty,lwd=legendlwd,
         col=legendcol,ncol=legendncol,cex=legendcex)    
  ### add a title
  title(main=main, cex.main=cex.main, xlab=xlab, ylab=ylab, cex.lab=1.5)
}
##########################################################################################
##========================================================================================
## END OF VISUALSEMIMARKOV FUNCTION
##========================================================================================
##########################################################################################

