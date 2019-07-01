
###================================================================
### FIT PARAMETRIC SEMI-MARKOV STATE ARRIVAL EXTENDED MODELS
###================================================================

weimodel_pps <-phreg (Surv(time,status)~ treat
                     ,dist="weibull",data=msmcancer3) ; weimodel_pps
aic<- -2*weimodel_pps$loglik[2]+2*nrow(weimodel_pps$var);aic
expmodel_pps <-phreg (Surv(time,status)~ treat
                     ,dist="weibull",shape=1,data=msmcancer3) ; expmodel_pps
aic<- -2*expmodel_pps$loglik[2]+2*nrow(expmodel_pps$var);aic
gommodel_pps <-phreg (Surv(time,status)~ treat
                     ,dist="gompertz",param="rate",data=msmcancer3) ; gommodel_pps
aic<- -2*gommodel_pps$loglik[2]+2*nrow(gommodel_pps$var);aic
loglmodel_pps <-aftreg (Surv(time,status)~ treat
                       ,dist="loglogistic",data=msmcancer3) ; loglmodel_pps
aic<- -2*loglmodel_pps$loglik[2]+2*nrow(loglmodel_pps$var);aic
lognmodel_pps <-aftreg (Surv(time,status)~ treat
                       ,dist="lognormal",data=msmcancer3) ; lognmodel_pps
aic<- -2*lognmodel_pps$loglik[2]+2*nrow(lognmodel_pps$var);aic
gengammodel_pps <-flexsurvreg (Surv(time,status)~ treat
                              ,dist="gengamma",data=msmcancer3) ; gengammodel_pps

gengammodel_pps$AIC


###================================================================
### CALCULATE PREDICTIONS
###================================================================

tt15<-seq(0,15,1/12)

#### WEIBULL
gamma_pps <- exp(weimodel_pps$coef[3])
lambda_FC_pps<- exp(-exp(weimodel_pps$coef[3])*weimodel_pps$coef[2]); lambda_FC_pps
lambda_RFC_pps<- exp(-exp(weimodel_pps$coef[3])*weimodel_pps$coef[2] + weimodel_pps$coef[1]); lambda_RFC_pps

wei_FC_pps15 <- exp(-lambda_FC_pps*tt15^gamma_pps);
wei_RFC_pps15 <- exp(-lambda_RFC_pps*tt15^gamma_pps);

#### EXPONENTIAL

lambda_FC_pps<- exp(-expmodel_pps$coef[2]); lambda_FC_pps
lambda_RFC_pps<- exp(-expmodel_pps$coef[2] + expmodel_pps$coef[1]); lambda_RFC_pps

exp_FC_pps15 <- exp(-lambda_FC_pps*tt15);
exp_RFC_pps15 <- exp(-lambda_RFC_pps*tt15);

#### GOMPERTZ
gamma_pps<- gommodel_pps$coef[2]
lambdaFC_pps<-exp(gommodel_pps$coef[3])
lambdaRFC_pps<-exp(gommodel_pps$coef[3]+gommodel_pps$coef[1])

gom_FC_pps15 <- exp((-lambdaFC_pps/gamma_pps)*(exp(gamma_pps*tt15)-1))
gom_RFC_pps15 <- exp((-lambdaRFC_pps/gamma_pps)*(exp(gamma_pps*tt15)-1))

### log-logistic 
gamma_pps<-exp(-loglmodel_pps$coeff[3])
lambdaFC_pps<-exp(-loglmodel_pps$coeff[2])
lambdaRFC_pps<-exp(-loglmodel_pps$coeff[2]+ loglmodel_pps$coeff[1])

logl_FC_pps15<- 1/(1+(lambdaFC_pps*tt15)^(1/gamma_pps))
logl_RFC_pps15<- 1/(1+(lambdaRFC_pps*tt15)^(1/gamma_pps))

### log normal
muFC<-lognmodel_pps$coeff[2]
muRFC<-lognmodel_pps$coeff[2]-lognmodel_pps$coeff[1]
sigma<-exp(-lognmodel_pps$coeff[3])

logn_FC_pps15<-1-pnorm((log(tt15)-muFC)/sigma)
logn_RFC_pps15<-1-pnorm((log(tt15)-muRFC)/sigma)

#### generalised gamma 
(kappa<- gengammodel_pps$res[3,1]) # check that kappa is positive
gamma_pps<- kappa^(-2)
muFC<-  gengammodel_pps$res[1,1] 
muRFC<-  gengammodel_pps$res[1,1] + gengammodel_pps$res[4,1]
sigma<- gengammodel_pps$res[2,1]
# zFC<-rep(0,181)
# zRFC<-rep(0,181)
zFC<-(log(tt15)-muFC)/sigma
zRFC<-(log(tt15)-muRFC)/sigma
uFC<-gamma_pps*exp(kappa*zFC)
uRFC<-gamma_pps*exp(kappa*zRFC)
gengam_FC_pps15<- 1-pgamma(uFC,gamma_pps)
gengam_RFC_pps15<- 1-pgamma(uRFC,gamma_pps)

########################################################################
### OBSERVED VS PREDICTED PLOTS - OVER TRIAL OBSERVATION PERIOD
#######################################################################

## =====================================================================
### GRAYSCALE FOR PUBLICATION
## =====================================================================

pdf('Figure5_pps_obs.pdf', width=8, height=8)

par(mar= c(5, 7, 4, 2)+0.1)

plot(survfit(Surv(time,status)~ 1, data=msmcancer3RFC), mark.time=F, ylim=c(0,1), xlim=c(0,2.5),conf.int=FALSE,  las=1,
     fun="event", xlab="Years since progression",
     ylab = "Probability of death after progression",
     font=2, font.axis=2, font.lab=2, bty="n", cex.axis=1.5, cex.lab=1.5, axes=FALSE, lwd=3, xaxs="i",yaxs="i") 
axis(1, at = 0:4, labels = c(0:4),lwd=2, lwd.ticks = 2, font=2, cex.axis=1.5)
axis(2, at = c(1, 0.8,0.6, 0.4, 0.2,0), labels = c(1,0.8,0.6, 0.4, 0.2,0),lwd=2, 
     lwd.ticks = 2, font=2, las=1, cex.axis=1.5)

polygon(c(survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$time, rev(survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$time)),
        c(1-survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$lower, rev(1-survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$upper)), 
        col=gray(0.9), border = NA)
lines(tt15,1-exp_RFC_pps15,type='l',col="black", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,4.2),lwd=3,lty="dotted")
lines(tt15,1-wei_RFC_pps15,type='l',col=gray(0.2), xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="51")
lines(tt15,1-gom_RFC_pps15,type='l',col=gray(0.3), xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="21")
lines(tt15,1-logn_RFC_pps15,type='l',col=gray(0.4), xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="twodash")
lines(tt15,1-logl_RFC_pps15,type='l',col=gray(0.5), xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="longdash")
lines(tt15,1-gengam_RFC_pps15,type='l',col="black", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="16")
lines(survfit(Surv(time,status)~ 1, data=msmcancer3RFC), mark.time=F, ylim=c(0,1), xlim=c(0,4),conf.int=FALSE,  las=1,
      fun="event", xlab="Years since progression",
      ylab = "Probability of death after progression",
      font=2, font.axis=2, font.lab=2, bty="n", cex.axis=1.5, cex.lab=1.5, axes=FALSE, lwd=3, xaxs="i",yaxs="i") 
box(lwd=2)
legend(0, 1, c("Kaplan-Meier", "95% CI for KM",  "Weibull", "log-logistic", "generalised gamma","Gompertz", "exponential","log normal"), 
       col=c("black", gray(0.9), gray(0.2), gray(0.5),"black",gray(0.3),"black",gray(0.4)), bty="n",
       lwd=c(2,10,2,2,2,2,2,2),cex=1.5,
       lty=c("solid","solid", "51", "longdash","16","21","dotted","twodash"))
dev.off()


########################################################################
### OBSERVED VS PREDICTED - EXTRAPOLATION TO 15 YEARS
########################################################################

## =====================================================================
### GRAYSCALE FOR PUBLICATION
## =====================================================================

par(mar= c(5, 4, 4, 2)+0.1)

pdf('Figure6_pps_ext.pdf', width=8, height=8)
par(mar= c(5, 7, 4, 2)+0.1)

plot(survfit(Surv(time,status)~ 1, data=msmcancer3RFC), mark.time=F, ylim=c(0,1), xlim=c(0,15),conf.int=FALSE,  las=1,
     fun="event", xlab="Years since progression",
     ylab = "Probability of death after progression",
     font=2, font.axis=2, font.lab=2, bty="n", cex.axis=1.5, cex.lab=1.5, axes=FALSE, lwd=3, xaxs="i",yaxs="i") 
axis(1, at = c(0, 2,4, 6, 8,10,12,14), lwd=2, lwd.ticks = 2, font=2, cex.axis=1.5)
axis(2, at = c(1, 0.8,0.6, 0.4, 0.2,0), labels = c(1,0.8,0.6, 0.4, 0.2,0),lwd=2, 
     lwd.ticks = 2, font=2, las=1, cex.axis=1.5)

polygon(c(survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$time, rev(survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$time)),
        c(1-survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$lower, rev(1-survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$upper)), 
        col=gray(0.9), border = NA)
lines(tt15,1-exp_RFC_pps15,type='l',col="black", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,4.2),lwd=3,lty="dotted")
lines(tt15,1-wei_RFC_pps15,type='l',col=gray(0.2), xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="51")
lines(tt15,1-gom_RFC_pps15,type='l',col=gray(0.3), xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="21")
lines(tt15,1-logn_RFC_pps15,type='l',col=gray(0.4), xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="twodash")
lines(tt15,1-logl_RFC_pps15,type='l',col=gray(0.5), xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="longdash")
lines(tt15,1-gengam_RFC_pps15,type='l',col="black", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="16")
lines(survfit(Surv(time,status)~ 1, data=msmcancer3RFC), mark.time=F, ylim=c(0,1), xlim=c(0,4),conf.int=FALSE,  las=1,
      fun="event", xlab="Years since progression",
      ylab = "Probability of death after progression",
      font=2, font.axis=2, font.lab=2, bty="n", cex.axis=1.5, cex.lab=1.5, axes=FALSE, lwd=3, xaxs="i",yaxs="i") 
box(lwd=2)
legend(7, 0.6, c("Kaplan-Meier", "95% CI for KM",  "Weibull", "Gompertz",  "exponential","generalised gamma","log-logistic", "log normal"  ), 
       col=c("black", gray(0.9), gray(0.2),gray(0.3),"black","black",gray(0.5), gray(0.4)), bty="n",
       lwd=c(3,10,2,2,2,2,2,2),cex=1.5,
       lty=c("solid","solid","51" ,"21","dotted", "16","longdash","twodash"))
dev.off()




########################################################################
### OBSERVED VS PREDICTED - OVER TRIAL OBSERVATION PERIOD
#######################################################################

## =====================================================================
### COLOUR PLOTS
## =====================================================================

par(mar= c(5, 7, 4, 2)+0.1)

# remove the '#' from the line below to show all four plots together

#par(mfrow=c(2,2))

plot(survfit(Surv(time,status)~ 1, data=msmcancer3RFC), mark.time=F, ylim=c(0,1), xlim=c(0,2.5),conf.int=FALSE,  las=1,
     fun="event", xlab="Years since progression",
     ylab = "Probability of death after progression",
     font=2, font.axis=2, font.lab=2, bty="n", cex.axis=1.5, cex.lab=1.5, axes=FALSE, lwd=3, xaxs="i",yaxs="i") 
axis(1, at = 0:4, labels = c(0:4),lwd=2, lwd.ticks = 2, font=2, cex.axis=1.5)
axis(2, at = c(1, 0.8,0.6, 0.4, 0.2,0), labels = c(1,0.8,0.6, 0.4, 0.2,0),lwd=2, 
     lwd.ticks = 2, font=2, las=1, cex.axis=1.5)

polygon(c(survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$time, rev(survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$time)),
        c(1-survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$lower, rev(1-survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$upper)), 
        col = "lightgrey", border = NA)
lines(tt15,1-exp_RFC_pps15,type='l',col="red", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,4.2),lwd=3,lty="dotted")
lines(tt15,1-wei_RFC_pps15,type='l',col="blue", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="51")
lines(tt15,1-gom_RFC_pps15,type='l',col="green", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="21")
lines(tt15,1-logn_RFC_pps15,type='l',col="purple", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="twodash")
lines(tt15,1-logl_RFC_pps15,type='l',col="orange", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="longdash")
lines(tt15,1-gengam_RFC_pps15,type='l',col="cyan", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="16")
lines(survfit(Surv(time,status)~ 1, data=msmcancer3RFC), mark.time=F, ylim=c(0,1), xlim=c(0,4),conf.int=FALSE,  las=1,
     fun="event", xlab="Years since progression",
     ylab = "Probability of death after progression",
     font=2, font.axis=2, font.lab=2, bty="n", cex.axis=1.5, cex.lab=1.5, axes=FALSE, lwd=3, xaxs="i",yaxs="i") 
title(main="(a) RFC:observed ", cex.main=1.5)
box(lwd=2)
legend(0, 0.8, c("Kaplan-Meier", "Gompertz","generalised gamma","Weibull","Exponential", "log-logistic",  "log normal"), 
       col=c("black", "green", "cyan", "blue", "red", "orange",  "purple" ), bty="n",lwd=2,cex=1.5,
       lty=c("solid","21","16","51", "dotted", "longdash", "twodash"))


plot(survfit(Surv(time,status)~ 1, data=msmcancer3FC), mark.time=F, ylim=c(0,1), xlim=c(0,2.5),conf.int=FALSE,  las=1,
     fun="event", xlab="Years since progression",
     ylab = "Probability of death after progression",
     font=2, font.axis=2, font.lab=2, bty="n", cex.axis=1.5, cex.lab=1.5, axes=FALSE, lwd=3, xaxs="i",yaxs="i") 
axis(1, at = 0:4, labels = c(0:4),lwd=2, lwd.ticks = 2, font=2, cex.axis=1.5)
axis(2, at = c(1, 0.8,0.6, 0.4, 0.2,0), labels = c(1,0.8,0.6, 0.4, 0.2,0),lwd=2,
     lwd.ticks = 2, font=2, las=1, cex.axis=1.5)

polygon(c(survfit(Surv(time,status)~ 1, data=msmcancer3FC)$time, rev(survfit(Surv(time,status)~ 1, data=msmcancer3FC)$time)),
        c(1-survfit(Surv(time,status)~ 1, data=msmcancer3FC)$lower, rev(1-survfit(Surv(time,status)~ 1, data=msmcancer3FC)$upper)), 
        col = "lightgrey", border = NA)
lines(tt15,1-exp_FC_pps15,type='l',col="red", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,4.1),lwd=3,lty="dotted")
lines(tt15,1-wei_FC_pps15,type='l',col="blue", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="51")
lines(tt15,1-gom_FC_pps15,type='l',col="green", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="21")
lines(tt15,1-logn_FC_pps15,type='l',col="purple", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="twodash")
lines(tt15,1-logl_FC_pps15,type='l',col="orange", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="longdash")
lines(tt15,1-gengam_FC_pps15,type='l',col="cyan", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="16")
lines(survfit(Surv(time,status)~ 1, data=msmcancer3FC), mark.time=F, ylim=c(0,1), xlim=c(0,4),conf.int=FALSE,  las=1,
      fun="event", xlab="Years since progression",
      ylab = "Probability of death after progression",
      font=2, font.axis=2, font.lab=2, bty="n", cex.axis=1.5, cex.lab=1.5, axes=FALSE, lwd=3, xaxs="i",yaxs="i") 
box(lwd=2)
title(main="(b) FC: observed ", cex.main=1.5)
legend(0, 0.8, c("Kaplan-Meier", "Gompertz","generalised gamma","Weibull","Exponential", "log-logistic",  "log normal"), 
       col=c("black", "green", "cyan", "blue", "red", "orange",  "purple" ), bty="n",lwd=2,cex=1.5,
       lty=c("solid","21","16","51", "dotted", "longdash", "twodash"))


########################################################################
### OBSERVED VS PREDICTED - EXTRAPOLATION TO 15 YEARS
########################################################################

## =====================================================================
### COLOUR PLOTS
## =====================================================================

par(mar= c(5, 7, 4, 2)+0.1)

# remove the '#' from the line below to show all four plots together

#par(mfrow=c(2,2))

plot(survfit(Surv(time,status)~ 1, data=msmcancer3RFC), mark.time=F, ylim=c(0,1), xlim=c(0,15),conf.int=FALSE,  las=1,
     fun="event", xlab="Years since progression",
     ylab = "Probability of death after progression",
     font=2, font.axis=2, font.lab=2, bty="n", cex.axis=1.5, cex.lab=1.5, axes=FALSE, lwd=3, xaxs="i",yaxs="i") 
axis(1, at = 0:15, lwd=2, lwd.ticks = 2, font=2, cex.axis=1.5)
axis(2, at = c(1, 0.8,0.6, 0.4, 0.2,0), labels = c(1,0.8,0.6, 0.4, 0.2,0),lwd=2, 
     lwd.ticks = 2, font=2, las=1, cex.axis=1.5)

polygon(c(survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$time, rev(survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$time)),
        c(1-survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$lower, rev(1-survfit(Surv(time,status)~ 1, data=msmcancer3RFC)$upper)), 
        col = "lightgrey", border = NA)
lines(tt15,1-exp_RFC_pps15,type='l',col="red", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,4.2),lwd=3,lty="dotted")
lines(tt15,1-wei_RFC_pps15,type='l',col="blue", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="51")
lines(tt15,1-gom_RFC_pps15,type='l',col="green", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="21")
lines(tt15,1-logn_RFC_pps15,type='l',col="purple", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="twodash")
lines(tt15,1-logl_RFC_pps15,type='l',col="orange", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="longdash")
lines(tt15,1-gengam_RFC_pps15,type='l',col="cyan", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="16")
lines(survfit(Surv(time,status)~ 1, data=msmcancer3RFC), mark.time=F, ylim=c(0,1), xlim=c(0,4),conf.int=FALSE,  las=1,
      fun="event", xlab="Years since progression",
      ylab = "Probability of death after progression",
      font=2, font.axis=2, font.lab=2, bty="n", cex.axis=1.5, cex.lab=1.5, axes=FALSE, lwd=3, xaxs="i",yaxs="i") 
title(main="(c) RFC: extrapolation ", cex.main=1.5)
box(lwd=2)
legend("bottomright", c("Kaplan-Meier", "Gompertz","generalised gamma","Weibull","Exponential", "log-logistic",  "log normal"), 
       col=c("black", "green", "cyan", "blue", "red", "orange",  "purple" ), bty="n",lwd=2,cex=1.5,
       lty=c("solid","21","16","51", "dotted", "longdash", "twodash"))


plot(survfit(Surv(time,status)~ 1, data=msmcancer3FC), mark.time=F, ylim=c(0,1), xlim=c(0,15),conf.int=FALSE,  las=1,
     fun="event", xlab="Years since progression",
     ylab = "Probability of death after progression",
     font=2, font.axis=2, font.lab=2, bty="n", cex.axis=1.5, cex.lab=1.5, axes=FALSE, lwd=3, xaxs="i",yaxs="i") 
axis(1, at = 0:15, labels = c(0:15),lwd=2, lwd.ticks = 2, font=2, cex.axis=1.5)
axis(2, at = c(1, 0.8,0.6, 0.4, 0.2,0), labels = c(1,0.8,0.6, 0.4, 0.2,0),lwd=2,
     lwd.ticks = 2, font=2, las=1, cex.axis=1.5)

polygon(c(survfit(Surv(time,status)~ 1, data=msmcancer3FC)$time, rev(survfit(Surv(time,status)~ 1, data=msmcancer3FC)$time)),
        c(1-survfit(Surv(time,status)~ 1, data=msmcancer3FC)$lower, rev(1-survfit(Surv(time,status)~ 1, data=msmcancer3FC)$upper)), 
        col = "lightgrey", border = NA)
lines(tt15,1-exp_FC_pps15,type='l',col="red", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,4.1),lwd=3,lty="dotted")
lines(tt15,1-wei_FC_pps15,type='l',col="blue", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="51")
lines(tt15,1-gom_FC_pps15,type='l',col="green", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="21")
lines(tt15,1-logn_FC_pps15,type='l',col="purple", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="twodash")
lines(tt15,1-logl_FC_pps15,type='l',col="orange", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="longdash")
lines(tt15,1-gengam_FC_pps15,type='l',col="cyan", xlab = "Years since start of study",
      ylab = "Survival probabilities",ylim=c(0,1),xlim=c(0,15),lwd=3,lty="16")
lines(survfit(Surv(time,status)~ 1, data=msmcancer3FC), mark.time=F, ylim=c(0,1), xlim=c(0,4),conf.int=FALSE,  las=1,
      fun="event", xlab="Years since progression",
      ylab = "Probability of death after progression",
      font=2, font.axis=2, font.lab=2, bty="n", cex.axis=1.5, cex.lab=1.5, axes=FALSE, lwd=3, xaxs="i",yaxs="i") 
box(lwd=2)
title(main="(b) FC: extrapolation ", cex.main=1.5)
legend("bottomright", c("Kaplan-Meier", "Gompertz","generalised gamma","Weibull","Exponential", "log-logistic",  "log normal"), 
       col=c("black", "green", "cyan", "blue", "red", "orange",  "purple" ), bty="n",lwd=2,cex=1.5,
       lty=c("solid","21","16","51", "dotted", "longdash", "twodash"))








