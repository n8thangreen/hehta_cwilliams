####--------------------------------------------------------------

RFC_CR<-cuminc(RFCdata$progdeath_ty, RFCdata$crstatus)
FC_CR<-cuminc(FCdata$progdeath_ty, FCdata$crstatus)
####--------------------------------------------------------------
### PROGRESSION-FREE -> DEATH OVER TRIAL OBSERVED PERIOD

### WITH PROG->DEATH FITTED AS GOMPERTZ
### WITH PROGRESSION-FREE -> PROGRESSION FITTED AS GOMPERTZ


#### FOR RFC 

####--------------------------------------------------------------

pdf('Appendix1_FigureA3_progfreetodeath_obs_RFC.pdf', width=8, height=8)
par(mar= c(5, 7, 4, 2)+0.1)

visualsemiMarkov(nobjects=6,
                 objects=rbind(gom1exp2gom3simRFC, gom1wei2gom3simRFC, gom1gom2gom3simRFC,
                               gom1gam2gom3simRFC, gom1logn2gom3simRFC, gom1logl2gom3simRFC),
                 tteach=rep(2429,6),
                 state=4, observed="cuminc",obsdata=RFC_CR,CI=FALSE,
                 lty=c("dotted", "51", "longdash","16", "twodash","21"), col=c("red" ,"blue"  ,"green" ,"cyan"  ,"purple","orange"),
                 xlim=c(0,4),ylim=c(0,0.1),lwd=3,
                 ylab = "", main= "",
                 legendpos=c(2.0, 0.03), legendcex=1.25,
                 legendcurves=c("exponential", "Weibull", "log-logistic", "generalised gamma","log normal", "Gompertz" ),
                 legendcol=c("red", "blue", "orange", "cyan", "purple", "green"),
                 legendbty=rep("n",6),legendlwd=rep(2,6),legendlty=c("dotted", "51", "longdash","16", "twodash","21"))

title(  ylab = "Probability of progression-free-> death without progression", cex.lab = 1.5,
        line = 4)
text(1.75,0.098, "progression-free -> progression fitted as Gompertz", cex=1.25, font=1)
text(1.5,0.090, "progression -> death fitted as Gompertz", cex=1.25, font=1)
text(3,0.037, "progression-free -> death", cex=1.25, font=2)
text(2.5,0.032, "fitted as:", cex=1.25, font=2)

# legend(0, 0.07, c("CR cumulative incidence estimate", "95% CI"),
#        col=c("black", "lightgrey"), 
#        bty="n",lwd=c(2, 10), cex=1.25)

legend(2.0, 0.03 , c("Gompertz","log normal", "generalised gamma","exponential", "Weibull", "log-logistic" ),
       col=c("green", "purple","cyan","red", "blue", "orange"  ),
       bty="n",lwd=2, lty=c("21", "twodash", "16","dotted", "51", "longdash" ),cex=1.25)

legend(0, 0.05, "observed", lwd=3, bty="n")
dev.off()

####--------------------------------------------------------------

####--------------------------------------------------------------
### PROGRESSION-FREE -> DEATH EXTRAPOLATED TO 15 YEARS

### WITH PROG->DEATH FITTED AS GOMPERTZ
### WITH PROGRESSION-FREE -> PROGRESSION FITTED AS GOMPERTZ


#### FOR RFC 

####--------------------------------------------------------------

pdf('Appendix1_FigureA4_progfreetodeath_ext_RFC.pdf', width=8, height=8)

par(mar= c(5, 7, 4, 2)+0.1)

visualsemiMarkov(nobjects=6,
                 objects=rbind(gom1exp2gom3simRFC, gom1wei2gom3simRFC, gom1gom2gom3simRFC,
                               gom1gam2gom3simRFC, gom1logn2gom3simRFC, gom1logl2gom3simRFC),
                 tteach=rep(2429,6),
                 state=4, observed="cuminc",obsdata=RFC_CR,CI=FALSE,
                 lty=c("dotted", "51", "longdash","16", "twodash","21"), col=c("red" ,"blue"  ,"green" ,"cyan"  ,"purple","orange"),
                 xlim=c(0,15),ylim=c(0,0.11),lwd=3,
                 ylab = "", main= "",
                 legendpos=c(6.0, 0.03), legendcex=1.25,
                 legendcurves=c("exponential", "Weibull", "log-logistic", "generalised gamma","log normal", "Gompertz" ),
                 legendcol=c("red", "blue", "orange", "cyan", "purple", "green"),
                 legendbty=rep("n",6),legendlwd=rep(2,6),legendlty=c("dotted", "51", "longdash","16", "twodash","21"))

title(  ylab = "Probability of progression-free-> death without progression", cex.lab = 1.5,
        line = 4)
text(6,0.108, "progression-free -> progression fitted as Gompertz", cex=1.25, font=1)
text(5,0.098, "progression -> death fitted as Gompertz", cex=1.25, font=1)
text(7.5,0.035, "progression-free -> death", cex=1.25, font=2)
text(6,0.031, "fitted as:", cex=1.25, font=2)

# legend(4, 0.07, c("CR cumulative incidence estimate", "95% CI"),
#        col=c("black", "lightgrey"), 
#        bty="n",lwd=c(2, 10), cex=1.25)

legend(6.0, 0.03 , c("log normal","Gompertz","exponential",  "generalised gamma","Weibull", "log-logistic" ),
       col=c("purple","green", "red","cyan", "blue", "orange"  ),
       bty="n",lwd=2, lty=c("twodash","21", "dotted", "16", "51", "longdash" ),cex=1.25)

legend(5, 0.05, "observed", lwd=3, bty="n")
dev.off()


####--------------------------------------------------------------
### PROGRESSION-FREE -> DEATH OVER TRIAL OBSERVED PERIOD

### WITH PROG->DEATH FITTED AS GOMPERTZ
### WITH PROGRESSION-FREE -> PROGRESSION FITTED AS GOMPERTZ


#### FOR FC 

####--------------------------------------------------------------

par(mar= c(5, 7, 4, 2)+0.1)

visualsemiMarkov(nobjects=6,
                 objects=rbind(gom1exp2gom3simFC, gom1wei2gom3simFC, gom1gom2gom3simFC,
                               gom1gam2gom3simFC, gom1logn2gom3simFC, gom1logl2gom3simFC),
                 tteach=rep(2429,6),
                 state=4, observed="cuminc",obsdata=FC_CR,
                 CRtime=FCdata$progdeath_ty, CRstatus=FCdata$crstatus,
                 lty=c("dotted", "51", "longdash","16", "twodash","21"), col=c("red" ,"blue"  ,"green" ,"cyan"  ,"purple","orange"),
                 xlim=c(0,4),ylim=c(0,0.1),lwd=3,
                 ylab = "", main= "",
                 legendpos=c(2.0, 0.03), legendcex=1.25,
                 legendcurves=c("exponential", "Weibull", "log-logistic", "generalised gamma","log normal", "Gompertz" ),
                 legendcol=c("red", "blue", "orange", "cyan", "purple", "green"),
                 legendbty=rep("n",6),legendlwd=rep(2,6),legendlty=c("dotted", "51", "longdash","16", "twodash","21"))

title(  ylab = "Probability of progression-free-> death without progression", cex.lab = 1.5,
        line = 4)
text(1.75,0.098, "progression-free -> progression fitted as Gompertz", cex=1.25, font=1)
text(1,0.090, "progression -> death fitted as Gompertz", cex=1.25, font=1)
text(3,0.037, "progression-free -> death", cex=1.25, font=2)
text(2.5,0.032, "fitted as:", cex=1.25, font=2)

legend(0, 0.082, c("CR cumulative incidence estimate", "95% CI"),
       col=c("black", "lightgrey"), 
       bty="n",lwd=c(2, 10), cex=1.25)

legend(2.0, 0.03 , c("exponential", "Weibull", "log-logistic", "generalised gamma","log normal", "Gompertz" ),
       col=c("red", "blue", "orange", "cyan", "purple", "green"),
       bty="n",lwd=2, lty=c("dotted", "51", "longdash","16", "twodash","21"),cex=1.25)


####--------------------------------------------------------------

####--------------------------------------------------------------
### PROGRESSION-FREE -> DEATH EXTRAPOLATED TO 15 YEARS

### WITH PROG->DEATH FITTED AS GOMPERTZ
### WITH PROGRESSION-FREE -> PROGRESSION FITTED AS GOMPERTZ


#### FOR FC 

####--------------------------------------------------------------

par(mar= c(5, 7, 4, 2)+0.1)

visualsemiMarkov(nobjects=6,
                 objects=rbind(gom1exp2gom3simFC, gom1wei2gom3simFC, gom1gom2gom3simFC,
                               gom1gam2gom3simFC, gom1logn2gom3simFC, gom1logl2gom3simFC),
                 tteach=rep(2429,6),
                 state=4, observed="cuminc",obsdata=FC_CR,
                 CRtime=FCdata$progdeath_ty, CRstatus=FCdata$crstatus,
                 lty=c("dotted", "51", "longdash","16", "twodash","21"), col=c("red" ,"blue"  ,"green" ,"cyan"  ,"purple","orange"),
                 xlim=c(0,15),ylim=c(0,0.1),lwd=3,
                 ylab = "", main= "",
                 legendpos=c(6.0, 0.03), legendcex=1.25,
                 legendcurves=c("exponential", "Weibull", "log-logistic", "generalised gamma","log normal", "Gompertz" ),
                 legendcol=c("red", "blue", "orange", "cyan", "purple", "green"),
                 legendbty=rep("n",6),legendlwd=rep(2,6),legendlty=c("dotted", "51", "longdash","16", "twodash","21"))

title(  ylab = "Probability of progression-free-> death without progression", cex.lab = 1.5,
        line = 4)
text(9,0.078, "progression-free -> progression fitted as Gompertz", cex=1.25, font=1)
text(7,0.073, "progression -> death fitted as Gompertz", cex=1.25, font=1)
text(7.5,0.035, "progression-free -> death", cex=1.25, font=2)
text(6,0.031, "fitted as:", cex=1.25, font=2)

legend(4, 0.05, c("CR cumulative incidence estimate", "95% CI"),
       col=c("black", "lightgrey"), 
       bty="n",lwd=c(2, 10), cex=1.25)

legend(6.0, 0.03 , c("exponential", "Weibull", "log-logistic", "generalised gamma","log normal", "Gompertz" ),
       col=c("red", "blue", "orange", "cyan", "purple", "green"),
       bty="n",lwd=2, lty=c("dotted", "51", "longdash","16", "twodash","21"),cex=1.25)


























