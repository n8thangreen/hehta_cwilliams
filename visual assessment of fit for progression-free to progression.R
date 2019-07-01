### the data for the observed proportion is not shared for the tutorial.
### it needs  the actual data rather than digitised data
### However it is shown in the plots in the appendix.


load("allcombos.RData")

####--------------------------------------------------------------

####--------------------------------------------------------------
### PROGRESSION-FREE -> PROGRESSION OVER TRIAL OBSERVATION PERIOD

### WITH PROG->DEATH FITTED AS GOMPERTZ

#### FOR RFC 

####--------------------------------------------------------------

length(gom1exp2gom3simRFC$time)
length(gam1exp2gom3simRFC$time)

pdf('Appendix1_FigureA1_timeinprog_obs_RFC.pdf', width=8, height=8)

par(mar= c(5, 7, 4, 2)+0.1)

visualsemiMarkov(nobjects=36, 
                 objects=rbind(exp1exp2gom3simRFC, exp1wei2gom3simRFC, exp1gom2gom3simRFC,
                               exp1gam2gom3simRFC, exp1logn2gom3simRFC, exp1logl2gom3simRFC,
                               wei1exp2gom3simRFC, wei1wei2gom3simRFC,wei1gom2gom3simRFC,
                               wei1gam2gom3simRFC, wei1logn2gom3simRFC,wei1logl2gom3simRFC,
                               gam1exp2gom3simRFC, gam1wei2gom3simRFC, gam1gom2gom3simRFC,
                               gam1gam2gom3simRFC, gam1logn2gom3simRFC, gam1logl2gom3simRFC,
                               logn1exp2gom3simRFC, logn1wei2gom3simRFC, logn1gom2gom3simRFC,
                               logn1gam2gom3simRFC, logn1logn2gom3simRFC, logn1logl2gom3simRFC,
                               logl1exp2gom3simRFC, logl1wei2gom3simRFC, logl1gom2gom3simRFC,
                               logl1gam2gom3simRFC, logl1logn2gom3simRFC, logl1logl2gom3simRFC,
                               gom1exp2gom3simRFC, gom1wei2gom3simRFC, gom1gom2gom3simRFC,
                               gom1gam2gom3simRFC, gom1logn2gom3simRFC, gom1logl2gom3simRFC),
                 tteach=c(rep(181,12), rep(577,6), rep(181,12),rep(2429,6)),
                 state=2, observed=NA,
                lty=rep(1,36),
                 col=c("darkred","red", "indianred3","indianred","indianred2","red4",
                       "blue","dodgerblue","deepskyblue","lightblue" ,"midnightblue","steelblue3",
                       "lightgrey", "grey","grey22","grey44","grey82", "lightslategrey",
                       "lightpink2","lightpink","pink","lightpink3","lightpink1","pink2",
                       "darkorange","orange3","darkorange2","orange","orange2","darkorange4",
                       "darkgreen","green","lawngreen","lightgreen","springgreen" ,"limegreen"),
                 xlim=c(0,4),ylim=c(0,0.6),lwd=3,
                 ylab = "Probability of being in progression",
                 main= "Probability of being in progression",
                legendpos=c(0.25, 0.475), legendcex=1.25, 
                legendcurves=c("Gompertz", "generalised gamma",  "Weibull", "log-logistic", "log normal", "exponential" ),
                legendcol=c( "green", "grey44",  "blue", "orange", "pink" ,"red"), 
                legendbty=rep("n",6),legendlwd=rep(10,6),legendlty=rep("dotted",6)
                )

text(1.5,0.58, "progression -> death fitted as Gompertz", cex=1.25, font=2)
text(1.75,0.5, "progression-free -> progression", cex=1.25, font=2)
text(0.75,0.48, "fitted as:", cex=1.25, font=2)
lines(RFC_prevCI5000$time, RFC_prevCI5000$prev, type="l", lwd=3,lty=1, col="black") # this line requires the actual trial data
legend(0,0.2, "observed proportion", lwd=3 , bty="n")

dev.off()


####--------------------------------------------------------------
### PROGRESSION-FREE -> PROGRESSION EXTRAPOLATED TO 15 YEARS

### WITH PROG->DEATH FITTED AS GOMPERTZ

#### FOR RFC 

####--------------------------------------------------------------

pdf('Appendix1_FigureA2_timeinprog_ext_RFC.pdf', width=8, height=8)

par(mar= c(5, 7, 4, 2)+0.1)

visualsemiMarkov(nobjects=36, 
                 objects=rbind(exp1exp2gom3simRFC, exp1wei2gom3simRFC, exp1gom2gom3simRFC,
                               exp1gam2gom3simRFC, exp1logn2gom3simRFC, exp1logl2gom3simRFC,
                               wei1exp2gom3simRFC, wei1wei2gom3simRFC,wei1gom2gom3simRFC,
                               wei1gam2gom3simRFC, wei1logn2gom3simRFC,wei1logl2gom3simRFC,
                               gam1exp2gom3simRFC, gam1wei2gom3simRFC, gam1gom2gom3simRFC,
                               gam1gam2gom3simRFC, gam1logn2gom3simRFC, gam1logl2gom3simRFC,
                               logn1exp2gom3simRFC, logn1wei2gom3simRFC, logn1gom2gom3simRFC,
                               logn1gam2gom3simRFC, logn1logn2gom3simRFC, logn1logl2gom3simRFC,
                               logl1exp2gom3simRFC, logl1wei2gom3simRFC, logl1gom2gom3simRFC,
                               logl1gam2gom3simRFC, logl1logn2gom3simRFC, logl1logl2gom3simRFC,
                               gom1exp2gom3simRFC, gom1wei2gom3simRFC, gom1gom2gom3simRFC,
                               gom1gam2gom3simRFC, gom1logn2gom3simRFC, gom1logl2gom3simRFC),
                 tteach=c(rep(181,12), rep(577,6), rep(181,12),rep(2429,6)),
                 state=2, observed=NA,
                 lty=rep(1,36),
                 col=c("darkred","red", "indianred3","indianred","indianred2","red4",
                       "blue","dodgerblue","deepskyblue","lightblue" ,"midnightblue","steelblue3",
                       "lightgrey", "grey","grey22","grey44","grey82", "lightslategrey",
                       "lightpink2","lightpink","pink","lightpink3","lightpink1","pink2",
                       "darkorange","orange3","darkorange2","orange","orange2","darkorange4",
                       "darkgreen","green","lawngreen","lightgreen","springgreen" ,"limegreen"),
                 xlim=c(0,15),ylim=c(0,0.6),lwd=3,
                 ylab = "Probability of being in progression",
                 main= "Probability of being in progression",
                 legendpos=c(8, 0.475), legendcex=1.25, 
                 legendcurves=c("Gompertz", "generalised gamma",  "Weibull", "log-logistic", "log normal", "exponential" ),
                 legendcol=c( "green", "grey44",  "blue", "orange", "pink" ,"red"), 
                 legendbty=rep("n",6),legendlwd=rep(10,6),legendlty=rep("dotted",6)
)

text(8,0.58, "progression -> death fitted as Gompertz ", cex=1.25, font=2)

text(10,0.5, "progression-free -> progression", cex=1.25, font=2)
text(9,0.48, "fitted as:", cex=1.25, font=2)

legend(8, 0.475 , c("Gompertz", "generalised gamma",  "Weibull", "log-logistic", "log normal", "exponential" ),
       col=c( "green", "grey44",  "blue", "orange", "pink" ,"red"), 
       bty="n",lwd=10, lty="dotted", cex=1.25)
lines(RFC_prevCI5000$time, RFC_prevCI5000$prev, type="l", lwd=3,lty=1, col="black") # this line requires the actual trial data
legend(5,0.1, "observed proportion", lwd=3 , bty="n")
dev.off()


####--------------------------------------------------------------

####--------------------------------------------------------------
### PROGRESSION-FREE -> PROGRESSION OVER TRIAL OBSERVATION PERIOD

### WITH PROG->DEATH FITTED AS GOMPERTZ

#### FOR FC 

####--------------------------------------------------------------

par(mar= c(5, 7, 4, 2)+0.1)

visualsemiMarkov(nobjects=36, 
                 objects=rbind(exp1exp2gom3simFC, exp1wei2gom3simFC, exp1gom2gom3simFC,
                               exp1gam2gom3simFC, exp1logn2gom3simFC, exp1logl2gom3simFC,
                               wei1exp2gom3simFC, wei1wei2gom3simFC,wei1gom2gom3simFC,
                               wei1gam2gom3simFC, wei1logn2gom3simFC,wei1logl2gom3simFC,
                               gam1exp2gom3simFC, gam1wei2gom3simFC, gam1gom2gom3simFC,
                               gam1gam2gom3simFC, gam1logn2gom3simFC, gam1logl2gom3simFC,
                               logn1exp2gom3simFC, logn1wei2gom3simFC, logn1gom2gom3simFC,
                               logn1gam2gom3simFC, logn1logn2gom3simFC, logn1logl2gom3simFC,
                               logl1exp2gom3simFC, logl1wei2gom3simFC, logl1gom2gom3simFC,
                               logl1gam2gom3simFC, logl1logn2gom3simFC, logl1logl2gom3simFC,
                               gom1exp2gom3simFC, gom1wei2gom3simFC, gom1gom2gom3simFC,
                               gom1gam2gom3simFC, gom1logn2gom3simFC, gom1logl2gom3simFC),
                 tteach=c(rep(181,12), rep(577,6), rep(181,12),rep(2429,6)),
                 state=2, observed=NA,
                 lty=rep(1,36),
                 col=c("darkred","red", "indianred3","indianred","indianred2","red4",
                       "blue","dodgerblue","deepskyblue","lightblue" ,"midnightblue","steelblue3",
                       "lightgrey", "grey","grey22","grey44","grey82", "lightslategrey",
                       "lightpink2","lightpink","pink","lightpink3","lightpink1","pink2",
                       "darkorange","orange3","darkorange2","orange","orange2","darkorange4",
                       "darkgreen","green","lawngreen","lightgreen","springgreen" ,"limegreen"),
                 xlim=c(0,4),ylim=c(0,0.6),lwd=3,
                 ylab = "Probability of being in progression",
                 main= "Probability of being in progression",
                 legendpos=c(0.25, 0.475), legendcex=1.25, 
                 legendcurves=c("Gompertz", "generalised gamma",  "Weibull", "log-logistic", "log normal", "exponential" ),
                 legendcol=c( "green", "grey44",  "blue", "orange", "pink" ,"red"), 
                 legendbty=rep("n",6),legendlwd=rep(10,6),legendlty=rep("dotted",6)
)

text(1.25,0.59, "progression -> death fitted as Gompertz", cex=1.25, font=1)
text(1.75,0.5, "progression-free -> progression", cex=1.25, font=2)
text(0.75,0.48, "fitted as:", cex=1.25, font=2)


####--------------------------------------------------------------
### PROGRESSION-FREE -> PROGRESSION EXTRAPOLATED TO 15 YEARS

### WITH PROG->DEATH FITTED AS GOMPERTZ
#### FOR FC 

####--------------------------------------------------------------

par(mar= c(5, 7, 4, 2)+0.1)

visualsemiMarkov(nobjects=36, 
                 objects=rbind(exp1exp2gom3simFC, exp1wei2gom3simFC, exp1gom2gom3simFC,
                               exp1gam2gom3simFC, exp1logn2gom3simFC, exp1logl2gom3simFC,
                               wei1exp2gom3simFC, wei1wei2gom3simFC,wei1gom2gom3simFC,
                               wei1gam2gom3simFC, wei1logn2gom3simFC,wei1logl2gom3simFC,
                               gam1exp2gom3simFC, gam1wei2gom3simFC, gam1gom2gom3simFC,
                               gam1gam2gom3simFC, gam1logn2gom3simFC, gam1logl2gom3simFC,
                               logn1exp2gom3simFC, logn1wei2gom3simFC, logn1gom2gom3simFC,
                               logn1gam2gom3simFC, logn1logn2gom3simFC, logn1logl2gom3simFC,
                               logl1exp2gom3simFC, logl1wei2gom3simFC, logl1gom2gom3simFC,
                               logl1gam2gom3simFC, logl1logn2gom3simFC, logl1logl2gom3simFC,
                               gom1exp2gom3simFC, gom1wei2gom3simFC, gom1gom2gom3simFC,
                               gom1gam2gom3simFC, gom1logn2gom3simFC, gom1logl2gom3simFC),
                 tteach=c(rep(181,12), rep(577,6), rep(181,12),rep(2429,6)),
                 state=2, observed=NA,
                 lty=rep(1,36),
                 col=c("darkred","red", "indianred3","indianred","indianred2","red4",
                       "blue","dodgerblue","deepskyblue","lightblue" ,"midnightblue","steelblue3",
                       "lightgrey", "grey","grey22","grey44","grey82", "lightslategrey",
                       "lightpink2","lightpink","pink","lightpink3","lightpink1","pink2",
                       "darkorange","orange3","darkorange2","orange","orange2","darkorange4",
                       "darkgreen","green","lawngreen","lightgreen","springgreen" ,"limegreen"),
                 xlim=c(0,15),ylim=c(0,0.6),lwd=3,
                 ylab = "Probability of being in progression",
                 main= "Probability of being in progression",
                 legendpos=c(8, 0.475), legendcex=1.25, 
                 legendcurves=c("Gompertz", "generalised gamma",  "Weibull", "log-logistic", "log normal", "exponential" ),
                 legendcol=c( "green", "grey44",  "blue", "orange", "pink" ,"red"), 
                 legendbty=rep("n",6),legendlwd=rep(10,6),legendlty=rep("dotted",6)
)

text(8,0.59, "progression -> death fitted as Gompertz ", cex=1.25, font=1)
text(10,0.5, "progression-free -> progression", cex=1.25, font=2)
text(9,0.48, "fitted as:", cex=1.25, font=2)

legend(8, 0.475 , c("Gompertz", "generalised gamma",  "Weibull", "log-logistic", "log normal", "exponential" ),
       col=c( "green", "grey44",  "blue", "orange", "pink" ,"red"), 
       bty="n",lwd=10, lty="dotted", cex=1.25)


