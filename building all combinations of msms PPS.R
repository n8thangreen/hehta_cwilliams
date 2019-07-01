
##########################################################################################
##========================================================================================
## CALCULATING STATE OCCUPANCY PROBABILITIES
##========================================================================================
##########################################################################################

### DUE TO COMPUTATIONAL DIFFICULTIES USING GOMPERTZ AND GENERALISED GAMMA FOR THE FIRST
### TRANSITION, MODELS THAT INCLUDE THIS USE MORE TIME POINTS. ALL OTHER MODELS USE THE 
### EQUIVALENT OF MONTHLY CYCLES.

### THIS PROGRAM STARTS THE PREDICTIONS FROM THE PROGRESSION STATE, AND THEREFORE IS 
### USEFUL FOR ASSESSING FITS FOR THE PROGRESSION-> DEATH TRANSITION

##############################
### RFC TREATMENT ARM
#############################

##############################
### exp3
#############################

ptm <- proc.time() #start counting the time it takes
##### exp2exp3
exp1exp2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                      coveval=rbind(1,1,1),  
                      dist=cbind("exp", "exp", "exp"),
                       trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 
(proc.time() - ptm) #stop counting the time it took (seconds)
(proc.time() - ptm)/60 #stop counting the time it took (minutes)
(proc.time() - ptm)/3600 #stop counting the time it took (hours)

wei1exp2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "exp", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1exp2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               dist=cbind("gom", "exp", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

gam1exp2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "exp", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1exp2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("logn", "exp", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1exp2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("logl", "exp", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### wei2exp3
exp1wei2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "wei", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1wei2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "wei", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1wei2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1), 
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               dist=cbind("gom", "wei", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1wei2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "wei", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1wei2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "wei", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1wei2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "wei", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gom2exp3
exp1gom2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "gom", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gom2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "gom", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gom2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               dist=cbind("gom", "gom", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gom2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gom", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gom2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "gom", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gom2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "gom", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gam2exp3
exp1gam2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "gam", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gam2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "gam", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gam2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               dist=cbind("gom", "gam", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gam2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gam", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gam2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "gam", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gam2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "gam", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##### logn2exp3
exp1logn2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "logn", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logn2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "logn", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logn2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               dist=cbind("gom", "logn", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logn2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logn", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logn2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "logn", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logn2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "logn", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### logl2exp3
exp1logl2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("exp", "logl", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logl2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("wei", "logl", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logl2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                dist=cbind("gom", "logl", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logl2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logl", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logl2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logn", "logl", "exp"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logl2exp3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logl", "logl", "exp"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 



##############################
### wei3
#############################

##### exp2wei3
exp1exp2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "exp", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1exp2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "exp", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1exp2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("gom", "exp", "wei"),
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1exp2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "exp", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1exp2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "exp", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1exp2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "exp", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### wei2wei3
exp1wei2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "wei", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1wei2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "wei", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1wei2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               dist=cbind("gom", "wei", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1wei2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "wei", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1wei2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "wei", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1wei2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "wei", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gom2wei3
exp1gom2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "gom", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gom2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "gom", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gom2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1), 
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               dist=cbind("gom", "gom", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gom2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gom", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gom2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "gom", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gom2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "gom", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gam2wei3
exp1gam2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "gam", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gam2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "gam", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gam2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               dist=cbind("gom", "gam", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gam2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gam", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gam2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "gam", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gam2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "gam", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##### logn2wei3
exp1logn2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("exp", "logn", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logn2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("wei", "logn", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logn2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                dist=cbind("gom", "logn", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logn2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logn", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logn2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logn", "logn", "wei"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logn2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logl", "logn", "wei"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### logl2wei3
exp1logl2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("exp", "logl", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logl2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("wei", "logl", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logl2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                dist=cbind("gom", "logl", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logl2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logl", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logl2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logn", "logl", "wei"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logl2wei3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logl", "logl", "wei"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##############################
### gom3
#############################

##### exp2gom3
exp1exp2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "exp", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1exp2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "exp", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1exp2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               dist=cbind("gom", "exp", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1exp2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "exp", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1exp2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "exp", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1exp2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "exp", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### wei2gom3
exp1wei2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "wei", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1wei2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "wei", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1wei2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               dist=cbind("gom", "wei", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1wei2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "wei", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1wei2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "wei", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1wei2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "wei", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gom2gom3
exp1gom2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "gom", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gom2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "gom", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gom2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               dist=cbind("gom", "gom", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gom2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gom", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gom2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "gom", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gom2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "gom", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gam2gom3
exp1gam2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "gam", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gam2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "gam", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gam2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1), 
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               dist=cbind("gom", "gam", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gam2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gam", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gam2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "gam", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gam2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "gam", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##### logn2gom3
exp1logn2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("exp", "logn", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logn2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("wei", "logn", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logn2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                dist=cbind("gom", "logn", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logn2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logn", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logn2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logn", "logn", "gom"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logn2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logl", "logn", "gom"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### logl2gom3
exp1logl2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("exp", "logl", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logl2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("wei", "logl", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logl2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                dist=cbind("gom", "logl", "gom"),
                                
                                
                                
                                
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logl2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logl", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logl2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logn", "logl", "gom"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logl2gom3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logl", "logl", "gom"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##############################
### gam3
#############################

##### exp2gam3
exp1exp2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "exp", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1exp2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "exp", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1exp2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("gom", "exp", "gam"),
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               
                               
                               
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1exp2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "exp", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1exp2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "exp", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1exp2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "exp", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### wei2gam3
exp1wei2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "wei", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1wei2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "wei", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1wei2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("gom", "wei", "gam"),
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               
                               
                               
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1wei2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "wei", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1wei2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "wei", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1wei2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "wei", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gom2gam3
exp1gom2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "gom", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gom2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "gom", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gom2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("gom", "gom", "gam"),
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               
                               
                               
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gom2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gom", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gom2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "gom", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gom2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "gom", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gam2gam3
exp1gam2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "gam", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gam2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "gam", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gam2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("gom", "gam", "gam"),
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                               
                               
                               
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gam2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gam", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gam2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "gam", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gam2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "gam", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##### logn2gam3
exp1logn2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("exp", "logn", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logn2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("wei", "logn", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logn2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("gom", "logn", "gam"),
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                
                                
                                
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logn2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logn", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logn2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logn", "logn", "gam"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logn2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logl", "logn", "gam"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### logl2gam3
exp1logl2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("exp", "logl", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logl2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("wei", "logl", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logl2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("gom", "logl", "gam"),
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                
                                
                                
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logl2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logl", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logl2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logn", "logl", "gam"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logl2gam3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logl", "logl", "gam"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##############################
### logl3
#############################

##### exp2logl3
exp1exp2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "exp", "logl"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1exp2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "exp", "logl"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1exp2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("gom", "exp", "logl"),
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                
                                
                                
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1exp2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "exp", "logl"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1exp2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "exp", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1exp2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "exp", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### wei2logl3
exp1wei2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "wei", "logl"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1wei2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "wei", "logl"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1wei2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("gom", "wei", "logl"),
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                
                                
                                
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1wei2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "wei", "logl"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1wei2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "wei", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1wei2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "wei", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gom2logl3
exp1gom2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "gom", "logl"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gom2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "gom", "logl"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gom2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("gom", "gom", "logl"),
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                
                                
                                
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gom2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gom", "logl"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gom2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "gom", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gom2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "gom", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gam2logl3
exp1gam2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("exp", "gam", "logl"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gam2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("wei", "gam", "logl"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gam2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               dist=cbind("gom", "gam", "logl"),
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)), 
                                
                                
                                
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gam2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(1,1,1),  
                               timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gam", "logl"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gam2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logn", "gam", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gam2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("logl", "gam", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##### logn2logl3
exp1logn2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("exp", "logn", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logn2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("wei", "logn", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logn2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("gom", "logn", "logl"),
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)), 
                                 
                                 
                                 
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logn2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logn", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logn2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logn", "logn", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logn2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logl", "logn", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### logl2logl3
exp1logl2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("exp", "logl", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logl2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("wei", "logl", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logl2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("gom", "logl", "logl"),
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)), 
                                 
                                 
                                 
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logl2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logl", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logl2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logn", "logl", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logl2logl3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logl", "logl", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##############################
### logn3
#############################

##### exp2logn3
exp1exp2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("exp", "exp", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1exp2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("wei", "exp", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1exp2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("gom", "exp", "logn"),
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                
                                
                                
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1exp2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "exp", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1exp2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logn", "exp", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1exp2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logl", "exp", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### wei2logn3
exp1wei2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("exp", "wei", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1wei2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("wei", "wei", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1wei2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("gom", "wei", "logn"),
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                
                                
                                
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1wei2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "wei", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1wei2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logn", "wei", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1wei2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logl", "wei", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gom2logn3
exp1gom2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("exp", "gom", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gom2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("wei", "gom", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gom2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("gom", "gom", "logn"),
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                
                                
                                
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gom2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gom", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gom2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logn", "gom", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gom2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logl", "gom", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gam2logn3
exp1gam2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("exp", "gam", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gam2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("wei", "gam", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gam2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                dist=cbind("gom", "gam", "logn"),
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                
                                
                                
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gam2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(1,1,1),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gam", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gam2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logn", "gam", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gam2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("logl", "gam", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##### logn2logn3
exp1logn2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("exp", "logn", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logn2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("wei", "logn", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logn2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("gom", "logn", "logn"),
                                 timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                 
                                 
                                 
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logn2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logn", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logn2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                  coveval=rbind(1,1,1),  
                                  dist=cbind("logn", "logn", "logn"),
                                   trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logn2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                  coveval=rbind(1,1,1),  
                                  dist=cbind("logl", "logn", "logn"),
                                   trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### logl2logn3
exp1logl2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("exp", "logl", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logl2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("wei", "logl", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logl2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 dist=cbind("gom", "logl", "logn"),
                                 timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),
                                 
                                 
                                 
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logl2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(1,1,1),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logl", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logl2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                  coveval=rbind(1,1,1),  
                                  dist=cbind("logn", "logl", "logn"),
                                   trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logl2logn3simRFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                  coveval=rbind(1,1,1),  
                                  dist=cbind("logl", "logl", "logn"),
                                   trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 



##############################
# ### FC TREATMENT ARM
# #############################
# 
# ##############################
# ### exp3
# #############################
# 
##### exp2exp3
exp1exp2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "exp", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1exp2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "exp", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1exp2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "exp", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1exp2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                        coveval=rbind(0,0,0), 
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)), 
                               dist=cbind("gam", "exp", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logn1exp2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "exp", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1exp2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "exp", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### wei2exp3
exp1wei2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "wei", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1wei2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "wei", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1wei2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "wei", "exp"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1wei2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "wei", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1wei2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "wei", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1wei2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "wei", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gom2exp3
exp1gom2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "gom", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gom2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "gom", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gom2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "gom", "exp"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gom2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gom", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gom2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "gom", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gom2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "gom", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gam2exp3
exp1gam2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "gam", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gam2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "gam", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gam2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "gam", "exp"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gam2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gam", "exp"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gam2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "gam", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gam2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "gam", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##### logn2exp3
exp1logn2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "logn", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logn2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "logn", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logn2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "logn", "exp"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logn2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logn", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logn2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "logn", "exp"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logn2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "logn", "exp"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### logl2exp3
exp1logl2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "logl", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logl2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "logl", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logl2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "logl", "exp"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logl2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logl", "exp"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logl2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "logl", "exp"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logl2exp3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "logl", "exp"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 



##############################
### wei3
#############################

##### exp2wei3
exp1exp2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "exp", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1exp2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "exp", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1exp2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "exp", "wei"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1exp2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "exp", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1exp2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "exp", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1exp2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "exp", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### wei2wei3
exp1wei2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "wei", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1wei2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "wei", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1wei2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "wei", "wei"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1wei2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "wei", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1wei2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "wei", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1wei2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "wei", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gom2wei3
exp1gom2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "gom", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gom2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "gom", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gom2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "gom", "wei"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gom2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gom", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gom2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "gom", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gom2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "gom", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gam2wei3
exp1gam2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "gam", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gam2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "gam", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gam2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "gam", "wei"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gam2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gam", "wei"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gam2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "gam", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gam2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "gam", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##### logn2wei3
exp1logn2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "logn", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logn2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "logn", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logn2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "logn", "wei"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logn2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logn", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logn2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "logn", "wei"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logn2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "logn", "wei"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### logl2wei3
exp1logl2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "logl", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logl2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "logl", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logl2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "logl", "wei"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logl2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logl", "wei"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logl2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "logl", "wei"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logl2wei3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "logl", "wei"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##############################
### gom3
#############################

##### exp2gom3
exp1exp2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "exp", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1exp2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "exp", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1exp2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "exp", "gom"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1exp2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "exp", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1exp2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "exp", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1exp2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "exp", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### wei2gom3
exp1wei2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "wei", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1wei2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "wei", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1wei2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "wei", "gom"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1wei2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "wei", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1wei2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "wei", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1wei2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "wei", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gom2gom3
exp1gom2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "gom", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gom2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "gom", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gom2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "gom", "gom"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gom2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gom", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gom2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "gom", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gom2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "gom", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gam2gom3
exp1gam2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "gam", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gam2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "gam", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gam2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "gam", "gom"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gam2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gam", "gom"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gam2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "gam", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gam2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "gam", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##### logn2gom3
exp1logn2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "logn", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logn2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "logn", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logn2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "logn", "gom"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logn2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logn", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logn2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "logn", "gom"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logn2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "logn", "gom"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### logl2gom3
exp1logl2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "logl", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logl2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "logl", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logl2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "logl", "gom"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logl2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logl", "gom"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logl2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "logl", "gom"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logl2gom3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "logl", "gom"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##############################
### gam3
#############################

##### exp2gam3
exp1exp2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "exp", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1exp2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "exp", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1exp2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "exp", "gam"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1exp2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "exp", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1exp2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "exp", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1exp2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "exp", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### wei2gam3
exp1wei2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "wei", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1wei2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "wei", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1wei2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "wei", "gam"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1wei2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "wei", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1wei2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "wei", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1wei2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "wei", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gom2gam3
exp1gom2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "gom", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gom2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "gom", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gom2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "gom", "gam"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gom2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gom", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gom2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "gom", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gom2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "gom", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gam2gam3
exp1gam2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("exp", "gam", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gam2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               dist=cbind("wei", "gam", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gam2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                               timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "gam", "gam"),
                              
                              
                              
                              
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gam2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                               coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gam", "gam"),
                                trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gam2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logn", "gam", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gam2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("logl", "gam", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##### logn2gam3
exp1logn2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "logn", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logn2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "logn", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logn2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "logn", "gam"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logn2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logn", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logn2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "logn", "gam"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logn2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "logn", "gam"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### logl2gam3
exp1logl2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "logl", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logl2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "logl", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logl2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "logl", "gam"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logl2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logl", "gam"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logl2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "logl", "gam"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logl2gam3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "logl", "gam"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##############################
### logl3
#############################

##### exp2logl3
exp1exp2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "exp", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1exp2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "exp", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1exp2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "exp", "logl"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1exp2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "exp", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1exp2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "exp", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1exp2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "exp", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### wei2logl3
exp1wei2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "wei", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1wei2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "wei", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1wei2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "wei", "logl"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1wei2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "wei", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1wei2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "wei", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1wei2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "wei", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gom2logl3
exp1gom2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "gom", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gom2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "gom", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gom2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "gom", "logl"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gom2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gom", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gom2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "gom", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gom2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "gom", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gam2logl3
exp1gam2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "gam", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gam2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "gam", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gam2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "gam", "logl"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gam2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gam", "logl"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gam2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "gam", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gam2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "gam", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##### logn2logl3
exp1logn2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("exp", "logn", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logn2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("wei", "logn", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logn2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "logn", "logl"),
                                
                                
                                
                                
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logn2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                  timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logn", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logn2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                  coveval=rbind(0,0,0),  
                                  dist=cbind("logn", "logn", "logl"),
                                   trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logn2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                  coveval=rbind(0,0,0),  
                                  dist=cbind("logl", "logn", "logl"),
                                   trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### logl2logl3
exp1logl2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("exp", "logl", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logl2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("wei", "logl", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logl2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "logl", "logl"),
                                
                                
                                
                                
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logl2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                  timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logl", "logl"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logl2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                  coveval=rbind(0,0,0),  
                                  dist=cbind("logn", "logl", "logl"),
                                   trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logl2logl3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                  coveval=rbind(0,0,0),  
                                  dist=cbind("logl", "logl", "logl"),
                                   trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##############################
### logn3
#############################

##### exp2logn3
exp1exp2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "exp", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1exp2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "exp", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1exp2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "exp", "logn"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1exp2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "exp", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1exp2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "exp", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1exp2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "exp", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### wei2logn3
exp1wei2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "wei", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1wei2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "wei", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1wei2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "wei", "logn"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1wei2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "wei", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1wei2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "wei", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1wei2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "wei", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gom2logn3
exp1gom2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "gom", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gom2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "gom", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gom2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "gom", "logn"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gom2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gom", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gom2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "gom", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gom2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "gom", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### gam2logn3
exp1gam2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("exp", "gam", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1gam2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                dist=cbind("wei", "gam", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1gam2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "gam", "logn"),
                               
                               
                               
                               
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1gam2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "gam", "logn"),
                                 trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1gam2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logn", "gam", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1gam2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("logl", "gam", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


##### logn2logn3

exp1logn2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("exp", "logn", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logn2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("wei", "logn", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logn2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "logn", "logn"),
                                
                                
                                
                                
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gam1logn2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                  timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logn", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logn2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                  coveval=rbind(0,0,0),  
                                  dist=cbind("logn", "logn", "logn"),
                                   trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logn2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                  coveval=rbind(0,0,0),  
                                  dist=cbind("logl", "logn", "logn"),
                                   trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

##### logl2logn3
exp1logl2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("exp", "logl", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

wei1logl2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 dist=cbind("wei", "logl", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


gom1logl2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                 timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), seq(12+1/600,15, 1/600)),            dist=cbind("gom", "logl", "logn"),
                          trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 



gam1logl2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                 coveval=rbind(0,0,0),  
                                  timeseq_ext=c(seq(49/12,144/12,1/12), seq(144/12+1/144, 15, 1/144)),  dist=cbind("gam", "logl", "logn"),
                                  trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 

logn1logl2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                  coveval=rbind(0,0,0),  
                                  dist=cbind("logn", "logl", "logn"),
                                   trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 


logl1logl2logn3simFC2<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                  coveval=rbind(0,0,0),  
                                  dist=cbind("logl", "logl", "logn"),
                                   trans=tmat2, M=5000, predinitial=FALSE, predfrom=2) 





















