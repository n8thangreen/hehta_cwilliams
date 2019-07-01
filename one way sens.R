####### ONE-WAY SENSITIVITY ANALYSES FOR BASE CASE gom1gam2gom3

###### 
#############################################################

#############################################################
####### 20 YEAR TIME HORIZON
#############################################################
### RFC 

ptm <- proc.time() #start counting the time it takes
gom1gam2gom3simRFCSA1<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                                 coveval=rbind(1,1,1),  
                                                 dist=cbind("gom", "gam", "gom"),
                                                 timeseq=seq(0,4,1/12),
                                                 timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), 
                                                 seq(12+1/600,20, 1/600)),
                                                  trans=tmat2, M=5000)
(proc.time() - ptm) #stop counting the time it took (seconds)
(proc.time() - ptm)/60 #stop counting the time it took (minutes)
(proc.time() - ptm)/3600 #stop counting the time it took (hours)
#### TOOK 4.8 MINS

### FC 

ptm <- proc.time() #start counting the time it takes
gom1gam2gom3simFCSA1<-semiMarkov(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                                 coveval=rbind(0,0,0),  
                                                 dist=cbind("gom", "gam", "gom"),
                                                 timeseq=seq(0,4,1/12),
                                                 timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), 
                                                   
                                                                           seq(12+1/600,20, 1/600)),
                                                  trans=tmat2, M=5000)
(proc.time() - ptm) #stop counting the time it took (seconds)
(proc.time() - ptm)/60 #stop counting the time it took (minutes)
(proc.time() - ptm)/3600 #stop counting the time it took (hours)
#### TOOK 4.8 MINS



#############################################################
####### TREATMENT EFFECT NO LONGER PERSISTING IN EXTRAPOLATION PERIOD
#############################################################
source("function semiMarkov_varyHR.R")

### RFC 

ptm <- proc.time() #start counting the time it takes
gom1gam2gom3simRFCSA2<-semiMarkov_varyHR(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                                   coveval=rbind(1,1,1),  
                                                   dist=cbind("gom", "gam", "gom"),
                                                   timeseq=seq(0,4,1/12),
                                                   timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), 
                                                                 seq(12+1/600,15, 1/600)),
                                                   varyHR=TRUE, HRtreat=c(1,1,1),adjust=c(0,0,14),
                                                    trans=tmat2, M=5000)
(proc.time() - ptm) #stop counting the time it took (seconds)
(proc.time() - ptm)/60 #stop counting the time it took (minutes)
(proc.time() - ptm)/3600 #stop counting the time it took (hours)


#which(diff(gom1gam2gom3simRFCSA2$Haz)>1 & diff(gom1gam2gom3simRFCSA2$trans)==0)

## identify timepoints where the difference in consecutive cumulative hazards <0
which(diff(gom1gam2gom3simRFCSA2$Haz)<0 & diff(gom1gam2gom3simRFCSA2$trans)==0)

gom1gam2gom3simRFCSA2[4900:4930,]
#### TOOK 3.3 MINS

### FC 

ptm <- proc.time() #start counting the time it takes
gom1gam2gom3simFCSA2<-semiMarkov_varyHR(ncovs=c(1,1,1),covs=rbind("treat", "treat", "treat"),
                                                  coveval=rbind(0,0,0),  
                                                  dist=cbind("gom", "gam", "gom"),
                                                  timeseq=seq(0,4,1/12),
                                                  timeseq_ext=c(seq(49/12,100/12,1/12), seq(100/12+1/144, 12, 1/144), 
                                                                seq(12+1/600,15, 1/600)),
                                                  varyHR=TRUE, HRtreat=c(1,1,1),adjust=c(0,0,0),
                                                   trans=tmat2, M=5000)
(proc.time() - ptm) #stop counting the time it took (seconds)
(proc.time() - ptm)/60 #stop counting the time it took (minutes)
(proc.time() - ptm)/3600 #stop counting the time it took (hours)

#### TOOK 3.3 MINS



##----------------------------------------------------------
####FUNCTION FOR ICERs 
##-----------------------------------------------------------

ICERs_allcombs<-function(objectssimFC=gom1objectssimFC, objectssimRFC=gom1objectssimRFC, rowentriesnames=gom1entriesnames,  ncom=36, ntt=863,
                        utility1=0.8, utility2=0.6, 
                        cost1aRFC=10113, cost1bRFC=1224, cost1cRFC=2776, cost1dRFC=1109, cost1eRFC=21,cost1fRFC=1109,
                        cost1aFC=0, cost1bFC=0, cost1cFC=2790, cost1dFC=1115, cost1eFC=22,cost1fFC=1115,
                        cost3RFC=0,cost3FC=0, cost4RFC=1,cost4FC=1)
{
  allcombFC<-vector("list", ncom)
  allcombmonFC<-vector("list", ncom)
  allcombRFC<-vector("list", ncom)
  allcombmonRFC<-vector("list", ncom)
  meanLY_PFSFC<-vector("list", ncom)
  meanLY_progFC<-vector("list", ncom)
  meanLY_PFSRFC<-vector("list", ncom)
  meanLY_progRFC<-vector("list", ncom)
  incQALYs<-vector("list", ncom)
  mean_cost_FC<-vector("list", ncom)
  mean_cost_RFC<-vector("list", ncom)
  
  for (i in 1:ncom){
    allcombFC[[i]]<-matrix(nrow=ntt, ncol=8)
    allcombFC[[i]][,1]<-objectssimFC[[i]]$time
    allcombFC[[i]][,2]<-objectssimFC[[i]]$pstate1
    allcombFC[[i]][,3]<-objectssimFC[[i]]$pstate2
    allcombFC[[i]][,4]<-allcombFC[[i]][,2]/(1+0.035)^allcombFC[[i]][,1]
    allcombFC[[i]][,5]<-allcombFC[[i]][,3]/(1+0.035)^allcombFC[[i]][,1]
    allcombFC[[i]][,6]<-allcombFC[[i]][,4]
    allcombFC[[i]][,7]<-allcombFC[[i]][,5]
    allcombFC[[i]][c(1:12),6]<-allcombFC[[i]][c(1:12),2]
    allcombFC[[i]][c(1:12),7]<-allcombFC[[i]][c(1:12),3]
    allcombFC[[i]][,8]<-is.wholenumber(allcombFC[[i]][,1]*12)
    allcombmonFC[[i]]<-subset(allcombFC[[i]], allcombFC[[i]][,8]==1)
    
    ######## AUCs discounted from 1yr onwards 
    meanLY_PFSFC[[i]]<-trapz(allcombmonFC[[i]][,1], allcombmonFC[[i]][,6]) # PFS FC 
    meanLY_progFC[[i]]<-trapz(allcombmonFC[[i]][,1], allcombmonFC[[i]][,7]) # prog FC
    
    allcombRFC[[i]]<-matrix(nrow=ntt, ncol=8)
    allcombRFC[[i]][,1]<-objectssimRFC[[i]]$time
    allcombRFC[[i]][,2]<-objectssimRFC[[i]]$pstate1
    allcombRFC[[i]][,3]<-objectssimRFC[[i]]$pstate2
    allcombRFC[[i]][,4]<-allcombRFC[[i]][,2]/(1+0.035)^allcombRFC[[i]][,1]
    allcombRFC[[i]][,5]<-allcombRFC[[i]][,3]/(1+0.035)^allcombRFC[[i]][,1]
    allcombRFC[[i]][,6]<-allcombRFC[[i]][,4]
    allcombRFC[[i]][,7]<-allcombRFC[[i]][,5]
    allcombRFC[[i]][c(1:12),6]<-allcombRFC[[i]][c(1:12),2]
    allcombRFC[[i]][c(1:12),7]<-allcombRFC[[i]][c(1:12),3]
    allcombRFC[[i]][,8]<-is.wholenumber(allcombRFC[[i]][,1]*12)
    allcombmonRFC[[i]]<-subset(allcombRFC[[i]], allcombRFC[[i]][,8]==1)
    
    
    ######## AUCs discounted from 1yr onwards 
    meanLY_PFSRFC[[i]]<-trapz(allcombmonRFC[[i]][,1], allcombmonRFC[[i]][,6]) # PFS RFC 
    meanLY_progRFC[[i]]<-trapz(allcombmonRFC[[i]][,1], allcombmonRFC[[i]][,7]) # prog RFC
    
    incQALYs[[i]]<-utility1*(meanLY_PFSRFC[[i]]- meanLY_PFSFC[[i]])+utility2*(meanLY_progRFC[[i]]- meanLY_progFC[[i]])
    
    ######## costs
    
    mean_cost_FC[[i]]<-cost3FC+cost1aFC+cost1bFC+cost1cFC+cost1dFC+cost1eFC+cost1fFC+cost4FC*28*12*meanLY_PFSFC[[i]]+360+507+cost4FC*84*12*meanLY_progFC[[i]]+(5179/((meanLY_progRFC[[i]] + meanLY_progFC[[i]])/2*12))*12*meanLY_progFC[[i]]
    mean_cost_RFC[[i]]<-cost3RFC+cost1aRFC+cost1bRFC+cost1cRFC+cost1dRFC+cost1eRFC+cost1fRFC+cost4RFC*28*12*meanLY_PFSRFC[[i]]+592+640+cost4RFC*84*12*meanLY_progRFC[[i]]+(5179/((meanLY_progRFC[[i]] + meanLY_progFC[[i]])/2*12))*12*meanLY_progRFC[[i]]  
  }
  #incQALYs2<-unlist(incQALYs)
  ICERstable<-cbind(unlist(incQALYs),unlist(mean_cost_RFC) -unlist(mean_cost_FC),
                    (unlist(mean_cost_RFC) -unlist(mean_cost_FC))/unlist(incQALYs)) 
  rownames(ICERstable)<-rowentriesnames
  colnames(ICERstable)<-c("incQALYs", "incCosts", "ICER")
  return(ICERstable)
}

#################################################
SAobjects20FC<-vector("list", 1)
SAobjects20FC<-list(gom1gam2gom3simFCSA1)

SAobjects20RFC<-vector("list", 1)
SAobjects20RFC<-list(gom1gam2gom3simRFCSA1)


SAobjectsdistFC<-vector("list", 1)
SAobjectsdistFC<-list(gom1gam2gom3simFCSA2)

SAobjectsdistRFC<-vector("list", 1)
SAobjectsdistRFC<-list(gom1gam2gom3simRFCSA2)

nrow(SAobjects20FC[[1]])
nrow(SAobjectsdistRFC[[1]])

ICERs_20<-ICERs_allcombs(objectssimFC=SAobjects20FC, objectssimRFC=SAobjects20RFC, rowentriesnames=NA,  ncom=1, ntt=5429)
ICERs_dist<-ICERs_allcombs(objectssimFC=SAobjectsdistFC, objectssimRFC=SAobjectsdistRFC, rowentriesnames=NA,  ncom=1, ntt=2429)

#============================================================
####### ONE WAY SENSITIVITY ANALYSES INVOLVING CHANGES TO 
####### UTILITIES
#============================================================
load(file="basecase15gom1gam2gom3.RData")
SAobjectsbasecaseFC<-vector("list", 1)
SAobjectsbasecaseFC<-list(gom1gam2gom3simFC)

SAobjectsbasecaseRFC<-vector("list", 1)
SAobjectsbasecaseRFC<-list(gom1gam2gom3simRFC)

ICERs_util1<-ICERs_allcombs(objectssimFC=SAobjectsbasecaseFC, objectssimRFC=SAobjectsbasecaseRFC, rowentriesnames=NA,  ncom=1, ntt=2429,
                            utility1=0.9, utility2=0.5)
ICERs_util2<-ICERs_allcombs(objectssimFC=SAobjectsbasecaseFC, objectssimRFC=SAobjectsbasecaseRFC, rowentriesnames=NA,  ncom=1, ntt=2429,
                            utility1=0.75, utility2=0.65)

#============================================================
####### ONE WAY SENSITIVITY ANALYSES INVOLVING CHANGES TO 
####### COSTS
#============================================================

ICERs_cost1<-ICERs_allcombs(objectssimFC=SAobjectsbasecaseFC, objectssimRFC=SAobjectsbasecaseRFC, rowentriesnames=NA,  ncom=1, ntt=2429,
                            cost1aRFC=8868, cost1bRFC=947, cost1cRFC=2449, cost1dRFC=2481, cost1eRFC=55,cost1fRFC=2474,
                            cost1aFC=0, cost1bFC=0, cost1cFC=2184, cost1dFC=2212, cost1eFC=53,cost1fFC=2353)
ICERs_cost2<-ICERs_allcombs(objectssimFC=SAobjectsbasecaseFC, objectssimRFC=SAobjectsbasecaseRFC, rowentriesnames=NA,  ncom=1, ntt=2429,
                            cost1aRFC=21150.61, cost1bRFC=0, cost1cRFC=0, cost1dRFC=0, cost1eRFC=0,cost1fRFC=0,
                            cost1aFC=10000, cost1bFC=0, cost1cFC=0, cost1dFC=0, cost1eFC=0,cost1fFC=0)

ICERs_cost3<-ICERs_allcombs(objectssimFC=SAobjectsbasecaseFC, objectssimRFC=SAobjectsbasecaseRFC, rowentriesnames=NA,  ncom=1, ntt=2429,
                            cost3RFC=227.46,cost3FC=144.31818 )
ICERs_cost4<-ICERs_allcombs(objectssimFC=SAobjectsbasecaseFC, objectssimRFC=SAobjectsbasecaseRFC, rowentriesnames=NA,  ncom=1, ntt=2429,
                            cost4RFC=1.5,cost4FC=1.5 )
ICERs_cost5<-ICERs_allcombs(objectssimFC=SAobjectsbasecaseFC, objectssimRFC=SAobjectsbasecaseRFC, rowentriesnames=NA,  ncom=1, ntt=2429,
                            cost4RFC=0.5,cost4FC=0.5 )
ICERs_cost6<-ICERs_allcombs(objectssimFC=SAobjectsbasecaseFC, objectssimRFC=SAobjectsbasecaseRFC, rowentriesnames=NA,  ncom=1, ntt=2429,
                            cost1aRFC=10113, cost1bRFC=2401, cost1cRFC=2776, cost1dRFC=1712, cost1eRFC=21,cost1fRFC=1712,
                            cost1aFC=0, cost1bFC=0, cost1cFC=2790, cost1dFC=1721, cost1eFC=22,cost1fFC=1721)
ICERs_cost7<-ICERs_allcombs(objectssimFC=SAobjectsbasecaseFC, objectssimRFC=SAobjectsbasecaseRFC, rowentriesnames=NA,  ncom=1, ntt=2429,
                            cost1aRFC=10113, cost1bRFC=436, cost1cRFC=2776, cost1dRFC=793, cost1eRFC=21,cost1fRFC=793,
                            cost1aFC=0, cost1bFC=0, cost1cFC=2790, cost1dFC=797, cost1eFC=22,cost1fFC=797)


ICERs_all<-rbind(ICERs_20,ICERs_dist,ICERs_util1,ICERs_util2,
                 ICERs_cost1,ICERs_cost2, ICERs_cost3, ICERs_cost4, ICERs_cost5, ICERs_cost6, ICERs_cost7)
ICERs_all


