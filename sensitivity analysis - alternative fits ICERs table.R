#########ALTERNATIVE FITS FOR MULTI-STATE MODELLING

load("allcombos.RData") 

is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

gom1entriesnames<-c( " gom1gom2gom3 "  ,
                     " gom1exp2exp3 "  ,
                     " gom1wei2exp3 "  ,
                     " gom1gom2exp3 "  ,
                     " gom1gam2exp3 "  ,
                     " gom1logn2exp3 "  ,
                     " gom1logl2exp3 "  ,
                     " gom1exp2wei3 "  ,
                     " gom1wei2wei3 "  ,
                     " gom1gom2wei3 "  ,
                     " gom1gam2wei3 "  ,
                     " gom1logn2wei3 "  ,
                     " gom1logl2wei3 "  ,
                     " gom1exp2gom3 "  ,
                     " gom1wei2gom3 "  ,
                     " gom1gam2gom3 "  ,
                     " gom1logn2gom3 "  ,
                     " gom1logl2gom3 "  ,
                     " gom1exp2gam3 "  ,
                     " gom1wei2gam3 "  ,
                     " gom1gom2gam3 "  ,
                     " gom1gam2gam3 "  ,
                     " gom1logn2gam3 "  ,
                     " gom1logl2gam3 "  ,
                     " gom1exp2logn3 "  ,
                     " gom1wei2logn3 "  ,
                     " gom1gom2logn3 "  ,
                     " gom1gam2logn3 "  ,
                     " gom1logn2logn3 "  ,
                     " gom1logl2logn3 "  ,
                     " gom1exp2logl3 "  ,
                     " gom1wei2logl3 "  ,
                     " gom1gom2logl3 "  ,
                     " gom1gam2logl3 "  ,
                     " gom1logn2logl3 "  ,
                     " gom1logl2logl3 ")


gom1objectssimFC<-vector("list", 36)
gom1objectssimFC<-list(gom1gom2gom3simFC ,
                       gom1exp2exp3simFC ,
                       gom1wei2exp3simFC ,
                       gom1gom2exp3simFC ,
                       gom1gam2exp3simFC ,
                       gom1logn2exp3simFC ,
                       gom1logl2exp3simFC ,
                       gom1exp2wei3simFC ,
                       gom1wei2wei3simFC ,
                       gom1gom2wei3simFC ,
                       gom1gam2wei3simFC ,
                       gom1logn2wei3simFC ,
                       gom1logl2wei3simFC ,
                       gom1exp2gom3simFC ,
                       gom1wei2gom3simFC ,
                       gom1gam2gom3simFC ,
                       gom1logn2gom3simFC ,
                       gom1logl2gom3simFC ,
                       gom1exp2gam3simFC ,
                       gom1wei2gam3simFC ,
                       gom1gom2gam3simFC ,
                       gom1gam2gam3simFC ,
                       gom1logn2gam3simFC ,
                       gom1logl2gam3simFC ,
                       gom1exp2logn3simFC ,
                       gom1wei2logn3simFC ,
                       gom1gom2logn3simFC ,
                       gom1gam2logn3simFC ,
                       gom1logn2logn3simFC ,
                       gom1logl2logn3simFC ,
                       gom1exp2logl3simFC ,
                       gom1wei2logl3simFC ,
                       gom1gom2logl3simFC ,
                       gom1gam2logl3simFC ,
                       gom1logn2logl3simFC ,
                       gom1logl2logl3simFC)

gom1objectssimRFC<-vector("list", 36)
gom1objectssimRFC<-list(gom1gom2gom3simRFC ,
                        gom1exp2exp3simRFC ,
                        gom1wei2exp3simRFC ,
                        gom1gom2exp3simRFC ,
                        gom1gam2exp3simRFC ,
                        gom1logn2exp3simRFC ,
                        gom1logl2exp3simRFC ,
                        gom1exp2wei3simRFC ,
                        gom1wei2wei3simRFC ,
                        gom1gom2wei3simRFC ,
                        gom1gam2wei3simRFC ,
                        gom1logn2wei3simRFC ,
                        gom1logl2wei3simRFC ,
                        gom1exp2gom3simRFC ,
                        gom1wei2gom3simRFC ,
                        gom1gam2gom3simRFC ,
                        gom1logn2gom3simRFC ,
                        gom1logl2gom3simRFC ,
                        gom1exp2gam3simRFC ,
                        gom1wei2gam3simRFC ,
                        gom1gom2gam3simRFC ,
                        gom1gam2gam3simRFC ,
                        gom1logn2gam3simRFC ,
                        gom1logl2gam3simRFC ,
                        gom1exp2logn3simRFC ,
                        gom1wei2logn3simRFC ,
                        gom1gom2logn3simRFC ,
                        gom1gam2logn3simRFC ,
                        gom1logn2logn3simRFC ,
                        gom1logl2logn3simRFC ,
                        gom1exp2logl3simRFC ,
                        gom1wei2logl3simRFC ,
                        gom1gom2logl3simRFC ,
                        gom1gam2logl3simRFC ,
                        gom1logn2logl3simRFC ,
                        gom1logl2logl3simRFC)

gam1entriesnames<-c( " gam1gom2gom3 "  ,
                     " gam1exp2exp3 "  ,
                     " gam1wei2exp3 "  ,
                     " gam1gom2exp3 "  ,
                     " gam1gam2exp3 "  ,
                     " gam1logn2exp3 "  ,
                     " gam1logl2exp3 "  ,
                     " gam1exp2wei3 "  ,
                     " gam1wei2wei3 "  ,
                     " gam1gom2wei3 "  ,
                     " gam1gam2wei3 "  ,
                     " gam1logn2wei3 "  ,
                     " gam1logl2wei3 "  ,
                     " gam1exp2gom3 "  ,
                     " gam1wei2gom3 "  ,
                     " gam1gam2gom3 "  ,
                     " gam1logn2gom3 "  ,
                     " gam1logl2gom3 "  ,
                     " gam1exp2gam3 "  ,
                     " gam1wei2gam3 "  ,
                     " gam1gom2gam3 "  ,
                     " gam1gam2gam3 "  ,
                     " gam1logn2gam3 "  ,
                     " gam1logl2gam3 "  ,
                     " gam1exp2logn3 "  ,
                     " gam1wei2logn3 "  ,
                     " gam1gom2logn3 "  ,
                     " gam1gam2logn3 "  ,
                     " gam1logn2logn3 "  ,
                     " gam1logl2logn3 "  ,
                     " gam1exp2logl3 "  ,
                     " gam1wei2logl3 "  ,
                     " gam1gom2logl3 "  ,
                     " gam1gam2logl3 "  ,
                     " gam1logn2logl3 "  ,
                     " gam1logl2logl3 ")


gam1objectssimFC<-vector("list", 36)
gam1objectssimFC<-list(gam1gom2gom3simFC ,
                               gam1exp2exp3simFC ,
                               gam1wei2exp3simFC ,
                               gam1gom2exp3simFC ,
                               gam1gam2exp3simFC ,
                               gam1logn2exp3simFC ,
                               gam1logl2exp3simFC ,
                               gam1exp2wei3simFC ,
                               gam1wei2wei3simFC ,
                               gam1gom2wei3simFC ,
                               gam1gam2wei3simFC ,
                               gam1logn2wei3simFC ,
                               gam1logl2wei3simFC ,
                               gam1exp2gom3simFC ,
                               gam1wei2gom3simFC ,
                               gam1gam2gom3simFC ,
                               gam1logn2gom3simFC ,
                               gam1logl2gom3simFC ,
                               gam1exp2gam3simFC ,
                               gam1wei2gam3simFC ,
                               gam1gom2gam3simFC ,
                               gam1gam2gam3simFC ,
                               gam1logn2gam3simFC ,
                               gam1logl2gam3simFC ,
                               gam1exp2logn3simFC ,
                               gam1wei2logn3simFC ,
                               gam1gom2logn3simFC ,
                               gam1gam2logn3simFC ,
                               gam1logn2logn3simFC ,
                               gam1logl2logn3simFC ,
                               gam1exp2logl3simFC ,
                               gam1wei2logl3simFC ,
                               gam1gom2logl3simFC ,
                               gam1gam2logl3simFC ,
                               gam1logn2logl3simFC ,
                               gam1logl2logl3simFC)

gam1objectssimRFC<-vector("list", 36)
gam1objectssimRFC<-list(gam1gom2gom3simRFC ,
                                gam1exp2exp3simRFC ,
                                gam1wei2exp3simRFC ,
                                gam1gom2exp3simRFC ,
                                gam1gam2exp3simRFC ,
                                gam1logn2exp3simRFC ,
                                gam1logl2exp3simRFC ,
                                gam1exp2wei3simRFC ,
                                gam1wei2wei3simRFC ,
                                gam1gom2wei3simRFC ,
                                gam1gam2wei3simRFC ,
                                gam1logn2wei3simRFC ,
                                gam1logl2wei3simRFC ,
                                gam1exp2gom3simRFC ,
                                gam1wei2gom3simRFC ,
                                gam1gam2gom3simRFC ,
                                gam1logn2gom3simRFC ,
                                gam1logl2gom3simRFC ,
                                gam1exp2gam3simRFC ,
                                gam1wei2gam3simRFC ,
                                gam1gom2gam3simRFC ,
                                gam1gam2gam3simRFC ,
                                gam1logn2gam3simRFC ,
                                gam1logl2gam3simRFC ,
                                gam1exp2logn3simRFC ,
                                gam1wei2logn3simRFC ,
                                gam1gom2logn3simRFC ,
                                gam1gam2logn3simRFC ,
                                gam1logn2logn3simRFC ,
                                gam1logl2logn3simRFC ,
                                gam1exp2logl3simRFC ,
                                gam1wei2logl3simRFC ,
                                gam1gom2logl3simRFC ,
                                gam1gam2logl3simRFC ,
                                gam1logn2logl3simRFC ,
                                gam1logl2logl3simRFC)

restentriesnames<-c( " exp1exp2exp3 "  ,
                     " wei1exp2exp3 "  ,
                     " logn1exp2exp3 "  ,
                     " logl1exp2exp3 "  ,
                     " exp1wei2exp3 "  ,
                     " wei1wei2exp3 "  ,
                     " logn1wei2exp3 "  ,
                     " logl1wei2exp3 "  ,
                     " exp1gom2exp3 "  ,
                     " wei1gom2exp3 "  ,
                     " logn1gom2exp3 "  ,
                     " logl1gom2exp3 "  ,
                     " exp1gam2exp3 "  ,
                     " wei1gam2exp3 "  ,
                     " logn1gam2exp3 "  ,
                     " logl1gam2exp3 "  ,
                     " exp1logn2exp3 "  ,
                     " wei1logn2exp3 "  ,
                     " logn1logn2exp3 "  ,
                     " logl1logn2exp3 "  ,
                     " exp1logl2exp3 "  ,
                     " wei1logl2exp3 "  ,
                     " logn1logl2exp3 "  ,
                     " logl1logl2exp3 "  ,
                     " exp1exp2wei3 "  ,
                     " wei1exp2wei3 "  ,
                     " logn1exp2wei3 "  ,
                     " logl1exp2wei3 "  ,
                     " exp1wei2wei3 "  ,
                     " wei1wei2wei3 "  ,
                     " logn1wei2wei3 "  ,
                     " logl1wei2wei3 "  ,
                     " exp1gom2wei3 "  ,
                     " wei1gom2wei3 "  ,
                     " logn1gom2wei3 "  ,
                     " logl1gom2wei3 "  ,
                     " exp1gam2wei3 "  ,
                     " wei1gam2wei3 "  ,
                     " logn1gam2wei3 "  ,
                     " logl1gam2wei3 "  ,
                     " exp1logn2wei3 "  ,
                     " wei1logn2wei3 " ,
                     " logn1logn2wei3 "  ,
                     " logl1logn2wei3 "  ,
                     " exp1logl2wei3 "  ,
                     " wei1logl2wei3 "  ,
                     " logn1logl2wei3 "  ,
                     " logl1logl2wei3 "  ,
                     " exp1exp2gom3 "  ,
                     " wei1exp2gom3 "  ,
                     " logn1exp2gom3 "  ,
                     " logl1exp2gom3 "  ,
                     " exp1wei2gom3 "  ,
                     " wei1wei2gom3 "  ,
                     " logn1wei2gom3 "  ,
                     " logl1wei2gom3 "  ,
                     " exp1gom2gom3 "  ,
                     " wei1gom2gom3 "  ,
                     " logn1gom2gom3 "  ,
                     " logl1gom2gom3 "  ,
                     " exp1gam2gom3 "  ,
                     " wei1gam2gom3 "  ,
                     " logn1gam2gom3 "  ,
                     " logl1gam2gom3 "  ,
                     " exp1logn2gom3 "  ,
                     " wei1logn2gom3 "  ,
                     " logn1logn2gom3 "  ,
                     " logl1logn2gom3 "  ,
                     " exp1logl2gom3 "  ,
                     " wei1logl2gom3 "  ,
                     " logn1logl2gom3 "  ,
                     " logl1logl2gom3 "  ,
                     " exp1exp2gam3 "  ,
                     " wei1exp2gam3 "  ,
                     " logn1exp2gam3 "  ,
                     " logl1exp2gam3 "  ,
                     " exp1wei2gam3 "  ,
                     " wei1wei2gam3 "  ,
                     " logn1wei2gam3 "  ,
                     " logl1wei2gam3 "  ,
                     " exp1gom2gam3 "  ,
                     " wei1gom2gam3 "  ,
                     " logn1gom2gam3 "  ,
                     " logl1gom2gam3 "  ,
                     " exp1gam2gam3 "  ,
                     " wei1gam2gam3 "  ,
                     " logn1gam2gam3 "  ,
                     " logl1gam2gam3 "  ,
                     " exp1logn2gam3 "  ,
                     " wei1logn2gam3 "  ,
                     " logn1logn2gam3 "  ,
                     " logl1logn2gam3 "  ,
                     " exp1logl2gam3 "  ,
                     " wei1logl2gam3 "  ,
                     " logn1logl2gam3 "  ,
                     " logl1logl2gam3 "  ,
                     " exp1exp2logn3 "  ,
                     " wei1exp2logn3 "  ,
                     " logn1exp2logn3 "  ,
                     " logl1exp2logn3 "  ,
                     " exp1wei2logn3 "  ,
                     " wei1wei2logn3 "  ,
                     " logn1wei2logn3 "  ,
                     " logl1wei2logn3 "  ,
                     " exp1gom2logn3 "  ,
                     " wei1gom2logn3 "  ,
                     " logn1gom2logn3 "  ,
                     " logl1gom2logn3 "  ,
                     " exp1gam2logn3 "  ,
                     " wei1gam2logn3 "  ,
                     " logn1gam2logn3 "  ,
                     " logl1gam2logn3 "  ,
                     " exp1logn2logn3 "  ,
                     " wei1logn2logn3 "  ,
                     " logn1logn2logn3 "  ,
                     " logl1logn2logn3 "  ,
                     " exp1logl2logn3 "  ,
                     " wei1logl2logn3 "  ,
                     " logn1logl2logn3 "  ,
                     " logl1logl2logn3 "  ,
                     " exp1exp2logl3 "  ,
                     " wei1exp2logl3 "  ,
                     " logn1exp2logl3 "  ,
                     " logl1exp2logl3 "  ,
                     " exp1wei2logl3 "  ,
                     " wei1wei2logl3 "  ,
                     " logn1wei2logl3 "  ,
                     " logl1wei2logl3 "  ,
                     " exp1gom2logl3 "  ,
                     " wei1gom2logl3 "  ,
                     " logn1gom2logl3 "  ,
                     " logl1gom2logl3 "  ,
                     " exp1gam2logl3 "  ,
                     " wei1gam2logl3 "  ,
                     " logn1gam2logl3 "  ,
                     " logl1gam2logl3 "  ,
                     " exp1logn2logl3 "  ,
                     " wei1logn2logl3 "  ,
                     " logn1logn2logl3 "  ,
                     " logl1logn2logl3 "  ,
                     " exp1logl2logl3 "  ,
                     " wei1logl2logl3 "  ,
                     " logn1logl2logl3 "  ,
                     " logl1logl2logl3 ")

restobjectssimFC<-vector("list", 144)
restobjectssimFC<-list(exp1exp2exp3simFC ,
                       wei1exp2exp3simFC ,
 
                       logn1exp2exp3simFC ,
                       logl1exp2exp3simFC ,
                       exp1wei2exp3simFC ,
                       wei1wei2exp3simFC ,
              
                       logn1wei2exp3simFC ,
                       logl1wei2exp3simFC ,
                       exp1gom2exp3simFC ,
                       wei1gom2exp3simFC ,
                
                       logn1gom2exp3simFC ,
                       logl1gom2exp3simFC ,
                       exp1gam2exp3simFC ,
                       wei1gam2exp3simFC ,
                 
                       logn1gam2exp3simFC ,
                       logl1gam2exp3simFC ,
                       exp1logn2exp3simFC ,
                       wei1logn2exp3simFC ,
                   
                       logn1logn2exp3simFC ,
                       logl1logn2exp3simFC ,
                       exp1logl2exp3simFC ,
                       wei1logl2exp3simFC ,
                    
                       logn1logl2exp3simFC ,
                       logl1logl2exp3simFC ,
                       exp1exp2wei3simFC ,
                       wei1exp2wei3simFC ,
                    
                       logn1exp2wei3simFC ,
                       logl1exp2wei3simFC ,
                       exp1wei2wei3simFC ,
                       wei1wei2wei3simFC ,
                    
                       logn1wei2wei3simFC ,
                       logl1wei2wei3simFC ,
                       exp1gom2wei3simFC ,
                       wei1gom2wei3simFC ,
                    
                       logn1gom2wei3simFC ,
                       logl1gom2wei3simFC ,
                       exp1gam2wei3simFC ,
                       wei1gam2wei3simFC ,
                     
                       logn1gam2wei3simFC ,
                       logl1gam2wei3simFC ,
                       exp1logn2wei3simFC ,
                       wei1logn2wei3simFC ,
                    
                       logn1logn2wei3simFC ,
                       logl1logn2wei3simFC ,
                       exp1logl2wei3simFC ,
                       wei1logl2wei3simFC ,
                    
                       logn1logl2wei3simFC ,
                       logl1logl2wei3simFC ,
                       exp1exp2gom3simFC ,
                       wei1exp2gom3simFC ,
                     
                       logn1exp2gom3simFC ,
                       logl1exp2gom3simFC ,
                       exp1wei2gom3simFC ,
                       wei1wei2gom3simFC ,
                    
                       logn1wei2gom3simFC ,
                       logl1wei2gom3simFC ,
                       exp1gom2gom3simFC ,
                       wei1gom2gom3simFC ,
                    
                       logn1gom2gom3simFC ,
                       logl1gom2gom3simFC ,
                       exp1gam2gom3simFC ,
                       wei1gam2gom3simFC ,
                      
                       logn1gam2gom3simFC ,
                       logl1gam2gom3simFC ,
                       exp1logn2gom3simFC ,
                       wei1logn2gom3simFC ,
                    
                       logn1logn2gom3simFC ,
                       logl1logn2gom3simFC ,
                       exp1logl2gom3simFC ,
                       wei1logl2gom3simFC ,
                     
                       logn1logl2gom3simFC ,
                       logl1logl2gom3simFC ,
                       exp1exp2gam3simFC ,
                       wei1exp2gam3simFC ,
                     
                       logn1exp2gam3simFC ,
                       logl1exp2gam3simFC ,
                       exp1wei2gam3simFC ,
                       wei1wei2gam3simFC ,
                     
                       logn1wei2gam3simFC ,
                       logl1wei2gam3simFC ,
                       exp1gom2gam3simFC ,
                       wei1gom2gam3simFC ,
                     
                       logn1gom2gam3simFC ,
                       logl1gom2gam3simFC ,
                       exp1gam2gam3simFC ,
                       wei1gam2gam3simFC ,
                     
                       logn1gam2gam3simFC ,
                       logl1gam2gam3simFC ,
                       exp1logn2gam3simFC ,
                       wei1logn2gam3simFC ,
                      
                       logn1logn2gam3simFC ,
                       logl1logn2gam3simFC ,
                       exp1logl2gam3simFC ,
                       wei1logl2gam3simFC ,
                    
                       logn1logl2gam3simFC ,
                       logl1logl2gam3simFC ,
                       exp1exp2logn3simFC ,
                       wei1exp2logn3simFC ,
                     
                       logn1exp2logn3simFC ,
                       logl1exp2logn3simFC ,
                       exp1wei2logn3simFC ,
                       wei1wei2logn3simFC ,
                     
                       logn1wei2logn3simFC ,
                       logl1wei2logn3simFC ,
                       exp1gom2logn3simFC ,
                       wei1gom2logn3simFC ,
                      
                       logn1gom2logn3simFC ,
                       logl1gom2logn3simFC ,
                       exp1gam2logn3simFC ,
                       wei1gam2logn3simFC ,
                      
                       logn1gam2logn3simFC ,
                       logl1gam2logn3simFC ,
                       exp1logn2logn3simFC ,
                       wei1logn2logn3simFC ,
                      
                       logn1logn2logn3simFC ,
                       logl1logn2logn3simFC ,
                       exp1logl2logn3simFC ,
                       wei1logl2logn3simFC ,
                     
                       logn1logl2logn3simFC ,
                       logl1logl2logn3simFC ,
                       exp1exp2logl3simFC ,
                       wei1exp2logl3simFC ,
                     
                       logn1exp2logl3simFC ,
                       logl1exp2logl3simFC ,
                       exp1wei2logl3simFC ,
                       wei1wei2logl3simFC ,
                   
                       logn1wei2logl3simFC ,
                       logl1wei2logl3simFC ,
                       exp1gom2logl3simFC ,
                       wei1gom2logl3simFC ,
                     
                       logn1gom2logl3simFC ,
                       logl1gom2logl3simFC ,
                       exp1gam2logl3simFC ,
                       wei1gam2logl3simFC ,
                     
                       logn1gam2logl3simFC ,
                       logl1gam2logl3simFC ,
                       exp1logn2logl3simFC ,
                       wei1logn2logl3simFC ,
                      
                       logn1logn2logl3simFC ,
                       logl1logn2logl3simFC ,
                       exp1logl2logl3simFC ,
                       wei1logl2logl3simFC ,
                     
                       logn1logl2logl3simFC ,
                       logl1logl2logl3simFC)


restobjectssimRFC<-vector("list", 144)
restobjectssimRFC<-list(exp1exp2exp3simRFC ,
                        wei1exp2exp3simRFC ,
                    
                        logn1exp2exp3simRFC ,
                        logl1exp2exp3simRFC ,
                        exp1wei2exp3simRFC ,
                        wei1wei2exp3simRFC ,
                     
                        logn1wei2exp3simRFC ,
                        logl1wei2exp3simRFC ,
                        exp1gom2exp3simRFC ,
                        wei1gom2exp3simRFC ,
                       
                        logn1gom2exp3simRFC ,
                        logl1gom2exp3simRFC ,
                        exp1gam2exp3simRFC ,
                        wei1gam2exp3simRFC ,
                      
                        logn1gam2exp3simRFC ,
                        logl1gam2exp3simRFC ,
                        exp1logn2exp3simRFC ,
                        wei1logn2exp3simRFC ,
                       
                        logn1logn2exp3simRFC ,
                        logl1logn2exp3simRFC ,
                        exp1logl2exp3simRFC ,
                        wei1logl2exp3simRFC ,
                       
                        logn1logl2exp3simRFC ,
                        logl1logl2exp3simRFC ,
                        exp1exp2wei3simRFC ,
                        wei1exp2wei3simRFC ,
                       
                        logn1exp2wei3simRFC ,
                        logl1exp2wei3simRFC ,
                        exp1wei2wei3simRFC ,
                        wei1wei2wei3simRFC ,
                       
                        logn1wei2wei3simRFC ,
                        logl1wei2wei3simRFC ,
                        exp1gom2wei3simRFC ,
                        wei1gom2wei3simRFC ,
                       
                        logn1gom2wei3simRFC ,
                        logl1gom2wei3simRFC ,
                        exp1gam2wei3simRFC ,
                        wei1gam2wei3simRFC ,
                      
                        logn1gam2wei3simRFC ,
                        logl1gam2wei3simRFC ,
                        exp1logn2wei3simRFC ,
                        wei1logn2wei3simRFC ,
                      
                        logn1logn2wei3simRFC ,
                        logl1logn2wei3simRFC ,
                        exp1logl2wei3simRFC ,
                        wei1logl2wei3simRFC ,
                      
                        logn1logl2wei3simRFC ,
                        logl1logl2wei3simRFC ,
                        exp1exp2gom3simRFC ,
                        wei1exp2gom3simRFC ,
                     
                        logn1exp2gom3simRFC ,
                        logl1exp2gom3simRFC ,
                        exp1wei2gom3simRFC ,
                        wei1wei2gom3simRFC ,
                      
                        logn1wei2gom3simRFC ,
                        logl1wei2gom3simRFC ,
                        exp1gom2gom3simRFC ,
                        wei1gom2gom3simRFC ,
                      
                        logn1gom2gom3simRFC ,
                        logl1gom2gom3simRFC ,
                        exp1gam2gom3simRFC ,
                        wei1gam2gom3simRFC ,
                       
                        logn1gam2gom3simRFC ,
                        logl1gam2gom3simRFC ,
                        exp1logn2gom3simRFC ,
                        wei1logn2gom3simRFC ,
                       
                        logn1logn2gom3simRFC ,
                        logl1logn2gom3simRFC ,
                        exp1logl2gom3simRFC ,
                        wei1logl2gom3simRFC ,
                       
                        logn1logl2gom3simRFC ,
                        logl1logl2gom3simRFC ,
                        exp1exp2gam3simRFC ,
                        wei1exp2gam3simRFC ,
                       
                        logn1exp2gam3simRFC ,
                        logl1exp2gam3simRFC ,
                        exp1wei2gam3simRFC ,
                        wei1wei2gam3simRFC ,
                       
                        logn1wei2gam3simRFC ,
                        logl1wei2gam3simRFC ,
                        exp1gom2gam3simRFC ,
                        wei1gom2gam3simRFC ,
                      
                        logn1gom2gam3simRFC ,
                        logl1gom2gam3simRFC ,
                        exp1gam2gam3simRFC ,
                        wei1gam2gam3simRFC ,
                      
                        logn1gam2gam3simRFC ,
                        logl1gam2gam3simRFC ,
                        exp1logn2gam3simRFC ,
                        wei1logn2gam3simRFC ,
                        
                        logn1logn2gam3simRFC ,
                        logl1logn2gam3simRFC ,
                        exp1logl2gam3simRFC ,
                        wei1logl2gam3simRFC ,
                       
                        logn1logl2gam3simRFC ,
                        logl1logl2gam3simRFC ,
                        exp1exp2logn3simRFC ,
                        wei1exp2logn3simRFC ,
                       
                        logn1exp2logn3simRFC ,
                        logl1exp2logn3simRFC ,
                        exp1wei2logn3simRFC ,
                        wei1wei2logn3simRFC ,
                       
                        logn1wei2logn3simRFC ,
                        logl1wei2logn3simRFC ,
                        exp1gom2logn3simRFC ,
                        wei1gom2logn3simRFC ,
                       
                        logn1gom2logn3simRFC ,
                        logl1gom2logn3simRFC ,
                        exp1gam2logn3simRFC ,
                        wei1gam2logn3simRFC ,
                       
                        logn1gam2logn3simRFC ,
                        logl1gam2logn3simRFC ,
                        exp1logn2logn3simRFC ,
                        wei1logn2logn3simRFC ,
                    
                        logn1logn2logn3simRFC ,
                        logl1logn2logn3simRFC ,
                        exp1logl2logn3simRFC ,
                        wei1logl2logn3simRFC ,
                     
                        logn1logl2logn3simRFC ,
                        logl1logl2logn3simRFC ,
                        exp1exp2logl3simRFC ,
                        wei1exp2logl3simRFC ,
                       
                        logn1exp2logl3simRFC ,
                        logl1exp2logl3simRFC ,
                        exp1wei2logl3simRFC ,
                        wei1wei2logl3simRFC ,
                       
                        logn1wei2logl3simRFC ,
                        logl1wei2logl3simRFC ,
                        exp1gom2logl3simRFC ,
                        wei1gom2logl3simRFC ,
                      
                        logn1gom2logl3simRFC ,
                        logl1gom2logl3simRFC ,
                        exp1gam2logl3simRFC ,
                        wei1gam2logl3simRFC ,
                       
                        logn1gam2logl3simRFC ,
                        logl1gam2logl3simRFC ,
                        exp1logn2logl3simRFC ,
                        wei1logn2logl3simRFC ,
                      
                        logn1logn2logl3simRFC ,
                        logl1logn2logl3simRFC ,
                        exp1logl2logl3simRFC ,
                        wei1logl2logl3simRFC ,
                       
                        logn1logl2logl3simRFC ,
                        logl1logl2logl3simRFC)
                        

##----------------------------------------------------------
####FUNCTION FOR ICERs 
##-----------------------------------------------------------

ICERs_allcombs<-function(objectssimFC=gom1objectssimFC, objectssimRFC=gom1objectssimRFC, rowentriesnames=gom1entriesnames,  ncom=36, ntt=863)
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
  
  incQALYs[[i]]<-0.8*(meanLY_PFSRFC[[i]]- meanLY_PFSFC[[i]])+0.6*(meanLY_progRFC[[i]]- meanLY_progFC[[i]])
  
  ######## costs
  
  mean_cost_FC[[i]]<-2790+1115+22+1115+28*12*meanLY_PFSFC[[i]]+360+507+84*12*meanLY_progFC[[i]]+(5179/((meanLY_progRFC[[i]] + meanLY_progFC[[i]])/2*12))*12*meanLY_progFC[[i]]
  mean_cost_RFC[[i]]<-10113+1224+2776+1109+21+1109+28*12*meanLY_PFSRFC[[i]]+592+640+84*12*meanLY_progRFC[[i]]+(5179/((meanLY_progRFC[[i]] + meanLY_progFC[[i]])/2*12))*12*meanLY_progRFC[[i]]  
}
#incQALYs2<-unlist(incQALYs)
ICERstable<-cbind(unlist(incQALYs),unlist(mean_cost_RFC) -unlist(mean_cost_FC),
            (unlist(mean_cost_RFC) -unlist(mean_cost_FC))/unlist(incQALYs)) 
rownames(ICERstable)<-rowentriesnames
colnames(ICERstable)<-c("incQALYs", "incCosts", "ICER")
return(ICERstable)
}

#################################################
nrow(gom1objectssimFC[[1]])
nrow(gam1objectssimFC[[1]])
nrow(restobjectssimFC[[1]])

ICERs_gom1<-ICERs_allcombs(objectssimFC=gom1objectssimFC, objectssimRFC=gom1objectssimRFC, rowentriesnames=gom1entriesnames,  ncom=36, ntt=2429)
ICERs_gam1<-ICERs_allcombs(objectssimFC=gam1objectssimFC, objectssimRFC=gam1objectssimRFC, rowentriesnames=gam1entriesnames,  ncom=36, ntt=577)
ICERs_rest<-ICERs_allcombs(objectssimFC=restobjectssimFC, objectssimRFC=restobjectssimRFC, rowentriesnames=restentriesnames,  ncom=144, ntt=181)

ICERs_all<-rbind(ICERs_gom1,ICERs_gam1, ICERs_rest )

ICERs_alls<- ICERs_all[order(ICERs_all[,3]),] 
ICERs_alls



  

