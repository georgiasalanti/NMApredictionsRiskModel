#################################  This function develops the Risk score (step 1) - Internal Method #########################################################
#############################################################################################################################################################

internalRisk.fun=function(dataset,outcome) {
  ##loading libraries
  library(glmnet)
  library(Hmisc)
  library(rms)
  
#########if  the outcome is 1 year MS relapse 
if(outcome=="RELAPSE1year") {  
  ################# LASSO preperation ##########################################
  ################################################################
  #drop STUDYID, USIBJID (no needed for the matrix), drop TRT01A (blinded to treatment)
  ##### drop RELAPSE2year - we need only the analysis foor 1 year
  
  todrop<-c("STUDYID","USUBJID","TRT01A","RELAPSE2year")
  X<-dataset[ , !(names(dataset) %in% todrop)]
  ### develop the matrix for the model
  X.both<-model.matrix(X$RELAPSE1year~.,data=X)
  #### delete NAs values (LASSO requirement)
  X.both<-na.omit(X.both)
  X<-na.omit(X)
  set.seed(1)
  #############################LASSO###########################################
  ########################################################################
  ###LASSO (alpha=1-default) with 10 cross validations
  cv.fit.both<-cv.glmnet(x=X.both,y=X$RELAPSE1year,family="binomial")
  ###coefficients of variables via LASSO
  cv.coef.both<-coef(cv.fit.both,s="lambda.1se")
  ###############RESULTS
  ### Non zero coefficients -> selected variables
  
  cv.pf.em.both<-rownames(cv.coef.both)[as.numeric(cv.coef.both)!=0]
  
  ###selected variables
  cv.pf.em.both
  ##plot lambda
  plot(cv.fit.both)
  ##plot coeff
  plot(coef(cv.fit.both))
  ###logistic model based on LASSO selected variables
  model<-lrm(RELAPSE1year~AGE+RLPS3YR+REGION+GDLESBL+SFPCSBL+SENSORBL,data=X)
  model
  ##check the assumption of linearity
  a<-lrm(RELAPSE1year~rcs(AGE,5)+rcs(RLPS3YR,5)+REGION+rcs(GDLESBL,5)+rcs(SFPCSBL,5)+SENSORBL,data=X)
  anova(a)
  b<-lrm(RELAPSE1year~rcs(AGE,4)+rcs(RLPS3YR,4)+REGION+rcs(GDLESBL,4)+rcs(SFPCSBL,4)+SENSORBL,data=X)
  anova(b)
  AIC(a)
  AIC(b) ###best choice of knots is knots=4
  
  ###final model (as anova does not indicate nonlinear relationship)
  model<-lrm(RELAPSE1year~AGE+RLPS3YR+REGION+GDLESBL+SFPCSBL+SENSORBL,data=X)
  model
  df<-model[["stats"]][["d.f."]]
  events<- nrow(X[which(X$RELAPSE1year==1),])
  ####return back
  cat("The selected variables are:" , cv.pf.em.both, fill = TRUE)
  EPV<-events/df
  cat("The EPV of the model is", EPV, fill=TRUE)
  if (EPV>10){cat("The EPV of the model is >10, as recommended.", fill=TRUE)}
  if (EPV<10){cat("The EPV of the model is <10, possibility of overfitting problem.", fill = TRUE)}
  
  cat("The final model is:",fill=TRUE)
  return(model)}

###### if the outcome is 2 years MS relapses
  if (outcome=="RELAPSE2year") {
    ################# LASSO preperation ##########################################
    ################################################################
    MSrelapse<-dataset[which(MSrelapse$STUDYID!="ADVANCE"),]
    dataset<-MSrelapse
    #drop STUDYID, USIBJID (no needed for the matrix), drop TRT01A (blinded to treatment)
    ##### drop RELAPSE1year - we need only the analysis for 2 years
    todrop<-c("STUDYID","USUBJID","TRT01A","RELAPSE1year")
    X<-dataset[ , !(names(dataset) %in% todrop)]
    ### develop the matrix for the model
    X.both<-model.matrix(X$RELAPSE2year~.,data=X)
    #### delete NAs values (LASSO requirement)
    X.both<-na.omit(X.both)
    X<-na.omit(X)
    set.seed(1)
    #############################LASSO###########################################
    ########################################################################
    ###LASSO (alpha=1-default) with 10 cross validations
    cv.fit.both<-cv.glmnet(x=X.both,y=X$RELAPSE2year,family="binomial")
    ###coefficients of variables via LASSO
    cv.coef.both<-coef(cv.fit.both,s="lambda.1se")
    ###############RESULTS
    ### Non zero coefficients -> selected variables
    
    cv.pf.em.both<-rownames(cv.coef.both)[as.numeric(cv.coef.both)!=0]
    
    ###selected variables
    cv.pf.em.both
    ##plot lambda
    plot(cv.fit.both)
    ##plot coeff
    plot(coef(cv.fit.both))
    ###logistic model based on LASSO selected variables
    model<-lrm(RELAPSE2year~AGE+WEIGHTBL+EDSSBL+RLPS3YR+TRELMOS+PRMSGR+REGION+NHPTDHBL+GDLESBL+SFPCSBL+SENSORBL+DISTWKBL,data=X)
    model 
    ##check the assumption of linearity
    ##best choice of knots and test for linearity
    a<-lrm(RELAPSE2year~rcs(AGE,5)+rcs(WEIGHTBL,5)+rcs(EDSSBL,5)+rcs(RLPS3YR,5)+rcs(TRELMOS,5)+PRMSGR+REGION+rcs(NHPTDHBL,5)+rcs(GDLESBL,5)+rcs(SFPCSBL,5)+SENSORBL+DISTWKBL,data=X)
    anova(a)##relapse3years
    b<-lrm(RELAPSE2year~rcs(AGE,4)+rcs(WEIGHTBL,4)+rcs(EDSSBL,4)+rcs(RLPS3YR,4)+rcs(TRELMOS,4)+PRMSGR+REGION+rcs(NHPTDHBL,4)+rcs(GDLESBL,4)+rcs(SFPCSBL,4)+SENSORBL+DISTWKBL,data=X)
    anova(b)##relapse3years non linear
    AIC(a)
    AIC(b) ###best choice of knots is knots=4
    
    finalmodel<-lrm(RELAPSE2year~AGE+WEIGHTBL+EDSSBL+rcs(RLPS3YR,4)+TRELMOS+PRMSGR+REGION+NHPTDHBL+GDLESBL+SFPCSBL+SENSORBL+DISTWKBL,data=X)
    
    df<-finalmodel[["stats"]][["d.f."]]
    events<- nrow(X[which(X$RELAPSE2year==1),])
    ####return back
    cat("The selected variables are:" , cv.pf.em.both, fill = TRUE)
    EPV<-events/df
    cat("The EPV of the model is", EPV, fill=TRUE)
    if (EPV>10){cat("The EPV of the model is >10, as recommended.", fill=TRUE)}
    if (EPV<10){cat("The EPV of the model is <10, possibility of overfitting problem.", fill = TRUE)}
    
    cat("The final model is:",fill=TRUE)
    return(finalmodel)
  }

  }