#################################  This function develops the Risk score (step 1) - Cross Internal Method #########################################################
#############################################################################################################################################################


CrossInternalRisk.fun=function(dataset, outcome){
#load libraries
library(glmnet)
library(Hmisc)
if (outcome=="RELAPSE1year") {
####################random half RCTs from studies###############################
set.seed(1)
Advance<-dataset[which(dataset$STUDYID=="ADVANCE"),]
Advance.risk<-Advance[sample(nrow(Advance), nrow(Advance)/2),]
todrop<-c("STUDYID","USUBJID","RELAPSE2year")
Advance.risk<-Advance.risk[ , !(names(Advance.risk) %in% todrop)]
set.seed(1)
Define<-dataset[which(dataset$STUDYID=="DEFINE"),]
Define.risk<-Define[sample(nrow(Define), nrow(Define)/2),]
Define.risk<-Define.risk[ , !(names(Define.risk) %in% todrop)]
set.seed(1)
Confirm<-dataset[which(dataset$STUDYID=="CONFIRM"),]
Confirm.risk<-Confirm[sample(nrow(Confirm), nrow(Confirm)/2),]
Confirm.risk<-Confirm.risk[ , !(names(Confirm.risk) %in% todrop)]
set.seed(1)
Affirm<-dataset[which(dataset$STUDYID=="AFFIRM"),]
Affirm.risk<-Affirm[sample(nrow(Affirm), nrow(Affirm)/2),]
Affirm.risk<-Affirm.risk[ , !(names(Affirm.risk) %in% todrop)]
set.seed(1)
Mscrg<-dataset[which(dataset$STUDYID=="MSCRG"),]
Mscrg.risk<-Mscrg[sample(nrow(Mscrg), nrow(Mscrg)/2),]
Mscrg.risk<-Mscrg.risk[ , !(names(Mscrg.risk) %in% todrop)]
##all half studies together
mrg<-rbind(Advance.risk,Define.risk,Confirm.risk,Affirm.risk,Mscrg.risk)
#####################LASSO preparation####################

###blinded to treatment so drop variable TRT01A
todrop<-c("TRT01A")
mrg.both<-mrg[ , !(names(mrg) %in% todrop)]
### delete NA values (LASSO requierement)
mrg.both<-na.omit(mrg.both)

##logistic regression based on selected variables from LASSO
model<-lrm(RELAPSE1year~AGE+REGION+RLPS3YR+GDLESBL+SFPCSBL, data=mrg.both)

##check the assumption of linearity
a<-lrm(RELAPSE1year~rcs(AGE,5)+REGION+rcs(RLPS3YR,5)+rcs(GDLESBL,5)+rcs(SFPCSBL,5), data=mrg.both)
anova(a) #all linear
b<-lrm(RELAPSE1year~rcs(AGE,4)+REGION+rcs(RLPS3YR,4)+rcs(GDLESBL,4)+rcs(SFPCSBL,4), data=mrg.both)
anova(b) #all linear
AIC(a)
AIC(b) ###best choice of knots is knots=4

###final model (as anova does not indicate nonlinear relationship)
finalmodel<-lrm(RELAPSE1year~AGE+REGION+RLPS3YR+GDLESBL+SFPCSBL, data=mrg.both)
finalmodel
df<-finalmodel[["stats"]][["d.f."]]
events<- nrow(mrg.both[which(mrg.both$RELAPSE1year==1),])
####return back
cat("The selected variables are:AGE, REGION, RLPS3YR, GDLESBL, SFPCSBL ", fill = TRUE)
EPV<-events/df
cat("The EPV of the model is", EPV, fill=TRUE)
if (EPV>10){cat("The EPV of the model is >10, as recommended.", fill=TRUE)}
if (EPV<10){cat("The EPV of the model is <10, possibility of overfitting problem.", fill = TRUE)}

cat("The final model is:",fill=TRUE)
return(finalmodel)
}

  if (outcome=="RELAPSE2year") { 
    ####################random half RCTs from studies###############################
    todrop<-c("STUDYID","USUBJID","RELAPSE1year")
    set.seed(1)
    Define<-dataset[which(dataset$STUDYID=="DEFINE"),]
    Define.risk<-Define[sample(nrow(Define), nrow(Define)/2),]
    Define.risk<-Define.risk[ , !(names(Define.risk) %in% todrop)]
    set.seed(1)
    Confirm<-dataset[which(dataset$STUDYID=="CONFIRM"),]
    Confirm.risk<-Confirm[sample(nrow(Confirm), nrow(Confirm)/2),]
    Confirm.risk<-Confirm.risk[ , !(names(Confirm.risk) %in% todrop)]
    set.seed(1)
    Affirm<-dataset[which(dataset$STUDYID=="AFFIRM"),]
    Affirm.risk<-Affirm[sample(nrow(Affirm), nrow(Affirm)/2),]
    Affirm.risk<-Affirm.risk[ , !(names(Affirm.risk) %in% todrop)]
    set.seed(1)
    Mscrg<-dataset[which(dataset$STUDYID=="MSCRG"),]
    Mscrg.risk<-Mscrg[sample(nrow(Mscrg), nrow(Mscrg)/2),]
    Mscrg.risk<-Mscrg.risk[ , !(names(Mscrg.risk) %in% todrop)]
    ##all half studies together
    mrg<-rbind(Define.risk,Confirm.risk,Affirm.risk,Mscrg.risk)
    #####################LASSO preparation####################
    
    ###blinded to treatment so drop variable TRT01A
    todrop<-c("TRT01A")
    mrg.both<-mrg[ , !(names(mrg) %in% todrop)]
    ### delete NA values (LASSO requierement)
    mrg.both<-na.omit(mrg.both)
    #### model matrix needed for LASSO
    half.matrix<-model.matrix(mrg.both$RELAPSE2year~.,data=mrg.both)
    half.matrix<-na.omit(half.matrix)
    
    ##logistic regression based on selected variables from LASSO
    model<-lrm(RELAPSE2year~REGION+RLPS3YR+GDLESBL+DISTWKBL+EDSSBL+SFPCSBL, data=mrg.both)
    
    ##check the assumption of linearity
    a<-lrm(RELAPSE2year~REGION+rcs(RLPS3YR,5)+rcs(GDLESBL,5)+DISTWKBL+rcs(EDSSBL,5)+rcs(SFPCSBL,5), data=mrg.both)
    anova(a) ##all linear
    b<-lrm(RELAPSE2year~REGION+rcs(RLPS3YR,4)+rcs(GDLESBL,4)+DISTWKBL+rcs(EDSSBL,4)+rcs(SFPCSBL,4), data=mrg.both)
    anova(b)##rlps3y non linear
    AIC(a)
    AIC(b) ###best choice of knots is knots=4
    
    ###final model (as anova does not indicate nonlinear relationship)
    finalmodel<-lrm(RELAPSE2year~REGION+rcs(RLPS3YR,4)+GDLESBL+DISTWKBL+EDSSBL+SFPCSBL, data=mrg.both)
    finalmodel
    df<-finalmodel[["stats"]][["d.f."]]
    events<- nrow(mrg.both[which(mrg.both$RELAPSE2year==1),])
    ####return back
    cat("The selected variables are: REGION, RLPS3Y, GDLESBL, DISTWKBL, EDSSBL, SFPCSBL ", fill = TRUE)
    EPV<-events/df
    cat("The EPV of the model is", EPV, fill=TRUE)
    if (EPV>10){cat("The EPV of the model is >10, as recommended.", fill=TRUE)}
    if (EPV<10){cat("The EPV of the model is <10, possibility of overfitting problem.", fill = TRUE)}
    
    cat("The final model is:",fill=TRUE)
    return(finalmodel)
    
  }
  

}