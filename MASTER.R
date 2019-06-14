############################################################
#         Master analysis for MS NMA Prediction MODEL 
############################################################
###load the github libraries
library(devtools)
install_github("htx-r/CleaningData",force=TRUE)
install_github("htx-r/RiskModelNMApredictions", force = TRUE)
library(RiskModelNMApredictions)
library(CleaningData)

###### Give your path of data
mydatapath="C:/Users/kc19o338/Desktop/HTx/data/IPD data from 6 Biogen trials"

### load data
cleanBIOGENtrials<-cleanBIOGENtrials.fun(mydatapath)
adsl01<-cleanBIOGENtrials$adsl01
### Select variables that I need and recode them in numerical values (e.g. Male=1, Female=0)
MSrelapse<-numericalDataRisk.fun(adsl01)
#########results of Internal risk score
##for 1 year relapses
model1<-internalRisk.fun(MSrelapse,"RELAPSE1year" )
model1
##for 2 year relapses
model2<-internalRisk.fun(MSrelapse,"RELAPSE2year" )
model2
###########      results of Cross Internal Risk score          #################################
###############################################################################################
#selection of variables based on n half RCTs (script named CrossSelection.Lasso)        #######################
##### because of the "for" loop it takes a lot of time and is only needed for 
##### the selection of variables in cross Internal Risk score
##### its results are integrated into Cross internal risk score
#results for 1 year MS relapses
CrossInternalRisk.fun(MSrelapse, "RELAPSE1year")
#results for 2 year MS relapses
CrossInternalRisk.fun(MSrelapse, "RELAPSE2year")
