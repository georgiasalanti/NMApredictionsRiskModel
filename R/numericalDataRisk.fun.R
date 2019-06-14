###########################################################################################################################
##############  function that selects variables of interest and ####################
##############         recodes some of factors        #############################
#####################################################################################################################################################

numericalDataRisk.fun=function(dataset){
  library(car)
  ####recode treatment names####
  dataset$TRT01A<-as.numeric(dataset$TRT01A)
  dataset$TRT01A<-recode(dataset$TRT01A, "1='Avonex'; 2='Dimethyl fumarate';3='Peginterferon Beta-1a'; 4='Peginterferon Beta-1a';5='Glatiramer acetate';6='Natalizumab';7='Natalizumab + Avonex';8='Placebo';9='Avonex'")
  dataset$TRT01A<-as.factor(dataset$TRT01A)
  ##################################HANDLING VARIABLES##############################################
  ##keep only needed variables
  keep<-c("STUDYID","USUBJID","AGE","SEX","RACE","TRT01A","HEIGHTBL","WEIGHTBL","EDSSBL","ONSYRS","DIAGYRS",
          "DOMIHAND","RLPS3YR","TRELMOS", "MCDBL","PRMSGR","REGION","T25FWABL","NHPTMBL","NHPTDHBL","NHPTNHBL",
          "PASATABL","MSFCBL","GDLESBL","T2VOLBL","T1VOLBL","BVZBL","VFT100BL","VFT25BL",
          "VFT125BL","SFPCSBL","SFMCSBL","VISUALBL","BRAINBL","PYRAMIBL","SENSORBL","BOWLBLBL",
          "CEREBRBL","DISTWKBL","T25FWP1","NHPTMP1","NHPTDHP1","NHPTNHP1","PASATP1","T25FWPC","NHPTMPC","NHPTDHPC",
          "NHPTNHPC","PASATPC","RELAPSE1year","RELAPSE2year")
  MSrelapse<-dataset[,keep]
  
  ##################################HANDLING VARIABLES############################
  ###########################RECODE VARIABLES and make them factors############################
  MSrelapse$SEX<-recode(MSrelapse$SEX, "'M'=1; 'F'=0")
  MSrelapse$RACE<-recode(MSrelapse$RACE, "'WHITE'=1; 'NON-WHITE'=0")
  MSrelapse$DOMIHAND<-recode(MSrelapse$DOMIHAND, "'Left'=1;'LEFT'=1;  'Right'=0; 'RIGHT'=0;" )
  MSrelapse$DOMIHAND[which(MSrelapse$DOMIHAND=="")]<-NA
  ###MCDBL categories instead of 1, 2, 3, 4 there are  1, 2, >=3
  MSrelapse$MCDBL[which(MSrelapse$MCDBL==4)]<-3
  MSrelapse$MCDBL<-as.factor(MSrelapse$MCDBL)
  MSrelapse$PRMSGR<-as.factor(MSrelapse$PRMSGR)
  MSrelapse$REGION<-recode(MSrelapse$REGION, "'Eastern Europe'=1;'India'=2;  'North America'=3; 'ROW'=4;'Western Europe'=5 " )
  ###VISUALBL categories instead of 0, 1, 2, 3, 4, 5, 6 there are 0, 1, 2, >=3
  MSrelapse$VISUALBL[which(MSrelapse$VISUALBL==4 | MSrelapse$VISUALBL==5 | MSrelapse$VISUALBL==6)]<-3
  MSrelapse$VISUALBL<-as.factor(MSrelapse$VISUALBL)
  ###BRAINBL categories instead of 0, 1, 2, 3, 4 there are 0, 1, >=2
  MSrelapse$BRAINBL[which(MSrelapse$BRAINBL==3 | MSrelapse$BRAINBL==4)]<-2
  MSrelapse$BRAINBL<-as.factor(MSrelapse$BRAINBL) 
  ###PYRAMIBL categories instead of 0, 1, 2, 3, 4, 5, 6 there are 0, 1, 2, >=3
  MSrelapse$PYRAMIBL[which(MSrelapse$PYRAMIBL==4 | MSrelapse$PYRAMIBL==5)]<-3
  MSrelapse$PYRAMIBL<-as.factor(MSrelapse$PYRAMIBL)

  ###SENSOR BL categories instead of 0, 1, 2, 3, 4, 5, 6 there are 0, 1, 2, >=3
  MSrelapse$SENSORBL[which(MSrelapse$SENSORBL==4 | MSrelapse$SENSORBL==5 | MSrelapse$SENSORBL==6)]<-3
  MSrelapse$SENSORBL<-as.factor(MSrelapse$SENSORBL)
  ###BOWLBLBLL categories instead of 0, 1, 2, 3, 4, 5, 6 there are 0, 1, 2, >=3
   MSrelapse$BOWLBLBL[which(MSrelapse$BOWLBLBL==4 | MSrelapse$BOWLBLBL==5 | MSrelapse$BOWLBLBL==6)]<-3
  MSrelapse$BOWLBLBL<-as.factor(MSrelapse$BOWLBLBL)
  ###CEREBRBL categories instead of 0 ,1, 2, 3, 4 there are 0, 1, >=2
   MSrelapse$CEREBRBL[which(MSrelapse$CEREBRBL==3 | MSrelapse$CEREBRBL==4)]<-2
  MSrelapse$CEREBRBL<-as.factor(MSrelapse$CEREBRBL)
  
  mylevels <- c("< 100 METRES",
                ">= 100 METRES AND < 200 METRES", 
                ">= 200 METRES AND < 300 METRES", 
                ">= 300 METRES AND < 500 METRES", 
                ">= 500 METRES")
  
  MSrelapse$DISTWKBL<- as.numeric(factor(MSrelapse$DISTWKBL, levels=mylevels))
  MSrelapse$DISTWKBL<-as.factor(MSrelapse$DISTWKBL)
 
  #remove Sentinel study - Not included in AD data - Combination of therapies 
  MSrelapse<-MSrelapse[which(MSrelapse$STUDYID!="SENTINEL"),]
  return(MSrelapse)
}
