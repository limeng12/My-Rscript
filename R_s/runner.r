###################################global variable and class#############################################
RPath<-"/Users/mengli/Documents/projects/splicingSNP/R/";

setwd(RPath);
source("BED.r");
source("Global.r");

###################################process HGMD and PSI data#############################################
#setwd(RPath);
#source("getHGMDIntronSNP.r");
#setwd(RPath);
#source("getPSINeutralSNPRun.r");
#setwd(RPath);
#source("testOn_bigger0_PSI.r");


####################read HGMD, YqoQi Neutral and PSI data,plot figure2 and figure3#######################
setwd(RPath);
source("readDataHGMDNeutral.r");
setwd(RPath)
source("figure2.r");
setwd(RPath)
source("figure3.r");


###################################feature selection and plot figure4####################################
setwd(RPath);
#source("featureSelection.r");


###################################plot figure5 and figure6##############################################

setwd(RPath);
source("figure5.r");

setwd(RPath);
source("testOnClinvar.R");

setwd(RPath);
source("testPsiOnClinvar.R");

setwd(RPath);
source("figure6.R");

