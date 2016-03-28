setwd("/home/limeng/splicingSNP/R");
source("BED.r");
source("getPSINeutralSNP.r");
library('BSgenome.Hsapiens.UCSC.hg19');

#setwd("/home/limeng/splicingSNP/1000genome/SplicingJunction20/PSIResultAF10/AF10");
#sortAndUniquePSI("allPsiResult.txt");
#setwd("/home/limeng/splicingSNP/1000genome/WholeExonIntron20/PSIResult/AF10/");
#sortAndUniquePSI("allPsiResultWholeExon20.txt");
#setwd("/home/limeng/splicingSNP/1000genome/WholeExonIntron20/PSIResult/AF5/");
#sortAndUniquePSI("allPsiResultWholeExon20AF5.txt");
#setwd("/home/limeng/splicingSNP/1000genome/SplicingJunction20/PSIResultAF10/AF10");
#getExonPositionForPSIFile("High dPSI unique.csv","/home/limeng/splicingSNP/annotation/refgene.bed");
#setwd("/home/limeng/splicingSNP/1000genome/WholeExonIntron20/PSIResult/AF5/");
#getExonPositionForPSIFile("High dPSI unique.csv","/home/limeng/splicingSNP/annotation/refgene.bed",5);


setwd(psiFilesDir);
psiFileName<-sortAndUniquePSI("allPsiResultWholeExon20AF5Species.txt");

#psiFileName="High dPSI unique.csv";

setwd(psiFilesDir);

getNeutralExonCentralPos(psiFileName,0);


#setwd("/home/limeng/splicingSNP/1000genome/WholeExonIntron20/PSIResult/AF10/")
#getNeutralExonCentralPos();

getYaoqiNeutral();

