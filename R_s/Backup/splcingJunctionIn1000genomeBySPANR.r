setwd(splicingJunction1000Genome);

sortAndUniquePSI("PSIpsiResultLinks1000genome.phase3.splicingjunction2.SNPs_only.species.AF005.vcf.forPSI.txt");

#data<-read.csv("High dPSI unique.csv");

getNeutralExonCentralPos("High dPSI unique.csv",20);




#psi<-data[,"dPSI"];
#hist(clinsig5PredictByPsi,main="Clinvar disease predict by SPANR",xlab="|dPSI|");
#hist(clinsig2PredictByPsi,main="Clinvar disease predict by SPANR",xlab="|dPSI|");
#hist(abs(psi),main="1000genome splicingJunction SNPs(MAF>5%) predict by SPANR",xlab="|dPSI|");

