setwd("E:\\limeng\\RBP_network\\r");

########combine Miso and FPKM value#############################################
source("combine_psi_expression.r",echo=FALSE);
################################################################################


#########combine Miso,FRPM and RBP. Linear regression between RBP and PSI#######
source("rbp_pos(exon_intron)_cor&psi_lm.r",echo=FALSE);
################################################################################


##############Gene different expression between two groups######################
source("RBP_gene_expression_heatmap.r",echo=FALSE);
################################################################################


########logistic regression between RBP binding and PSI#########################
source("RBP_binding_trend_logit_AUC_test.r",echo=FALSE);
################################################################################


###########RBP binding specificity and RBPs' function###########################
source("RBP_binding_trend_pos.r",echo=FALSE)
################################################################################


########Validation of RBP interaction using bioGrid. Cal each RBP's digree######
source("Biogrid_cor_comp.r",echo=FALSE)
################################################################################


####################RBPs' binding trend and motif###############################
source("pwm_seqLogo.r",echo=FALSE)
################################################################################


###############miso's PSIs quality control######################################
source("miso_psi_validation.r",echo=FALSE)
################################################################################


####################Quality control of RNA-seq and Clip-seq#####################
source("CLIP_MISO_FPKM_qc.r",echo=FALSE)
################################################################################





