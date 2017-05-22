SELECT distinct gene_product.symbol	//从gene_product数据表中提取symbol字段，且去重复distinct(删除重复的行)
FROM					
gene_product				//从gene_product表中提取的每一行数据data1，与以下7个表中的数据满足相应的条件
INNER JOIN species ON (gene_product.species_id = species.id)	//与表species 的每一行数据的字段id与data1中id字段相等
INNER JOIN dbxref ON (gene_product.dbxref_id = dbxref.id)	//与表dbxref的每一行数据的字段id与data1的id字段相等
INNER JOIN association ON (gene_product.id = association.gene_product_id)	//与表association的gene_product_id与data1中的id相等
INNER JOIN evidence ON (association.id=evidence.association_id)			//与表evidence的association_id与data1中的id相等
INNER JOIN graph_path ON (association.term_id = graph_path.term2_id)		//表graph_path 的term2_id与表association 中的term_id相等相等
INNER JOIN term ON (graph_path.term1_id = term.id)				//表term的id字段与表graph_path.term1_id 相等
INNER JOIN term AS associated_term ON (graph_path.term2_id = associated_term.id)//表term的id与graph_path.term2_id相等
WHERE	//将满足上述的所有数据进行条件过滤
term.acc ='GO:0004713' AND		//表term.acc =='GO:0004713'
species.common_name ='human' AND	//且species.common_name =='human'
dbxref.xref_dbname='UniProtKB' AND	//且dbxref.xref_dbname='UniProtKB'来自UniProtKB数据库的项目
gene_product.symbol in			//且gene_product.symbol字段包含在以下括号中的字符串中
('SCPEP1','MANBA','MRPL22','RB1','CTHRC1','DNAJC6','CAST','TRERF1','SBK3','MEF2C','FES','CTSC','HTR7','VPRBP','KCNE3','NRP2','MFSD11','ZNF462','TDRP','
IK','SERPINE2','GPR68','NCF1','PER3','L3MBTL3','ZNF28','LOC399815','SMEK2','ABCA3','LOC727751','PRKAA2','ABCA10','LOC100288974','
PIK3CG','KLHDC7A','GPR161','PLEKHO1','APMAP','GBP5','MMP1','HBEGF','ARPC4','KIAA1211','C1RL','SAMD13','ZNF100','SLC45A3','CLK4','
MIR6732','AIFM3','TMEM206','KIAA0368','PCTP','MIR936','ADAM28','BTK','ZNF280B','LHX3','IGFALS','PHACTR1','CERK','CDH24','EXOC2','
ATP1A1','C15orf37','CCDC142','STEAP2','B4GALNT1','PDE7A','FAM117A','CXCL3','PDZD2','DCAF7','COMP','RIMS4','CHRNA4','KCNH2','LCN12','
C6orf141','IL36G','ARPC4-TTLL3','GPR125','UIMC1','BCAP31','ZNF561-AS1','PPAP2A','GADD45B','RUNX2','ZNF136','MAP2K2','ANAPC16','ATG7','SIRT6','FRMD5','BRWD1','WDR11','PTK7','PAX5','CTNNA1','ABHD2','
ZNF529','TRAPPC11','RAB5C','PCM1','TEKT2','PTPRZ1','KBTBD12','ADARB2','IGSF9B','GLB1L2','FBN3','OLFM1','DSC1','SLC6A14','CD1D','
CELF2','SATB1','SLC35D2','FAM35DP','DERL1','ELK4','CRK','PSEN1','GGA3','ENTPD6','TM9SF4','MIR4761','SPNS2','PTPRF','METTL10','MAP3K11','
UBC','CLSTN1','GALNT9','TINAGL1','ANKFY1','GANAB','FAH','SRC','RABGEF1','AKAP10','ENGASE','RINT1','ELF2','RDX','SP140','SLC16A14','PRDM16','
PDE3A','EDA2R','TMSB15B','HPS5','LMCD1','KIAA1467','CTDSP2','ULK4','RNMTL1','FLJ44313','CDC25B','MKLN1','NOTCH2NL','ZNF816','FARS2','RANBP9','
ARMCX4','MFSD2A','FOXA2','PWWP2A','RALGPS1','RUSC1','IKZF4','RANBP3','AGPAT9','PHYKPL','ULBP2','LOC100129973','CNOT8','MEN1','CEP44','ATPAF1','
PRPF4B','PSIP1','AP2A2','SLC22A20','HIST1H2BN','ZDHHC5','SLC41A1','C11orf57','TSACC','EIF4EBP3','DDX11L2','TNFSF10','TMTC2');
