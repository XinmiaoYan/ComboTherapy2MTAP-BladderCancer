library(Seurat)
library(dplyr)
library(Matrix)
library(harmony)
library(patchwork)
library(ggplot2)
library(devtools)
library(CelliD)
library(tidyverse)
library(ggpubr)
library(clustree)
library(tidyverse)
library(readxl)
library(openxlsx)
library(colorspace)
library(RColorBrewer)
library(viridis)
library(scales)
library(ggmap)
library(DoubletFinder)
library(sctransform)
library(circlize)
library("RColorBrewer")
library(glmGamPoi)
library(openxlsx)
rm(list = ls())

setwd('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/02_result/')
set.seed(999)

dir = '15_Source_data'
dir.create(dir)
setwd(dir)

dir = '02_FigS1-2'
dir.create(dir)
setwd(dir)

dir = '01_CD4T'
dir.create(dir)
setwd(dir)

load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/01_data/01_obj_cd4_7847_v1.Rdata')
obj=obj_cd4
table(obj$level3_cell)

cell_rank = c('CD4_c0_T_ANXA1','CD4_c1_T_SERINC5','CD4_c2_Treg_IL32','CD4_c3_CTL_GZMA',
'CD4_c4_Treg_IKZF2','CD4_c5_Tfh_TOX')
obj$level3_cell = factor(obj$level3_cell, levels = cell_rank)

my_color = c('#F27794','#B3DE69','#B9E5EE','#EEED7F','#e1c9fc','#BD56A2','#FBCB1F') 
p1 = DimPlot(obj, label = F, reduction = 'umap', group.by = 'level3_cell',cols = my_color)
ggsave("01_cell.pdf",plot=p1, width = 5, height = 3)

annotation_gene = c('EEF1B2','ANXA1','NOSIP','CCR7','TPT1','FAU','NACA','EEF1G','UBA52','TOMM7','SNHG29','PFDN5',#c0
'CD69','MACF1','SERINC5','RIPOR2','ANK3','SORL1','NKTR','AHNAK','ATM','PAG1','PRKY', #c1
'FOXP3','TNFRSF4','IL32','BATF','TNFRSF18','TIGIT', #c2 treg
"GZMA", "GZMB",'GZMH',"GZMK",'PRF1','GNLY','NKG7','IFNG','CCL4','CCL5','CXCR3',# c3 ctl
'IKZF2','IL2RA','CTLA4','TNFRSF9','TNFRSF1B', 'CCR4','IL10RA', 'CCR8', 'IL2RB','CXCR6','ICOS', #CD4_c4_Treg_IKZF2
'CXCL13', 'PASK','TOX','TOX2', 'PDCD1','BTLA', 'ITM2A','LIMS1','NR3C1' #CD4_c5_Tfh
)

d1 = DotPlot(obj, features = unique(annotation_gene), cluster.idents = F) +
  RotatedAxis() + scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "BuPu")) + 
  theme(axis.title.x = element_blank(),strip.text.x = element_text(angle=90, size = 10)) +
  theme(strip.background =element_rect(fill="#F0F8FF"),legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=6),legend.text = element_text(size=6)) +
  FontSize(x.text = 9, main = 12)+labs(x=NULL, y=NULL)
ggsave("02_dp_clusters.pdf",plot=d1, width = 12, height = 2.2)
setwd('../')

dir = '02_CD8T'
dir.create(dir)
setwd(dir)

load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/01_data/00_obj_cd8_7168.Rdata')
obj=obj_cd8

cell_rank = c('CD8_c0_Tcm','CD8_c1_T_KLRG1','CD8_c2_Tex','CD8_c3_Teff',
'CD8_c4_Trm','CD8_c5_MAIT')
obj$level3_cell = factor(obj$level3_cell, levels = cell_rank)

Idents(obj) <- "level3_cell"
my_color = c('#B9E5EE','#4DB7D4','#EEED7F','#dedc00','#F781BF',"#E7298A") 
p1 = DimPlot(obj, label = F, reduction = 'umap', group.by = 'level3_cell',cols = my_color)
ggsave("01_cell.pdf",plot=p1, width = 5, height = 3)

annotation_gene = c('GZMK','ITM2C','CXCR4','GPR183','CCR7','SELL','CMC1', #c0_Tcm
                    'NEAT1','SYNE1','MACF1','SYNE2','ATM','KLRG1','SORL1','RIPOR2','MIAT','AOAH', #c1_T
                    'CXCL13','CCL3','PDCD1','LAG3','HAVCR2','CTLA4','TNFRSF9','TNFRSF18','TOX','TIGIT','IFNG',#'CD38','EOMES', #c2_Tex
                    "GZMA","GZMB",'GZMH',"GZMK",'GZMM','PRF1', 'GNLY','NKG7',"FCGR3A",'FGFBP2','CX3CR1','KLF2', #c3_Teff
                    'CD69', 'ITGAE','ZNF683', 'ITGA1', 'CXCR6','SRGAP3','ENTPD1',#c4_Trm
                    'TRAV1-2', 'SLC4A10', 'KLRB1','IL18R1','KLRB1','RORC' #c5_MAIT
)
d1 = DotPlot(obj, features = unique(annotation_gene), cluster.idents = F) +
  RotatedAxis() + scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "BuPu")) + 
  theme(axis.title.x = element_blank(),strip.text.x = element_text(angle=90, size = 10)) +
  theme(strip.background =element_rect(fill="#F0F8FF"),legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=6),legend.text = element_text(size=6)) +
  FontSize(x.text = 9, main = 12)+labs(x=NULL, y=NULL)
ggsave("02_dp_clusters.pdf",plot=d1, width = 11, height = 2.2)

setwd('../')


dir = '03_NK'
dir.create(dir)
setwd(dir)

load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/01_data/00_obj_nk_1228.Rdata')
obj = obj_nk
table(obj$level3_cell)

cell_rank = c('NK_c0_CD16_GNLY','NK_c1_CD56_XCL2','NK_c2_CD16_ITGAL')
obj$level3_cell = factor(obj$level3_cell, levels = cell_rank)

Idents(obj) <- "level3_cell"
my_color = c('#1BB6AFFF',"#FDDCC9",  "#CE1261") 
p1 = DimPlot(obj, label = F, reduction = 'umap', group.by = 'level3_cell',cols = my_color, pt = 0.8)
ggsave("01_cell.pdf",plot=p1, width = 5, height = 3)

annotation_gene = c('FCGR3A','FGFBP2','GZMB','GZMH','GZMM','GNLY','TBX21', #c0
'NCAM1','GZMK','XCL1','XCL2','HSPA1A','ZNF683','CXCR4',
'KLRG1','ITGAL','IKZF3')#, #c1)
d1 = DotPlot(obj, features = unique(annotation_gene), cluster.idents = F) +
  RotatedAxis() + scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "BuPu")) + 
  theme(axis.title.x = element_blank(),strip.text.x = element_text(angle=90, size = 10)) +
  theme(strip.background =element_rect(fill="#F0F8FF"),legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=6),legend.text = element_text(size=6)) +
  FontSize(x.text = 9, main = 12)+labs(x=NULL, y=NULL)
ggsave("02_dp_clusters.pdf",plot=d1, width = 6, height = 1.5)

setwd('../')

dir = '04_B'
dir.create(dir)
setwd(dir)

load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/01_data/01_obj_b_5240_v1.Rdata')
obj=obj_b

table(obj$level3_cell)

cell_rank = names(table(obj$level3_cell))
obj$level3_cell = factor(obj$level3_cell, levels = cell_rank)

Idents(obj) <- "level3_cell"
my_color = c('#CB78A6','#D35F00','#F7EC44','#FCB93E','#009D73')
p1 = DimPlot(obj, label = F, reduction = 'umap', group.by = 'level3_cell',cols = my_color, pt = 0.5)
ggsave("01_cell.pdf",plot=p1, width = 5, height = 3)

annotation_gene = c("CD79A","CD79B","BANK1","SELL", "FCER2", "IGHD","IGHM", #c0
"CD27", "JCHAIN","MZB1", "SDC1", "TNFRSF17", "XBP1","IGHG1", "IGHG2", "IGHG3","IGHG4","IGHA1", "IGHA2",#c1
"CD24", "BCL6", "IL4R", "BACH2","CXCR5", "FCRLA", "AICDA", "BCL2A1",#c3
"MKI67", "TUBB", "STMN1","TYMS", "IRF4", "PRDM1" #c4
)

d1 = DotPlot(obj, features = unique(annotation_gene), cluster.idents = F) +
  RotatedAxis() + scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "BuPu")) + 
  theme(axis.title.x = element_blank(),strip.text.x = element_text(angle=90, size = 10)) +
  theme(strip.background =element_rect(fill="#F0F8FF"),legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=6),legend.text = element_text(size=6)) +
  FontSize(x.text = 9, main = 12)+labs(x=NULL, y=NULL)
ggsave("02_dp_clusters.pdf",plot=d1, width = 9, height = 1.9)

setwd('../')

dir = '05_myeloid'
dir.create(dir)
setwd(dir)

load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/02_result/11_check_myeloid_250424/00_obj.RData')
obj
table(obj$level4_cell)

cell_rank = c('Mac_c0_APOE_CD163','Neu_c1_IL1B_CSF3R','Mac_c2_C1QC_MKI67','cDC2_c3_FCER1A',
'Mono_c4_CD16','pDCs_c5_GZMB','Mast_c6_TPSAB1','cDC3_c7_FSCN1','cDC1_c8_CLEC9A')
obj$level4_cell = factor(obj$level4_cell, levels = cell_rank)


Idents(obj) <- "level4_cell"
my_color = c('#B6D98C','#BD56A2','#1FA759','#FBCB1F','#FFE2DE','#F0C1ED','#1879B6','#9B2463','#86CEFF')
p1 = DimPlot(obj, label = F, reduction = 'umap', group.by = 'level4_cell',
cols = my_color, order = rev(cell_rank), pt = 0.45)
ggsave("01_cell.pdf",plot=p1, width = 5, height = 3)

annotation_gene = c('C1QA','C1QB','C1QC','APOE','CCL18',"MSR1", "MRC1","CD163",'MERTK', 'FPR3', 'TREM2', 'MS4A4A', 
                    'SLCO2B1', 'PLTP', 'STAB1', 'SNX6', 'CD68','CXCL9',#C0 - macro
                    "S100A8","S100A9","FCN1",'S100A12','CSF3R',"CD14",'CYP1B1',"NLRP3","EREG",'CD36', 
                    'IL1B',#C1- monocyte-cd14
                    'MKI67','STMN1','TUBB',#C2 PRO
                    "FCER1A","CD1C","CD1E",'HLA-DQA1',#C3 cdc2
                    "FCGR3A","LST1","LILRB2",'KLF2',  #C4 monocyte 
                    "GZMB","LILRA4","LRRC26","IL3RA","CLEC4C","SCT", 'CXCR4',#C5 pdc
                    "MS4A2","CPA3","TPSAB1","TPSB2", #C6- mast
                    'FSCN1','CCR7','LAMP3','CCL19', 'CCL17','NFKB2', 'GPR132',  #C7- cdc3
                    "CLEC9A","BATF3","XCR1" #C8 cdc1
                    )

d1 = DotPlot(obj, features = unique(annotation_gene), cluster.idents = F) +
  RotatedAxis() + scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "BuPu")) + 
  theme(axis.title.x = element_blank(),strip.text.x = element_text(angle=90, size = 10)) +
  theme(strip.background =element_rect(fill="#F0F8FF"),legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=6),legend.text = element_text(size=6)) +
  FontSize(x.text = 9, main = 12)+labs(x=NULL, y=NULL)
ggsave("02_dp_clusters.pdf",plot=d1, width = 14, height = 2.5)

setwd('../')

dir = '06_Fib'
dir.create(dir)
setwd(dir)

load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/01_data/00_obj_fi_4512.Rdata')
obj=obj_fi

table(obj$level3_cell)

cell_rank = c('myCAF_c0_ASPN','myCAF_c1_SERPINE1','myCAF_c2_ACTA2','iCAF_c3_C3',
'apCAF_c4_CD74','prolifCAF_c5_MKI67')
obj$level3_cell = factor(obj$level3_cell, levels = cell_rank)


Idents(obj) <- "level3_cell"
my_color = c('#1B9E77','#B2DF8A','#E3BE00','#FB9A99','#E7298A',
            '#910241','#00CDD1')
p1 = DimPlot(obj, label = F, reduction = 'umap', group.by = 'level3_cell',
cols = my_color, order = rev(cell_rank), pt = 0.45)
ggsave("01_cell.pdf",plot=p1, width = 5, height = 3)

annotation_gene = c("ASPN",'COL12A1',"COL11A1", 'COL1A1','COL10A1','ITGA11', 'SDC1', 'COL3A1','POSTN','MMP14', #c0 myCAF
'SERPINE1','TIMP1','CCL2','MT2A','VEGFA','CXCL10','CCL2', 'ANGPTL4','GAPDH','SLC2A3','LOX','IGFBP3',  #c1,CAF-SERPINE1
'MCAM','A2M','CRIP1','SPARCL1','ACTA2','COL4A1','COL4A2','MYH11','NOTCH3', 'MYLK', #c2, CAF-
'PLA2G2A','CFD','IGF1','APOD','C3','C7', 'EFEMP1','PI16','CLU','CXCL12','CCDC80','PLAC9','CXCL14','ADH1B',#c3 iCAF
'CD74','PDPN','IL7R','MMP11','VCAM1','MMP13','HLA-B','BST2','IL32','SOD2','CCL11','IFI6','B2M','HLA-A','IFI30',#c4 apCAF
'MKI67','STMN1','TUBB','TYMS','CENPF','H2AFZ','HIST1H4C','TOP2A','ASPM','TPX2','RRM2'  #c5 proliferatingCAF
)

d1 = DotPlot(obj, features = unique(annotation_gene), cluster.idents = F) +
  RotatedAxis() + scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "BuPu")) + 
  theme(axis.title.x = element_blank(),strip.text.x = element_text(angle=90, size = 10)) +
  theme(strip.background =element_rect(fill="#F0F8FF"),legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=6),legend.text = element_text(size=6)) +
  FontSize(x.text = 9, main = 12)+labs(x=NULL, y=NULL)
ggsave("02_dp_clusters.pdf",plot=d1, width = 14, height = 2.5)

setwd('../')


dir = '07_End'
dir.create(dir)
setwd(dir)

load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/01_data/01_obj_endo_518_v1.Rdata')
obj=obj_endo

table(obj$level3_cell)

cell_rank = names(table(obj$level3_cell))
obj$level3_cell = factor(obj$level3_cell, levels = cell_rank)

Idents(obj) <- "level3_cell"
my_color = c('#F4931A','#E41E27','#FDE0DD','#B3DE69','#BD56A2','#B9E5EE')
p1 = DimPlot(obj, label = F, reduction = 'umap', group.by = 'level3_cell',
cols = my_color, order = rev(cell_rank), pt = 0.8)
ggsave("01_cell.pdf",plot=p1, width = 5, height = 3)

annotation_gene = c('PECAM1','CDH5','CLDN5','ERG', #pan-endothelial c0_PanECS_CDH5_CLDN5
'ACKR1','VWF','NR2F2','SELP','SELE', #c1_Venous_ACKR1_SELP
'RGCC','SPARC','EDN1', 'CD93', 'PLVAP',#c2_GeneralCapillary_RGCC_SPARC
'ZNF385D','VCAM1','EBF1','EBF3','MEOX1','MEOX2','IGFBP7', #c3_Venous_ACKR1_ZNF385D
'PROX1','FLT4','PDPN','CCL21','SEMA3A','SEMA3D','TBX1','NR2F1', #c4_LymphaticECS_PROX1_PDPN
'S100A4',#Capillary 
'EFNB2','SOX17','BMX','SEMA3G','HEY1','LTBP4','FBLN5','GJA5','GJA4' #c5_Arterial_EFNB2_GJA5
)

d1 = DotPlot(obj, features = unique(annotation_gene), cluster.idents = F) +
  RotatedAxis() + scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "BuPu")) + 
  theme(axis.title.x = element_blank(),strip.text.x = element_text(angle=90, size = 10)) +
  theme(strip.background =element_rect(fill="#F0F8FF"),legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=6),legend.text = element_text(size=6)) +
  FontSize(x.text = 9, main = 12)+labs(x=NULL, y=NULL)
ggsave("02_dp_clusters.pdf",plot=d1, width = 11, height = 2.3)

setwd('../')


