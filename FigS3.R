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

dir = '01_All'
dir.create(dir)
setwd(dir)

load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/01_data/00_obj_MTAP_neg_49161.Rdata')
obj = obj_MTAP_neg

Idents(obj) <- "level1_cell"
levels(obj) <- c("T&NK", "B","Plasma", 'Myeloid','Mast', "pDCs","Epithelial","Fibroblasts", "Endothelial")
annotation_gene = c("PTPRC", "CD3D", "CD3E", "GZMA", "NKG7",
                    "CD19", "MS4A1", "CD79A","JCHAIN","XBP1","MZB1",
                    "FCGR3A", "CD14", "LYZ", 
                    "C1QB","C1QC","CD68","TPSAB1", "TPSB2","LILRA4",'CLEC4C',
                    "EPCAM","KRT8", "COL1A1","COL1A2","COL3A1", "FAP", "PECAM1", "VWF")

d1 = DotPlot(obj, features = unique(annotation_gene), cluster.idents = F) +
  RotatedAxis() +  scale_color_gradient2(low = "grey", mid = "white", high = "#FB8604") + 
  theme(axis.title.x = element_blank(),strip.text.x = element_text(angle=90, size = 10)) +
  theme(strip.background =element_rect(fill="#F0F8FF"),legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=6),legend.text = element_text(size=6)) +
  FontSize(x.text = 9, main = 12)+labs(x=NULL, y=NULL)
ggsave("01_FigS3a.pdf",plot=d1, width = 8, height = 2.5)





