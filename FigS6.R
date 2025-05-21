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

dir = '06_FigS6'
dir.create(dir)
setwd(dir)

load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/01_data/00_obj_MTAP_neg_49161.Rdata')
obj = obj_MTAP_neg

Idents(obj) <- "level1_cell"
levels(obj) <- c("Epithelial","T&NK", "B","Plasma", 'Myeloid','Mast', "pDCs","Fibroblasts", "Endothelial")
my_color = c('#9ccb86', #epi
             "#e9e29c", #T
             "#C75DAB", "#D691C1", #b
             '#FFAE19', #my
             '#cc0000',# mast
             '#ffcccc', #pdc
             "#A7D3D4",# fib
             '#009B9E', #endo
             "#A7D3D4", '#E4E4E4',"#009B9E","#42B7B9" )
pt=0.5
p1 = DimPlot(obj, label = F, reduction = 'umap',group.by = 'Timepoint', pt.size = pt,cols = my_color)
ggsave("01_FigS6a.pdf",plot=p1, width = 7, height = 5)




