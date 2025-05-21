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



#FigS6c
load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/01_data/00_obj_MTAP_neg_49161.Rdata')
obj_MTAP_neg

dt = as.data.frame(cbind(table(obj$Timepoint,obj$level1_cell_new)))
dt$`B&Plasma` = rowSums(dt[,c('B','Plasma')])
dt$`Myeloids` = rowSums(dt[,c('Myeloid','Mast','pDCs')])
# dt = dt[,cell_order]
# dt = dt[,c(LymphoidCells,MyeloidCells, StromalCells, TumorCells)]
output_df = dt[, c('T&NK','B&Plasma','Myeloids','Fibroblasts','Endothelial','Epithelial')]
# write.xlsx(output_df, file = "01_Data2ROE_l1.xlsx", rowNames = TRUE)
# data = read.csv('01_sample_cells.csv',row.names = 1)

M = output_df
Xsq <- chisq.test(M,correct=F)
ROE <- Xsq$observed/Xsq$expected
# write.xlsx(ROE, file = "02_Data_ROE_l1.xlsx", rowNames = TRUE)

library(ComplexHeatmap)
bk <- c(seq(0,0.9,by=0.01),seq(1,2,by=0.01))
# Define custom legend range and colors
library(circlize)
breaks = seq(0, 2, 0.5)
custom_colors <- colorRamp2(
  breaks = seq(0, 2, 0.5),            # Define custom range for the legend
  colors = colorRampPalette(colors = c("white", "pink", "#CE1261"))(length(breaks))
)

label_matrix <- matrix("", nrow = nrow(ROE), ncol = ncol(ROE))
label_matrix[ROE > 1 & ROE <= 1.2] <- "*"
label_matrix[ROE > 1.2 & ROE <= 1.4] <- "**"
label_matrix[ROE > 1.4] <- "***"
label_matrix = t(label_matrix)


pdf(file = 'FigS6c.pdf',width = 3, height = 3)
Heatmap(
  t(ROE),
  name = "Ro/e",  
  col = custom_colors,
  cluster_rows = F,
  cluster_columns = FALSE,
  #row_split = cell_categories,  # Add row breaks by categories
  #row_order = original_row_order, 
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(label_matrix[i, j], x, y, gp = gpar(fontsize = 12, col = "black"))
  }
)
dev.off()


#### FigS6d
load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/02_result/11_check_myeloid_250424/00_obj.RData')
obj
table(obj$level4_cell)

# cell_rank = c('Mac_c0_APOE_CD163','Neu_c1_IL1B_CSF3R','Mac_c2_C1QC_MKI67','cDC2_c3_FCER1A',
# 'Mono_c4_CD16','pDCs_c5_GZMB','Mast_c6_TPSAB1','cDC3_c7_FSCN1','cDC1_c8_CLEC9A')
# obj$level4_cell = factor(obj$level4_cell, levels = cell_rank)

dt = as.data.frame(cbind(table(obj$Timepoint,obj$level4_cell)))
# dt = dt[,cell_order]
# dt = dt[,c(LymphoidCells,MyeloidCells, StromalCells, TumorCells)]
output_df = dt
# write.xlsx(output_df, file = "04_Data2ROE_myeloid.xlsx", rowNames = TRUE)
# data = read.csv('01_sample_cells.csv',row.names = 1)

M = output_df
Xsq <- chisq.test(M,correct=F)
ROE <- Xsq$observed/Xsq$expected

library(ComplexHeatmap)
bk <- c(seq(0,0.9,by=0.01),seq(1,2,by=0.01))

# Define custom legend range and colors
library(circlize)
breaks = seq(0, 2, 0.5)
custom_colors <- colorRamp2(
  breaks = seq(0, 2, 0.5),            # Define custom range for the legend
  colors = colorRampPalette(colors = c("white", "pink", "#CE1261"))(length(breaks))
)

label_matrix <- matrix("", nrow = nrow(ROE), ncol = ncol(ROE))
label_matrix[ROE > 1 & ROE <= 1.2] <- "*"
label_matrix[ROE > 1.2 & ROE <= 1.4] <- "**"
label_matrix[ROE > 1.4] <- "***"
label_matrix = t(label_matrix)

pdf(file = 'FigS6d.pdf',width = 3.4, height = 3.7)
Heatmap(
  t(ROE),
  name = "Ro/e",  
  col = custom_colors,
  cluster_rows = F,
  cluster_columns = FALSE,
  #row_split = cell_categories,  # Add row breaks by categories
  #row_order = original_row_order, 
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(label_matrix[i, j], x, y, gp = gpar(fontsize = 12, col = "black"))
  }
)
dev.off()




