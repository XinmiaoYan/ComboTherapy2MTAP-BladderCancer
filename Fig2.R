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
p1 = DimPlot(obj, label = F, reduction = 'umap',pt.size = pt,cols = my_color)
ggsave("01_Fig2b.pdf",plot=p1, width = 7, height = 5)


#FigS2f-g
load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/02_result/11_check_myeloid_250424/00_obj.RData')
obj

cell_rank = c('Mac_c0_APOE_CD163','Neu_c1_IL1B_CSF3R','Mac_c2_C1QC_MKI67','cDC2_c3_FCER1A',
'Mono_c4_CD16','pDCs_c5_GZMB','Mast_c6_TPSAB1','cDC3_c7_FSCN1','cDC1_c8_CLEC9A')
obj$level4_cell = factor(obj$level4_cell, levels = cell_rank)

A2A_BR_M2_genes = c('ADORA2A', 'ADORA2B','NR4A1','NR4A2','MRC1','APOE','SPP1','C1QA','C1QB','C1QC')
Antigen_Presentation = c('CD74','HLA-DQA1','HLA-DQB1',"HLA-DRA","HLA-DRB1", "HLA-DOA","HLA-DPA1", "HLA-DPB1","HLA-DMA", "HLA-DMB") # PMID: 30215098


# Create annotation dataframe
annotation_gene <- data.frame(
  gene = c(A2A_BR_M2_genes,
           Antigen_Presentation),
  type_f = c(rep('A2A/BR-M2_genes', length(A2A_BR_M2_genes)),
             rep('Antigen_Presentation', length(Antigen_Presentation))
))
head(annotation_gene)
rank = c('A2A/BR-M2_genes','Antigen_Presentation')
annotation_gene$type_f = factor(annotation_gene$type_f, levels = rank)
head(annotation_gene)

obj@meta.data = mutate(obj@meta.data,
                      BestResponse_Group =case_when(BestResponse_Class %in% c('PR','SD')~ 'PR/SD',
                                              BestResponse_Class == 'PD'~ 'PD',
                                                    TRUE ~ NA_character_))

# calculate the module score
MPs = list(A2A_BR_M2_genes=A2A_BR_M2_genes,
           Antigen_Presentation_HLAII = Antigen_Presentation
           )

for(i in seq(1,length(names(MPs)),1)){
    # name = i
    dir_name = names(MPs[i])
    obj = AddModuleScore(obj,features = list(MPs[[i]]), name = dir_name)
    colnames(obj@meta.data) <- gsub(x = colnames(obj@meta.data),
                                    pattern = paste0(dir_name,1),
                                    replacement = dir_name
                                    )
}

# save(obj, file = '00_obj.Rdata')


#plot a box plot to show the MPs
dir = '01_BestResponse'
dir.create(dir)
setwd(dir)

for (cell in unique(obj$level4_cell)){
  
  obj1 = subset(obj, subset=level4_cell == cell & Timepoint == 'Baseline')

  dir.create(cell)
  setwd(cell)

  # 01 dotplot
  d1 <- DotPlot(obj1, group.by = 'BestResponse_Group', 
  features = annotation_gene$gene, cluster.idents = FALSE,scale = T) +
  RotatedAxis() + scale_color_gradientn(colors = RColorBrewer::brewer.pal(4, "Oranges"))+
  facet_grid(. ~ annotation_gene$type_f, scales = 'free', space = 'free') +
  theme(
    axis.title.x = element_blank(),
    strip.text.x = element_text(angle = 0, size = 10), # Strip text at 0 degrees (horizontal)
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10), # X-axis text at 90 degrees (vertical)
    axis.text.y = element_text(size = 10),
    strip.background = element_rect(fill = "#F0F8FF"),
    legend.key.size = unit(0.3, 'cm'), 
    legend.title = element_text(size = 9),
    panel.spacing.y = unit(3, "lines"), 
    panel.background = element_rect(fill = "white"), # Change panel background to white
    plot.background = element_rect(fill = "white"),  # Change plot background to white
    panel.grid.major = element_line(color = "lightgray", size = 0.15), # Change major grid lines to gray
    panel.grid.minor = element_line(color = "lightgray", size = 0.15), # Change minor grid lines to gray
    panel.border = element_rect(color = "darkgray", fill = NA, size = 1) 
  ) +
  FontSize(x.text = 9, main = 12) +
  xlab(NULL) + ylab(NULL)
  ggsave("01_dp_CB.pdf",plot=d1, width = 7, height = 1.8)

#   dot_data <- d1$data
#   # Preview the data
#   head(dot_data)
#   # Write to a tab-separated file
#   write.xlsx(dot_data, file = "02_dotplot_data.xlsx", rowNames = TRUE)


  #boxplot
  dt = obj1@meta.data[,c('level4_cell','Timepoint','BestResponse','Clinical Benefit','BestResponse_Group',names(MPs))]

  for (i in seq_along(MPs)) {
      name = names(MPs)[i] 

      p = ggplot(dt, aes_string(x = "BestResponse_Group", y = name, fill = "BestResponse_Group")) +
          geom_boxplot(colour = '#86817c', outlier.size = 0, outlier.stroke = 0) +
          theme(legend.position = "none") +
          scale_fill_manual(values = c('#FFD966','#a8adb4')) +
          stat_compare_means(method = "wilcox.test", aes(label = paste0("'P = ', ..p.format.."))) +
          theme(
              axis.title.x = element_blank(),
              strip.text.x = element_text(angle = 0, size = 10),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
              axis.text.y = element_text(angle = 90, vjust = 0.5, size = 10),
              strip.background = element_rect(fill = "#F0F8FF"),
              legend.key.size = unit(0.3, 'cm'),
              legend.title = element_text(size = 9),
              panel.spacing.y = unit(3, "lines"),
              panel.background = element_rect(fill = "white"),
              plot.background = element_rect(fill = "white"),
              panel.grid.major = element_line(color = "lightgray", size = 0.15),
              panel.grid.minor = element_line(color = "lightgray", size = 0.15),
              panel.border = element_rect(color = "darkgray", fill = NA, size = 1)
          )

      ggsave(paste0('03_', name, "_boxplot.pdf"), plot = p, width = 2, height = 3)
  }

#   write.xlsx(dt, file = "03_boxplot_data.xlsx", rowNames = TRUE)

  setwd('../')
}