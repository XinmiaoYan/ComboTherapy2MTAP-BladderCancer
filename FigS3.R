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

dir = '03_FigS3'
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

setwd('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/02_result/15_Source_data/Figs')
set.seed(999)

dir = 'Fig_myeloids_rank_adenosine'
dir.create(dir)
setwd(dir)

load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/02_result/11_check_myeloid_250424/00_obj.RData')
obj

cell_rank = c('Mac_c0_APOE_CD163','Neu_c1_IL1B_CSF3R','Mac_c2_C1QC_MKI67','cDC2_c3_FCER1A',
'Mono_c4_CD16','pDCs_c5_GZMB','Mast_c6_TPSAB1','cDC3_c7_FSCN1','cDC1_c8_CLEC9A')
obj$level4_cell = factor(obj$level4_cell, levels = cell_rank)

my_color = c("#C75DAB", 
             "#e9e29c", 
             "#42B7B9", 
             '#eeb479',
             '#9ccb86','#39b185','#2887a1',
             '#D3D3D3','#E4E4E4','#A1A1A1')

p1 = DimPlot(obj, label = F, reduction = 'umap', group.by = 'level4_cell',cols = my_color)
ggsave("01_cell.pdf",plot=p1, width = 5, height = 3)

A2A_BR_M2_genes = c('ADORA2A', 'ADORA2B','NR4A1','NR4A2','MRC1','APOE','SPP1','C1QA','C1QB','C1QC')
Antigen_Presentation = c('CD74',
#'TAP1','TAP2','TAPBP',
'HLA-DQA1','HLA-DQB1',"HLA-DRA","HLA-DRB1", "HLA-DOA",
"HLA-DPA1", "HLA-DPB1","HLA-DMA", "HLA-DMB"
) # PMID: 30215098


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

obj$Phagocytosis = NULL
obj$Angiogenesis = NULL
obj$Checkpoint = NULL

for(i in seq(1,length(names(MPs)),1)){
    # name = i
    dir_name = names(MPs[i])
    obj = AddModuleScore(obj,features = list(MPs[[i]]), name = dir_name)
    colnames(obj@meta.data) <- gsub(x = colnames(obj@meta.data),
                                    pattern = paste0(dir_name,1),
                                    replacement = dir_name
                                    )
}

save(obj, file = '00_obj.Rdata')


# 
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

  dot_data <- d1$data
  # Preview the data
  head(dot_data)
  # Write to a tab-separated file
  write.xlsx(dot_data, file = "02_dotplot_data.xlsx", rowNames = TRUE)


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

  write.xlsx(dt, file = "03_boxplot_data.xlsx", rowNames = TRUE)

  setwd('../')
}

setwd('../')



### rank
dir = '02_ClinicalBenefit_rank'
dir.create(dir)
setwd(dir)

# # M2
# load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/02_result/03_level3/02_Fig3_baseline/00_object_l3_20375.Rdata')
# object  #20375 l3
load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/02_result/02_level1/02_Fig3_baseline/00_object_20375.Rdata')
object #l1
table(object$level1_cell)

obj_baseline = subset(object, subset = Timepoint == 'Baseline')
MAC2_baseline = subset(obj, subset = Timepoint == 'Baseline')

patient_list <- NULL  # 创建一个空的列表
for (i in 1:11) {
  element <- paste("Pt-", i, sep = "") 
  patient_list = c(patient_list,element) 
}
print(patient_list) 

clinical_df = object@meta.data[c('PatientID_v1','Biopsy site','Gender','Prior IO','BestResponse','Mixed Response','Clinical Benefit','PFS (months)')]
clinical_info <- clinical_df[!duplicated(clinical_df$PatientID_v1), ]
rownames(clinical_info)=NULL
clinical_info <- column_to_rownames(clinical_info, var = "PatientID_v1")
clinical_info = clinical_info[c('Biopsy site','Gender','Prior IO','BestResponse','Mixed Response','Clinical Benefit','PFS (months)')]
clinical_info = clinical_info[patient_list,]

# colnames(clinical_info) = c("BiopsySite", "Gender",  "PriorIO", "BestResponse", "MixedResponse","ClinicalBenefit") 
annotation_colors <- list(
  `BiopsySite` = c(`Lymph node` ='#deeaee', Bone='#b1cbbb', Lung='#eea29a',`Pelvic met` = '#c94c4c'),
  Gender = c(Male = "#afc9e9", Female = "#ffb3d9"),
  `PrioIO` = c(No = 'grey', Yes = '#82b74b'),
  `MixedResponse` = c(No = 'grey', Yes = '#F5DF4D'),
  BestResponse = c(PR = '#e06000', SD = '#ffab00',PD = '#D4E4F7'),
  `ClinicalBenefit` = c(CB = '#FFD966',  `Non-CB`   = '#a8adb4')
)


dir = '01_All'
dir.create(dir)
setwd(dir)


object

cell_number = as.data.frame(cbind(table(object$PatientID_v1, object$level1_cell)))
cell_number = cell_number[patient_list,]

target_cell = as.data.frame(cbind(table(obj$PatientID_v1, obj$level4_cell)))
target_cell = target_cell[patient_list,]
cell_compisition = (target_cell)/rowSums(cell_number)
cell_compisition = cell_compisition[,c('Mac_c0_APOE_CD163','Mac_c2_C1QC_MKI67')]


cell_compisition
# clinical_info = clinical_info[c('Biopsy site','Gender','Prior IO','BestResponse','Mixed Response','Clinical Benefit','PFS (months)')]
# colnames(clinical_info) = c('BiopsySite','Gender','PrioIO','BestResponse','MixedResponse','ClinicalBenefit','PFS(months)')
clinical_info$sample = rownames(clinical_info)
identical(clinical_info$sample , rownames(cell_compisition))
data = cbind(clinical_info,cell_compisition)

dt = data
dt = arrange(dt, desc(CellProportion))
dt$Patient <- factor(dt$Patient, levels = rownames(dt))

df = dt
write.xlsx(df, file = "02_data.xlsx", rowNames = TRUE)

for (i in colnames(df)[1:6]){
p = ggplot(df, aes(x=Patient, y=CellProportion, color = get(i))) +
  geom_segment(aes(x=Patient, xend=Patient, y=0, yend=CellProportion)) +
  geom_point(size=5,fill=alpha(0.9),alpha=0.9)  +
  theme_light() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 9)
  ) +ggtitle('Mac_c0&c2')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1),
            plot.title = element_text(hjust = 0.5))+
  labs(x = NULL, y = "Proportion in all cells(%)") +
  scale_color_manual(values = annotation_colors[[i]])+
  labs(color = i) 
ggsave(plot =p, file = paste( i,'.pdf', sep=''), width = 4, height=4)
}


setwd('../')


# Non-epithelial
dir = '02_Non_Epi'
dir.create(dir)
setwd(dir)

object = subset(object, subset =level1_cell =='Epithelial',invert = T )
cell_number = as.data.frame(cbind(table(object$PatientID_v1, object$level1_cell)))
cell_number = cell_number[patient_list,]

target_cell = as.data.frame(cbind(table(obj$PatientID_v1, obj$level4_cell)))
target_cell = target_cell[patient_list,]

cell_compisition = (target_cell)/rowSums(cell_number)
cell_compisition = cell_compisition[,c('Mac_c0_APOE_CD163','Mac_c2_C1QC_MKI67')]


cell_compisition
# clinical_info = clinical_info[c('Biopsy site','Gender','Prior IO','BestResponse','Mixed Response','Clinical Benefit','PFS (months)')]
# colnames(clinical_info) = c('BiopsySite','Gender','PrioIO','BestResponse','MixedResponse','ClinicalBenefit','PFS(months)')
clinical_info$sample = rownames(clinical_info)
identical(clinical_info$sample , rownames(cell_compisition))
data = cbind(clinical_info,cell_compisition)

dt = data
dt = arrange(dt, desc(CellProportion))
dt$Patient <- factor(dt$Patient, levels = rownames(dt))

df = dt
write.xlsx(df, file = "02_data.xlsx", rowNames = TRUE)


cell_compisition = data.frame(
  Patient = names(rowSums(target_cell)/rowSums(cell_number)),
  CellProportion = rowSums(target_cell)/rowSums(cell_number)
)


cell_compisition
identical(clinical_info$sample , rownames(cell_compisition))
data = cbind(clinical_info,cell_compisition)

dt = data
dt = arrange(dt, desc(CellProportion))
dt$Patient <- factor(dt$Patient, levels = rownames(dt))

df = dt
for (i in colnames(df)[1:6]){
p = ggplot(df, aes(x=Patient, y=CellProportion, color = get(i))) +
  geom_segment(aes(x=Patient, xend=Patient, y=0, yend=CellProportion)) +
  geom_point(size=5,fill=alpha(0.9),alpha=0.9)  +
  theme_light() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 9)
  ) +ggtitle('Mac_c0&c2')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1),
            plot.title = element_text(hjust = 0.5))+
  labs(x = NULL, y = "Proportion in Non-epithelial cells(%)") +
  scale_color_manual(values = annotation_colors[[i]])+
  labs(color = i) 
ggsave(plot =p, file = paste( i,'.pdf', sep=''), width = 4, height=4)
}

setwd('../')


# Non-epithelial
dir = '03_Immune'
dir.create(dir)
setwd(dir)

object = subset(object, subset =level1_cell %in% c('Endothelial','Fibroblasts'),invert = T )
cell_number = as.data.frame(cbind(table(object$PatientID_v1, object$level1_cell)))
cell_number = cell_number[patient_list,]

target_cell = as.data.frame(cbind(table(obj$PatientID_v1, obj$level4_cell)))
target_cell = target_cell[patient_list,]
cell_compisition = (target_cell)/rowSums(cell_number)
cell_compisition = cell_compisition[,c('Mac_c0_APOE_CD163','Mac_c2_C1QC_MKI67')]


cell_compisition = data.frame(
  Patient = names(rowSums(target_cell)/rowSums(cell_number)),
  CellProportion = rowSums(target_cell)/rowSums(cell_number)
)

cell_compisition
# clinical_info = clinical_info[c('Biopsy site','Gender','Prior IO','BestResponse','Mixed Response','Clinical Benefit','PFS (months)')]
# colnames(clinical_info) = c('BiopsySite','Gender','PrioIO','BestResponse','MixedResponse','ClinicalBenefit','PFS(months)')
clinical_info$sample = rownames(clinical_info)
identical(clinical_info$sample , rownames(cell_compisition))
data = cbind(clinical_info,cell_compisition)

dt = data
dt = arrange(dt, desc(CellProportion))
dt$Patient <- factor(dt$Patient, levels = rownames(dt))

df = dt
write.xlsx(df, file = "02_data.xlsx", rowNames = TRUE)



cell_compisition
identical(clinical_info$sample , rownames(cell_compisition))
data = cbind(clinical_info,cell_compisition)

dt = data
dt = arrange(dt, desc(CellProportion))
dt$Patient <- factor(dt$Patient, levels = rownames(dt))

df = dt
for (i in colnames(df)[1:6]){
p = ggplot(df, aes(x=Patient, y=CellProportion, color = get(i))) +
  geom_segment(aes(x=Patient, xend=Patient, y=0, yend=CellProportion)) +
  geom_point(size=5,fill=alpha(0.9),alpha=0.9)  +
  theme_light() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 9)
  ) +ggtitle('Mac_c0&c2')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1),
            plot.title = element_text(hjust = 0.5))+
  labs(x = NULL, y = "Proportion in Immune cells(%)") +
  scale_color_manual(values = annotation_colors[[i]])+
  labs(color = i) 
ggsave(plot =p, file = paste( i,'.pdf', sep=''), width = 4, height=4)
}

setwd('../')


## FigS3c-d
# # M2
# load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/02_result/03_level3/02_Fig3_baseline/00_object_l3_20375.Rdata')
# object  #20375 l3
load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/02_result/02_level1/02_Fig3_baseline/00_object_20375.Rdata')
object #l1
table(object$level1_cell)

obj_baseline = subset(object, subset = Timepoint == 'Baseline')
obj = subset(obj, subset = Timepoint == 'Baseline')

patient_list <- NULL  # 创建一个空的列表
for (i in 1:11) {
  element <- paste("Pt-", i, sep = "") 
  patient_list = c(patient_list,element) 
}
print(patient_list) 

clinical_df = object@meta.data[c('PatientID_v1','Biopsy site','Gender','Prior IO','BestResponse','Mixed Response','Clinical Benefit','PFS (months)')]
clinical_info <- clinical_df[!duplicated(clinical_df$PatientID_v1), ]
rownames(clinical_info)=NULL
clinical_info <- column_to_rownames(clinical_info, var = "PatientID_v1")
clinical_info = clinical_info[c('Biopsy site','Gender','Prior IO','BestResponse','Mixed Response','Clinical Benefit','PFS (months)')]
colnames(clinical_info) = c('BiopsySite','Gender','PrioIO','BestResponse','MixedResponse','ClinicalBenefit','PFS(months)')
clinical_info = clinical_info[patient_list,]

# colnames(clinical_info) = c("BiopsySite", "Gender",  "PriorIO", "BestResponse", "MixedResponse","ClinicalBenefit") 
annotation_colors <- list(
  `BiopsySite` = c(`Lymph node` ='#deeaee', Bone='#b1cbbb', Lung='#eea29a',`Pelvic met` = '#c94c4c'),
  Gender = c(Male = "#afc9e9", Female = "#ffb3d9"),
  `PrioIO` = c(No = 'grey', Yes = '#82b74b'),
  `MixedResponse` = c(No = 'grey', Yes = '#F5DF4D'),
  BestResponse = c(PR = '#e06000', SD = '#ffab00',PD = '#D4E4F7'),
  `ClinicalBenefit` = c(CB = '#FFD966',  `Non-CB`   = '#a8adb4')
)


dir = '01_All'
dir.create(dir)
setwd(dir)

object
cell_number = as.data.frame(cbind(table(object$PatientID_v1, object$level1_cell)))
cell_number = cell_number[patient_list,]

target_cell = as.data.frame(cbind(table(obj$PatientID_v1, obj$level4_cell)))
target_cell = target_cell[patient_list,]
cell_compisition = (target_cell)/rowSums(cell_number)
cell_compisition = cell_compisition[,c('Mac_c0_APOE_CD163','Mac_c2_C1QC_MKI67')]


cell_compisition
# clinical_info = clinical_info[c('Biopsy site','Gender','Prior IO','BestResponse','Mixed Response','Clinical Benefit','PFS (months)')]
# colnames(clinical_info) = c('BiopsySite','Gender','PrioIO','BestResponse','MixedResponse','ClinicalBenefit','PFS(months)')
clinical_info$sample = rownames(clinical_info)
identical(clinical_info$sample , rownames(cell_compisition))


for(i in colnames(cell_compisition)){

  dt1 = cell_compisition[i]
  colnames(dt1) = 'CellProportion'

  
  data = cbind(clinical_info,dt1)
  dt = data
  dt = arrange(dt, desc(CellProportion), desc(sample))
  dt$Patient = dt$sample
  dt$Patient <- factor(dt$Patient, levels = rownames(dt))

  df = dt
  cell = i

  for (j in colnames(df)[1:6]){
  p = ggplot(df, aes(x=Patient, y=CellProportion, color = get(j))) +
    geom_segment(aes(x=Patient, xend=Patient, y=0, yend=CellProportion)) +
    geom_point(size=5,fill=alpha(0.9),alpha=0.9)  +
    theme_light() + 
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      text = element_text(size = 10)
    ) +ggtitle(cell)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1),
              plot.title = element_text(hjust = 0.5))+
    labs(x = NULL, y = "Proportion in all cells(%)") +
    scale_color_manual(values = annotation_colors[[j]])+
    labs(color = j) 
  ggsave(plot =p, file = paste(cell, '_',j,'.pdf', sep=''), width = 4.5, height=4)
  }

}

setwd('../')


# Non-epithelial
dir = '02_Non_Epi'
dir.create(dir)
setwd(dir)

object = subset(object, subset =level1_cell =='Epithelial',invert = T )
cell_number = as.data.frame(cbind(table(object$PatientID_v1, object$level1_cell)))
cell_number = cell_number[patient_list,]

target_cell = as.data.frame(cbind(table(obj$PatientID_v1, obj$level4_cell)))
target_cell = target_cell[patient_list,]
cell_compisition = (target_cell)/rowSums(cell_number)
cell_compisition = cell_compisition[,c('Mac_c0_APOE_CD163','Mac_c2_C1QC_MKI67')]


cell_compisition
# clinical_info = clinical_info[c('Biopsy site','Gender','Prior IO','BestResponse','Mixed Response','Clinical Benefit','PFS (months)')]
# colnames(clinical_info) = c('BiopsySite','Gender','PrioIO','BestResponse','MixedResponse','ClinicalBenefit','PFS(months)')
clinical_info$sample = rownames(clinical_info)
identical(clinical_info$sample , rownames(cell_compisition))


for(i in colnames(cell_compisition)){

  dt1 = cell_compisition[i]
  colnames(dt1) = 'CellProportion'

  
  data = cbind(clinical_info,dt1)
  dt = data
  dt = arrange(dt, desc(CellProportion), desc(sample))
  dt$Patient = dt$sample
  dt$Patient <- factor(dt$Patient, levels = rownames(dt))

  df = dt
  cell = i

  for (j in colnames(df)[1:6]){
  p = ggplot(df, aes(x=Patient, y=CellProportion, color = get(j))) +
    geom_segment(aes(x=Patient, xend=Patient, y=0, yend=CellProportion)) +
    geom_point(size=5,fill=alpha(0.9),alpha=0.9)  +
    theme_light() + 
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      text = element_text(size = 10)
    ) +ggtitle(cell)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1),
              plot.title = element_text(hjust = 0.5))+
    labs(x = NULL, y = "Proportion in Non-epithelial cells(%)") +
    scale_color_manual(values = annotation_colors[[j]])+
    labs(color = j) 
  ggsave(plot =p, file = paste(cell, '_',j,'.pdf', sep=''), width = 4.5, height=4)
  }

}


setwd('../')


# Non-epithelial
dir = '03_Immune'
dir.create(dir)
setwd(dir)

object = subset(object, subset =level1_cell %in% c('Endothelial','Fibroblasts'),invert = T )
cell_number = as.data.frame(cbind(table(object$PatientID_v1, object$level1_cell)))
cell_number = cell_number[patient_list,]

target_cell = as.data.frame(cbind(table(obj$PatientID_v1, obj$level4_cell)))
target_cell = target_cell[patient_list,]
cell_compisition = (target_cell)/rowSums(cell_number)
cell_compisition = cell_compisition[,c('Mac_c0_APOE_CD163','Mac_c2_C1QC_MKI67')]


cell_compisition
# clinical_info = clinical_info[c('Biopsy site','Gender','Prior IO','BestResponse','Mixed Response','Clinical Benefit','PFS (months)')]
# colnames(clinical_info) = c('BiopsySite','Gender','PrioIO','BestResponse','MixedResponse','ClinicalBenefit','PFS(months)')
clinical_info$sample = rownames(clinical_info)
identical(clinical_info$sample , rownames(cell_compisition))


for(i in colnames(cell_compisition)){

  dt1 = cell_compisition[i]
  colnames(dt1) = 'CellProportion'

  
  data = cbind(clinical_info,dt1)
  dt = data
  dt = arrange(dt, desc(CellProportion), desc(sample))
  dt$Patient = dt$sample
  dt$Patient <- factor(dt$Patient, levels = rownames(dt))

  df = dt
  cell = i

  for (j in colnames(df)[1:6]){
  p = ggplot(df, aes(x=Patient, y=CellProportion, color = get(j))) +
    geom_segment(aes(x=Patient, xend=Patient, y=0, yend=CellProportion)) +
    geom_point(size=5,fill=alpha(0.9),alpha=0.9)  +
    theme_light() + 
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      text = element_text(size = 10)
    ) +ggtitle(cell)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust= 1),
              plot.title = element_text(hjust = 0.5))+
    labs(x = NULL, y = "Proportion in Immune cells(%)") +
    scale_color_manual(values = annotation_colors[[j]])+
    labs(color = j) 
  ggsave(plot =p, file = paste(cell, '_',j,'.pdf', sep=''), width = 4.5, height=4)
  }

}

setwd('../')

#FigS3b
# cell composition l1 & l3



#FigS3e

load('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/02_result/11_check_myeloid_250424/00_obj.RData')
obj
# cell_rank = c('Mac_c0_APOE_CD163','Neu_c1_IL1B_CSF3R','Mac_c2_C1QC_MKI67','cDC2_c3_FCER1A',
# 'Mono_c4_CD16','pDCs_c5_GZMB','Mast_c6_TPSAB1','cDC3_c7_FSCN1','cDC1_c8_CLEC9A')
# obj$level4_cell = factor(obj$level4_cell, levels = cell_rank)

M1=c('TNF','CXCL9','CXCL10','CD86','NOS2','IL1A','IL1B','IL12A','IL12B','IL6')
M2=c('CD163','MRC1','FN1','APOE','ARG1','MRC1','IL10','SPP1','CCL9','CCL13','CCL18','CCL15','C1QA',
'C1QB','TGFB1','TGFB2','TGFB3','VEGFA','VEGFB','VEGFC','VEGFD')

## add module score#
sc1 = subset(obj, subset = level4_cell %in% c('Mac_c0_APOE_CD163','Mac_c2_C1QC_MKI67')
                           & Timepoint == 'Baseline')
MPs = list(M1_like=M1,
           M2_like = M2)
           
for(i in seq(1,length(names(MPs)),1)){
    dir_name = names(MPs[i])
    sc1 = AddModuleScore(sc1,features = list(MPs[[i]]), name = dir_name)
    colnames(sc1@meta.data) <- gsub(x = colnames(sc1@meta.data),
                                    pattern = paste0(dir_name,1),
                                    replacement = dir_name
                                    )
}

Modules = c('M1_like','M2_like')
d1 = DotPlot(sc1, features = Modules,cluster.idents = F, group.by='level4_cell', scale = F) +
  RotatedAxis() + scale_color_viridis() + 
  theme(axis.title.x = element_blank(),strip.text.x = element_text(angle=90, size = 10)) +
  theme(strip.background =element_rect(fill="#F0F8FF"),legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=6),legend.text = element_text(size=6)) +
  FontSize(x.text = 9, main = 12)+ylab(NULL)+xlab(NULL)
ggsave("FigS3e.pdf",plot=d1, width = 3.5, height = 1.2)


#FigS3f-g
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

  dot_data <- d1$data
  # Preview the data
  head(dot_data)
  # Write to a tab-separated file
  write.xlsx(dot_data, file = "02_dotplot_data.xlsx", rowNames = TRUE)


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

  write.xlsx(dt, file = "03_boxplot_data.xlsx", rowNames = TRUE)

  setwd('../')
}







