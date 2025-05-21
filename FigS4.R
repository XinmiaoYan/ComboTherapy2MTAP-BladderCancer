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

setwd('/rsrch6/home/genomic_med/lwang22_lab/Xinmiao/02_MTAP_clinical_trial_remove_pt19/03_result/03_overview_data_230626/03_result_230731/02_result/15_Source_data/Figs')
set.seed(999)

dir = 'FigS4'
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
ggsave("01_FigS4a.pdf",plot=p1, width = 5, height = 3)


# Extract UMAP coordinates
umap_coords <- Embeddings(obj, reduction = "umap")
metadata <- obj@meta.data
table(metadata$seurat_clusters, metadata$seurat_clusters_new)
identical(rownames(umap_coords), rownames(metadata))
# Combine UMAP and metadata
output_df <- cbind(UMAP_1 = umap_coords[, 1],
                   UMAP_2 = umap_coords[, 2],
                   metadata)

head(output_df)
# Write to a tab-separated file
write.xlsx(output_df, file = "01_umap_with_metadata.xlsx", rowNames = TRUE)

# check adenosine between CB vs NCB and bestresponse.
Adenosine_Receptor = c('ADORA1','ADORA2A', 'ADORA2B','ADORA3')
Pro_Inflammatory = c('NFKB1','RELA','IL1B','IL15','IL18','CD40','PTGS2','CXCL9','CXCL10')
Anti_Inflammatory = c('IL10','TGFB1','NR4A1','NR4A2','MRC1','APOE','SPP1','C1QA','C1QC')
Angiogenesis=c('VEGFA','VEGFB','HIF1A',"CXCL8", "CXCL2",'MMP2','MMP9','TYMP','VCAN','CD44',
'FYN','E2F3','ITGAV','CXCR4','PTK2','CCND2','EZH2')
Phagocytosis=c('FCGR1A','FCGR2A','FCGR2B','FCGR3A','MARCO','MSR1','CD163','CD68',
'MERTK','C1QB','CD36','AXL')
Antigen_Presentation = c('B2M','CALR','CANX','CD4','CD74','CIITA','CREB1','CTSB','CTSS',
'PSMB1','PSMB2','PSMB3','PSMB4','PSMB5','PSMB6','PSMB7','PSMB8','PSMB9','PSMB10',
'TAP1','TAP2','TAPBP','PDIA3',
"HLA-A","HLA-B","HLA-C","HLA-E","HLA-F",
'HLA-DQA1','HLA-DQB1',"HLA-DRA","HLA-DRB1","HLA-DRB5", "HLA-DOA",
"HLA-DPA1", "HLA-DPB1","HLA-DMA", "HLA-DMB"
) # PMID: 30215098
Checkpoint = c('CD274','CD80','CD86','VSIR','LGALS9','HAVCR2','NECTIN2',
'SIGLEC10','SIRPA','TREM2','CD276','LAIR1','IDO1','BTLA')
Cytokine = c('LTB','TNFSF10','TNF','CD70','TNFSF4','IL16','TNFSF14','IL1RN','OSM',
'TXLNA','TNFSF13','TNFSF13B','TNFSF12')
Chemokine =	c('CCL5','CCL18','CCL2','CCL17','CCL19','CCL22','CXCL16','CCL3','CCL4')

# Create annotation dataframe
annotation_gene <- data.frame(
  gene = c(Adenosine_Receptor,Pro_Inflammatory, Anti_Inflammatory, 
           Angiogenesis, Phagocytosis,
           Antigen_Presentation, Checkpoint,Cytokine,Chemokine),
  type_f = c(rep('Adenosine_Receptor', length(Adenosine_Receptor)),
             rep('Pro_Inflammatory', length(Pro_Inflammatory)),
             rep('Anti_Inflammatory', length(Anti_Inflammatory)),
             rep('Angiogenesis', length(Angiogenesis)),
             rep('Phagocytosis', length(Phagocytosis)),
             rep('Antigen_Presentation', length(Antigen_Presentation)),
             rep('Checkpoint', length(Checkpoint)),
             rep('Cytokine', length(Cytokine)),
             rep('Chemokine', length(Chemokine))
))
head(annotation_gene)
rank = c('Adenosine_Receptor','Pro_Inflammatory','Anti_Inflammatory','Angiogenesis','Phagocytosis',
         'Antigen_Presentation', 'Checkpoint', 'Cytokine','Chemokine')
annotation_gene$type_f = factor(annotation_gene$type_f, levels = rank)
head(annotation_gene)

# calculate the module score
MPs = list(Pro_Inflammatory=Pro_Inflammatory,
           Anti_Inflammatory = Anti_Inflammatory,
           Angiogenesis = Angiogenesis,
           Phagocytosis = Phagocytosis,
           Antigen_Presentation = Antigen_Presentation,
           Checkpoint = Checkpoint,
           Cytokine = Cytokine,
           Chemokine = Chemokine
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



dt = obj@meta.data[,c('level4_cell','Timepoint','BestResponse','Clinical Benefit','BestResponse_Group',names(MPs))]

for (i in seq_along(MPs)) {
    name = names(MPs)[i] 

    p = ggplot(dt, aes_string(x = "level4_cell", y = name, fill = "level4_cell")) +
        geom_boxplot(colour = '#86817c', outlier.size = 0, outlier.stroke = 0) +
        theme(legend.position = "none") +
        scale_fill_manual(values = my_color) +
        stat_compare_means(method = "kruskal.test", aes(label = ..p.format..)) +
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

    ggsave(paste0('02_FigS4b_', name, "_boxplot.pdf"), plot = p, width = 4, height = 3.5)
}

write.xlsx(dt, file = "02_FigS4b.xlsx", rowNames = TRUE)

setwd('../')