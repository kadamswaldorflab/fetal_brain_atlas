library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(Seurat)
library(scCustomize)

#####

setwd("/gscratch/kawaldorflab/tykim/R_files/fetal_brain")
set.seed(123)

fb_seurat_v25 <- readRDS("/mmfs1/gscratch/kawaldorflab/jcorn427/fet_brain/fb_seurat_v25.RDS")
fb_seurat_v25 <- NormalizeData(fb_seurat_v25, assay = "RNA")
DefaultAssay(fb_seurat_v25) <- "RNA"

######

# Astro Markers Dot plot
new_x <- c("vRG", "oRG", "GPC", "AS2", "AS1", "AS0")
astro_subset <- subset(fb_seurat, idents = new_x)
Idents(fb_seurat) <- factor(Idents(fb_seurat), levels = new_x)
astro_subset$cell_type_v10 <- factor(
  x = astro_subset$cell_type_v10,
  levels = new_x
)
gradient <- c("blue", "white", "red")

astro_plot <- DotPlot_scCustom(astro_subset, features = c("GFAP", "NFIA", "SLC1A2", "SLC1A3", "SPARCL1",
                                                          "TIMP3", "APOE", "S100B", "VIM", "MT2A"), 
                               group.by = "cell_type_v10", flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

astro_plot +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1)) +
  scale_size_continuous(breaks = c(0, 20, 40, 60),
                        limits = c(0, 100)) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# Oligodendrocyte lineage Dot plot
x_axis = c("OPC0", "OPC1", "OPC2", "COP", "BCAS1+OL", "OL1", "OL2")
fb_subset <- subset(fb_seurat, idents = x_axis)
Idents(fb_subset) <- factor(Idents(fb_subset), levels = x_axis)
fb_subset$cell_type_v4 <- factor(
  x = fb_subset$cell_type_v4,
  levels = x_axis
)

gradient <- c("blue", "white", "red")

olig_plot <- DotPlot_scCustom(fb_subset, features = c("SOX10", "PDGFRA", "OLIG2", "PCDH17", "ENSMMUG00000056728", "PCDH15", "MMP16",
                                                      "CA10", "BCAS1", "ENPP6", "MAL", "MOG", "PLP1", "MBP"), group.by = "cell_type_v4", 
                              flip_axes = TRUE, remove_axis_titles = FALSE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster") 

olig_plot +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1)) +
  scale_size_continuous(breaks = c(0, 20, 40, 60), 
                        limits = c(0, 100)) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


# IPC Neural Lineage Dot plot (Subset)
excluded <- c("M1-MG", "M2-MG", "aMG", "Astro", "OPC2",
              "PC", "EC", "OL1", "EP", "BCAS1+OL",
              "COP", "OL2", "CD3+T", "EB", "OPC1", "OPC0", "Astro1", "Astro2", "B1 NSC")
filtered <- subset(fb_seurat_v20, subset = !(cell_type_v19 %in% excluded))

xaxis = c("IPC", "DL-ExN", "UL-ExN")
filtered <- subset(filtered, idents = xaxis)
Idents(filtered) <- factor(Idents(filtered), levels = xaxis)
filtered$cell_type_v19 <- factor(
  x = filtered$cell_type_v19,
  levels = xaxis
)
gradient <- c("blue", "white", "red")

# GLAST = SLC1A3, MPP5 = PALS1, CTGF = CCN2
new_features <- c( "PROX1", 'SP9', "DLX5", "DLX1", "FABP7", "GAD1", "GAD2", "NES", "PAX6", "PTPRZ1",
                   "SOX6", "SLC1A3", "VIM", "HOPX", "TNC", "CALB2", "TAC3", 'CXCR4', "EOMES", "TBR1", 
                   "SATB1", "SOX5", "POU3F2", "NR2F1", "NR2F2", 'NRP1', "CRYM", "SLC17A7", "SLC17A6", "SLC17A8", 
                   "NEUROD2", "NEUROD6", "TLE4", "SATB2", "UNC5D", "GRIA2", "CUX1", "CUX2",  "RORB", "GRIN2B", 
                   'NRP2', "NRGN"
)

p1 <- DotPlot_scCustom(filtered, features = new_features, group.by = "cell_type_v19", 
                       flip_axes = TRUE, x_lab_rotate = TRUE, remove_axis_titles = FALSE, 
                       assay = "RNA", colors_use = gradient) + labs(x = "Gene", y = "Cluster")
p1 +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1)) +
  scale_size_continuous(breaks = c(0, 20, 40, 60),
                        limits = c(0, 100)) +
  scale_color_gradientn(
    colors = gradient,
    breaks = c(-1, 0, 1)
  ) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
p1


# IPC Neural Lineage Dot plot (Entire Seurat Object)
ipc_plot <- DotPlot_scCustom(fb_seurat_v25, 
                             features = c("DCX", "NEUROD1", "EOMES", "TBR1", "SOX5", "POU3F2", "NR2F1", 'NRP1', "CRYM", "SLC17A6", 
                                          "NEUROD2", "TLE4", "SATB2", "UNC5D", "GRIA2", "CUX1", "CUX2", "RORB", "GRIN2B", "SLC17A7", 
                                          "CAMK2A", "NRGN"
                                          ), 
                             group.by = "cell_type_v25", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

ipc_plot_filtered <- ipc_plot + 
  scale_y_discrete(limits = c("ExN0", "ExN1", "ExN2", "ExN3"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
ipc_plot_filtered



# Oligodendrocyte lineage Dot plot (Subset)
x_axis = c("GPC", "OPC0", "OPC1", "OPC2", "COP", "OL1", "OL2")
fb_subset <- subset(fb_seurat_v20, idents = x_axis)
Idents(fb_subset) <- factor(Idents(fb_subset), levels = x_axis)
fb_subset$cell_type_v19 <- factor(
  x = fb_subset$cell_type_v19,
  levels = x_axis
)

gradient <- c("blue", "white", "red")

olig_plot <- DotPlot_scCustom(fb_subset, 
                              features = c("SOX10", "PDGFRA", "OLIG2", "PCDH17", "PCDH15", "MMP16", "CA10", "ENSMMUG00000056728", "BCAS1", "ENPP6", 
                                           "MAL", "MOG", "PLP1", "MBP"), 
                              group.by = "cell_type_v19", flip_axes = TRUE, remove_axis_titles = FALSE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster") 

olig_plot +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1)) +
  scale_size_continuous(breaks = c(0, 20, 40, 60),
                        limits = c(0, 100)) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


# Oligodendrocyte lineage Dot plot (Entire Seurat Object)
oligo_plot <- DotPlot_scCustom(fb_seurat_v25, features = c("SOX10", "PDGFRA", "OLIG2", "PCDH17", "PCDH15", "MMP16", "CA10", "ENSMMUG00000056728", "BCAS1", "ENPP6", 
                                                           "MAL", "MOG", "PLP1", "MBP"), 
                               group.by = "cell_type_v25", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

oligo_plot_filtered <- oligo_plot + 
  scale_y_discrete(limits = c("OPC0", "OPC1", "OPC2", "COP", "OL1", "OL2"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
oligo_plot_filtered


# Astrocyte Lineage Dot plot (Subset)
new_x <- c("GPC", "AS0", "AS1", "AS2")
astro_subset <- subset(fb_seurat_v20, idents = new_x)
Idents(fb_seurat_v20) <- factor(Idents(fb_seurat_v20), levels = new_x)
astro_subset$cell_type_v19 <- factor(
  x = astro_subset$cell_type_v19,
  levels = new_x
)
gradient <- c("blue", "white", "red")


astro_plot <- DotPlot_scCustom(astro_subset, features = c("MT2A", "VIM", "S100B", "APOE", "TIMP3", "SPARCL1", "SLC1A3", "SLC1A2", "NFIA", "GFAP"), 
                               group.by = "cell_type_v19", flip_axes = TRUE, remove_axis_titles = FALSE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

astro_plot  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1)) +
  scale_size_continuous(breaks = c(0, 20, 40, 60),
                        limits = c(0, 100)) +
  scale_color_gradientn(
    colors = gradient,
    breaks = c(-1, 0, 1)
  ) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)) 


# Astrocyte Lineage Dot plot (Entire Seurat Object)
astro_plot <- DotPlot_scCustom(fb_seurat_v25, features = c("MT2A", "VIM", "S100B", "APOE", "TIMP3", "SPARCL1", "SLC1A3", "SLC1A2", "NFIA", "GFAP"), 
                               group.by = "cell_type_v25", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

astro_plot_filtered <- astro_plot + 
  scale_y_discrete(limits = c("AS0", "AS1", "AS2"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
astro_plot_filtered


# Miscellaneous Dot Plot (Subset)
newx2 <- c("EC", "PC", "EB", "TRM", "EP")
misc_sub <- subset(fb_seurat_v20, idents = newx2)
Idents(misc_sub) <- factor(Idents(misc_sub), levels = newx2)
misc_sub$cell_type_v19 <- factor(
  x = misc_sub$cell_type_v19,
  levels = newx2
)
gradient <- c("blue", "white", "red")

misc_plot <- DotPlot_scCustom(misc_sub, features = c("MKI67", "PECAM1", "VWF", "CD34", "CD248", "MYOF", "ABCC9", "GJA4", "HBM", "HBG1", 
                                                     "HBG2", "HBF", "ALAS2", "GATA1", "KLF1", "HBA2", "CD3G", "CD3D", "CD69", "ICOS", 
                                                     "CD4", "CD8A", "FOXJ1", "SOX2", "CD133", "S100B", "TPPP3", "CCDC153", "GFAP", "CD163", 
                                                     "CD68", "AIF1", "TPI1", "P2RY12", "IRF8", "IL1B", "CD80", "CD86", "IL1", "ITGAM", 
                                                     "CD40", "IL6", "MAMUDRB1", "LDA"), 
                              group.by = "cell_type_v19", remove_axis_titles = FALSE, flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

misc_plot  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1)) +
  scale_size_continuous(breaks = c(0, 20, 40, 60),
                        limits = c(0, 100)) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


# Miscellaneous Dot Plot (Entire Seurat Object)
miscel_plot <- DotPlot_scCustom(fb_seurat_v25, 
                                features = c("MKI67", "PECAM1", "VWF", "CD34", "CD248", "MYOF", "ABCC9", "GJA4", "HBA", "HBM", 
                                             "HBG1", "HBG2", "HBF", "ALAS2", "ENSMMUG00000044429", "ENSMMUG00000041831", "HBD", "HBG1", "HBG2", "FECH", 
                                             "BLVRB", "SLC25A37", "GATA1", "KLF1", "ANK1", "HBA2",  "CD3G", "CD3D", "CD69", "ICOS", 
                                             "CD4", "CD8A", "FOXJ1", "SOX2", "CD133", "S100B", "TPPP3", "CCDC153", "GFAP", "CD68", 
                                             "CD163", "AIF1", "C1QB", "C1QC",  "SALL1", "CX3CR1", "P2RY12", "MRC1", "CSF1R",  "GPR34", 
                                             "CD80", "CD86", "IL1B"
                                             ), 
                                group.by = "cell_type_v25", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

miscel_plot_filtered <- miscel_plot + 
  scale_y_discrete(limits = c("EC", "PC", "EB0", "EB1", "CD4+T", "EP", "MG"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
miscel_plot_filtered


# Immune Dot Plot (Entire Seurat Object)
immune_plot <- DotPlot_scCustom(fb_seurat_v25, 
                                features = c("PECAM1", "VWF", "HBA2", "CD3G", "CD3D", "CD69", "ICOS", "CD4", "CD8A", "CD163", 
                                             "CD68", "AIF1", "TPI1", "P2RY12", "IRF8", "IL1B", "CD80", "CD86", "IL1", "ITGAM", 
                                             "CD40", "IL6", "MAMUDRB1", "LDA", "C1QB", "C1QC", "SALL1", "CX3CR1", "MRC1", "CSF1R", 
                                             "TMEM119", "LYVE1", "GPR34", "TREM2", "APOE", "CTSD"), 
                                group.by = "cell_type_v25", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

immune_plot_filtered <- immune_plot + 
  scale_y_discrete(limits = c("CD4+T", "Mono", "MG"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
immune_plot_filtered


# Microglia Dot plot (Subsetting)
micro_x <- c("Mono", "MG1", "MG2",  "MG3")
micro_sub <- subset(fb_seurat_v20, idents = micro_x)
Idents(micro_sub) <- factor(Idents(micro_sub), levels = micro_x)
micro_sub$cell_type_v19 <- factor(
  x = micro_sub$cell_type_v19,
  levels = micro_x
)
gradient <- c("blue", "white", "red")

micro_plot <- DotPlot_scCustom(micro_sub, 
                               features = c( "PTPRC", "CD14", "FCGR3", "CD68", "CD163", "AIF1", "C1QB", "C1QC", "SALL1", "CX3CR1", 
                                             "ITGAL", "MERTK", "STAB1", "P2RY12", "MRC1", "ITGAM", "CSF1R", "CCR2", "TMEM119", "LYVE1", 
                                             "GPR34", "CD80", "CD86", "IL1B", "IL6", "TREM2", "APOE", "CTSD", "ARG1", "CHI3L1"), 
                               group.by = "cell_type_v19", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

micro_plot  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1)) +
  scale_size_continuous(breaks = c(0, 20, 40, 60),
                        limits = c(0, 100)) +
  scale_color_gradientn(
    colors = gradient,
    breaks = c(-1, 0, 1)
  ) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


# Microglia Dot plot (Entire Seurat Object)
micro_plot <- DotPlot_scCustom(fb_seurat_v25, 
                               features = c("CD68", "CD163", "AIF1", "C1QB", "C1QC", "SALL1", "CX3CR1", "P2RY12", "MRC1", "CSF1R",  
                                            "GPR34", "CD80", "CD86", "IL1B", "CTSD"), 
                               group.by = "cell_type_v25", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

micro_plot_filtered <- micro_plot + 
  scale_y_discrete(limits = c("MG"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
micro_plot_filtered


# CGE & MGE Neural Lineage Dot plot (Subset)
nl_x <- c("vRG", "oRG", "CGE-P", "iCGE-IN", "iCGE-VIP/CR-IN", "MGE-P", "MGE-INP")
nl_sub <- subset(fb_seurat_v20, idents = nl_x)
Idents(nl_sub) <- factor(Idents(nl_sub), levels = nl_x)
nl_sub$cell_type_v19 <- factor(
  x = nl_sub$cell_type_v19,
  levels = nl_x
)
gradient <- c("blue", "white", "red")

nl_plot <- DotPlot_scCustom(nl_sub, 
                            features = c("DLX1", "DLX5", "GAD1", "GAD2", "NR2F1", "NR2F2", "CCK", "SP8", "NRP2", "PROX1", 
                                         "CALB2", "KCNC1", "CHRNA2", "VIP", "NKX2-1", "SOX2", "LHX6", "DLX2", "SOX6", "SP9", 
                                         "NPY", "ERBB4", "GRIK3", "KCNC2", "SST", "PVALB"), 
                            group.by = "cell_type_v19", flip_axes = TRUE, remove_axis_titles = FALSE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

nl_plot  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1)) +
  scale_size_continuous(breaks = c(0, 20, 40, 60),
                        limits = c(0, 100)) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


# Integrated vRG, oRG Genes to CGE & MGE Neural Lineage Dot plot (Subset)
nl_x <- c("vRG", "oRG", "CGE0", "CGE1", "CGE2", "MGE0", "MGE1")
nl_sub <- subset(fb_seurat_v20, idents = nl_x)
Idents(nl_sub) <- factor(Idents(nl_sub), levels = nl_x)
nl_sub$cell_type_v19 <- factor(
  x = nl_sub$cell_type_v19,
  levels = nl_x
)
gradient <- c("blue", "white", "red")

nl_plot <- DotPlot_scCustom(nl_sub, 
                            features = c("MKI67", "CDK1", "TOP2A", "PCNA", "ASPM", "ANLN", "HOPX", "FAM107A", "TNC", "HES5", 
                                         "NR2F1", "NR2F2", "DLX1", "DLX5", "GAD1", "GAD2", "TUBB3", "SP8", "NRP2", "PROX1", 
                                         "CALB2", "KCNC1", "CHRNA2", "VIP", "NKX2-1", "SOX2", "DLX2", "SOX6", "SP9", "NPY", 
                                         "LHX6",  "ERBB4", "GRIK3", "KCNC2", "SST", "PVALB"), 
                            group.by = "cell_type_v19", flip_axes = TRUE, remove_axis_titles = FALSE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

nl_plot  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1)) +
  scale_size_continuous(breaks = c(0, 20, 40, 60),
                        limits = c(0, 100)) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


# Integrated vRG, oRG Genes to CGE & MGE Neural Lineage Dot plot (Entire Seurat Object)
cge_mge_plot <- DotPlot_scCustom(fb_seurat_v25, 
                                 features = c("MKI67", "CDK1", "TOP2A", "PCNA", "ASPM", "ANLN", "HOPX", "FAM107A", "TNC", "HES5", 
                                              "NR2F1", "NR2F2", "DLX1", "DLX5", "GAD1", "GAD2", "TUBB3", "SP8", "NRP2", "PROX1", 
                                              "CALB2", "KCNC1", "CHRNA2", "VIP", "NKX2-1", "SOX2", "DLX2", "SOX6", "SP9", "NPY", 
                                              "LHX6",  "ERBB4", "GRIK3", "KCNC2", "SST", "PVALB"), 
                                 group.by = "cell_type_v25", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

cge_mge_plot_filtered <- cge_mge_plot + 
  scale_y_discrete(limits = c("vRG", "oRG", "CGE0", "CGE1", "CGE2", "MGE0", "MGE1"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
cge_mge_plot_filtered


###

# Astro Subcluster & UMAP
levels(fb_seurat_v17)
clustersToSubset <- c("vRG", "oRG", "GPC", "AS0", "AS1", "AS2")
data.sub <- subset(fb_seurat_v17, idents = clustersToSubset)

DimPlot_scCustom(data.sub, reduction = "umap.harmony")

# Run sctransform
data.sub <- SCTransform(object=data.sub, vst.flavor="v2",  conserve.memory=T, return.only.var.genes=F)

# Run PCA and UMAP
DefaultAssay(data.sub) <- "SCT"
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:20)

data.sub <- FindNeighbors(data.sub, dims = 1:20)
data.sub <- FindClusters(data.sub, resolution = 0.25)

p1 <- DimPlot_scCustom(data.sub)
p2 <- DimPlot_scCustom(data.sub, group.by = "cell_type_v16")
p1+p2

p2
# saveRDS(data.sub, file = "astro_sub_v2.RDS")


# Neuron Subcluster & UMAP
levels(fb_seurat_v17)
clustersToSubset <- c("CGE-P", "iCGE-IN", "iCGE-VIP/CR-IN", "MGE-P", "MGE-INP")
data.sub <- subset(fb_seurat_v17, idents = clustersToSubset)

DimPlot_scCustom(data.sub, reduction = "umap.harmony")

# Run sctransform
data.sub <- SCTransform(object=data.sub, vst.flavor="v2",  conserve.memory=T, return.only.var.genes=F)

# Run PCA and UMAP
DefaultAssay(data.sub) <- "SCT"
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:20)

data.sub <- FindNeighbors(data.sub, dims = 1:20)
data.sub <- FindClusters(data.sub, resolution = 0.25)

p1 <- DimPlot_scCustom(data.sub)
p2 <- DimPlot_scCustom(data.sub, group.by = "cell_type_v16")
p1+p2

p2


# Oligo Subcluster & UMAP
levels(fb_seurat_v17)
clustersToSubset <- c("vRG", "oRG", "GPC", "OPC0", "OPC1", "OPC2", "COP", "OL1", "OL2")
data.sub <- subset(fb_seurat_v17, idents = clustersToSubset)

DimPlot_scCustom(data.sub, reduction = "umap.harmony")

# Run sctransform
data.sub <- SCTransform(object=data.sub, vst.flavor="v2",  conserve.memory=T, return.only.var.genes=F)

# Run PCA and UMAP
DefaultAssay(data.sub) <- "SCT"
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:30)

data.sub <- FindNeighbors(data.sub, dims = 1:30)
data.sub <- FindClusters(data.sub, resolution = 0.25)

p1 <- DimPlot_scCustom(data.sub)
p2 <- DimPlot_scCustom(data.sub, group.by = "cell_type_v16")
p1+p2

p2


# Microglia Subcluster & UMAP
levels(fb_seurat_v20)
clustersToSubset <- c("Mono", "MG1", "MG2", "BAM")
data.sub <- subset(fb_seurat_v20, idents = clustersToSubset)

DimPlot_scCustom(data.sub, reduction = "umap.harmony")

# Run sctransform
data.sub <- SCTransform(object=data.sub, vst.flavor="v2",  conserve.memory=T, return.only.var.genes=F)

# Run PCA and UMAP
DefaultAssay(data.sub) <- "SCT"
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:30)

data.sub <- FindNeighbors(data.sub, dims = 1:30)
data.sub <- FindClusters(data.sub, resolution = 0.25)

p1 <- DimPlot_scCustom(data.sub)
p2 <- DimPlot_scCustom(data.sub, group.by = "cell_type_v16")
p1+p2

p2


# Misc Subcluster & UMAP
levels(fb_seurat_v17)
clustersToSubset <- c("EC", "PC", "EB", "CD3+T", "EP", "AIF1++ MG", "MG")
data.sub <- subset(fb_seurat_v17, idents = clustersToSubset)

DimPlot_scCustom(data.sub, reduction = "umap.harmony")

# Run sctransform
data.sub <- SCTransform(object=data.sub, vst.flavor="v2",  conserve.memory=T, return.only.var.genes=F)

# Run PCA and UMAP
DefaultAssay(data.sub) <- "SCT"
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:30)

data.sub <- FindNeighbors(data.sub, dims = 1:30)
data.sub <- FindClusters(data.sub, resolution = 0.25)


p1 <- DimPlot_scCustom(data.sub)
p2 <- DimPlot_scCustom(data.sub, group.by = "cell_type_v16")
p1+p2

p2


# Microglial Sub-Cluster UMAP
levels(micro_sub_v3)
clustersToSubset <- c("MG1", "MG2", "Mono", "BAM")
data.sub <- subset(micro_sub_v3, idents = clustersToSubset)

DimPlot_scCustom(data.sub, reduction = "umap.harmony")

# Run sctransform
data.sub <- SCTransform(object=data.sub, vst.flavor="v2",  conserve.memory=T, return.only.var.genes=F)

# Run PCA and UMAP
DefaultAssay(data.sub) <- "SCT"
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:30)

data.sub <- FindNeighbors(data.sub, dims = 1:30)
data.sub <- FindClusters(data.sub, resolution = 0.25)

p1 <- DimPlot_scCustom(data.sub)
p2 <- DimPlot_scCustom(data.sub, group.by = "cell_type_v18")
p1+p2

p2


# CGE and MGE Neural Progenitors UMAP
neuron_sub <- readRDS("/mmfs1/gscratch/kawaldorflab/jcorn427/fet_brain/neuron_sub_v6.RDS")
neuron_sub <- NormalizeData(neuron_sub, assay = "RNA")
DefaultAssay(neuron_sub) <- "RNA"

levels(neuron_sub)
Idents(neuron_sub) <- "cell_type_v18"
clustersToSubset <- c("vRG", "oRG", "CGE-P", "iCGE-IN", "iCGE-VIP/CR-IN", "MGE-INP", "MGE-P")
data.sub <- subset(neuron_sub, idents = clustersToSubset)

DimPlot_scCustom(data.sub, reduction = "umap.harmony")

# Run sctransform
data.sub <- SCTransform(object=data.sub, vst.flavor="v2",  conserve.memory=T, return.only.var.genes=F)

# Run PCA and UMAP
DefaultAssay(data.sub) <- "SCT"
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:30)

data.sub <- FindNeighbors(data.sub, dims = 1:30)
data.sub <- FindClusters(data.sub, resolution = 0.25)

p1 <- DimPlot_scCustom(data.sub)
p2 <- DimPlot_scCustom(data.sub, group.by = "cell_type_v18")
p1+p2

p2


# CGE and MGE Neural Progenitors UMAP (Entire Seurat Object)
levels(fb_seurat_v21)
Idents(fb_seurat_v21) <- "cell_type_v20"
clustersToSubset <- c("vRG", "oRG", "CGE0", "CGE1", "CGE2", "MGE0", "MGE1")
data.sub <- subset(fb_seurat_v21, idents = clustersToSubset)

DimPlot_scCustom(data.sub, reduction = "umap.harmony")

DefaultAssay(data.sub) <- "SCT"
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:30)

data.sub <- FindNeighbors(data.sub, dims = 1:30)
data.sub <- FindClusters(data.sub, resolution = 0.25)

p1 <- DimPlot_scCustom(data.sub)
p2 <- DimPlot_scCustom(data.sub, group.by = "cell_type_v20")
p1+p2

p2


# Highest expressed genes
neuron_sub <- PrepSCTFindMarkers(neuron_sub)

all_markers <- FindAllMarkers(object = neuron_sub) %>%
  Add_Pct_Diff()

top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7, named_vector = FALSE,
                                   make_unique = TRUE, rank_by = "avg_log2FC")
gradient <- c("blue", "white", "red")
plots <- DotPlot_scCustom(seurat_object = neuron_sub, features = top_markers, flip = TRUE, remove_axis_titles = FALSE, x_lab_rotate = TRUE, colors_use = gradient) + labs(x = "Gene", y = "Cluster")
#plots[[1]]
#plots <- Clustered_DotPlot(seurat_object = neuron_sub, features = top_markers, flip = TRUE, x_lab_rotate = 90, k = 20)

plots +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1)) +
  scale_size_continuous(breaks = c(0, 20, 40, 60),
                        limits = c(0, 100)) +
  scale_color_gradientn(
    colors = gradient,
    breaks = c(-1, 0, 1)
  ) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)) 

top_20_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 20, data_frame = TRUE, rank_by = "avg_log2FC")
top_20_markers[top_20_markers$cluster == "iCGE-VIP/CR-IN",]


# Oligodendrocyte + Astrocyte Lineage
levels(oligo_astro)
# "GPC", "OPC0", "OPC1", "OPC2", "COP", "OL1", "OL2", "AS0", "AS1", "AS2"
clustersToSubset <- c()
data.sub <- subset(oligo_astro, idents = clustersToSubset)

DimPlot_scCustom(data.sub, reduction = "umap.harmony")

# Run sctransform
data.sub <- SCTransform(object=data.sub, vst.flavor="v2",  conserve.memory=T, return.only.var.genes=F)

# Run PCA and UMAP
DefaultAssay(data.sub) <- "SCT"
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:40)

data.sub <- FindNeighbors(data.sub, dims = 1:40)
data.sub <- FindClusters(data.sub, resolution = 0.25)

p1 <- DimPlot_scCustom(data.sub)
p2 <- DimPlot_scCustom(data.sub, group.by = "cell_type_v18")
p1+p2

p2

oligo_astro <- readRDS("/mmfs1/gscratch/kawaldorflab/jcorn427/fet_brain/Oligo_Astro_sub.RDS")
oligo_astro <- NormalizeData(oligo_astro, assay = "RNA")
DefaultAssay(oligo_astro) <- "RNA"


# Split Oligo and Astro into 2 UMAPS with GPC in both
levels(oligo_astro)
# "GPC", "OPC0", "OPC1", "OPC2", "COP", "OL1", "OL2", "AS0", "AS1", "AS2"
Idents(oligo_astro) <- "cell_type_v18"
clustersToSubset <- c("GPC", "AS0", "AS1", "AS2")
data.sub <- subset(oligo_astro, idents = clustersToSubset)

DimPlot_scCustom(data.sub, reduction = "umap.harmony")

# Run sctransform
data.sub <- SCTransform(object=data.sub, vst.flavor="v2",  conserve.memory=T, return.only.var.genes=F)

# Run PCA and UMAP
DefaultAssay(data.sub) <- "SCT"
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:30)

data.sub <- FindNeighbors(data.sub, dims = 1:30)
data.sub <- FindClusters(data.sub, resolution = 0.25)

p1 <- DimPlot_scCustom(data.sub)
p2 <- DimPlot_scCustom(data.sub, group.by = "cell_type_v18")
p1+p2

p2


###

# CGE & MGE
neuro_sub <- readRDS("/mmfs1/gscratch/kawaldorflab/jcorn427/fet_brain/neuro_sub.RDS")
neuro_sub <- NormalizeData(neuro_sub, assay = "RNA")
DefaultAssay(neuro_sub) <- "RNA"

# Oligodendrocyte
oligo_sub <- readRDS("/mmfs1/gscratch/kawaldorflab/jcorn427/fet_brain/oligo_sub.RDS")
oligo_sub <- NormalizeData(oligo_sub, assay = "RNA")
DefaultAssay(oligo_sub) <- "RNA"

# Astrocyte
astro_sub <- readRDS("/mmfs1/gscratch/kawaldorflab/jcorn427/fet_brain/astro_sub.RDS")
astro_sub <- NormalizeData(astro_sub, assay = "RNA")
DefaultAssay(astro_sub) <- "RNA"

# Microglia
micro_sub <- readRDS("/mmfs1/gscratch/kawaldorflab/jcorn427/fet_brain/micro_sub.RDS")
micro_sub <- NormalizeData(micro_sub, assay = "RNA")
DefaultAssay(micro_sub) <- "RNA"

###

### CGE MGE B
# Plot genes of interest
Idents(neuro_sub) <- "cell_type_updated"
subset_to_plot <- subset(neuro_sub, idents = "other", invert = TRUE) # replace neur_sub with subset of interest
umap_coords <- as.data.frame(subset_to_plot@reductions$umap.harmony@cell.embeddings)

# Ensure the genes exist in your object before proceeding
genes_to_plot <- c("MKI67", "NR2F2", "NRP2", "KCNC2") # replace with your genes of interest
if (!all(genes_to_plot %in% rownames(fb_seurat_v25))) {
  stop("One or more genes not found in the Seurat object.")
}

# Fetch expression data. `FetchData` is the safest way to get expression.
expression_data <- FetchData(object = subset_to_plot, vars = genes_to_plot)

# Data for the background layer (all cells)
background_data <- umap_coords

# Prepare data for the shuffled foreground layers (expressing cells only)
plot_data <- cbind(umap_coords, expression_data)

long_data <- plot_data %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  # Filter for cells that express the gene (expression > 0)
  filter(expression > 0)

# Randomly shuffle the rows of the combined data frame
shuffled_data <- long_data %>%
  sample_frac(1)

ggplot() +
  # 1. Plot all cells as a light grey background layer
  geom_point(data = background_data,
             aes(x = umapharmony_1, y = umapharmony_2),
             color = "lightgrey",
             size = 0.5,
             alpha = 0.5) +
  
  # 2. Plot the shuffled, expressing cells on top
  geom_point(data = shuffled_data,
             aes(x = umapharmony_1, y = umapharmony_2, color = gene),
             alpha = 0.5,
             size = 0.8) +
  
  # Set manual colors and customize the theme, also, remember to replace genes with those you are interested in
  scale_color_manual(values = c(
    "MKI67" = "purple",
    "NR2F2" = "orange",
    "NRP2" = "blue",
    "KCNC2" = "green"
  )) +
  # Use guides() to override the size of the points in the legend
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "MKI67, NR2F2, NRP2, KCNC2",
       color = "Gene Expression") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Remove grid lines
    axis.line = element_line(color = "black"), # Add axis lines back in
    # Increase legend title and text size
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)
  )


### CGE MGE C
# Plot genes of interest
Idents(neuro_sub) <- "cell_type_updated"
subset_to_plot <- subset(neuro_sub, idents = "other", invert = TRUE) # replace neur_sub with subset of interest
umap_coords <- as.data.frame(subset_to_plot@reductions$umap.harmony@cell.embeddings)

# Ensure the genes exist in your object before proceeding
genes_to_plot <- c("HOPX", "NKX2-1", "LHX6") # replace with your genes of interest
if (!all(genes_to_plot %in% rownames(fb_seurat_v25))) {
  stop("One or more genes not found in the Seurat object.")
}

# Fetch expression data. `FetchData` is the safest way to get expression.
expression_data <- FetchData(object = subset_to_plot, vars = genes_to_plot)

# Data for the background layer (all cells)
background_data <- umap_coords

# Prepare data for the shuffled foreground layers (expressing cells only)
plot_data <- cbind(umap_coords, expression_data)

long_data <- plot_data %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  # Filter for cells that express the gene (expression > 0)
  filter(expression > 0)

# Randomly shuffle the rows of the combined data frame
shuffled_data <- long_data %>%
  sample_frac(1)

ggplot() +
  # 1. Plot all cells as a light grey background layer
  geom_point(data = background_data,
             aes(x = umapharmony_1, y = umapharmony_2),
             color = "lightgrey",
             size = 0.5,
             alpha = 0.5) +
  
  # 2. Plot the shuffled, expressing cells on top
  geom_point(data = shuffled_data,
             aes(x = umapharmony_1, y = umapharmony_2, color = gene),
             alpha = 0.5,
             size = 0.8) +
  
  # Set manual colors and customize the theme, also, remember to replace genes with those you are interested in
  scale_color_manual(values = c(
    "HOPX" = "blue",
    "NKX2-1" = "magenta",
    "LHX6" = "green"
  )) +
  # Use guides() to override the size of the points in the legend
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "HOPX, NKX2-1, LHX6",
       color = "Gene Expression") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Remove grid lines
    axis.line = element_line(color = "black"), # Add axis lines back in
    # Increase legend title and text size
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)
  )


### Oligo B
# Plot genes of interest
Idents(oligo_sub) <- "cell_type_updated"
subset_to_plot <- subset(oligo_sub, idents = "other", invert = TRUE) # replace neur_sub with subset of interest
umap_coords <- as.data.frame(subset_to_plot@reductions$umap.harmony@cell.embeddings)

# Ensure the genes exist in your object before proceeding
genes_to_plot <- c("PDGFRA", "MAL", "CA10", "ENPP6") # replace with your genes of interest
if (!all(genes_to_plot %in% rownames(fb_seurat_v25))) {
  stop("One or more genes not found in the Seurat object.")
}

# Fetch expression data. `FetchData` is the safest way to get expression.
expression_data <- FetchData(object = subset_to_plot, vars = genes_to_plot)

# Data for the background layer (all cells)
background_data <- umap_coords

# Prepare data for the shuffled foreground layers (expressing cells only)
plot_data <- cbind(umap_coords, expression_data)

long_data <- plot_data %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  # Filter for cells that express the gene (expression > 0)
  filter(expression > 0)

# Randomly shuffle the rows of the combined data frame
shuffled_data <- long_data %>%
  sample_frac(1)

ggplot() +
  # 1. Plot all cells as a light grey background layer
  geom_point(data = background_data,
             aes(x = umapharmony_1, y = umapharmony_2),
             color = "lightgrey",
             size = 0.5,
             alpha = 0.5) +
  
  # 2. Plot the shuffled, expressing cells on top
  geom_point(data = shuffled_data,
             aes(x = umapharmony_1, y = umapharmony_2, color = gene),
             alpha = 0.5,
             size = 0.8) +
  
  # Set manual colors and customize the theme, also, remember to replace genes with those you are interested in
  scale_color_manual(values = c(
    "PDGFRA" = "blue",
    "MAL" = "red",
    "CA10" = "orange"
  )) +
  # Use guides() to override the size of the points in the legend
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "PDGFRA, MAL, CA10, ENPP6",
       color = "Gene Expression") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Remove grid lines
    axis.line = element_line(color = "black"), # Add axis lines back in
    # Increase legend title and text size
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)
  )


### Oligo C
# Plot genes of interest
Idents(oligo_sub) <- "cell_type_updated"
subset_to_plot <- subset(oligo_sub, idents = "other", invert = TRUE) # replace neur_sub with subset of interest
umap_coords <- as.data.frame(subset_to_plot@reductions$umap.harmony@cell.embeddings)

# Ensure the genes exist in your object before proceeding
genes_to_plot <- c("SOX10", "PCDH15", "BCAS1", "PLP1") # replace with your genes of interest
if (!all(genes_to_plot %in% rownames(fb_seurat_v25))) {
  stop("One or more genes not found in the Seurat object.")
}

# Fetch expression data. `FetchData` is the safest way to get expression.
expression_data <- FetchData(object = subset_to_plot, vars = genes_to_plot)

# Data for the background layer (all cells)
background_data <- umap_coords

# Prepare data for the shuffled foreground layers (expressing cells only)
plot_data <- cbind(umap_coords, expression_data)

long_data <- plot_data %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  # Filter for cells that express the gene (expression > 0)
  filter(expression > 0)

# Randomly shuffle the rows of the combined data frame
shuffled_data <- long_data %>%
  sample_frac(1)

ggplot() +
  # 1. Plot all cells as a light grey background layer
  geom_point(data = background_data,
             aes(x = umapharmony_1, y = umapharmony_2),
             color = "lightgrey",
             size = 0.5,
             alpha = 0.5) +
  
  # 2. Plot the shuffled, expressing cells on top
  geom_point(data = shuffled_data,
             aes(x = umapharmony_1, y = umapharmony_2, color = gene),
             alpha = 0.5,
             size = 0.8) +
  
  # Set manual colors and customize the theme, also, remember to replace genes with those you are interested in
  scale_color_manual(values = c(
    "SOX10" = "purple",
    "PCDH15" = "orange",
    "BCAS1" = "blue",
    "PLP1" = "green"
  )) +
  # Use guides() to override the size of the points in the legend
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "SOX10, PCDH15, BCAS1, PLP1",
       color = "Gene Expression") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Remove grid lines
    axis.line = element_line(color = "black"), # Add axis lines back in
    # Increase legend title and text size
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)
  )


### Oligo
# Plot genes of interest
Idents(oligo_sub) <- "cell_type_updated"
subset_to_plot <- subset(oligo_sub, idents = "other", invert = TRUE) # replace neur_sub with subset of interest
umap_coords <- as.data.frame(subset_to_plot@reductions$umap.harmony@cell.embeddings)

# Ensure the genes exist in your object before proceeding
genes_to_plot <- c("ENPP6", "MBP", "PCDH17") # replace with your genes of interest
if (!all(genes_to_plot %in% rownames(fb_seurat_v25))) {
  stop("One or more genes not found in the Seurat object.")
}

# Fetch expression data. `FetchData` is the safest way to get expression.
expression_data <- FetchData(object = subset_to_plot, vars = genes_to_plot)

# Data for the background layer (all cells)
background_data <- umap_coords

# Prepare data for the shuffled foreground layers (expressing cells only)
plot_data <- cbind(umap_coords, expression_data)

long_data <- plot_data %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  # Filter for cells that express the gene (expression > 0)
  filter(expression > 0)

# Randomly shuffle the rows of the combined data frame
shuffled_data <- long_data %>%
  sample_frac(1)

ggplot() +
  # 1. Plot all cells as a light grey background layer
  geom_point(data = background_data,
             aes(x = umapharmony_1, y = umapharmony_2),
             color = "lightgrey",
             size = 0.5,
             alpha = 0.5) +
  
  # 2. Plot the shuffled, expressing cells on top
  geom_point(data = shuffled_data,
             aes(x = umapharmony_1, y = umapharmony_2, color = gene),
             alpha = 0.5,
             size = 0.8) +
  
  # Set manual colors and customize the theme, also, remember to replace genes with those you are interested in
  scale_color_manual(values = c(
    "ENPP6" = "red",
    "MBP" = "orange",
    "PCDH17" = "blue"
  )) +
  # Use guides() to override the size of the points in the legend
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "ENPP6, MBP, PCDH17",
       color = "Gene Expression") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Remove grid lines
    axis.line = element_line(color = "black"), # Add axis lines back in
    # Increase legend title and text size
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)
  )


### Astro B
# Plot genes of interest
Idents(astro_sub) <- "cell_type_updated"
subset_to_plot <- subset(astro_sub, idents = "other", invert = TRUE) # replace neur_sub with subset of interest
umap_coords <- as.data.frame(subset_to_plot@reductions$umap.harmony@cell.embeddings)

# Ensure the genes exist in your object before proceeding
genes_to_plot <- c("S100B", "GFAP") # replace with your genes of interest
if (!all(genes_to_plot %in% rownames(fb_seurat_v25))) {
  stop("One or more genes not found in the Seurat object.")
}

# Fetch expression data. `FetchData` is the safest way to get expression.
expression_data <- FetchData(object = subset_to_plot, vars = genes_to_plot)

# Data for the background layer (all cells)
background_data <- umap_coords

# Prepare data for the shuffled foreground layers (expressing cells only)
plot_data <- cbind(umap_coords, expression_data)

long_data <- plot_data %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  # Filter for cells that express the gene (expression > 0)
  filter(expression > 0)

# Randomly shuffle the rows of the combined data frame
shuffled_data <- long_data %>%
  sample_frac(1)

ggplot() +
  # 1. Plot all cells as a light grey background layer
  geom_point(data = background_data,
             aes(x = umapharmony_1, y = umapharmony_2),
             color = "lightgrey",
             size = 0.5,
             alpha = 0.5) +
  
  # 2. Plot the shuffled, expressing cells on top
  geom_point(data = shuffled_data,
             aes(x = umapharmony_1, y = umapharmony_2, color = gene),
             alpha = 0.5,
             size = 0.8) +
  
  # Set manual colors and customize the theme, also, remember to replace genes with those you are interested in
  scale_color_manual(values = c(
    "GFAP" = "orange",
    "S100B" = "blue"
  )) +
  # Use guides() to override the size of the points in the legend
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "S100B, GFAP",
       color = "Gene Expression") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Remove grid lines
    axis.line = element_line(color = "black"), # Add axis lines back in
    # Increase legend title and text size
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)
  )


### Astro C
# Plot genes of interest
Idents(astro_sub) <- "cell_type_updated"
subset_to_plot <- subset(astro_sub, idents = "other", invert = TRUE) # replace neur_sub with subset of interest
umap_coords <- as.data.frame(subset_to_plot@reductions$umap.harmony@cell.embeddings)

# Ensure the genes exist in your object before proceeding
genes_to_plot <- c("S100B", "APOE", "TIMP3", "SPARCL1", "NFIA") # replace with your genes of interest
if (!all(genes_to_plot %in% rownames(fb_seurat_v25))) {
  stop("One or more genes not found in the Seurat object.")
}

# Fetch expression data. `FetchData` is the safest way to get expression.
expression_data <- FetchData(object = subset_to_plot, vars = genes_to_plot)

# Data for the background layer (all cells)
background_data <- umap_coords

# Prepare data for the shuffled foreground layers (expressing cells only)
plot_data <- cbind(umap_coords, expression_data)

long_data <- plot_data %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  # Filter for cells that express the gene (expression > 0)
  filter(expression > 0)

# Randomly shuffle the rows of the combined data frame
shuffled_data <- long_data %>%
  sample_frac(1)

ggplot() +
  # 1. Plot all cells as a light grey background layer
  geom_point(data = background_data,
             aes(x = umapharmony_1, y = umapharmony_2),
             color = "lightgrey",
             size = 0.5,
             alpha = 0.5) +
  
  # 2. Plot the shuffled, expressing cells on top
  geom_point(data = shuffled_data,
             aes(x = umapharmony_1, y = umapharmony_2, color = gene),
             alpha = 0.5,
             size = 0.8) +
  
  # Set manual colors and customize the theme, also, remember to replace genes with those you are interested in
  scale_color_manual(values = c(
    "S100B" = "purple",
    "APOE" = "orange",
    "TIMP3" = "blue",
    "SPARCL1" = "green",
    "NFIA" = "red"
  )) +
  # Use guides() to override the size of the points in the legend
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "S100B, APOE, TIMP3, SPARCL1, NFIA",
       color = "Gene Expression") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Remove grid lines
    axis.line = element_line(color = "black"), # Add axis lines back in
    # Increase legend title and text size
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)
  )


### Micro B
# Plot genes of interest
Idents(micro_sub) <- "cell_type_updated"
subset_to_plot <- subset(micro_sub, idents = "other", invert = TRUE) # replace neur_sub with subset of interest
umap_coords <- as.data.frame(subset_to_plot@reductions$umap.harmony@cell.embeddings)

# Ensure the genes exist in your object before proceeding
genes_to_plot <- c("CD86", "SALL1") # replace with your genes of interest
if (!all(genes_to_plot %in% rownames(fb_seurat_v25))) {
  stop("One or more genes not found in the Seurat object.")
}

# Fetch expression data. `FetchData` is the safest way to get expression.
expression_data <- FetchData(object = subset_to_plot, vars = genes_to_plot)

# Data for the background layer (all cells)
background_data <- umap_coords

# Prepare data for the shuffled foreground layers (expressing cells only)
plot_data <- cbind(umap_coords, expression_data)

long_data <- plot_data %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  # Filter for cells that express the gene (expression > 0)
  filter(expression > 0)

# Randomly shuffle the rows of the combined data frame
shuffled_data <- long_data %>%
  sample_frac(1)

ggplot() +
  # 1. Plot all cells as a light grey background layer
  geom_point(data = background_data,
             aes(x = umapharmony_1, y = umapharmony_2),
             color = "lightgrey",
             size = 0.5,
             alpha = 0.5) +
  
  # 2. Plot the shuffled, expressing cells on top
  geom_point(data = shuffled_data,
             aes(x = umapharmony_1, y = umapharmony_2, color = gene),
             alpha = 0.5,
             size = 0.8) +
  
  # Set manual colors and customize the theme, also, remember to replace genes with those you are interested in
  scale_color_manual(values = c(
    "CD86" = "blue",
    "SALL1" = "orange"
  )) +
  # Use guides() to override the size of the points in the legend
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "CD86, SALL1",
       color = "Gene Expression") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Remove grid lines
    axis.line = element_line(color = "black"), # Add axis lines back in
    # Increase legend title and text size
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)
  )


### Micro C
# Plot genes of interest
Idents(micro_sub) <- "cell_type_updated"
subset_to_plot <- subset(micro_sub, idents = "other", invert = TRUE) # replace neur_sub with subset of interest
umap_coords <- as.data.frame(subset_to_plot@reductions$umap.harmony@cell.embeddings)

# Ensure the genes exist in your object before proceeding
genes_to_plot <- c("CX3CR1", "CD86") # replace with your genes of interest
if (!all(genes_to_plot %in% rownames(fb_seurat_v25))) {
  stop("One or more genes not found in the Seurat object.")
}

# Fetch expression data. `FetchData` is the safest way to get expression.
expression_data <- FetchData(object = subset_to_plot, vars = genes_to_plot)

# Data for the background layer (all cells)
background_data <- umap_coords

# Prepare data for the shuffled foreground layers (expressing cells only)
plot_data <- cbind(umap_coords, expression_data)

long_data <- plot_data %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  # Filter for cells that express the gene (expression > 0)
  filter(expression > 0)

# Randomly shuffle the rows of the combined data frame
shuffled_data <- long_data %>%
  sample_frac(1)

ggplot() +
  # 1. Plot all cells as a light grey background layer
  geom_point(data = background_data,
             aes(x = umapharmony_1, y = umapharmony_2),
             color = "lightgrey",
             size = 0.5,
             alpha = 0.5) +
  
  # 2. Plot the shuffled, expressing cells on top
  geom_point(data = shuffled_data,
             aes(x = umapharmony_1, y = umapharmony_2, color = gene),
             alpha = 0.5,
             size = 0.8) +
  
  # Set manual colors and customize the theme, also, remember to replace genes with those you are interested in
  scale_color_manual(values = c(
    "CX3CR1" = "blue",
    "CD86" = "orange"
  )) +
  # Use guides() to override the size of the points in the legend
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "CX3CR1, CD86",
       color = "Gene Expression") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Remove grid lines
    axis.line = element_line(color = "black"), # Add axis lines back in
    # Increase legend title and text size
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)
  )



### ExN3 ###
FeaturePlot(fb_seurat_pre, features = "nFeature_RNA", 
            reduction = "umap", pt.size = 1) + 
  ggtitle("Number of Detected Genes per Cell") 


FeaturePlot(fb_seurat_pre, features = "nCount_RNA", 
            reduction = "umap", pt.size = 1) +
  ggtitle("Total UMIs per Cell") +
  scale_color_gradient(
    low = "blue",
    high = "red",
    limits = c(0, 10000)
  )

FeatureScatter(fb_seurat_v25, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(fb_seurat_v25, features = "nFeature_RNA", group.by = "cell_type_v25", pt.size = 0.1)


meta %>%
  group_by(cell_type_v25) %>%
  summarize(mean_nFeature = mean(nFeature_RNA),
            median_nFeature = median(nFeature_RNA),
            n_cells = n())


# Wilcoxon rank-sum test
meta <- fb_seurat_v25@meta.data
wilcox.test(nFeature_RNA ~ cell_type_v25, data = subset(meta, cell_type_v25 %in% c("ExN2", "ExN3")))
meta$group <- ifelse(meta$cell_type_v25 == "ExN3", "ExN3", "Other")
wilcox.test(nCount_RNA ~ group, data = meta)



all_markers <- FindAllMarkers(object = fb_seurat_v25) %>%
  Add_Pct_Diff()

top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7, named_vector = FALSE,
                                   make_unique = TRUE, rank_by = "avg_log2FC")

plots <- Clustered_DotPlot(seurat_object = fb_seurat_v25, features = top_markers, flip = TRUE, x_lab_rotate = 90)
plots[[1]]
plots <- Clustered_DotPlot(seurat_object = fb_seurat_v25, features = top_markers, flip = TRUE, x_lab_rotate = 90, k = 20)

top_20_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 20, data_frame = TRUE, rank_by = "avg_log2FC")
write.csv(top_20_markers, file = "top_20_markers.csv")



########## pre-scDEED ##########
fb_seurat_pre <- readRDS("/mmfs1/gscratch/kawaldorflab/jcorn427/fet_brain/fb_03/fb_seurat_pre_scDEED.RDS")
fb_seurat_pre <- NormalizeData(fb_seurat_pre, assay = "RNA")
DefaultAssay(fb_seurat_pre) <- "RNA"


# CGE-MGE 
cge_mge <- DotPlot_scCustom(fb_seurat_pre, 
                            features = c("MKI67", "CDK1", "TOP2A", "PCNA", "ASPM", "ANLN", "HOPX", "FAM107A", "TNC", "HES5", 
                                         "NR2F1", "NR2F2", "DLX1", "DLX5", "GAD1", "GAD2", "TUBB3", "SP8", "NRP2", "PROX1", 
                                         "CALB2", "KCNC1", "CHRNA2", "VIP", "NKX2-1", "SOX2", "DLX2", "SOX6", "SP9", "NPY", 
                                         "LHX6",  "ERBB4", "GRIK3", "KCNC2", "SST", "PVALB"), 
                                 group.by = "seurat_clusters", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")
cge_mge


# Microglia
micro <- DotPlot_scCustom(fb_seurat_pre, 
                          features = c("CD68", "CD163", "AIF1", "C1QB", "C1QC", "SALL1", "CX3CR1", "P2RY12", "MRC1", "CSF1R", 
                                       "GPR34", "CD80", "CD86", "IL1B", "CTSD", "TMEM119", "TREM2", "CST7", "LPL", "SPP1", 
                                       "GPNMB", "CD74", "STAT1", "IRF1", "MKI67", "TOP2A"), 
                               group.by = "seurat_clusters", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

micro_filt <- micro +
  scale_y_discrete(limits = c("20"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
micro_filt


# Miscellaneous
miscel <- DotPlot_scCustom(fb_seurat_pre, 
                           features = c("MKI67", "PECAM1", "VWF", "CD34", "CD248", "MYOF", "ABCC9", "GJA4", "ENSMMUG00000041831", "HBD", 
                                        "HBG1", "HBG2", "FECH", "BLVRB", "SLC25A37", "GATA1", "KLF1", "ANK1", "HBA2", "CD3G", 
                                        "CD3D", "CD69", "ICOS", "CD4", "CD8A", "FOXJ1", "SOX2", "CD133", "S100B", "TPPP3", 
                                        "CCDC153", "GFAP", "CD68", "CD163", "AIF1", "C1QB", "C1QC", "SALL1", "CX3CR1", "P2RY12", 
                                        "MRC1", "CSF1R", "GPR34", "CD80", "CD86", "IL1B", "COL1A1", "COL1A2", "COL3A1", "DCN", 
                                        "LUM", "FMOD", "FBLN1", "FN1", "VIM"),
                           group.by = "seurat_clusters", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")
miscel


# Astro
astro <- DotPlot_scCustom(fb_seurat_pre, 
                          features = c("MT2A", "VIM", "S100B", "APOE", "TIMP3", "SPARCL1", "SLC1A3", "SLC1A2", "NFIA", "GFAP"), 
                               group.by = "seurat_clusters", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

astro


# Oligo
oligo <- DotPlot_scCustom(fb_seurat_pre, 
                          features = c("SOX10", "PDGFRA", "OLIG2", "PCDH17", "PCDH15", "MMP16", "CA10", "ENSMMUG00000056728", "BCAS1", "ENPP6", 
                                       "MAL", "MOG", "PLP1", "MBP"), 
                               group.by = "seurat_clusters", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")
oligo


# IPC
ipc <- DotPlot_scCustom(fb_seurat_pre,
                        features = c("DCX", "NEUROD1", "EOMES", "TBR1", "SOX5", "POU3F2", "NR2F1", "NRP1", "CRYM", "SLC17A6", 
                                     "NEUROD2", "TLE4", "SATB2", "UNC5D", "GRIA2", "CUX1", "CUX2", "RORB", "GRIN2B", "SLC17A7", 
                                     "CAMK2A", "NRGN", "SATB1"), 
                        group.by = "seurat_clusters", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")
ipc


###
levels(fb_seurat_pre)
clustersToSubset <- c("20")
data.sub <- subset(fb_seurat_pre, idents = clustersToSubset)

DimPlot_scCustom(data.sub, reduction = "umap")

#Run sctransform
data.sub <- SCTransform(object=data.sub, vst.flavor="v2",  conserve.memory=T, return.only.var.genes=F)

# Run PCA and UMAP
DefaultAssay(data.sub) <- "SCT"
data.sub <- RunPCA(data.sub)
ElbowPlot(data.sub, ndims = 50)
data.sub <- RunUMAP(data.sub, dims = 1:20)

data.sub <- FindNeighbors(data.sub, dims = 1:20)
data.sub <- FindClusters(data.sub, resolution = 0.25)

p1 <- DimPlot_scCustom(data.sub)
p2 <- DimPlot_scCustom(data.sub, group.by = "seurat_clusters")
p1+p2

p2



########## Final Seurat Object ##########
fb_mrg<- readRDS("/mmfs1/gscratch/kawaldorflab/jcorn427/fet_brain/fb_03/fb_merged.RDS")
fb_mrg <- NormalizeData(fb_mrg, assay = "RNA")
DefaultAssay(fb_mrg) <- "RNA"

head(fb_mrg)
unique(fb_mrg$cell_type_v3)
table(fb_mrg$cell_type_v3)

# Convert to character first
fb_mrg$cell_type_v3 <- as.character(fb_mrg$cell_type_v3)

# Merge MG1 and MG2
fb_mrg$cell_type_v3[fb_mrg$cell_type_v3 %in% c("MG1", "MG2")] <- "MG1"

# (Optional) Convert back to factor if you prefer
fb_mrg$cell_type_v3 <- factor(fb_mrg$cell_type_v3)

# Check results
unique(fb_mrg$cell_type_v3)
table(fb_mrg$cell_type_v3)

saveRDS(fb_mrg, file = "/mmfs1/gscratch/kawaldorflab/tykim/R_files/fetal_brain/fb_merged.RDS")

fb_mrg<- readRDS("/mmfs1/gscratch/kawaldorflab/tykim/R_files/fetal_brain/fb_merged.RDS")
fb_mrg <- NormalizeData(fb_mrg, assay = "RNA")
DefaultAssay(fb_mrg) <- "RNA"

Idents(fb_mrg) <- fb_mrg$cell_type_v3


# CGE + MGE Dot Plot
lge_dot <- DotPlot_scCustom(fb_mrg, 
                            features = c("MEIS2", "FOXP1", "FOXP2", "ASCL1", "DLL1", "DLL3", "EGFR", "DLX1", "DLX5", "GAD1", 
                                         "GAD2", "PROX1", "CALB2", "KCNC1", "CHRNA2", "VIP", "NKX2-1", "SOX2", "LHX6", "DLX2", 
                                         "SOX6", "SP9", "NPY", "ERBB4", "GRIK3", "KCNC2", "SST", "PVALB"), 
                          group.by = "cell_type_v3", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

lge_dot_filtered <- lge_dot +
  scale_y_discrete(limits = c("CGE0", "CGE1", "MGE"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
lge_dot_filtered


# IPC Dot Plot
ipc_dot <- DotPlot_scCustom(fb_mrg, features = c("EOMES", "TBR1"), 
                            group.by = "cell_type_v3", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

ipc_dot_filtered <- ipc_dot +
  scale_y_discrete(limits = c("IPC", "UP/I"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
ipc_dot_filtered


# ExN Dot Plot
exn_dot <- DotPlot_scCustom(fb_mrg, 
                            features = c("DCX", "NEUROD1", "EOMES", "TBR1", "SOX5", "POU3F2", "NR2F1", "SLC17A6", "NEUROD2", "SATB2", 
                                         "UNC5D", "GRIA2", "CUX1", "CUX2", "GRIN2B", "SLC17A7", "CAMK2A", "NRGN", "NKX2-1", "LHX6", 
                                         "SOX6", "DLX1", "DLX5", "DLX6", "GAD1", "GAD2", "PROX1", "CALB2", "HTR3A", "NRSN1", 
                                         "POU2F2", "RELN" ,"VIP", "CCK", "SST", "NPY", "PVALB"), 
                             group.by = "cell_type_v3", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

exn_dot_filtered <- exn_dot +
  scale_y_discrete(limits = c("IPC", "ExN0","ExN1", "CGE0", "CGE1", "MGE"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
exn_dot_filtered


# T Dot Plot
t_dot <- DotPlot_scCustom(fb_mrg, features = c("CD3G", "CD3D", "CD3E", "MKI67", "PECAM1", "CD69", "ICOS", "CD4", "CD8A"), 
                          group.by = "cell_type_v3", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

t_dot_filtered <- t_dot +
  scale_y_discrete(limits = c("T", "UP/I"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
t_dot_filtered


# B Dot Plot
b_dot <- DotPlot_scCustom(fb_mrg, features = c("CD19", "MS4A1", "IGHM"), 
                             group.by = "cell_type_v3", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

b_dot_filtered <- b_dot +
  scale_y_discrete(limits = c("B", "UP/I"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
b_dot_filtered


# ChP Dot Plot
ChP_dot <- DotPlot_scCustom(fb_mrg, features = c("TTR", "CLDN1", "CLDN3", "CLDN4",  "CLDN5"), 
                          group.by = "cell_type_v3", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

ChP_dot_filtered <- ChP_dot +
  scale_y_discrete(limits = c("ChP", "UP/I"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
ChP_dot_filtered


# Miscel Dot
fib_dot <- DotPlot_scCustom(fb_mrg, 
                            features = c("COL1A2", "VIM", "DCN", "FMOD", "PDGFRB", "MYOF", "TAGLN", "ABCC9", "PECAM1", "VWF", 
                                         "CD34", "FLT1", "SOX2", "CD133", "CCDC153", "TTR", "CLDN3", "CLDN4", "FOXJ1", "RFX2", 
                                         "RFX3", "DNALI1", "TUBA4B", "CETN2", "AQP4"), 
                            group.by = "cell_type_v3", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

fib_dot_filtered <- fib_dot +
  scale_y_discrete(limits = c("FIB", "PC", "EC", "ChP", "Ep"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
fib_dot_filtered


# EC Dot Plot
ec_dot <- DotPlot_scCustom(fb_mrg, 
                           features = c("PECAM1", "VWF", "CD34", "MKI67", "CD248", "MYOF", "ABCC9", "GJA4", "HBA", "HBM", 
                                        "HBG1", "HBG2", "HBF", "ALAS2"), 
                            group.by = "cell_type_v3", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

ec_dot_filtered <- ec_dot +
  scale_y_discrete(limits = c("EC", "UP/I"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
ec_dot_filtered


# Astro + vRG +oRG
vrg_dot <- DotPlot_scCustom(fb_mrg, 
                            features = c("MKI67", "TOP2A", "PCNA", "ASPM", "CENPF", "CDK1", "ANLN", "HOPX", "FAM107A", "TNC",
                                         "MT2A", "VIM", "S100B", "APOE", "TIMP3", "SPARCL1", "SLC1A3", "SLC1A2", "NFIA", "GFAP"), 
                            group.by = "cell_type_v3", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

vrg_dot_filtered <- vrg_dot +
  scale_y_discrete(limits = c("vRG", "UP/I", "oRG-G", "AS0", "AS1", "AS2"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
vrg_dot_filtered


# Immune Cells
mg_dot <- DotPlot_scCustom(fb_mrg, 
                           features = c("CD3G", "CD3D", "ICOS", "CD4", "CD68", "LYZ", "FCN1", "AIF1", "IL1B", "CD86", 
                                        "CD14", "SPI1", "CCR2", "CD80", "P2RY12", "CX3CR1", "TREM2", "SALL1", "MRC1", "CSF1R", 
                                        "HBA2"), 
                            group.by = "cell_type_v3", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

mg_dot_filtered <- mg_dot +
  scale_y_discrete(limits = c("T", "Mono", "MG0", "MG1", "BAM"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
mg_dot_filtered


# PC Dot Plot
pc_dot <- DotPlot_scCustom(fb_mrg, 
                           features = c("PDGFRB", "ACTA2", "MYH11", "CNN1", "TAGLN", "ABCC9", "MKI67", "PECAM1", "VWF", "CD34", 
                                        "CD248", "MYOF", "GJA4", "HBA", "TPPP3"), 
                            group.by = "cell_type_v3", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

pc_dot_filtered <- pc_dot +
  scale_y_discrete(limits = c("PC", "UP/I"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
pc_dot_filtered


# AS1 + AS2 Dot Plot
as_plot <- DotPlot_scCustom(fb_mrg, 
                            features = c("MT2A", "VIM", "S100B", "APOE", "TIMP3", "SPARCL1", "SLC1A3", "SLC1A2", "NFIA", "GFAP"), 
                               group.by = "cell_type_v3", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

as_plot_filtered <- as_plot + 
  scale_y_discrete(limits = c("AS0", "AS1", "AS2"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
as_plot_filtered


# Oligo
oligo_plot <- DotPlot_scCustom(fb_mrg, 
                               features = c("SOX10", "PDGFRA", "OLIG2", "PCDH17", "PCDH15", "MMP16", "CA10", "ENSMMUG00000056728", "BCAS1", "ENPP6", 
                                            "MAL", "MOG", "PLP1", "MBP"), 
                               group.by = "cell_type_v3", remove_axis_titles = FALSE,  flip_axes = TRUE, x_lab_rotate = TRUE, colors_use = gradient) + 
  labs(x = "Gene", y = "Cluster")

oligo_plot_filtered <- oligo_plot + 
  scale_y_discrete(limits = c("UP/I", "OPC", "COP1", "COP2", "OL1", "OL2"))  +
  guides(size = guide_legend(title = "Percent expressed", order = 2),
         color = guide_colorbar(title = "Avg. Expression (scaled)", order = 1))
oligo_plot_filtered



### Top Markers
fb_seurat_pre <- PrepSCTFindMarkers(fb_seurat_pre)

fb_seurat_pre <- JoinLayers(fb_seurat_pre)

all_markers <- FindAllMarkers(object = fb_seurat_pre) %>%
  Add_Pct_Diff()

top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7, named_vector = FALSE,
                                   make_unique = TRUE, rank_by = "avg_log2FC")

plots <- Clustered_DotPlot(seurat_object = fb_seurat_pre, features = top_markers, flip = TRUE, x_lab_rotate = 90)
plots[[1]]
plots <- Clustered_DotPlot(seurat_object = fb_seurat_pre, features = top_markers, flip = TRUE, x_lab_rotate = 90, k = 20)



########## UMAP with Opacity ########## 
umap_coords <- as.data.frame(fb_mrg@reductions$umap@cell.embeddings)
genes_to_plot <- c("VIM", "HOPX", "FABP7")

if (!all(genes_to_plot %in% rownames(fb_mrg))) {
  stop("One or more genes not found in the Seurat object.")
}

expression_data <- FetchData(object = fb_mrg, vars = genes_to_plot)
background_data <- umap_coords

plot_data <- cbind(umap_coords, expression_data)

long_data <- plot_data %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  filter(expression > 0)

shuffled_data <- long_data %>%
  sample_frac(1)

ggplot() +
  # 1. Plot all cells as a light grey background layer
  geom_point(data = background_data,
             aes(x = umap_1, y = umap_2),
             color = "lightgrey",
             size = 0.5,
             alpha = 0.5) +
  
  # 2. Plot the shuffled, expressing cells on top
  geom_point(data = shuffled_data,
             aes(x = umap_1, y = umap_2, color = gene, alpha = expression),
             size = 0.8) +
  
  # Set manual colors for the 'color' aesthetic (gene identity)
  scale_color_manual(values = c(
                                "VIM" = "orange",
                                "HOPX" = "blue", 
                                "FABP7" = "green")) +
  
  # Control the alpha scale. Use guide_legend() to show alpha variations correctly.
  scale_alpha_continuous(
    range = c(0.1, 0.8),
    name = "Expression Level (Opacity)",
    # Change from guide_colorbar() to guide_legend()
    guide = guide_legend(title.position = "top", title.hjust = 0.5) 
  ) + 
  
  # Use guides() to override the *color* legend's appearance
  guides(
    color = guide_legend(title = "Gene Identity", override.aes = list(size = 4, alpha = 0.8)),
    # Use guides() to further customize the alpha legend if needed
    alpha = guide_legend(override.aes = list(size = 4, color = "black"))
  ) +
  labs(title = "TNC") +
  theme_minimal()
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
  
FeaturePlot(fb_mrg, features = c("EGFR", "VIM", "HOPX", "FABP7", "ID2", "ID3", "AQP4", "SLC1A3"))


FeatureScatter(fb_mrg, feature1 = "EGFR", feature2 = "ID2")

FeaturePlot(
  fb_mrg,
  features = c("EGFR", "ID2"),
  blend = TRUE,
  pt.size = 0.5
)


FeatureScatter(fb_mrg, feature1 = "EGFR", feature2 = "AQP4")

FeaturePlot(
  fb_mrg,
  features = c("EGFR", "AQP4"),
  blend = TRUE,
  pt.size = 0.5
)


FeatureScatter(fb_mrg, feature1 = "EGFR", feature2 = "GFAP")

FeaturePlot(
  fb_mrg,
  features = c("EGFR", "GFAP"),
  blend = TRUE,
  pt.size = 0.5
)


FeatureScatter(fb_mrg, feature1 = "FABP7", feature2 = "GFAP")

FeaturePlot(
  fb_mrg,
  features = c("FABP7", "GFAP"),
  blend = TRUE,
  pt.size = 0.5
)


##### Final UMAP #####
DimPlot_scCustom(fb_mrg, reduction = "umap", label = TRUE, label.box = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, 
                 colors_use = DiscretePalette_scCustomize(num_colors = 29, palette = "varibow", shuffle_pal = FALSE), color_seed = 5)


##### UMAP with Lasso & CellSelector #####
Idents(fb_mrg) <- "cell_type"
levels(fb_mrg)
clustersToSubset <- c("AS0", "AS1", "AS2", "vRG", "UP/I", "oRG-G")
data.sub <- subset(fb_mrg, subset = cell_type %in% clustersToSubset)

DimPlot_scCustom(data.sub, reduction = "umap", label.box = TRUE)


# Generate the UMAP plot
plot <- DimPlot(object = data.sub, reduction = "umap")

# Pass the plot to CellSelector()
# The R console will prompt you to select points on the plot
cells.located <- CellSelector(plot = plot)
# cells.located2 <- CellSelector(plot = plot)
# Create a new identity class for the selected cells
data.sub <- SetIdent(data.sub, cells = cells.located, value = 'selected_cells')
# data.sub <- SetIdent(data.sub, cells = cells.located2, value = 'selected_cells2')

data.sub <- RenameIdents(
  object = data.sub,
  'selected_cells' = "CombinedGroup"
)
# Verify the selection with a new plot
DimPlot(data.sub, reduction = "umap", group.by = 'ident')

data.sub$lasso <- Idents(data.sub)

saveRDS(data.sub, 'data.sub.RDS')

combined_group_subset <- subset(data.sub, idents = "CombinedGroup")

DimPlot(data.sub, reduction = "umap", group.by = 'lasso')

# Correctly create the new metadata column using as.character()
data.sub$cell_type_updated <- ifelse(
  test = data.sub$lasso == "CombinedGroup",
  yes = as.character(data.sub$cell_type),
  no = "other"
)

data.sub$cell_type_updated <- as.factor(data.sub$cell_type_updated)

# Check the new assignments
table(data.sub$cell_type_updated)

# Now, plot with the new, correctly labeled column
DimPlot_scCustom(data.sub, group.by = "cell_type_updated", reduction = "umap")

Idents(data.sub) <- "cell_type_updated"
subset_to_plot <- subset(data.sub, idents = "other", invert = TRUE)

DimPlot_scCustom(subset_to_plot, group.by = "cell_type_updated", reduction = "umap", label.box = TRUE, label = TRUE)

DimPlot_scCustom(subset_to_plot, reduction = "umap", label = TRUE, label.box = TRUE, pt.size = 0.5, repel = TRUE, shuffle = TRUE, 
                 colors_use = DiscretePalette_scCustomize(num_colors = 29, palette = "polychrome", shuffle_pal = TRUE), color_seed = 5)
