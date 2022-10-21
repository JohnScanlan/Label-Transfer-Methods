library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(SingleR)
library(dplyr)
library(SingleR)
library(scRNAseq)
library(ExperimentHub)
library(scuttle)

colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")

epi <- LoadH5Seurat("C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\MiceWeightLossCohort\\MOUSE RNA-seq\\preprocessed_epi.h5seurat")

plots <- FeaturePlot(epi, c('Cd4', 'Cd8a'), pt.size=1, 
                     blend=T, label = TRUE, combine = TRUE, 
                     cols = c("lightgrey", "red", "blue"))

plots[[3]]

FeaturePlot(epi, c('Zbtb32'), pt.size = 1)

DimPlot(epi, label = T)

plots
head(epi)

#DotPlots
features <- c('Cd8a','Cd4', 'Cd19', 'Ccr7', 'Il7r', 'Cd3d', 'Cd3g','Itgam', 'Ncr1', 'Cd44', 'Tcf7')
DotPlot(epi, features = features) + RotatedAxis() #All Cells

top1 <- c('S100a9', 'Cd79a', 'Lyz2', 'Igkc', 'Ccl5', 'Gzma', 'Iglc1', 'Trdc', 'Ighd', 'Tcf7', 'Cxcr6', 'Cd209a', 'Cd4',
          'Siglech', 'Xcl1', 'Aif1', 'Stmn1', 'Cxcr2', 'Camp', 'Ace', 'Apol7c')
DotPlot(epi, features = top1) + RotatedAxis() #All Cells

t_markers <- c('Cd8a', 'Cd4', 'Ccr7', 'Tcf7')
DotPlot(epi, features = t_markers) + RotatedAxis() #All Cells

b <- c('Zbtb32')
DotPlot(epi, features = features) + RotatedAxis() #B Cells

#Finding Markers
epi.markers <- FindAllMarkers(epi)
View(epi.markers)

#Labelling
new.cluster.ids <- c('NK', '1', 'yd_T', '3', '4', '5', '6', '7', '8', 
                     'Bs', '10', '11', '12', '13', '14','15', '16', '17', '18', '19')
names(new.cluster.ids) <- levels(epi)
epi <- RenameIdents(epi, new.cluster.ids)

DimPlot(epi,reduction = "umap", pt.size = 1, label.size = 5, label=T) + NoLegend()

#comparing cluster groups
m.markers <- FindMarkers(epi, ident.1 = '17', ident.2 = c('18'), min.pct = 0.25)
View(m.markers)

#Method 1 - Label Transfer
eh <- ExperimentHub()

query(eh, 'TabulaMurisData')

eh[['EH1617']]

epi_ref <- eh[['EH1617']]
epi_ref <- epi_ref[,!is.na(epi_ref$cell_ontology_class)]


unique(epi_ref$cell_ontology_class)


epi_ref <- logNormCounts(epi_ref)

results <- SingleR(test = as.SingleCellExperiment(epi), ref = epi_ref, labels = epi_ref$cell_ontology_class)

epi$singlr_labels

epi$singlr_labels <- results$labels
DimPlot(epi, reduction = 'umap', group.by = 'single_r', label = T, label.size = 5) + NoLegend()

#Method 2 - Label Transfer
ref <- celldex::MouseRNAseqData()

unique(ref$label.fine)#more fine labels
unique(ref$label.main)#more general labels

results <- SingleR(test = as.SingleCellExperiment(epi), ref = ref, labels = ref$label.fine)
epi$celldex <- results$labels

DimPlot(epi, reduction = 'umap', group.by = 'celldex', label = T, label.size = 5) + NoLegend()
