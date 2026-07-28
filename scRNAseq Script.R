#This script contains the code used to analyze mRNA from the MOE of two mouse lines: WT and OR10G4
#Seurat is used to process the count matrices and harmony is used to integrate the two datasets.
#Numerous blocks of code are exploratory or do not produce final publication figures.
#Using R V4.4.0

#Section I: Full Data set Prep work----

#Packages to install and load; some packages are on BioConductor and require BiocManager::install("package") to install
library(Seurat) #Base package
library(SeuratObject) #Required/loaded by Seurat
library(sctransform) #An alternative normalization method 
library(tidyverse) #Loads a bunch of tidy packages
library(magrittr) #pipe
library(harmony) #Integrate datasets
library(Rcpp) #Required/loaded by harmony
library(clustree) #Allows analysis of clustering resolution
library(patchwork) #Glue some figures together
library(ggrepel) #Some figures use this with text labels
library(ggpubr) #Mostly used for compare_means() function
library(scDblFinder) #For doublets
library(remotes) #Used to install a non-CRAN or BioC package (on Github)
#remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk) #Extension of base Seurat functions to interface with HDF5 files; Unused 
library(ggridges) #To make a ridgeplot
library(reticulate) #Part of GEP cNMF work

BiocManager::install("slingshot") #Block for Pseudotime Analysis
suppressPackageStartupMessages({
  library(plotly)
  options(rgl.printRglwidget = TRUE)
  library(Matrix)
  library(sparseMatrixStats)
  library(slingshot)
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("tradeSeq")
  library(tradeSeq)
  library(patchwork)
})

#Useful filter function
`%notin%` <- Negate(`%in%`)

#Ouput Folder for Figures and related
if (!dir.exists("Output2")) {dir.create("Output2")}

#Loading the 10X datasets. Make sure to set the correct directory.
WT.data <- Read10X(data.dir = "../scRNAseq Figures/Final_Approach/WT_raw_feature_bc_matrix/")
OR10G4.data <- Read10X(data.dir = "../scRNAseq Figures/Final_Approach/OR10G4_raw_feature_bc_matrix/")

#Initialize the Seurat objects with the raw (non-normalized data). 
#Individual Olfrs are expressed in fewer than 1% of OSNs, so we accept a one cell minimum.
#A 1500-feature minimum per cell plus additional filters based on Counts and Mitochondrial counts (Prefiltered results not displayed)
WT <- CreateSeuratObject(counts = WT.data, project = "Feinstein-PF-12681", min.cells = 1, min.features = 1500) 
WT <- PercentageFeatureSet(WT, pattern = "^mt-", col.name = "percent.mt")
WT <- subset(WT, subset = percent.mt <= 10 & nCount_RNA <= 35000)
WT[["Genotype"]] <- "WT"
OR10G4 <- CreateSeuratObject(counts = OR10G4.data, project = "Feinstein-PF-12681", min.cells = 1, min.features = 1500) 
OR10G4 <- PercentageFeatureSet(OR10G4, pattern = "^mt-", col.name = "percent.mt")
OR10G4 <- subset(OR10G4, subset = percent.mt <= 10 & nCount_RNA <= 35000)
OR10G4[["Genotype"]] <- "OR10G4"

#Merging both Seurat Objects. 
#Due to the potential impact of Olfrs and OR10G4 to impact later analyses, those genes will later be extracted, stored, and removed from Merge.
Merge <- merge(WT, y = c(OR10G4), add.cell.ids = c("WT", "OR10G4"), project = "Feinstein-PF-12681")
Merge <- PercentageFeatureSet(Merge, features = "10G4", col.name = "percent.10G4")
Merge <- PercentageFeatureSet(Merge, pattern = "^Olfr", col.name = "percent.OR")

#A Metadata plot, plus TG and ORs
#01_Merge_Metadata_10G4_Olfr
VlnPlot(Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.10G4", "percent.OR"), 
        ncol = 5, group.by = "Genotype") &
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5))
ggsave("Output2/01_Merge_Metadata_10G4_Olfr.png", width = 25, height = 15, units = "cm")

#SCTransform, this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures(). 
#Regressing out percent.mt feature. Also Feature and Count metadata due to potential impact of sequencing depth
Merge <- SCTransform(Merge, vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA"), verbose = FALSE)

#In order to check for doublets, PCs are first generated.
set.seed(123456)
Merge <- RunPCA(Merge, verbose = FALSE)

#Doublet Analysis followed by rerunning SCTransform and PCA and continuing to Harmony and clustering
#Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(Merge, assay = "SCT")
sceR <- scDblFinder(sce)

#Add results back to Seurat object
Merge$doublet_class <- colData(sceR)$scDblFinder.class
Merge$doublet_score <- colData(sceR)$scDblFinder.score

#Due to random nature of analysis, repeating to check for consistent labeling
sceR <- scDblFinder(sce)
Merge$doublet_class2 <- colData(sceR)$scDblFinder.class
Merge$doublet_score2 <- colData(sceR)$scDblFinder.score
sceR <- scDblFinder(sce)
Merge$doublet_class3 <- colData(sceR)$scDblFinder.class
Merge$doublet_score3 <- colData(sceR)$scDblFinder.score

Merge$doubletFinal <- ifelse((Merge$doublet_class == "doublet") + (Merge$doublet_class2 == "doublet") + (Merge$doublet_class3 == "doublet") > 1,
                             "doublet", "singlet")
#02_Merge_Doublets
DimPlot(Merge, split.by = "doubletFinal", group.by = "Genotype") + 
  ggtitle("PC Plot Showing Doublets. Final Consensus: ≥2 of 3 doublet tests")

ggsave("Output2/02_Merge_Doublets.png", width = 20, height = 15, units = "cm")

#In this iteration of the doublets script, 105 cells were doublets in at least 2/3 doublet tests
SingletCells <- Cells(subset(Merge, subset = doubletFinal == "singlet"))

#Repeating all code from Merge creation to PCA, keeping only singlets; Merge2 to keep Merge object intact for review
Merge2 <- merge(WT, y = c(OR10G4), add.cell.ids = c("WT", "OR10G4"), project = "Feinstein-PF-12681")
Merge2 <- subset(Merge2, cells = SingletCells)
Merge2 <- PercentageFeatureSet(Merge2, features = "10G4", col.name = "percent.10G4")
Merge2 <- PercentageFeatureSet(Merge2, pattern = "^Olfr", col.name = "percent.OR")
Merge2 <- SCTransform(Merge2, vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA"), verbose = FALSE)
set.seed(123456)
Merge2 <- RunPCA(Merge2, verbose = FALSE)

#Running RunHarmony followed by standard Clustering process using the "harmony" reduction as opposed to the usual "pca". 
#Then, a comparison of how the points are better mixed post-harmony.
H_Merge <- RunHarmony(Merge2, 
                      group.by.vars = c("Genotype"), 
                      reduction = "pca", 
                      assay.use = "SCT", 
                      reduction.save = "harmony")

#03_H_Merge_BeforeAfter_RunHarmony
p1 <- DimPlot(H_Merge, reduction = "pca", group.by = "Genotype", pt.size = 0.5) + 
  NoLegend() + 
  labs(title = "Before RunHarmony()") +
  theme(axis.title = element_blank())
p2 <- DimPlot(H_Merge, reduction = "harmony", group.by = "Genotype", pt.size = 0.5) + 
  NoLegend() + 
  labs(title = "After RunHarmony()") +
  theme(axis.title = element_blank())

p1 + p2
ggsave("Output2/03_H_Merge_BeforeAfter_RunHarmony.png", width = 25, height = 20, units = "cm")

#04_H_Merge_ElbowPlot_selectDims
ElbowPlot(H_Merge, ndims = 50)
ggsave("Output2/04_H_Merge_ElbowPlot_selectDims.png", width = 25, height = 18, units = "cm")

H_Merge <- RunUMAP(H_Merge, reduction = "harmony", assay = "SCT", dims = 1:12) 
H_Merge <- FindNeighbors(object = H_Merge, reduction = "harmony")
H_Merge <- FindClusters(H_Merge, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))

#04a_H_Merge_Dimsplit_Genotype
DimPlot(H_Merge, split.by = "Genotype", group.by = "Genotype") & NoLegend() & NoAxes() 
ggsave("Output2/04a_H_Merge_Dimsplit_Genotype.png", width = 30, height = 30, units = "cm")

#Reviewing the cluster changes over resolution and making a selection
#05_H_Merge_Clustree
clustree(H_Merge, prefix = "SCT_snn_res.")
ggsave("Output2/05_H_Merge_Clustree_SF3E.png", width = 30, height = 20, units = "cm")

#06_H_Merge_ResolutionCompare_UMAP
Idents(H_Merge) <- "SCT_snn_res.0.2"
p1 <- DimPlot(H_Merge, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.0.2")
Idents(H_Merge) <- "SCT_snn_res.0.4"
p2 <- DimPlot(H_Merge, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.0.4")
Idents(H_Merge) <- "SCT_snn_res.0.6"
p3 <- DimPlot(H_Merge, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.0.6")
Idents(H_Merge) <- "SCT_snn_res.0.8"
p4 <- DimPlot(H_Merge, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.0.8")
Idents(H_Merge) <- "SCT_snn_res.1"
p5 <- DimPlot(H_Merge, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.1")
Idents(H_Merge) <- "SCT_snn_res.1.2"
p6 <- DimPlot(H_Merge, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.1.2")
Idents(H_Merge) <- "SCT_snn_res.1.4"
p7 <- DimPlot(H_Merge, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.1.4")

(p1 | p2 | p3 | p4 | p5 | p6 | p7) & NoAxes()
ggsave("Output2/06_H_Merge_ResolutionCompare_UMAP_SF3D.png", width = 50, height = 10, units = "cm")

#Final Decision about resolution: 1.4 resolution: 27 clusters. Default might already be this resolution, but to be sure...
Idents(H_Merge) <- "SCT_snn_res.1.4"

#Clean up unnecessary objects to keep the session running smoothly. Use at your own discretion
rm(OR10G4, OR10G4.data, WT, WT.data, sce, sceR, Merge, Merge2)

#Creating metadata for mitotic cells
s.genes <- str_to_title(cc.genes.updated.2019$s.genes)
g2m.genes <- str_to_title(cc.genes.updated.2019$g2m.genes)
H_Merge <- CellCycleScoring(H_Merge, s.features = s.genes, g2m.features = g2m.genes)

#Creating log1p(UMI/Total UMI * 10000) transformation for all Genes in each cell for the RNA assay to allow for appropriate visualization
DefaultAssay(H_Merge) <- "RNA"

H_Merge <- NormalizeData(H_Merge)
H_Merge <- FindVariableFeatures(H_Merge)
H_Merge <- ScaleData(H_Merge, vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA"))

DefaultAssay(H_Merge) <- "SCT"

#Some overview of the metadata qualities of each cluster. 10G4 cells clearly have less genes and counts.
#07_H_Merge_Metadata_TG_Olfr_byCluster
VlnPlot(H_Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.10G4", "percent.OR"), ncol = 5, split.by = "Genotype")
ggsave("Output2/07_H_Merge_Metadata_TG_Olfr_byCluster.png", width = 50, height = 10, units = "cm")


#07a_H_Merge_HouseKeeping_Genotype_Cluster Are housekeeping genes similar between genotypes?
VlnPlot(H_Merge, features = c("Actb", "Gapdh", "Ywhaz", "Ubc", "Ppia"), split.by = "Genotype", ncol = 5, assay = "RNA")
ggsave("Output2/07a_H_Merge_HouseKeeping_Genotype_Cluster.png", width = 50, height = 10, units = "cm")

#The next step is to Identify what each Cluster represents before selecting OSN Lineage cells and rerunning the entire pipeline from the top, basically.

#Section II: Cluster Labeling----

#Running a strict DEG function to speed this thing up a little. Plotted Top Markers with DotPlot
#Top5 is used in further cluster exploration, not explicitly shown here.
H_Merge <- PrepSCTFindMarkers(H_Merge) #Required for FindMarkers functions with current setup.
MarkersV1 <- FindAllMarkers(H_Merge, logfc.threshold = 0.5, min.pct = 0.4, only.pos = TRUE)
Top5 <- MarkersV1 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

##Looking at the Markers of OSN Lineage, from GBC to mOSN
#08_H_Merge_LineageOSN_Markers_VlnUMAP
p1 <- VlnPlot(H_Merge, features = c("G2M.Score","S.Score", "Ascl1", "Neurog1", "Neurod1", "Tex15", "Lhx2", 
                                    "Myt1l", "Gnas","Gap43", "Omp", "Gnal", "Adcy3", "Stoml3"), 
              ncol = 7, assay = "RNA") & 
  NoAxes() &
  theme(axis.text.x = element_text(size = 6, angle = 0, vjust = 0)) 
p2 <- FeaturePlot(H_Merge, features = c("G2M.Score","S.Score", "Ascl1", "Neurog1", "Neurod1", "Tex15", "Lhx2", 
                                        "Myt1l", "Gnas","Gap43", "Omp", "Gnal", "Adcy3", "Stoml3"), ncol = 7) & 
  NoAxes() & NoLegend()

p1/p2
ggsave("Output2/08_H_Merge_LineageOSN_Markers_VlnUMAP_SF4A.png", width = 70, height = 40, units = "cm")

#Based on above, previous exploration/knowledge, and published figs:
#Cluster 17: GBC/INP 
#Cluster 1: iOSN1/iOSN2
#Cluster 2/19: iOSN3
#Cluster 0/3/4/5/8/10/11: mOSN

##Looking at other Epithelial cells
#09_H_Merge_Epithelial_Markers_VlnUMAP
p1 <- VlnPlot(H_Merge, features = c("Trp63", "Krt5", "Krt14", "Ascl3", "Cftr", "Sox9",
                                    "Trpm5", "Cyp2g1", "Cyp1a2", "Aqp5"), 
              ncol = 5, assay = "RNA") & 
  NoAxes() &
  theme(axis.text.x = element_text(size = 6, angle = 0, vjust = 0)) 
p2 <- FeaturePlot(H_Merge, features = c("Trp63", "Krt5", "Krt14", "Ascl3", "Cftr", "Sox9",
                                        "Trpm5", "Cyp2g1", "Cyp1a2", "Aqp5"), 
                  ncol = 5) & 
  NoAxes() & NoLegend()

p1/p2
ggsave("Output2/09_H_Merge_Epithelial_Markers_VlnUMAP_SF4B.png", width = 70, height = 40, units = "cm")

#My best interpretation is the following: 
#Cluster 16: MV_Brush
#Cluster 21: MV_Ionocyte 
#Cluster 12: Bowman's Gland Cells (Aqp5+)
#Cluster 12/20: Sustentacular
#Cluster 13: HBC

##Looking at Hematopoietic cells
#10_H_Merge_HematopoeticMarkers_VlnUMAP
p1 <- VlnPlot(H_Merge, features = c("G2M.Score", "Igkc", "Il7r", "Nkg7", "Ccr2", "S100a4", 
                                    "Cd14", "C1qa", "Ms4a7", "Ngp", "Cxcr2", "Gypa"), ncol = 6, 
              assay = "RNA") & 
  NoAxes() &
  theme(axis.text.x = element_text(size = 6, angle = 0, vjust = 0))  
p2 <- FeaturePlot(H_Merge, features = c("G2M.Score", "Igkc", "Il7r", "Nkg7", "Ccr2", "S100a4", 
                                        "Cd14", "C1qa", "Ms4a7", "Ngp", "Cxcr2", "Gypa"), ncol = 6) & 
  NoAxes() & NoLegend()

p1/p2
ggsave("Output2/10_H_Merge_HematopoeticMarkers_VlnUMAP_SF4C.png", width = 70, height = 40, units = "cm")

#Based on my understanding and lots of research: 
#cluster 14: B cells/T cells/NK cells
#Cluster 7: S100a4 Myeloid 
#Cluster 15: Ms4a7 Myeloid
#Cluster 23: Pre-Neutrophil
#Cluster 18: Maturing Neutrophils
#Cluster 22: Aged Mature Neutrophils
#Cluster 24: RBC

##Looking at connective tissue support cells
#11_H_Merge_LaminaSupportMarkers_VlnUMAP
p1 <- VlnPlot(H_Merge, features = c("S100b", "Acta2", "Cnn1", "Col1a1", "Col1a2", "Vtn", "Thy1", "Eng", "Pecam1", "Cd34"), ncol = 5,
              assay = "RNA") & 
  NoAxes() &
  theme(axis.text.x = element_text(size = 6, angle = 0, vjust = 0)) 
p2 <- FeaturePlot(H_Merge, features = c("S100b", "Acta2", "Cnn1", "Col1a1", "Col1a2", "Vtn", "Thy1", "Eng", "Pecam1", "Cd34"), ncol = 5) & 
  NoAxes() & NoLegend()

p1/p2
ggsave("Output2/11_H_Merge_LaminaSupportMarkers_VlnUMAP_SF4D.png", width = 70, height = 40, units = "cm")

#The remaining 2 clusters are a mixtures of OECs, Pericytes, Fibroblasts, Smooth muscle, Stromal cells
#Cluster 6: Endothelial/Stromal
#Cluster 9: OEC/Smooth Muscle/Myofibroblast/Fibroblasts/Pericytes/Stromal

#I will now assign cluster labels to the harmonized seurat, saving the original numeric clusters
H_Merge[["Backup_Cluster_Labels"]] <- Idents(object = H_Merge)
new.cluster.ids <- c("mOSN", "iOSN1/2", "iOSN3", "mOSN", "mOSN", "mOSN",
                     "Endo/Str", "Myeloid_S100a4", "mOSN", "OEC/SM/Fibro/Peri", "mOSN", 
                     "mOSN", "BG/Sus", "HBC", "B/T/NK Cell", "Myeloid_Ms4a7", "MV_Brush",
                     "GBC/INP", "iNeutrophil", "iOSN3", "Sus", "MV_Ionocyte",
                     "mNeutrophil",  "Pre-Neutrophil", "RBC")
names(new.cluster.ids) <- levels(H_Merge)
H_Merge <- RenameIdents(H_Merge, new.cluster.ids)

#12_H_Merge_Labeled_UMAP
DimPlot(H_Merge, reduction = "umap", label = TRUE, pt.size = 0.3, label.box = TRUE, 
        repel = TRUE, label.size = 2.5) + theme_bw() + NoLegend() + NoAxes()

ggsave("Output2/12_H_Merge_Labeled_UMAP_MF3A.png", width = 30, height = 30, units = "cm")

#14_H_Merge_Labeled_Vln_Markers
levels(H_Merge) <- rev(c("GBC/INP", "iOSN1/2", "iOSN3", "mOSN",
                         "HBC", "Sus", "BG/Sus", "MV_Brush", "MV_Ionocyte", "OEC/SM/Fibro/Peri", "Endo/Str",
                         "Myeloid_S100a4", "Myeloid_Ms4a7", "B/T/NK Cell", "Pre-Neutrophil", 
                         "iNeutrophil", "mNeutrophil", "RBC"))

p1 <- VlnPlot(H_Merge, features = c("Neurod1", "Lhx2", "Gap43", "Omp", "Adcy3",
                                              "Krt5", "Aqp5", "Cyp2g1", "Trpm5", "Ascl3", "Acta2", "Col1a1", "Pecam1",
                                              "S100a4", "Ms4a7", "Igkc", "Il7r", "Nkg7", "Trem1", "Ngp", "Gypa"), 
              ncol = 21, pt.size = 0) & 
  coord_flip() & 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

p2 <- VlnPlot(H_Merge, features = "Ascl1", ncol = 1, pt.size = 0) + 
  coord_flip() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

p2 + p1 + plot_layout(widths = c(0.2, 6.0), ncol = 2)
ggsave("Output2/14_H_Merge_Labeled_Vln_Markers_SF4F.png", width = 70, height = 25, units = "cm")

#Redoing the DEG table with labeled clusters and slightly different selection
#15_H_Merge_Labeled_DotPlot_Top5v2DEG_SF4E
MarkersV2 <- FindAllMarkers(H_Merge, logfc.threshold = 0.5, min.pct = 0.4, only.pos = TRUE)

Top5V2 <- MarkersV2 %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(order_by = avg_log2FC, n = 5)

DotPlot(H_Merge, features = unique(Top5V2$gene)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, vjust = 0.5)) + 
  theme(axis.title = element_blank())

ggsave("Output2/15_H_Merge_Labeled_DotPlot_Top5v2DEG_SF4E.png", width = 70, height = 20, units = "cm")

#13_H_Merge_Labeled_Vln_OR10G4_Mbnl1
p1 <- VlnPlot(H_Merge, features = "Omp", assay = "RNA", split.by = "Genotype") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")
p2 <- VlnPlot(H_Merge, features = "10G4", assay = "RNA", split.by = "Genotype") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")
p3 <- VlnPlot(H_Merge, features = "Mbnl1", assay = "RNA", split.by = "Genotype") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

p1/p2/p3
ggsave("Output2/13_H_Merge_Labeled_Vln_OR10G4_Mbnl1_SF3F.png", width = 40, height = 20, units = "cm")

#Section III: OSN Lineage Analysis----
#Generating an OSN Lineage subset: subset Merge data and rerun setup prior to additional analysis
OSNL <- subset(H_Merge, idents = c("GBC/INP", "iOSN1/2", "iOSN3", "mOSN"))

#Keep only genes that are expressed (count ≥ 1) in at least one cell. Does this work? Not for Mature_OSNs
DefaultAssay(OSNL) <- "RNA"
OSNL <- OSNL[rowSums(OSNL[["RNA"]]$counts > 0) > 0, ]

#Manually rebuilding the Seurat object with some metadata and only RNA counts
idents.df = data.frame("orig.ident" = OSNL$orig.ident, "nCount_RNA" = OSNL$nCount_RNA, 
                       "nFeature_RNA" = OSNL$nFeature_RNA, "percent.mt" = OSNL$percent.mt,
                       "Genotype" = OSNL$Genotype, "percent.10G4" = OSNL$percent.10G4,
                       "percent.OR" = OSNL$percent.OR, "S.Score" = OSNL$S.Score, 
                       "G2M.Score" = OSNL$G2M.Score, "Phase" = OSNL$Phase)

OSNL = CreateSeuratObject(counts = LayerData(OSNL, assay = "RNA", layer = "counts"), project = "Export", meta.data = idents.df)

Idents(OSNL) = "Genotype"

#Saving the original count data for Olfrs and 10G4, then removing those genes from the downstream analysis. 
#Previous runs suggest even using regression with the percent metadata, the genes impact clustering.
OSNL_Matrix <- as.matrix(GetAssayData(object = OSNL, assay = "RNA", layer = "counts"))
OSNL_Tibble <- as_tibble(rownames = "Symbol", OSNL_Matrix)
OSNL_nCount_RNA <- OSNL@meta.data %>%
  select(nCount_RNA) %>%
  rownames_to_column("Cell")

Excised_OR_genes <- OSNL_Tibble %>% 
  filter(str_detect(Symbol, "^Olfr|10G4")) %>%
  pivot_longer(cols= -1, names_to = "Cell")  %>%
  left_join(OSNL_nCount_RNA, by = "Cell") %>%
  mutate(Normalized_Count = log1p(value/nCount_RNA * 10000))

OR_gene_vector <- OSNL_Tibble %>% 
  filter(str_detect(Symbol, "^Olfr|10G4")) %>%
  pull(Symbol)

#Saving OR10G4 RNA count data as metadata
OSNL$RNA_NormCount_10G4 <- Excised_OR_genes %>% 
  select(-c(nCount_RNA, value)) %>%
  pivot_wider(names_from = "Symbol", values_from = "Normalized_Count") %>%
  pull(`10G4`)

#Removing Olfr and 10G4 from features
counts <- GetAssayData(OSNL, assay = "RNA", layer = "counts")
counts <- counts[-(which(rownames(counts) %in% OR_gene_vector)),]
OSNL <- subset(OSNL, features = rownames(counts))

#Creating log1p(UMI/Total UMI * 10000) transformation for all Genes in each cell for the RNA assay to allow for appropriate visualization
DefaultAssay(OSNL) <- "RNA"
OSNL <- NormalizeData(OSNL)

#SCT and PCA
OSNL <- SCTransform(OSNL, vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA"), verbose = FALSE)
set.seed(123456)
OSNL <- RunPCA(OSNL, verbose = FALSE)

#Running RunHarmony followed by standard Clustering process using the "harmony" reduction as opposed to the usual "pca". 
#Then, a comparison of how the points are better mixed post-harmony.
OSNL <- RunHarmony(OSNL, 
                      group.by.vars = c("Genotype"), 
                      reduction = "pca", 
                      assay.use = "SCT", 
                      reduction.save = "harmony")

OSNL <- RunUMAP(OSNL, reduction = "harmony", assay = "SCT", dims = 1:12) 
OSNL <- FindNeighbors(OSNL, reduction = "harmony")
OSNL <- FindClusters(OSNL, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0))

#Decision about resolution:
#This resolution might oversplit some clusters, renaming clusters later will allow regrouping if necessary.
Idents(OSNL) <- "SCT_snn_res.2"

#In between these two steps, the clusters were reviewed for clarity and identity (not shown). See below.

#Upon reviewing the clusters and their identities (not shown here explicitly, but based on previously established markers)... 
#Two clusters of cells will be removed due to being non-OSN Lineage cell or containing an excess of such cells
#The previous steps for OSNL will be repeated
Cells_for_OSNL <- Cells(subset(OSNL, idents = c("0",  "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", 
                                                "11", "12", "13", "14", "15", "17")))
OSNL <- subset(H_Merge, cells = Cells_for_OSNL)

#Keep only genes that are expressed (count ≥ 1) in at least one cell
DefaultAssay(OSNL) <- "RNA"
OSNL <- OSNL[rowSums(OSNL[["RNA"]]$counts > 0) > 0, ]

#Manually rebuilding the Seurat object with some metadata and only RNA counts
idents.df = data.frame("orig.ident" = OSNL$orig.ident, "nCount_RNA" = OSNL$nCount_RNA, 
                       "nFeature_RNA" = OSNL$nFeature_RNA, "percent.mt" = OSNL$percent.mt,
                       "Genotype" = OSNL$Genotype, "percent.10G4" = OSNL$percent.10G4,
                       "percent.OR" = OSNL$percent.OR, "S.Score" = OSNL$S.Score, 
                       "G2M.Score" = OSNL$G2M.Score, "Phase" = OSNL$Phase)

OSNL = CreateSeuratObject(counts = LayerData(OSNL, assay = "RNA", layer = "counts"), project = "Export", meta.data = idents.df)

Idents(OSNL) = "Genotype"

#Saving the original count data for Olfrs and 10G4, then removing those genes from the downstream analysis. 
#Previous runs suggest even using regression with the percent metadata, the genes impact clustering.
OSNL_Matrix <- as.matrix(GetAssayData(object = OSNL, assay = "RNA", layer = "counts"))
OSNL_Tibble <- as_tibble(rownames = "Symbol", OSNL_Matrix)
OSNL_nCount_RNA <- OSNL@meta.data %>%
  select(nCount_RNA) %>%
  rownames_to_column("Cell")

Excised_OR_genes <- OSNL_Tibble %>% 
  filter(str_detect(Symbol, "^Olfr|10G4")) %>%
  pivot_longer(cols= -1, names_to = "Cell")  %>%
  left_join(OSNL_nCount_RNA, by = "Cell") %>%
  mutate(Normalized_Count = log1p(value/nCount_RNA * 10000))

OR_gene_vector <- OSNL_Tibble %>% 
  filter(str_detect(Symbol, "^Olfr|10G4")) %>%
  pull(Symbol)

#Saving OR10G4 RNA count data as metadata
OSNL$RNA_NormCount_10G4 <- Excised_OR_genes %>% 
  select(-c(nCount_RNA, value)) %>%
  pivot_wider(names_from = "Symbol", values_from = "Normalized_Count") %>%
  pull(`10G4`)

#Removing Olfr and 10G4 from features
counts <- GetAssayData(OSNL, assay = "RNA", layer = "counts")
counts <- counts[-(which(rownames(counts) %in% OR_gene_vector)),]
OSNL <- subset(OSNL, features = rownames(counts))

#17_OSNL_Metadata_TG_Olfr_byGenotype
VlnPlot(OSNL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.10G4", "percent.OR"), ncol = 5, group.by = "Genotype") &
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave("Output2/17_OSNL_Metadata_TG_Olfr_byGenotype.png", width = 30, height = 15, units = "cm")

#Creating log1p(UMI/Total UMI * 10000) transformation for all Genes in each cell for the RNA assay to allow for appropriate visualization
DefaultAssay(OSNL) <- "RNA"
OSNL <- NormalizeData(OSNL)

#SCT and PCA
OSNL <- SCTransform(OSNL, vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA"), verbose = FALSE)
set.seed(123456)
OSNL <- RunPCA(OSNL, verbose = FALSE)

#Running RunHarmony followed by standard Clustering process using the "harmony" reduction as opposed to the usual "pca". 
#Then, a comparison of how the points are better mixed post-harmony.
OSNL <- RunHarmony(OSNL, 
                   group.by.vars = c("Genotype"), 
                   reduction = "pca", 
                   assay.use = "SCT", 
                   reduction.save = "harmony")

#18_OSNL_BeforeAfter_RunHarmony
p1 <- DimPlot(OSNL, reduction = "pca", group.by = "Genotype", pt.size = 0.5) + 
  NoLegend() + 
  labs(title = "Before RunHarmony()") +
  theme(axis.title = element_blank())
p2 <- DimPlot(OSNL, reduction = "harmony", group.by = "Genotype", pt.size = 0.5) + 
  NoLegend() + 
  labs(title = "After RunHarmony()") +
  theme(axis.title = element_blank())

p1 + p2
ggsave("Output2/18_OSNL_BeforeAfter_RunHarmony_SF5A.png", width = 30, height = 20, units = "cm")

#19_OSNL_ElbowPlot_selectDimsCleaner
ElbowPlot(OSNL, ndims = 50) #20 seems fine

ggsave("Output2/19_OSNL_ElbowPlot_selectDimsCleaner.png", width = 25, height = 18, units = "cm")

OSNL <- RunUMAP(OSNL, reduction = "harmony", assay = "SCT", dims = 1:13) 
OSNL <- FindNeighbors(OSNL, reduction = "harmony")
OSNL <- FindClusters(OSNL, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0))

#19B_OSNL_UMAP_GenotypeSplit
DimPlot(OSNL, split.by = "Genotype") & theme_bw() & NoAxes() & NoLegend()
ggsave("Output2/19B_OSNL_UMAP_GenotypeSplit_SF5B.png", width = 30, height = 15, units = "cm")

#Reviewing the cluster changes over resolution and making a selection
#20_OSNL_ClustreeCleaner
clustree(OSNL, prefix = "SCT_snn_res.")

ggsave("Output2/20_OSNL_ClustreeCleaner_SF5C.png", width = 30, height = 30, units = "cm")

Idents(OSNL) <- "SCT_snn_res.0.2"
p1 <- DimPlot(OSNL, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.0.2")
Idents(OSNL) <- "SCT_snn_res.0.4"
p2 <- DimPlot(OSNL, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.0.4")
Idents(OSNL) <- "SCT_snn_res.0.6"
p3 <- DimPlot(OSNL, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.0.6")
Idents(OSNL) <- "SCT_snn_res.0.8"
p4 <- DimPlot(OSNL, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.0.8")
Idents(OSNL) <- "SCT_snn_res.1"
p5 <- DimPlot(OSNL, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.1")
Idents(OSNL) <- "SCT_snn_res.1.2"
p6 <- DimPlot(OSNL, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.1.2")
Idents(OSNL) <- "SCT_snn_res.1.4"
p7 <- DimPlot(OSNL, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.1.4")
Idents(OSNL) <- "SCT_snn_res.1.6"
p8 <- DimPlot(OSNL, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.1.6")
Idents(OSNL) <- "SCT_snn_res.1.8"
p9 <- DimPlot(OSNL, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.1.8")
Idents(OSNL) <- "SCT_snn_res.2"
p10 <- DimPlot(OSNL, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "SCT_snn_res.2")

#21_OSNL_ResolutionCompare_UMAPCleaner
(p1 | p2 | p3 | p4 | p5 )/( p6 | p7 | p8 | p9 | p10) + 
  plot_annotation(theme = theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))) & 
  NoAxes()

ggsave("Output2/21_OSNL_ResolutionCompare_UMAPCleaner_SF5D.png", width = 50, height = 20, units = "cm")

#Final Decision about resolution: 2 resolution: 20 clusters
#This resolution might oversplit some clusters, renaming clusters later will allow regrouping if necessary.
Idents(OSNL) <- "SCT_snn_res.2"

#Some overview of the metadata qualities of each cluster. 10G4 cells clearly have less genes and counts.

#22_OSNL_Metadata_TG_Olfr_byClusterV2
VlnPlot(OSNL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.10G4", "percent.OR"), ncol = 5, split.by = "Genotype") &
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave("Output2/22_OSNL_Metadata_TG_Olfr_byClusterV2.png", width = 50, height = 10, units = "cm")

##Looking at the Markers of OSN Lineage, from GBC to mOSN
#23_OSNL_Known_Markers_VlnUMAPLabeled To run after labeling clusters
p1 <- VlnPlot(OSNL, features = c("G2M.Score", "S.Score", "Hes6", "Top2a", "Kit", "Ascl1", "Neurog1", "Neurod1", "Lhx2", "Gng8", 
                                 "Myt1l", "Tex15", "Gap43", "Omp", "Gnal", "Gnb1", "Gng13", "Adcy3", "Cnga2", "Stoml3"),
              ncol = 10, assay = "RNA", group.by = "Lineage") &
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))
p2 <- FeaturePlot(OSNL, features = c("G2M.Score", "S.Score", "Hes6", "Top2a", "Kit", "Ascl1", "Neurog1", "Neurod1", "Lhx2", "Gng8", 
                                     "Myt1l", "Tex15", "Gap43", "Omp", "Gnal", "Gnb1", "Gng13", "Adcy3", "Cnga2", "Stoml3"),
                  ncol = 10, order = TRUE) &
  theme_bw() &
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) &
  NoAxes() & 
  NoLegend()

p1/p2
ggsave("Output2/23_OSNL_Known_Markers_VlnUMAPLabeled_SF6A.png", width = 60, height = 25, units = "cm")

#Subsetting OSNL to remove cluster 19 due to its unusual features.
OSNL <- subset(OSNL, idents = "19", invert = TRUE)

#Adding some presumptive, overly specific Cluster labels
OSNL[["Clusters_res2.0"]] <- Idents(OSNL)
new.cluster.ids <- c("mOSN_Ventral_G4bias_A", "mOSN_Dorsal", "mOSN_Mixed_Cd55_ORbias",
                     "mOSN_Dorsal_Calb2_ORbias", "iOSN3", "iOSN3",
                     "mOSN_Ventral_G4bias_A", "mOSN_Ventral_G4bias_A", "mOSN_Ventral_ORbias",
                     "mOSN_Ventral_G4bias_B", "iOSN3", "iOSN1",
                     "mOSN_Ventral_Dlg2_ORbias", "iOSN2", "iOSN3",
                     "GBC/INP", "iOSN2", "mOSN_Ventral_G4bias_A", "mOSN_Dorsal_Calb2",
                     "mOSN_Ventral")
names(new.cluster.ids) <- levels(OSNL)
OSNL <- RenameIdents(OSNL, new.cluster.ids)

levels(OSNL) <- c("mOSN_Ventral_G4bias_A", "mOSN_Ventral_G4bias_B", "mOSN_Ventral_ORbias",
                  "mOSN_Ventral_Dlg2_ORbias", "mOSN_Ventral",
                  "mOSN_Mixed_Cd55_ORbias",
                  "mOSN_Dorsal", "mOSN_Dorsal_Calb2_ORbias", "mOSN_Dorsal_Calb2",
                  "iOSN3", "iOSN2", "iOSN1", "GBC/INP")

OSNL[["Lineage"]] <- Idents(OSNL)

#25_OSNL_Labeled_UMAP
DimPlot(OSNL, reduction = "umap", label = TRUE, repel = TRUE, label.box = TRUE) + 
  theme_bw() + 
  NoLegend() + 
  NoAxes()

ggsave("Output2/25_OSNL_Labeled_UMAP_MF3B.png", width = 25, height = 25, units = "cm")

#Rerunning the DEG to properly order the DotPlot
MarkerOSNL2 <- FindAllMarkers(OSNL, logfc.threshold = 0.5, min.pct = 0.4, only.pos = TRUE)

Top5_OSNL2 <- MarkerOSNL2 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 7)

#26_OSNL_Labeled_DotPlot_Top5DEG
DotPlot(OSNL, features = unique(Top5_OSNL2$gene)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(axis.title = element_blank())

ggsave("Output2/26_OSNL_Labeled_DotPlot_Top7DEG_SF6C.png", width = 60, height = 12, units = "cm")

#28_OSNL_Labeled_Vln_Markers
p1 <- VlnPlot(OSNL, features = c("Ascl1", "Neurog1", "Lhx2", "Myt1l", "Gap43", "Omp", "Stoml3", 
                                        "Acsm4", "Nfix", "Calb2", "Cd55", "Dlg2", "percent.10G4", "percent.OR"), ncol = 14, pt.size = 0, assay = "RNA") & 
  coord_flip() & 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

p2 <- VlnPlot(OSNL, features = "G2M.Score", ncol = 1, pt.size = 0) + 
  coord_flip() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

p2 + p1 + plot_layout(widths = c(0.4, 6.0), ncol = 2)

ggsave("Output2/28_OSNL_Labeled_Vln_Markers_MF3C.png", width = 60, height = 15, units = "cm")

#29_OSNL_OtherMarkers_VlnUMAP_Labeled
#A review of the unlabeled Top5 dotplot reveals Dorsal/Ventral Markers that split the presumptive mOSN
p1 <- VlnPlot(OSNL, features = c("Nfix", "Ncam2", "Acsm4", "Nqo1", "Calb2", "Cd55", "Cd36", "Dlg2", "percent.10G4", "percent.OR"), 
              ncol = 10, assay = "RNA") &
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))
p2 <- FeaturePlot(OSNL, features = c("Nfix", "Ncam2", "Acsm4", "Nqo1", "Calb2", "Cd55", "Cd36", "Dlg2", "percent.10G4", "percent.OR"),
                  ncol = 10, order = TRUE) &
  theme_bw() &
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) &
  NoAxes() & 
  NoLegend()

p1/p2
ggsave("Output2/29_OSNL_OtherMarkers_VlnUMAP_Labeled_SF6B.png", width = 70, height = 15, units = "cm")

#29x_OSNL_Labeled_Vln_OR10G4_Mbnl1
p1 <- VlnPlot(OSNL, features = "Omp", assay = "RNA", split.by = "Genotype") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")
p2 <- VlnPlot(OSNL, features = "percent.10G4", split.by = "Genotype") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")
p3 <- VlnPlot(OSNL, features = "Mbnl1", assay = "RNA", split.by = "Genotype") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

p4 <- VlnPlot(OSNL, features = "percent.mt", split.by = "Genotype") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")
p5 <- VlnPlot(OSNL, features = "nCount_RNA", split.by = "Genotype") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")
p6 <- VlnPlot(OSNL, features = "nFeature_RNA", split.by = "Genotype") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

(p1 | p4) / (p2 | p5) / (p3 | p6) 

ggsave("Output2/29x_H_OSNL_Labeled_Vln_OR10G4_Mbnl1_SF5E.png", width = 70, height = 20, units = "cm")

#A second, simpler clustering scheme below
Idents(OSNL) <- "Lineage"
new.cluster.ids <- c("mOSN", "mOSN", "mOSN", "mOSN", "mOSN", "mOSN", "mOSN", "mOSN", "mOSN",
                     "iOSN3", "iOSN2", "iOSN1", "GBC/INP")
names(new.cluster.ids) <- levels(OSNL)
OSNL <- RenameIdents(OSNL, new.cluster.ids)
levels(OSNL) <- c("mOSN", "iOSN3", "iOSN2", "iOSN1", "GBC/INP")
OSNL[["Simple_Clusters"]] <- Idents(OSNL)

#30A_OSNL_UMAP_Simple_Clusters
DimPlot(OSNL, reduction = "umap", label = TRUE, repel = TRUE, label.box = TRUE, label.size = 7) + 
  NoLegend() +
  NoAxes()

ggsave("Output2/30A_OSNL_UMAP_Simple_Clusters.png", width = 10, height = 10, units = "cm")

#30B_OSNL_UMAP_Markers_Simple_Clusters
p1 <- VlnPlot(OSNL, features = c("Ascl1", "Neurog1", "Lhx2", "Myt1l", "Gap43", "Stoml3"), ncol = 6, assay = "RNA") & 
  theme(axis.title = element_blank())
p2 <- FeaturePlot(OSNL, features = c("Ascl1", "Neurog1", "Lhx2", "Myt1l", "Gap43", "Stoml3"), ncol = 6, order = TRUE) & 
  theme_bw() &
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) &
  NoLegend() & 
  NoAxes()

p3 <- VlnPlot(OSNL, features = c("percent.10G4", "percent.OR"), split.by = "Genotype") & 
  theme(axis.title = element_blank())
p4 <- FeaturePlot(OSNL, features = c("percent.10G4", "percent.OR"), ncol = 2, order = TRUE) & 
  theme_bw() &
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) &
  NoLegend() & 
  NoAxes()

design <- "ABCCCCCCC
           DDEEEEEEE"
p3 + p1 + p4 + p2 +
  plot_layout(design = design)

ggsave("Output2/30B_OSNL_UMAP_Markers_Simple_Clusters.png", width = 60, height = 20, units = "cm")

#Section IV: Threshold subsection----
OSN_Metadata <- OSNL@meta.data %>%
  select(Genotype, Lineage, Simple_Clusters) %>%
  rownames_to_column("Cell")

OR_Counts <- Excised_OR_genes %>%
  right_join(OSN_Metadata, by = "Cell") %>%
  select(Cell, Genotype, Lineage, Simple_Clusters, nCount_RNA, everything()) %>%
  dplyr::rename(Count = value)

#40_OSN_Lineage_ORGenesvsNon-ORgenes_MaxCount
OSNL_Matrix <- as.matrix(GetAssayData(object = OSNL, assay = "RNA", layer = "counts"))
OSNL_Tibble <- as_tibble(rownames = "Symbol", OSNL_Matrix)

OR_CellSymbolCount <- OR_Counts %>%
  select(Cell, Symbol, Count)

OR_CellOther <- OR_Counts %>%
  select(-c(Symbol, Count, Normalized_Count)) %>%
  unique()

Another_intermediateV2 <- OSNL_Tibble %>% 
  pivot_longer(cols = -Symbol, names_to = "Cell", values_to = "Count") %>%
  bind_rows(OR_CellSymbolCount) %>%
  left_join(OR_CellOther, by = "Cell") %>%
  filter(Simple_Clusters == "mOSN" & Count > 0) %>%
  mutate(Ncount = Count/nCount_RNA * 10000) %>%
  group_by(Symbol) %>% 
  summarize(MaxCount = max(Ncount), MedCount = median(Ncount)) %>% 
  arrange(desc(MedCount)) %>% 
  mutate(GeneType = case_when(str_detect(Symbol, "Olfr|10G4") ~ "OR",
                              TRUE ~ "Non-OR")) %>% 
  mutate(GeneType2 = case_when(str_detect(Symbol, "Olfr") ~ "OR",
                               Symbol == "10G4" ~ "OR10G4",
                              TRUE ~ "Non-OR")) %>% 
  mutate(rank_all_max = rank(-MaxCount, ties.method = "min"),
         rank_all_med = rank(-MedCount, ties.method = "min")) %>%
  group_by(GeneType) %>%
  mutate(rank_OR_max = if_else(GeneType %in% c("OR","OR10G4"), rank(-MaxCount, ties.method = "min"), NA_real_),
         rank_OR_med = if_else(GeneType %in% c("OR", "OR10G4"), rank(-MedCount, ties.method = "min"), NA_real_)) %>%
  ungroup() %>% 
  mutate(Symbol = fct_reorder(Symbol, MedCount, .desc = TRUE)) 

Another_intermediateV2 %>% 
  ggplot(aes(Symbol, log1p(MedCount), fill = GeneType2)) + 
  geom_col(width = 1, alpha = 0.8) +
  geom_col(data = filter(Another_intermediateV2, Symbol == "10G4"),
    aes(Symbol, log1p(MedCount), fill = NULL), color = "red", fill = "red") +
  labs(x = "All Expressed Genes", y = "Log Transformed Median non-Zero Expression", fill = "Gene Category") +
  scale_fill_manual(values = c("OR" = "steelblue", "OR10G4" = "red", "Non-OR" = "green3")) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.5))

ggsave("Output2/40_OSN_Lineage_ORGenesvsNon-ORgenes_MaxCount_SF7B.png", width = 50, height = 10, units = "cm")

#4x_OSN_Lineage_mOSN_NormalizedCount_byStatus_ComparisonCountsAlt
OR_Counts %>%
  filter(Simple_Clusters == "mOSN") %>%
  filter(!(Genotype == "WT" & Symbol == "10G4")) %>%
  group_by(Cell) %>% 
  filter(any(Normalized_Count > 0) & Normalized_Count > 0 | (!any(Normalized_Count > 0) & row_number() == 1)) %>%
  mutate(Kind = ifelse(sum(Normalized_Count) == 0, "No OR", Symbol)) %>%
  ungroup() %>%
  mutate(Kind2 = case_when(Kind == "No OR" ~ "No OR",
                           Symbol == "10G4" ~ "OR10G4 Expression",
                           TRUE ~ "Endogenous OR Expression")) %>%
  filter(Kind2 != "No OR") %>%
  group_by(Cell) %>%
  mutate(Status1 = ifelse(n() == 1, "Single OR", "Multiple OR"),
         Status2 = case_when("OR10G4 Expression" %notin% Kind2 ~ "Endogenous OR Only",
                             "OR10G4 Expression" %in% Kind2 & Status1 == "Single OR" ~ "OR10G4 Only",
                             TRUE ~ "Endogenous OR/OR10G4")) %>%
  ungroup() %>%
  mutate(Status3 = paste0("Genotype: ", Genotype, "\n", Status1, " mOSN\n", Status2, "\n", Kind2)) %>%
  mutate(Group = factor(Status3)) %>%
  mutate(Group = fct_relevel(Group, levels(Group)[c(7, 6, 4, 1, 2, 3, 5)])) %>%
  mutate(Combo1 = paste(Genotype, Kind2)) %>%
  ggplot(aes(Group, Normalized_Count)) +
  geom_violin(aes(fill = Combo1), alpha = 0.6) +
  geom_jitter(alpha = 0.5, width = 0.15) +
  labs(title = "OR10G4 does not operate on the same Bimodal Expression Pattern of Endogenous ORs",
       x = "Genotype, mOSN type, OR kinds in mOSN, Expression source", 
       y = "Log Transformed Normalized OR Counts") + 
  scale_fill_discrete(guide = "none") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) 

ggsave("Output2/4x_OSN_Lineage_mOSN_NormalizedCount_byStatus_ComparisonCountsAlt_SF7A.png", width = 40, height = 15, units = "cm")

OSNL_Labeled_OR_Counts <- OR_Counts %>% 
  filter(!(Genotype == "WT" & Symbol == "10G4")) %>%
  group_by(Cell) %>% 
  filter(any(Normalized_Count > 0) & Normalized_Count > 0 | (!any(Normalized_Count > 0) & row_number() == 1)) %>%
  filter(Simple_Clusters != "GBC/INP") %>% 
  mutate(PercentTotalOR = round((exp(Normalized_Count) - 1) /sum((exp(Normalized_Count) - 1) + 0.0001) * 100, 4)) %>%
  ungroup() %>%
  mutate(Symbol_copy = Symbol) %>%
  mutate(Symbol = ifelse(Symbol == "10G4", "OR10G4 Expressed", "Non-G4 OR Expressed")) %>%
  mutate(Symbol = ifelse(Count == 0, "No OR", Symbol)) %>%
  group_by(Cell) %>%
  mutate(Status = case_when(n() > 1 & "OR10G4 Expressed" %in% Symbol ~ "OR10G4/EndoOR",
                            n() > 1 & "OR10G4 Expressed" %notin% Symbol ~ "EndoOR/EndoOR",
                            TRUE ~ Symbol)) %>%
  mutate(CellMean = mean(Normalized_Count)) %>%
  ungroup() %>%
  arrange(Status, desc(CellMean)) %>%
  mutate(Cell = factor(Cell, levels = unique(Cell))) %>%
  mutate(Count_Cat = case_when(Normalized_Count < 1 ~ "Sub1",
                               1 <= Normalized_Count & Normalized_Count < 2 ~ "Sub2",
                               2 <= Normalized_Count & Normalized_Count < 3 ~ "Sub3", 
                               3 <= Normalized_Count & Normalized_Count < 4 ~ "Sub4",
                               4 <= Normalized_Count & Normalized_Count < 5 ~ "Sub5",
                               5 <= Normalized_Count & Normalized_Count < 6 ~ "Sub6",
                               TRUE ~ "Above6")) %>%
  mutate(OR_Expression_Status = case_when(Status == "EndoOR/EndoOR" ~ "Non-S_OSN Endo_OR_Only Endo_OR_Exp",
                                          Status == "No OR" ~ "No_OR",
                                          Status == "Non-G4 OR Expressed" ~ "S_OSN Endo_OR_Exp",
                                          Status == "OR10G4 Expressed" ~ "S_OSN 10G4_Exp",
                                          Status == "OR10G4/EndoOR" & Symbol == "OR10G4 Expressed" ~ "Non_S-OSN Mixed_OR 10G4_Exp",
                                          TRUE ~ "Non_S-OSN Mixed_OR Endo_OR_Exp")) %>%
  group_by(Cell) %>%
  mutate(OR_Rank = dense_rank(dplyr::desc(Normalized_Count))) %>%
  mutate(OR_Expression_Status = factor(OR_Expression_Status, 
                                       levels = c("S_OSN Endo_OR_Exp", 
                                                  "Non-S_OSN Endo_OR_Only Endo_OR_Exp",
                                                  "Non_S-OSN Mixed_OR Endo_OR_Exp", 
                                                  "Non_S-OSN Mixed_OR 10G4_Exp",
                                                  "S_OSN 10G4_Exp", "No_OR")),
         OR_Rank = factor(OR_Rank)) %>%
  filter(OR_Expression_Status != "No_OR") %>%
  ungroup() %>%
  mutate(New_Status = case_when(OR_Expression_Status == "S_OSN Endo_OR_Exp" ~ "Single OR OSN\nEndogenous OR Counts", 
                                OR_Expression_Status == "Non-S_OSN Endo_OR_Only Endo_OR_Exp" ~ "Multiple OR OSN\nOnly Endogenous ORs\nEndogenous OR Counts",
                                OR_Expression_Status == "Non_S-OSN Mixed_OR Endo_OR_Exp" ~ "Multiple OR OSN\nEndogenous OR & OR10G4\nEndogenous OR Counts",
                                OR_Expression_Status == "Non_S-OSN Mixed_OR 10G4_Exp" ~ "Multiple OR OSN\nEndogenous OR & OR10G4\nOR10G4 Counts", 
                                OR_Expression_Status == "S_OSN 10G4_Exp" ~ "Single OR OSN\nOR10G4 Counts")) %>%
  mutate(OR_Rankv2 = ifelse(OR_Rank == 1, "Top", "Not Top"))

#46_OSN_Lineage_New_OR_Categories_PercentTotal
celltype_levels <- c(
  "High Multiple",
  "High/Intermediate",
  "High",
  "Intermediate Multiple",
  "Intermediate",
  "Low")

#See the cutoff values to be used in High/Intermediate/Low distinction.
OSNL_Labeled_OR_Counts %>% 
  group_by(OR_Rankv2, New_Status, Simple_Clusters) %>% 
  summarize(MeanNC = mean(Normalized_Count), MedNC = median(Normalized_Count), SDNC = sd(Normalized_Count)) %>% 
  mutate(OneSD = ifelse(OR_Rankv2 == "Top", MeanNC - SDNC, MeanNC + SDNC)) %>% 
  select(OR_Rankv2, New_Status, Simple_Clusters, OneSD) %>%
  pivot_wider(names_from = OR_Rankv2, values_from = OneSD) %>%
  filter(Simple_Clusters == "mOSN")

Intermediate_OSNL_Labeled_OR_Counts <- OSNL_Labeled_OR_Counts %>% 
  group_by(OR_Rankv2, New_Status, Simple_Clusters) %>% 
  summarize(MeanNC = mean(Normalized_Count), MedNC = median(Normalized_Count), SDNC = sd(Normalized_Count)) %>% 
  mutate(OneSD = ifelse(OR_Rankv2 == "Top", MeanNC - SDNC, MeanNC + SDNC)) %>% 
  select(OR_Rankv2, New_Status, Simple_Clusters, OneSD) %>%
  pivot_wider(names_from = OR_Rankv2, values_from = OneSD) %>%
  filter(Simple_Clusters == "mOSN") %>%
  select(-Simple_Clusters) %>%
  right_join(OSNL_Labeled_OR_Counts, by = c("New_Status")) %>%
  mutate(CountCat3 = case_when(Normalized_Count >= Top ~ "High", 
                             Normalized_Count <= 1.04 ~ "Low",
                             TRUE ~ "Intermediate")) %>%
  group_by(Cell, Simple_Clusters, Genotype) %>%
  summarise(
    categories = list(CountCat3),
    n_high = sum(CountCat3 == "High"),
    n_trans = sum(CountCat3 == "Intermediate"),
    n_low = sum(CountCat3 == "Low"),
    n_total = n()) %>%
  mutate(Cell_Type = case_when(
      n_high >= 2 ~ "High Multiple",
      n_high >= 1 & n_trans >= 1 ~ "High/Intermediate",
      n_high >= 1 & n_low >= 1 ~ "High/Low",
      n_high == 1 & n_total == 1 ~ "High Single",
      n_trans >= 2 ~ "Intermediate Multiple",
      n_trans >= 1 & n_low >= 1 ~ "Intermediate/Low",
      n_trans == 1 & n_total == 1 ~ "Intermediate Single",
      n_low >= 2 ~ "Low Multiple",
      n_low == 1 & n_total == 1 ~ "Low Single",
      TRUE ~ "Unknown")) %>%
  mutate(Cell_Type = case_when(Cell_Type %in% c("High/Low", "High Single") ~ "High",
                               Cell_Type %in% c("Intermediate/Low", "Intermediate Single") ~ "Intermediate",
                               Cell_Type %in% c("Low Multiple", "Low Single") ~ "Low",
                               TRUE ~ Cell_Type))

No_OR_OSNs <- OR_Counts %>% 
  group_by(Cell) %>% 
  filter(any(Normalized_Count > 0) & Normalized_Count > 0 | (!any(Normalized_Count > 0) & row_number() == 1)) %>%
  filter(Simple_Clusters != "GBC/INP") %>% 
  mutate(PercentTotalOR = round((exp(Normalized_Count) - 1) /sum((exp(Normalized_Count) - 1) + 0.0001) * 100, 4)) %>%
  ungroup() %>%
  mutate(Symbol_copy = Symbol) %>%
  mutate(Symbol = ifelse(Symbol == "10G4", "OR10G4 Expressed", "Non-G4 OR Expressed")) %>%
  mutate(Symbol = ifelse(Count == 0, "No OR", Symbol)) %>%
  filter(Symbol == "No OR")

No_OR_OSNs %>%
  select(Cell, Genotype, Simple_Clusters) %>%
  bind_cols(tibble(Cell_Type = "Low", n_total = 0)) %>%
  bind_rows(Intermediate_OSNL_Labeled_OR_Counts) %>%
  group_by(Cell, Simple_Clusters, Cell_Type, Genotype) %>%
  summarise(n = n()) %>%
  group_by(Simple_Clusters, Genotype) %>%
  mutate(percent = n / sum(n) * 100) %>%
  mutate(Cell_Type = factor(Cell_Type, levels = celltype_levels)) %>%
  ggplot(aes(x = Genotype, y = percent, fill = Cell_Type)) +
  geom_col(position = "fill") + 
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Genotype", y = "Proportion of Cells", fill = "Expressed OR Category") +
  facet_wrap(~Simple_Clusters, nrow = 1)+
  theme_bw() +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1))

ggsave("Output2/46_OSN_Lineage_New_OR_Categories_PercentTotal_MF3E.png", width = 25, height = 20, units = "cm")

#47_OSN_Lineage_New_OR_Categories_PercentTotal_G4ExpSplit
Intermediate_Again_splittingbyG4Exp <- OSNL_Labeled_OR_Counts %>% 
  group_by(OR_Rankv2, New_Status, Simple_Clusters) %>% 
  summarize(MeanNC = mean(Normalized_Count), MedNC = median(Normalized_Count), SDNC = sd(Normalized_Count)) %>% 
  mutate(OneSD = ifelse(OR_Rankv2 == "Top", MeanNC - SDNC, MeanNC + SDNC)) %>% 
  select(OR_Rankv2, New_Status, Simple_Clusters, OneSD) %>%
  pivot_wider(names_from = OR_Rankv2, values_from = OneSD) %>%
  filter(Simple_Clusters == "mOSN") %>%
  select(-Simple_Clusters) %>%
  right_join(OSNL_Labeled_OR_Counts, by = c("New_Status")) %>%
  mutate(CountCat3 = case_when(Normalized_Count >= Top ~ "High", 
                               Normalized_Count <= 1.04 ~ "Low",
                               TRUE ~ "Intermediate")) %>% 
  ungroup() %>% 
  mutate(G4NonLow = Symbol_copy == "10G4" & CountCat3 != "Low") %>% 
  group_by(Cell) %>%
  mutate(G4Exp = ifelse("10G4" %in% Symbol_copy & TRUE %in% G4NonLow, "OR10G4 Expressed", "No OR10G4 Expressed")) %>%
  ungroup() %>%
  group_by(Cell, Simple_Clusters, Genotype, G4Exp) %>%
  summarise(
    categories = list(CountCat3),
    n_high = sum(CountCat3 == "High"),
    n_trans = sum(CountCat3 == "Intermediate"),
    n_low = sum(CountCat3 == "Low"),
    n_total = n()) %>%
  mutate(Cell_Type = case_when(
    n_high >= 2 ~ "High Multiple",
    n_high >= 1 & n_trans >= 1 ~ "High/Intermediate",
    n_high >= 1 & n_low >= 1 ~ "High/Low",
    n_high == 1 & n_total == 1 ~ "High Single",
    n_trans >= 2 ~ "Intermediate Multiple",
    n_trans >= 1 & n_low >= 1 ~ "Intermediate/Low",
    n_trans == 1 & n_total == 1 ~ "Intermediate Single",
    n_low >= 2 ~ "Low Multiple",
    n_low == 1 & n_total == 1 ~ "Low Single",
    TRUE ~ "Unknown")) %>%
  mutate(Cell_Type = case_when(Cell_Type %in% c("High/Low", "High Single") ~ "High",
                               Cell_Type %in% c("Intermediate/Low", "Intermediate Single") ~ "Intermediate",
                               Cell_Type %in% c("Low Multiple", "Low Single") ~ "Low",
                               TRUE ~ Cell_Type))

No_OR_OSNs %>%
  select(Cell, Genotype, Simple_Clusters) %>%
  bind_cols(tibble(Cell_Type = "Low", n_total = 0, G4Exp = "No OR10G4 Expressed")) %>%
  bind_rows(Intermediate_Again_splittingbyG4Exp) %>%
  group_by(Cell, Simple_Clusters, Cell_Type, Genotype, G4Exp) %>%
  summarise(n = n()) %>%
  group_by(Simple_Clusters, G4Exp) %>%
  mutate(percent = n / sum(n) * 100) %>%
  mutate(Cell_Type = factor(Cell_Type, levels = celltype_levels)) %>%
  ggplot(aes(x = Simple_Clusters, y = percent, fill = Cell_Type)) +
  geom_col(position = "fill") + 
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Genotype", y = "Proportion of Cells", fill = "Expressed OR Category") +
  facet_wrap(G4Exp~Genotype, drop = TRUE)+
  theme_bw() +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1))

ggsave("Output2/47_OSN_Lineage_New_OR_Categories_PercentTotal_G4ExpSplit_SF7D.png", width = 50, height = 15, units = "cm")

#Table3_OR_Count_Categories
OR_Count_Categorization <- OSNL_Labeled_OR_Counts %>% 
   group_by(OR_Rankv2, New_Status, Simple_Clusters) %>% 
   summarize(MeanNC = mean(Normalized_Count), MedNC = median(Normalized_Count), SDNC = sd(Normalized_Count)) %>% 
   mutate(OneSD = ifelse(OR_Rankv2 == "Top", MeanNC - SDNC, MeanNC + SDNC)) %>% 
   select(OR_Rankv2, New_Status, Simple_Clusters, OneSD) %>%
   pivot_wider(names_from = OR_Rankv2, values_from = OneSD) %>%
   filter(Simple_Clusters == "mOSN") %>%
   select(-Simple_Clusters) %>%
   right_join(OSNL_Labeled_OR_Counts, by = c("New_Status")) %>%
   mutate(CountCat3 = case_when(Normalized_Count >= Top ~ "High", Normalized_Count <= 1.04 ~ "Low", TRUE ~ "Intermediate")) %>% 
   ungroup() %>% 
   group_by(Cell) %>%
   mutate(G4Detect = "10G4" %in% Symbol_copy) %>%
   ungroup() %>%
   group_by(Cell, Simple_Clusters, Genotype, G4Detect) %>%
   mutate(
       categories = list(CountCat3),
       n_high = sum(CountCat3 == "High"),
       n_trans = sum(CountCat3 == "Intermediate"),
       n_low = sum(CountCat3 == "Low"),
       n_total = n()) %>%
   mutate(Cell_Type = case_when(
       n_high >= 2 ~ "High Multiple",
       n_high >= 1 & n_trans >= 1 ~ "High/Intermediate",
       n_high >= 1 & n_low >= 1 ~ "High/Low",
       n_high == 1 & n_total == 1 ~ "High Single",
       n_trans >= 2 ~ "Intermediate Multiple",
       n_trans >= 1 & n_low >= 1 ~ "Intermediate/Low",
       n_trans == 1 & n_total == 1 ~ "Intermediate Single",
       n_low >= 2 ~ "Low Multiple",
       n_low == 1 & n_total == 1 ~ "Low Single",
       TRUE ~ "Unknown")) %>%
   mutate(Cell_OR_Cat = case_when(Cell_Type %in% c("High/Low", "High Single") ~ "High",
                                  Cell_Type %in% c("Intermediate/Low", "Intermediate Single") ~ "Intermediate",
                                  Cell_Type %in% c("Low Multiple", "Low Single") ~ "Low",
                                  TRUE ~ Cell_Type)) %>% 
  ungroup() %>% 
  mutate(MultipleDetect = n_total > 1) %>% 
  dplyr::count(Genotype, Symbol, Simple_Clusters, MultipleDetect, CountCat3, G4Detect, Cell_OR_Cat) %>% 
  filter(Simple_Clusters == "mOSN") %>% 
  pivot_wider(names_from = Cell_OR_Cat, values_from = n, values_fill = 0, names_prefix = "mOSN: ") %>% 
  dplyr::rename(OR_Group = Symbol, 
         OR_Count_Category = CountCat3, 
         `mOSN: OR10G4 Detected?` = G4Detect,
         `mOSN: Multiple ORs Detected?` = MultipleDetect) %>% 
  mutate(OR_Group = ifelse(OR_Group == "OR10G4 Expressed", "OR10G4", "Endogenous OR")) %>% 
  select(Genotype, `mOSN: Multiple ORs Detected?`, `mOSN: OR10G4 Detected?`, OR_Group, OR_Count_Category, `mOSN: High Multiple`, 
         `mOSN: High/Intermediate`, `mOSN: High`, `mOSN: Intermediate Multiple`, `mOSN: Intermediate`, `mOSN: Low`)

OR_Count_Categorization2 <- OR_Count_Categorization %>%
  select(1:3) %>%
  unique() %>%
  bind_cols(tibble(`Total mOSN Number` = c(30, 14, 414, 528, 152, 231))) %>%
  right_join(OR_Count_Categorization, by = colnames(OR_Count_Categorization[,1:3]))

write.csv(OR_Count_Categorization2, file = "Output2/Table3_OR_Count_Categories.csv")
#Table Edited after writing to file

#Section V: Pseudotime work----
#Cluster relabeling using a simplified alphabetical approach GBC/INP -> A, mOSN -> J
Idents(OSNL) <- "Clusters_res2.0"
new.cluster.ids <- c("J", "J", "J", "J", "F", "G", "I", "J", "J", "J", "H", "B", "J", "C", "E", "A", "D", "J", "J", "J")
names(new.cluster.ids) <- levels(OSNL)
OSNL <- RenameIdents(OSNL, new.cluster.ids)
levels(OSNL) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
OSNL[["Simple_Lineagev2"]] <- Idents(OSNL)

# Calculate cluster centroids (for plotting the labels later)
mm <- sparse.model.matrix(~ 0 + factor(OSNL$Simple_Lineagev2))
colnames(mm) <- levels(factor(OSNL$Simple_Lineagev2))
centroids2d <- as.matrix(t(t(OSNL@reductions$pca@cell.embeddings[,1:3]) %*% mm) / Matrix::colSums(mm))

# Define lineage ends
ENDS <- c("J")

set.seed(1)
lineages <- as.SlingshotDataSet(getLineages(
  data           = OSNL@reductions$pca@cell.embeddings[,1:3],
  clusterLabels  = OSNL$Simple_Lineagev2,
  dist.method    = "mnn", 
  end.clus       = ENDS, 
  start.clus     = "A"
)) # define where to START the trajectories

# Define curves
curves <- as.SlingshotDataSet(getCurves(
  data          = lineages,
  thresh        = 15,
  stretch       = 10,
  allow.breaks  = F,
  approx_points = 100,
  extend = 'n'
))

pseudotime <- slingPseudotime(curves, na = FALSE)
cellWeights <- slingCurveWeights(curves)

Cluster_Labels <- tibble(Cell = names(OSNL$Simple_Lineagev2), Cluster = OSNL$Simple_Lineagev2)
UMAP <- OSNL@reductions$umap@cell.embeddings %>%
  as_tibble(rownames = "Cell")
Pseudotime <- as_tibble(pseudotime, rownames = "Cell")
LineageWeights <- as_tibble(cellWeights, rownames = "Cell")

PT_table <- Cluster_Labels %>%
  left_join(UMAP, by = "Cell") %>%
  left_join(Pseudotime, by = "Cell") %>%
  left_join(LineageWeights, by = "Cell") %>%
  pivot_longer(
    cols = starts_with("Lineage"),
    names_to = c("Lineage", "Type"),
    names_pattern = "Lineage(\\d+)\\.(x|y)") %>%
  pivot_wider(
    names_from = Type,
    values_from = value) %>%
  rename(PT_pos = x, Wt = y)

OSNL_RNA_Tcounts <- GetAssayData(object = OSNL, assay = "RNA", layer = "data") %>%
  t() %>%
  as_tibble(rownames = "Cell")

OSN_Metadata_forPT <- OSNL@meta.data %>%
  select(Genotype, Lineage, Simple_Clusters, Simple_Lineagev2, percent.10G4, percent.OR) %>%
  rownames_to_column("Cell")

set.seed(123456)
PT_table2 <- PT_table %>% 
  group_by(Cluster) %>% 
  mutate(MedianPos = median(PT_pos), SD = sd(PT_pos)) %>% 
  ungroup() %>% 
  mutate(Outlier = abs(MedianPos - PT_pos) > MedianPos) %>% 
  rowwise() %>% 
  mutate(PT_pos = ifelse(Outlier == TRUE, runif(1, 0, 2) * MedianPos, PT_pos)) #Corrected several values in Cluster A that were 
                                                                               #a product of PCA-based positioning. Eg. Clearly 
                                                                               #a GBC/INP, but within mOSN timing. 
                                                                               #Random Pseudotime within Cluster A expectations. 

PT_table_Combo <- PT_table2 %>%
  left_join(OSNL_RNA_Tcounts, by = "Cell") %>%
  left_join(OSN_Metadata_forPT, by = "Cell")

my_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "gold", "#0072B2",
  "#D55E00", "#6C79A7", "#999999", "#800000", "#66CC00", "#000080", "red", "black"
)

#48_OSNL_Pseudotime_slingshot_SmallGeneSubset
PT_table_Combo %>%
  select(Cell, PT_pos, Cluster, Simple_Clusters, Genotype, Gng8, Gnal, Gnas, Omp, Mbnl1, Gap43) %>%
  left_join(Alt_OR_forPT, by = "Cell") %>%
  pivot_longer(cols = -c(Cell, PT_pos, Cluster, Simple_Clusters, Genotype), names_to = "Gene", values_to = "LogNormCount") %>%
  group_by(Cluster) %>% 
  mutate(MedianPT = median(PT_pos)) %>% 
  group_by(Genotype, Cluster, Gene) %>% 
  mutate(MedianLNC = median(LogNormCount)) %>%
  mutate(Gene = factor(Gene, levels = c("Gnas", "Gng8", "Gap43", "EndogenousOR","Gnal", "Omp", "Mbnl1", "OR10G4"))) %>%
  ggplot() + 
  geom_point(aes(PT_pos, LogNormCount, color = Simple_Clusters), size = 3, alpha = 0.9) +
#  geom_line(aes(MedianPT, MedianLNC, color = Genotype, group = Genotype), linewidth = 2) +
#  geom_boxplot(aes(MedianPT, LogNormCount, group = interaction(Cluster, Genotype), fill = Genotype), width = 5, color = "black",
#               position = position_dodge(width = 6)) +
  geom_smooth(aes(PT_pos, LogNormCount, color = Genotype), method = "gam", 
              formula = y ~ s(x, k = 60, bs = "cs"), se = FALSE, linewidth = 2) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~Gene, ncol = 4, scales = "free_y") +
  labs(x = "Pseudotime", y = "Log Normalized Count", color = "") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 15), 
        legend.position = "top",
        legend.box.background = element_rect(color = "black", size = 0.8),
        strip.text = element_text(size = 14, face = "bold")) +
  guides(color = guide_legend(nrow = 1))

ggsave("Output2/48_OSNL_Pseudotime_slingshot_SmallGeneSubset_SF7C.png", width = 65, height = 20, units = "cm")

#49_OSNL_Pseudotime_slingshot_ORGenesubset
PT_table_Combo %>%
  select(Cell, PT_pos, Cluster, Simple_Clusters,  Genotype, Gng8, Gnal, Gnas, Omp, Mbnl1, Gap43) %>%
  left_join(Alt_OR_forPT, by = "Cell") %>%
  pivot_longer(cols = -c(Cell, PT_pos, Cluster, Simple_Clusters, Genotype), names_to = "Gene", values_to = "LogNormCount") %>% 
  filter(Gene %in% c("EndogenousOR", "OR10G4")) %>% 
  filter(!(Genotype == "WT" & Gene == "OR10G4")) %>% 
  mutate(Category = paste0(Genotype, ": ", Gene)) %>% 
  #filter(Cell != "OR10G4_TCCACCACACTTTAGG-1") %>%
  ggplot(aes(PT_pos, LogNormCount)) + 
  geom_point(aes(fill = Category), size = 2, alpha = 0.3, shape = 21, stroke = 0.5, color = "black") + 
  geom_smooth(aes(color = Category, group = Category), method = "gam", 
              formula = y ~ s(x, k = 60, bs = "cs"), se = FALSE, linewidth = 3) +
  #scale_fill_manual(values = my_colors) +
  labs(x = "Pseudotime", y = "Log Normalized Count", color = "Count Source:", fill = "Count Source:") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 15), 
        legend.position = "top",
        legend.box.background = element_rect(color = "black", size = 0.8),
        strip.text = element_text(size = 14, face = "bold")) +
  guides(color = guide_legend(nrow = 1))

ggsave("Output2/49_OSNL_Pseudotime_slingshot_ORGenesubset_MF3D.png", width = 30, height = 25, units = "cm")

#Section VI: UMAP Distance Analysis----
library(data.table)
library(Matrix)
library(Seurat)

#Running pipeline for mOSNs only with a GEP followup
Idents(OSNL) <- "Simple_Clusters"
Mature_OSNs <- subset(OSNL, idents = c("mOSN")) #Make initial subset

#Now it works... SCT was throwing things off
Mature_OSNs@assays$SCT <- NULL

#Keep only genes that are expressed (count ≥ 1) in at least one cell. ID those genes
DefaultAssay(Mature_OSNs) <- "RNA"
GenesNonZero <- rowSums(Mature_OSNs[["RNA"]]$counts > 0) > 0
GenesKeep <- tibble(Gene = names(GenesNonZero), Keep = GenesNonZero) %>% 
  filter(Keep == TRUE) %>% 
  pull(Gene) %>%
  as.vector()

counts <- GetAssayData(Mature_OSNs, assay = "RNA", layer = "counts")
counts <- counts[(which(rownames(counts) %in% GenesKeep)),]
Mature_OSNs <- subset(Mature_OSNs, features = rownames(counts)) 


#Manually rebuilding the Seurat object with some metadata and only RNA counts
idents.df <- data.frame("orig.ident" = Mature_OSNs$orig.ident, "nCount_RNA" = Mature_OSNs$nCount_RNA, 
                       "nFeature_RNA" = Mature_OSNs$nFeature_RNA, "percent.mt" = Mature_OSNs$percent.mt,
                       "Genotype" = Mature_OSNs$Genotype, "percent.10G4" = Mature_OSNs$percent.10G4,
                       "percent.OR" = Mature_OSNs$percent.OR, "S.Score" = Mature_OSNs$S.Score, 
                       "G2M.Score" = Mature_OSNs$G2M.Score, "Phase" = Mature_OSNs$Phase)

Mature_OSNs <- CreateSeuratObject(counts = LayerData(Mature_OSNs, assay = "RNA", layer = "counts"), 
                                 project = "Export", meta.data = idents.df)

Idents(Mature_OSNs) <- "Genotype"

#No need to save OR genes because already done so for OSNL, but I will remove them. 
#Previous runs suggest even using regression with the percent metadata, the genes impact clustering.
#Removing Olfr and 10G4 from features, but it looks like my original OSNL input already had them removed... Not working with a fresh H_Merge
counts <- GetAssayData(Mature_OSNs, assay = "RNA", layer = "counts")
counts <- counts[-(which(rownames(counts) %in% OR_gene_vector)),]
Mature_OSNs <- subset(Mature_OSNs, features = rownames(counts))

#Creating log1p(UMI/Total UMI * 10000) transformation for all Genes in each cell for the RNA assay to allow for appropriate visualization
DefaultAssay(Mature_OSNs) <- "RNA"
Mature_OSNs <- NormalizeData(Mature_OSNs)

#SCT and PCA
Mature_OSNs <- SCTransform(Mature_OSNs, vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA"), verbose = FALSE)
set.seed(123456)
Mature_OSNs <- RunPCA(Mature_OSNs, verbose = FALSE)

#Running RunHarmony followed by standard Clustering process using the "harmony" reduction as opposed to the usual "pca".
Mature_OSNs <- RunHarmony(Mature_OSNs, 
                   group.by.vars = c("Genotype"), 
                   reduction = "pca", 
                   assay.use = "SCT", 
                   reduction.save = "harmony")

Mature_OSNs <- RunUMAP(Mature_OSNs, reduction = "harmony", assay = "SCT", dims = 1:10) 
Mature_OSNs <- FindNeighbors(Mature_OSNs, reduction = "harmony")
Mature_OSNs <- FindClusters(Mature_OSNs, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 5.0))

#Final Decision about resolution: More!
Idents(Mature_OSNs) <- "seurat_clusters" #Max Resolution

#50_mOSNs_Single_OR_Clustering
Background <- Embeddings(Mature_OSNs, reduction = "umap") %>%
  as_tibble(rownames = "Cell")

Mature_OSNs_Meta <- Mature_OSNs@meta.data %>%
  select(seurat_clusters) %>%
  rownames_to_column("Cell") %>%
  rename(Cluster = seurat_clusters)

Map_for_ORs_mOSNs <- OSNL_Labeled_OR_Counts %>% 
  group_by(OR_Rankv2, New_Status, Simple_Clusters) %>% 
  summarize(MeanNC = mean(Normalized_Count), MedNC = median(Normalized_Count), SDNC = sd(Normalized_Count)) %>% 
  mutate(OneSD = ifelse(OR_Rankv2 == "Top", MeanNC - SDNC, MeanNC + SDNC)) %>% 
  select(OR_Rankv2, New_Status, Simple_Clusters, OneSD) %>%
  pivot_wider(names_from = OR_Rankv2, values_from = OneSD) %>%
  filter(Simple_Clusters == "mOSN") %>%
  select(-Simple_Clusters) %>%
  right_join(OSNL_Labeled_OR_Counts, by = c("New_Status")) %>%
  mutate(CountCat3 = case_when(Normalized_Count >= Top ~ "High", Normalized_Count <= 1.04 ~ "Low", TRUE ~ "Intermediate")) %>% 
  ungroup() %>% 
  group_by(Cell) %>%
  mutate(G4Detect = "10G4" %in% Symbol_copy) %>%
  ungroup() %>%
   group_by(Cell, Simple_Clusters, Genotype, G4Detect) %>%
   mutate(
       categories = list(CountCat3),
       n_high = sum(CountCat3 == "High"),
       n_trans = sum(CountCat3 == "Intermediate"),
       n_low = sum(CountCat3 == "Low"),
       n_total = n()) %>%
   mutate(Cell_Type = case_when(
       n_high >= 2 ~ "High Multiple",
       n_high >= 1 & n_trans >= 1 ~ "High/Intermediate",
       n_high >= 1 & n_low >= 1 ~ "High/Low",
       n_high == 1 & n_total == 1 ~ "High Single",
       n_trans >= 2 ~ "Intermediate Multiple",
       n_trans >= 1 & n_low >= 1 ~ "Intermediate/Low",
       n_trans == 1 & n_total == 1 ~ "Intermediate Single",
       n_low >= 2 ~ "Low Multiple",
       n_low == 1 & n_total == 1 ~ "Low Single",
       TRUE ~ "Unknown")) %>%
   mutate(Cell_OR_Cat = case_when(Cell_Type %in% c("High/Low", "High Single") ~ "High",
                                  Cell_Type %in% c("Intermediate/Low", "Intermediate Single") ~ "Intermediate",
                                  Cell_Type %in% c("Low Multiple", "Low Single") ~ "Low",
                                  TRUE ~ Cell_Type)) %>% 
  ungroup() %>% 
  mutate(MultipleDetect = n_total > 1) %>%
  filter(Simple_Clusters == "mOSN") %>%
  mutate(G4High = ifelse(Symbol_copy == "10G4" & CountCat3 == "High", TRUE, FALSE),
         G4Int = ifelse(Symbol_copy == "10G4" & CountCat3 == "Intermediate", TRUE, FALSE)) %>%
  group_by(Cell) %>%
  mutate(G4HighCell = ifelse(sum(G4High) == 1, TRUE, FALSE),
         G4IntCell = ifelse(sum(G4Int) == 1, TRUE, FALSE)) %>%
  ungroup() %>%
  left_join(Background, by = "Cell") %>%
  left_join(Mature_OSNs_Meta, by = "Cell")

Input_Map1 <- Map_for_ORs_mOSNs %>%
  filter(CountCat3 != "Low") %>%
  filter(Cell_OR_Cat %in% c("High", "Intermediate")) %>%
  group_by(Symbol_copy) %>%
  mutate(n_symbol = n()) %>%
  ungroup() %>%
  filter(n_symbol > 2) %>%
  arrange(desc(n_symbol), Symbol_copy) %>%
  mutate(Symbol_copy = factor(Symbol_copy, levels = unique(Symbol_copy))) %>%
  select(Cell, Symbol_copy) %>%
  left_join(Background, by = "Cell") %>%
  left_join(Mature_OSNs_Meta, by = "Cell")

ggplot() +
  geom_point(data = Background, aes(umap_1, umap_2), color = "grey94") +
  geom_point(data = Input_Map1, aes(umap_1, umap_2, color = Cluster), size = 3) +
  facet_wrap(~Symbol_copy, ncol = 8) +
  labs(title = "UMAP Distribution of Singular Mature OSNs expressing specific ORs in 3+ cells",
      subtitle = "Colors represent 25 seurat clusters; PCA and UMAP calculated without using OR genes; Non-Low count for expression") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

ggsave("Output2/50_mOSNs_Single_OR_Clustering_SF8A.png", width = 40, height = 20, units = "cm")

#51_mOSNs_Dimplot_byGenotype
p1 <- DimPlot(Mature_OSNs, group.by = "Genotype") + 
  theme_bw() +
  NoAxes() + 
  theme(plot.title = element_blank(), 
        legend.position = "bottom")

p2 <- DimPlot(Mature_OSNs, reduction = "umap", label = TRUE, repel = TRUE, label.box = TRUE) +
  theme_bw() + 
  NoLegend() + 
  NoAxes()

p1 + p2
ggsave("Output2/51_mOSNs_Dimplot_byGenotype_Clusters_v2_SF8B.png", width = 20, height = 11, units = "cm")

#52_mOSNs_OR10G4exp_byCluster
Background2 <- Map_for_ORs_mOSNs %>%
  select(Cell, umap_1, umap_2, Cluster) %>%
  unique()

Map_for_ORs_mOSNs %>% 
  filter(G4HighCell + G4IntCell > 0) %>% 
  ggplot()  + 
  geom_point(data = Background2, aes(umap_1, umap_2, group = Cluster), color = "grey", size = 3) + 
  geom_point(aes(umap_1, umap_2, color = Cell_OR_Cat), size = 3, alpha = 0.85) +
  labs(color = "OR10G4-positive mOSN Expression Status:") +
  facet_wrap(~Cluster, scales = "free", ncol = 9) +
  scale_color_manual(values = my_colors) + #borrowing colors from Pseudotime
  theme_bw() +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "top")

ggsave("Output2/52_mOSNs_OR10G4exp_byCluster_SF8C.png", width = 50, height = 20, units = "cm")

#53_mOSNs_MultipleOR_Distance_toClusters
Map2 <- Map_for_ORs_mOSNs %>% 
  mutate(Gene_2 = ifelse(Symbol_copy == "10G4" & CountCat3 != "Low" & Cell_OR_Cat %in% c("High", "Intermediate"), 
                         paste0(Symbol_copy, "_", Cluster), Symbol_copy))

Map_Gene_Status_G4 <- Map2 %>%
  filter(CountCat3 != "Low") %>%
  mutate(Cat_simple =  ifelse(Cell_OR_Cat %in% c("High", "Intermediate"), "Singular", "Multiple"),
         G4Present = G4HighCell + G4IntCell == 1) %>%
  group_by(Gene_2, Cat_simple, G4Present) %>%
  summarize(umap_1avg = mean(umap_1),
            umap_1sd = sd(umap_1),
            umap_2avg = mean(umap_2),
            umap_2sd = sd(umap_2),
            n_cells = n()) %>%
  rownames_to_column("row_num")

MGSG4_dist <- Map_Gene_Status_G4 %>%
  ungroup() %>%
  select(row_num, umap_1avg, umap_2avg) %>%
  column_to_rownames("row_num") %>%
  as.matrix() %>%
  dist() %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("row_num1") %>%
  pivot_longer(cols = -row_num1, names_to = "row_num2", values_to = "Distance") %>%
  filter(row_num1 < row_num2)

MGSG4_mapped <- MGSG4_dist %>%
  left_join(Map_Gene_Status_G4, by = c("row_num1" = "row_num")) %>%
  left_join(Map_Gene_Status_G4, by = c("row_num2" = "row_num"))

Distance_to_10G4 <- MGSG4_mapped %>% 
  filter(str_detect(Gene_2.x, "10G4_") | str_detect(Gene_2.y, "10G4_")) %>% 
  filter(Gene_2.x %notin% c("10G4_3", "10G4_5", "10G4_6", "10G4_8", "10G4_9", "10G4_11", "10G4_14",
                            "10G4_20", "10G4_21", "10G4_23", "10G4_24"),
         Gene_2.y %notin% c("10G4_3", "10G4_5", "10G4_6", "10G4_8", "10G4_9", "10G4_11", "10G4_14",
                            "10G4_20", "10G4_21", "10G4_23", "10G4_24")) %>%
  mutate(Gene_2.x = ifelse(str_detect(Gene_2.x, "10G4_"), "10G4_s", Gene_2.x), 
         Gene_2.y = ifelse(str_detect(Gene_2.y, "10G4_"), "10G4_s", Gene_2.y)) %>%
  mutate(Gene_x = ifelse(str_detect(Gene_2.y, "10G4_"), Gene_2.y, Gene_2.x),
         Gene_y = ifelse(str_detect(Gene_2.y, "10G4_"), Gene_2.x, Gene_2.y),
         Cat_x = ifelse(str_detect(Gene_2.y, "10G4_"), Cat_simple.y, Cat_simple.x),
         Cat_y = ifelse(str_detect(Gene_2.y, "10G4_"), Cat_simple.x, Cat_simple.y),
         G4_x = ifelse(str_detect(Gene_2.y, "10G4_"), G4Present.y, G4Present.x),
         G4_y = ifelse(str_detect(Gene_2.y, "10G4_"), G4Present.x, G4Present.y)) %>%
  group_by(Gene_x, Cat_x, G4_x, Gene_y, Cat_y, G4_y) %>% 
  summarize(min(Distance))

Distance_to_10G4 %>% 
  ungroup() %>% 
  select(4:7) %>% 
  mutate(G4_y = ifelse(G4_y == TRUE, "G4", "X"), Cat_G4 = paste(Cat_y, G4_y, sep = "_")) %>% 
  select(-Cat_y, -G4_y) %>% 
  pivot_wider(names_from = Cat_G4, values_from = `min(Distance)`) %>% 
  filter(!is.na(Singular_X)) %>% 
  select(-Singular_G4) %>% 
  filter(!is.na(Multiple_G4) | !is.na(Multiple_X)) %>% 
  mutate(EndoOR_Shift = Multiple_X/Singular_X, G4_Shift = Multiple_G4/Singular_X) %>% 
  filter(!is.na(G4_Shift)) %>% 
  arrange(G4_Shift) %>% 
  mutate(Gene_y = factor(Gene_y, levels = Gene_y)) %>% 
  ggplot(aes(Gene_y, log2(G4_Shift))) + 
  geom_point(aes(color = Singular_X, size = Singular_X)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Log2 Ratio of Co-expressed/Singular-OR Distance to nearest OR10G4 cluster", 
       title = "A subset of ORs coexpressed with OR10G4 correspond with a positional shift towards singular OR10G4 clusters",
       color = "Singular-OR mOSN Distance to nearest OR10G4 Cluster",
       size = "Singular-OR mOSN Distance to nearest OR10G4 Cluster") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "top") +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 3) +
  guides(color= guide_legend(), size=guide_legend())

ggsave("Output2/53_mOSNs_MultipleOR_Distance_toClusters_SF8D.png", width = 40, height = 15, units = "cm")

#Section VII: The actual GEP analysis----
#54_mOSNs_GEP_UMAP
counts <- Mature_OSNs@assays$SCT$counts
barcodes <- colnames(counts)
gene_names <- rownames(counts)


#Identify two directories, somewhat excessive, but otherwise more code editing
data_dir = '../GEP_mOSN2/'

if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

filtered_dir = '../GEP_mOSN2/data/'

if (!dir.exists(filtered_dir)) {
  dir.create(filtered_dir, recursive = TRUE)
}

# Output counts matrix
writeMM(counts, paste0(filtered_dir, 'matrix.mtx'))

# Output cell barcodes
barcodes <- colnames(counts)
write.table(as.data.frame(barcodes), paste0(filtered_dir, 'barcodes.tsv'),
            col.names = FALSE, row.names = FALSE, sep = "\t")

# Output feature names
gene_names <- rownames(counts)
features <- data.frame("gene_id" = gene_names,"gene_name" = gene_names,type = "Gene Expression")
write.table(as.data.frame(features), sep = "\t", paste0(filtered_dir, 'genes.tsv'),
            col.names = FALSE, row.names = FALSE)

#Adjust to your local conda install
use_condaenv("cnmf_env",
             conda = "C:/Users/Feinstein Lab/anaconda3/Scripts/conda.exe",
             required = TRUE)

Sys.setenv(PATH = paste(
  "C:/Users/Feinstein Lab/anaconda3/envs/cnmf_env/Scripts",
  Sys.getenv("PATH"),
  sep=";"
))

#Does a thing with the matrix data and some graphs
runname = "mOSN_cNMF"
cmd = paste("cnmf prepare --output-dir", data_dir,
            "--name", runname,
            "-c", paste0(filtered_dir, 'matrix.mtx'),
            "--max-nmf-iter 2000", 
            "-k 4 8 10 12 18 24 --n-iter 20", sep=" ")
print(cmd)
system(cmd)

cmd = paste("cnmf factorize --output-dir", data_dir,
            "--name", runname,
            "--worker-index 0 --total-workers 1", sep=" ")
print(cmd)
system(cmd)

cmd = paste("cnmf combine --output-dir", data_dir,
            "--name", runname, sep=" ")
print(cmd)
system(cmd)

cmd = paste("cnmf k_selection_plot --output-dir", data_dir,
            "--name", runname, sep=" ")
print(cmd)
system(cmd)

#Edit the components number based on selection plot from the output directory.
cmd = paste("cnmf consensus --output-dir", data_dir,
            "--name", runname,
            '--components', 24,
            '--local-density-threshold', 0.1,
            '--show-clustering', sep=" ")
print(cmd)
system(cmd)

#Does more things with things. Requires editing when choosing different K.
usage_file <- paste(data_dir[1:length(data_dir)], runname, paste(runname, "usages", "k_24.dt_0_1", 'consensus', 'txt', sep="."), sep="/")
spectra_score_file <- paste(data_dir[1:length(data_dir)], runname, paste(runname, "gene_spectra_score", "k_24.dt_0_1", 'txt', sep="."), sep="/")
spectra_tpm_file <- paste(data_dir[1:length(data_dir)], runname, paste(runname, "gene_spectra_tpm", "k_24.dt_0_1", 'txt', sep="."), sep="/")

usage <- read.table(usage_file, sep='\t', row.names=1, header=TRUE)
spectra_score <- read.table(spectra_score_file, sep='\t', row.names=1, header=TRUE)
spectra_tpm <- read.table(spectra_tpm_file, sep='\t', row.names=1, header=TRUE)
head(usage)

usage_norm <- as.data.frame(t(apply(usage, 1, function(x) x / sum(x))))
usage_norm_24 <- usage_norm

new_metadata <- merge(Mature_OSNs@meta.data, usage_norm_24, by = "row.names", all.x = TRUE)

new_metadata <- new_metadata %>% select(-1) #Use if Row.names is duplicated

#GEP names based on evaluation of top scoring Genes in comparison with published OSN GEPs and GOEA interpretation. 
#Subject to inaccuracy and might change
new_metadata <- new_metadata %>%
  dplyr::rename(Terminal_Identity = X1,
                Mito_Riboprotein = X2,
                RNA_Processing = X3,
                OxPhos_Histone = X4,
                Posterior = X5,
                Elevated_Homeostasis = X6,
                Dorsal = X7,
                Mito_Demand = X8,
                Volatility = X9,
                Ventral = X10,
                Signaling_Primed = X11,
                UPR_ER_Stress = X12,
                Anterior = X13,
                Low_Activity = X14, 
                High_Activity = X15,
                ProSynapse = X16,
                Immature = X17,
                Cd36_Subset = X18,
                Trpc2_TypeB_Gucy = X23)
rownames(new_metadata) <- new_metadata$Row.names
Mature_OSNs@meta.data <- new_metadata

GEPorder <- c("Dorsal", "Ventral", "Anterior", "Posterior", "High_Activity", "Low_Activity", "Signaling_Primed",
              "Elevated_Homeostasis", "Volatility", "UPR_ER_Stress",
              "Terminal_Identity", "Mito_Riboprotein", "RNA_Processing", 
              "OxPhos_Histone", "Immature", "Mito_Demand", "ProSynapse", "Cd36_Subset")

FeaturePlot(Mature_OSNs, features = c(GEPorder, 
                                           "Nfix", "Acsm4", "S100a5", "Casp3", "Cidea", "percent.10G4", "Gap43",
                                      "percent.mt", "nFeature_RNA"), 
                 combine=T, ncol = 9, order = TRUE) &
  NoAxes() &
  NoLegend()

ggsave("Output2/54_mOSNs_GEP_UMAP_SF9A.png", width = 50, height = 20, units = "cm")

#55_mOSNs_GEP_DVScore
GEP_mOSNs <- Mature_OSNs@meta.data %>%
  select(GEPorder) %>%
  rownames_to_column("Cell") %>%
  right_join(Map2, by = "Cell") %>%
  filter(CountCat3 != "Low") %>%
  filter(Cell_OR_Cat %in% c("High", "Intermediate")) %>%
  group_by(Gene_2) %>%
  summarize(across(Dorsal:Cd36_Subset, \(x) mean(x, na.rm = TRUE))) 
  
Mature_OSNs@meta.data %>%
  select(GEPorder) %>%
  rownames_to_column("Cell") %>%
  right_join(Map2, by = "Cell") %>%
  filter(CountCat3 != "Low") %>%
  filter(Cell_OR_Cat %in% c("High", "Intermediate")) %>%
  mutate(DV_Score = Dorsal-Ventral) %>%
 # filter(!str_detect(Gene_2, "10G4")) %>%
  group_by(Gene_2) %>%
  mutate(GeneMean = mean(DV_Score), GeneSD = sd(DV_Score)) %>%
  ungroup() %>%
  mutate(GeneCat = ifelse(str_detect(Gene_2, "10G4"), "Singular OR10G4", "Singular OR")) %>%
  arrange(GeneMean) %>%
  mutate(Gene_2 = factor(Gene_2, levels = unique(Gene_2))) %>%
  ggplot(aes(Gene_2, DV_Score)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_point(aes(color = GeneCat)) +
  geom_errorbar(aes(ymin = GeneMean - GeneSD, ymax = GeneMean + GeneSD), width = 0.2) +
  labs(y = "DV Score (Dorsal GEP - Ventral GEP)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        axis.title.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.9),
        legend.title = element_blank())

ggsave("Output2/55_mOSNs_GEP_DVScoreV2_SF9F.png", width = 60, height = 10, units = "cm")

#56_mOSNs_GEP_Correlation

mOSN_Matrix <- as.matrix(GetAssayData(object = Mature_OSNs, assay = "RNA", layer = "data"))
mOSN_Tibble <- as_tibble(rownames = "Symbol", mOSN_Matrix) %>%
  pivot_longer(cols = -Symbol, names_to = "Cell", values_to = "NC") %>%
  pivot_wider(names_from = Symbol, values_from = "NC")
  
mOSN_T_Select <- mOSN_Tibble %>%
  select(Cell, Nfix, Acsm4, Gap43, Cidea, Casp3)

mOSN_T_Select <- Mature_OSNs@meta.data %>%
  select(-1) %>%
  rownames_to_column("Cell") %>%
  select(Cell, all_of(GEPorder), percent.10G4) %>%
  left_join(mOSN_T_Select, by = "Cell")

GEP_cor <- cor(dplyr::select(mOSN_T_Select, -Cell), use = "pairwise.complete.obs")

library(corrplot)
png("Output2/56_mOSNs_GEP_Correlation_SF9B.png", width = 2700, height = 2700, res = 300)
corrplot(GEP_cor, method = 'circle', type = 'lower', addCoef.col ='black', number.cex = 0.6, diag=FALSE, mar = c(0,0,0,0),
         tl.col = 'black', tl.srt = 45, tl.cex = 0.6)
dev.off()

#57A_mOSNs_DV_Score_Cell_Cluster_compare
Mature_OSNs@meta.data %>%
  select(GEPorder) %>%
  rownames_to_column("Cell") %>%
  mutate(DV_Score = Dorsal - Ventral) %>%
  arrange(DV_Score) %>%
  mutate(Cell = factor(Cell, levels = unique(Cell))) %>%
  ggplot(aes(Cell, DV_Score)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_point() +
  labs(y = "DV Score (Dorsal GEP - Ventral GEP)", x = "Individual mOSNs") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave("Output2/57A_mOSNs_DV_Score_Cell_Cluster_compare_SF9C.png", width = 10, height = 7.5, units = "cm")

Mature_OSNs@meta.data %>%
  select(GEPorder, seurat_clusters) %>%
  rownames_to_column("Cell") %>%
  mutate(DV_Score = Dorsal - Ventral) %>%
  group_by(seurat_clusters) %>%
  mutate(ClusterMean = mean(DV_Score), ClusterSD = sd(DV_Score)) %>%
  ungroup() %>%
  arrange(ClusterMean) %>%
  mutate(seurat_clusters = factor(seurat_clusters, levels = unique(seurat_clusters))) %>%
  ggplot(aes(seurat_clusters, DV_Score)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_point() +
  geom_errorbar(aes(ymin = ClusterMean - ClusterSD, ymax = ClusterMean + ClusterSD), width = 0.3) +
  labs(y = "DV Score (Dorsal GEP - Ventral GEP)", x = "Seurat Clusters") +
  theme_bw()
  
ggsave("Output2/57B_mOSNs_DV_Score_Cell_Cluster_compare_SF9D.png", width = 10, height = 7.5, units = "cm")

MO_G4Score <- Mature_OSNs@meta.data %>%
  select(GEPorder) %>%
  rownames_to_column("Cell") %>%
  right_join(Map2, by = "Cell") %>%
 # filter(Cell_OR_Cat %in% c("High", "High/Intermediate", "Multiple High")) %>%
  filter(CountCat3 != "Low") %>%
  filter(!(!str_detect(Symbol_copy, "10G4") & CountCat3 == "Low")) %>%
  filter(!(!str_detect(Symbol_copy, "10G4") & CountCat3 == "Intermediate")) %>%
  select(Cell, Elevated_Homeostasis, Volatility, UPR_ER_Stress, Normalized_Count, CountCat3, Cell_OR_Cat, Symbol_copy) %>%
  mutate(Gene = ifelse(str_detect(Symbol_copy, "10G4"), "OR10G4", "OR")) %>%
  select(-Symbol_copy, -CountCat3) %>%
  mutate(Normalized_Count = exp(Normalized_Count) - 1) %>%
  pivot_wider(names_from = c(Gene), values_from = Normalized_Count, values_fill = 0) %>%
  mutate(Tension_Score = Elevated_Homeostasis + Volatility + UPR_ER_Stress) %>%
  select(-c(Elevated_Homeostasis, Volatility, UPR_ER_Stress), -Cell_OR_Cat) %>%
  mutate(Category = case_when(OR == 0 ~ "Only OR10G4",
                              OR10G4 == 0 ~ "Only OR",
                              TRUE ~ "OR/OR10G4")) %>%
  mutate(Category = factor(Category, levels = c("Only OR", "OR/OR10G4", "Only OR10G4"))) 

compare_means(Tension_Score ~ Category, data = MO_G4Score, method = "t.test")

MO_G4Score %>%
  ggplot(aes(Category, Tension_Score)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0) +
  geom_signif(
    comparisons = list(c("Only OR", "OR/OR10G4"), c("Only OR", "Only OR10G4"), c("OR/OR10G4", "Only OR10G4")),
    map_signif_level = FALSE, textsize = 5, y_position = c(0.625, 0.725, 0.675)) +
  labs(y = "Tension Score (Elevated + Volatility + Stress)") +
  theme_bw() +
  theme(axis.title.x = element_blank())

ggsave("Output2/59_mOSNs_GEP_G4Score_Compare_SF9H.png", width = 15, height = 15, units = "cm")

#60_mOSNs_GEP_G4Score_OR_Compare
MO_GEP_G4Score <- Mature_OSNs@meta.data %>%
  select(GEPorder) %>%
  rownames_to_column("Cell") %>%
  right_join(Map2, by = "Cell") %>%
  # filter(Cell_OR_Cat %in% c("High", "High/Intermediate", "Multiple High")) %>%
  filter(CountCat3 != "Low") %>%
  mutate(Tension_Score = Elevated_Homeostasis + Volatility + UPR_ER_Stress) %>%
  # filter(str_detect(Gene_2, "10G4")) %>%
  group_by(Gene_2) %>%
  mutate(GeneMean = mean(Tension_Score), GeneSD = sd(Tension_Score)) %>%
  ungroup() %>%
  mutate(G4Presence = ifelse(G4HighCell + G4IntCell > 0, "G4P", "NoG4")) %>%
  mutate(Cat_Simple = ifelse(Cell_OR_Cat %in% c("High", "Intermediate"), "Singular", "Multiple")) %>%
  group_by(Gene_2, Cat_Simple, G4Presence) %>%
  summarize(MG4 = mean(Tension_Score)) %>%
  pivot_wider(names_from = Cat_Simple, values_from = MG4) %>%
  pivot_wider(names_from = G4Presence, values_from = c(Multiple, Singular)) %>%
  filter(!is.na(Multiple_G4P) & !is.na(Singular_NoG4)) %>%
  select(Gene_2, Multiple_G4P, Singular_NoG4) %>%
  ungroup() %>%
  arrange(desc(Singular_NoG4)) %>%
  mutate(Gene_2 = factor(Gene_2, levels = unique(Gene_2))) %>%
  pivot_longer(cols = -Gene_2, names_to = "mOSN_Category", values_to = "Tension_Score") %>%
  mutate(mOSN_Category = ifelse(mOSN_Category == "Singular_NoG4", "Singular OR", "OR/OR10G4")) 

MO_GEP_G4Score %>%
  ggplot(aes(Gene_2, Tension_Score)) +
  geom_point(aes(color = mOSN_Category), size = 3) +
  geom_line(aes(group = Gene_2)) +
  labs(y = "Tension_Score (Elevated + Volatility + Stress)", color = "mOSN Expresses:") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top")

ggsave("Output2/60_mOSNs_GEP_G4Score_OR_Compare_SF9I.png", width = 40, height = 15, units = "cm")
  
#58_mOSNs_GEP_G4Score  
Mature_OSNs@meta.data %>%
  select(GEPorder) %>%
  rownames_to_column("Cell") %>%
  right_join(Map2, by = "Cell") %>%
  filter(CountCat3 != "Low") %>%
  filter(Cell_OR_Cat %in% c("High", "Intermediate")) %>%
  mutate(Tension_Score = Elevated_Homeostasis + Volatility + UPR_ER_Stress) %>%
  # filter(!str_detect(Gene_2, "10G4")) %>%
  group_by(Gene_2) %>%
  mutate(GeneMean = mean(Tension_Score), GeneSD = sd(Tension_Score)) %>%
  ungroup() %>%
  mutate(GeneCat = ifelse(str_detect(Gene_2, "10G4"), "Singular OR10G4", "Singular OR")) %>%
  arrange(GeneMean) %>%
  mutate(Gene_2 = factor(Gene_2, levels = unique(Gene_2))) %>%
  ggplot(aes(Gene_2, Tension_Score)) +
  geom_point(aes(color = GeneCat)) +
  geom_errorbar(aes(ymin = GeneMean - GeneSD, ymax = GeneMean + GeneSD), width = 0.2) +
  labs(y = "Tension Score (Elevated + Volatility + Stress)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        axis.title.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.9),
        legend.title = element_blank())

ggsave("Output2/58_mOSNs_GEP_G4Score_SF9G.png", width = 70, height = 10, units = "cm")

MO_GEP_DVScore <- Mature_OSNs@meta.data %>%
  dplyr::select(GEPorder) %>%
  rownames_to_column("Cell") %>%
  right_join(Map2, by = "Cell") %>%
  # filter(Cell_OR_Cat %in% c("High", "High/Intermediate", "Multiple High")) %>%
  dplyr::filter(CountCat3 != "Low") %>%
  mutate(DV_Score = Dorsal - Ventral) %>%
  # filter(str_detect(Gene_2, "10G4")) %>%
  dplyr::group_by(Gene_2) %>%
  mutate(GeneMean = mean(DV_Score), GeneSD = sd(DV_Score)) %>%
  ungroup() %>%
  mutate(G4Presence = ifelse(G4HighCell + G4IntCell > 0, "G4P", "NoG4")) %>%
  mutate(Cat_Simple = ifelse(Cell_OR_Cat %in% c("High", "Intermediate"), "Singular", "Multiple")) %>%
  dplyr::group_by(Gene_2, Cat_Simple, G4Presence) %>%
  summarize(MDV = mean(DV_Score)) %>%
  pivot_wider(names_from = Cat_Simple, values_from = MDV) %>%
  pivot_wider(names_from = G4Presence, values_from = c(Multiple, Singular)) %>%
  dplyr::filter(!is.na(Multiple_G4P) & !is.na(Singular_NoG4)) %>%
  dplyr::select(Gene_2, Multiple_G4P, Singular_NoG4) %>%
  ungroup() %>%
  dplyr::arrange(desc(Singular_NoG4)) %>%
  mutate(Gene_2 = factor(Gene_2, levels = unique(Gene_2))) %>%
  pivot_longer(cols = -Gene_2, names_to = "mOSN_Category", values_to = "DV_Score") %>%
  mutate(mOSN_Category = ifelse(mOSN_Category == "Singular_NoG4", "Singular OR", "OR/OR10G4")) 

#61_mOSNs_GEP_DVScore_OR_Compare
MO_GEP_DVScore %>%
  ggplot(aes(Gene_2, DV_Score)) +
  geom_point(aes(color = mOSN_Category), size = 3) +
  geom_line(aes(group = Gene_2)) +
  labs(y = "DV Score (Dorsal GEP - Ventral GEP)", color = "mOSN Expresses:") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top")

ggsave("Output2/61_mOSNs_GEP_DVScore_OR_Compare_SF9F.png", width = 30, height = 15, units = "cm")


