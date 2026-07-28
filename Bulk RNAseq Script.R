# This is a script to generate DEG Figures and Tables for the 10G4 Paper.

# Set your working directory to the "DEG 10G4" folder, the folder with all the imports.

#01 - Required Packages-----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
 
BiocManager::install("BiocFileCache")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("vsn")
BiocManager::install("biomaRt")

library(DESeq2)
library(apeglm)
library(vsn)
library(biomaRt)

if (!require(pheatmap)) install.packages('pheatmap')
library(pheatmap)

if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

if (!require(patchwork)) install.packages('patchwork')
library(patchwork)

if (!require(ggrepel)) install.packages('ggrepel')
library(ggrepel)

if (!require(RColorBrewer)) install.packages('RColorBrewer')
library(RColorBrewer)

if (!require(scales)) install.packages('scales')
library(scales)  

if (!require(broom)) install.packages('broom')
library(broom)  

library(ggh4x)
library(viridis)

`%notin%` <- Negate(`%in%`)

#02 - Import the raw data, prepare the dataset, run DESeq() and results() functions---- 
count_table <- read_delim("Inputs/raw_counts_genes.clean.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

#Extract the gene ids to use as rownames instead
genes <- count_table$gene_id

#Select specific sample columns and remove gene_id column. I'm ignoring the OB samples. Make the counts a matrix and add rownames
cts <- count_table[, c(2,5,6,8,3,4,7)]
cts <- as.matrix(cts)
rownames(cts) <- genes

#Making the coldata object with the condition of the samples, samples as rownames, and condition as a factor
coldata <- data.frame(condition = c("WT", "WT", "WT", "WT", "OR10G4", "OR10G4", "OR10G4"))
rownames(coldata) <- colnames(cts)
coldata$condition <- as.factor(coldata$condition)

#A check on matching the rownames and colnames, aka the samples. Should be true.
all(rownames(coldata) == colnames(cts))

#The rest of the process is kinda simple enough. The (experimental) design here is simple with only one feature to select from the coldata table: condition
dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = coldata, 
                              design = ~ condition)

#We will apply a small pre-filtering step to the data, which in this case will take dds from 49k to about 30k genes. 
#It only removes genes where the total from all samples is 10 or fewer counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Here the untreated condition is established as the reference
dds$condition <- relevel(dds$condition, ref = "WT")

#run the DESeq function, the results function, a LFC shrinkage function, and a variance-reducing function.
dds <- DESeq(dds)

res <- results(dds, alpha=0.05)

#There are alot of missing p values or padj values. Exploration confirms them to be the genes with low counts or genes where one sample of a condition is an outlier.
summary(res)

Missing <- as.data.frame(res) %>% 
  rownames_to_column("GeneID") %>%
  filter(is.na(padj))

MissingCounts <- count_table %>%
  rename(GeneID = gene_id) %>%
  inner_join(Missing, by = "GeneID")

#03 - Transforming data for visualization and QA/QC plots NO SAVED FILES----
#Shrinking the LFC values based on relative SE values for visualization
resLFC <- lfcShrink(dds, coef="condition_OR10G4_vs_WT", type="apeglm", res = res) 

#Transforming counts for visualization purposes
vsd <- vst(dds, blind=FALSE) 
ntd <- normTransform(dds) 

#Quality-testing plots; not for publication
#Does the Variance change with the mean?
meanSdPlot(assay(ntd)) 
meanSdPlot(assay(vsd)) 

#What's the relationship between L2FC and normalized counts and p? How does the shrunk LFC compare?
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

#Looks at the norm counts of the gene with the smallest padj: OR10G4. NOTE: WT have counts, which must represent noise. Baseline for what counts to ignore.
plotCounts(dds, gene=which.min(res$padj), intgroup="condition") 

#The majority of the highest count genes do not show condition-specific differences
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("condition","sizeFactor")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)

#Clear clustering evident with distance measure
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)

#PC1 shows 62% variance and very cleanly splits WT and OR10G4
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#This plot and summary show that there are not many outliers in the dataset and the "Cooks" values are similar for each sample
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
summary(res)

#Dispersion plot shows an ideal result: dispersion decreases with increasing mean
plotDispEsts(dds)

#04 - Annotation with biomaRt----
#What sets are available?
listMarts()
listEnsemblArchives()

#Set up connection to ensembl database. Use an archived version to access GRCm39, version 105, not 116
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://dec2021.archive.ensembl.org")

listDatasets(ensembl) %>% 
  filter(str_detect(description, "Mouse"))

#Specify a data set to use
ensembl <- useDataset("mmusculus_gene_ensembl", mart=ensembl)
listFilters(ensembl) %>% 
  filter(str_detect(name, "ensembl"))

#Set the filter type and values
ourFilterType <- "ensembl_gene_id"
filterValues <- rownames(res)

#Check the available "attributes" - things you can retrieve
listAttributes(ensembl) %>% 
  head(50)

#Set the list of attributes
attributeNames <- c('ensembl_gene_id', 'external_gene_name', 'description', 'chromosome_name', 'start_position', 'end_position')

#Run the query. This may time out or get an error. If you are savvy, try a mirror or figure out a manual method.
annot <- getBM(attributes = attributeNames, filters = ourFilterType, values = filterValues, mart = ensembl)

#What's missing post-annotation? Lots with the current version, so I went back and picked 105, not 116
print(filterValues[(filterValues %notin% annot$ensembl_gene_id)])

Missing <- filterValues[(filterValues %notin% annot$ensembl_gene_id)]

#With Chromosome Zero, I need to add 10G4 to the annot table manually.
First <- data.frame('ensembl_gene_id' = "PFMUSG0001", 
                    'external_gene_name' = "OR10G4", 
                    'description' = "10G4 Transgene", 
                    'chromosome_name' = "3", 
                    'start_position' = 60436556, 
                    'end_position' = 60487531)


#Weird Gene that doesn't have a home... Easier to nullify and remove
MissingOne <- data.frame('ensembl_gene_id' = "ENSMUSG00000106198",
                     'external_gene_name' = NA,
                     'description' = NA,
                     'chromosome_name' = NA,
                     'start_position' = NA,
                     'end_position' = NA)

annot <- bind_rows(annot, MissingOne)
annot <- annot %>% arrange(ensembl_gene_id)
annot <- bind_rows(First, annot)

annot_backup <- annot

#Prepare the annotation
all(rownames(dds) == annot$ensembl_gene_id)
colnames(annot) <- c("GeneID", "Symbol", "Description", "Chromosome", "Start", "End")

#05 - Combine the annotations with results and create some subsets for export----
res_annot <- as.data.frame(res) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(annot, "GeneID") %>%
  arrange(pvalue)
write_csv(res_annot, file="Output/DESeq10G4_30kgenes.csv")

OlfrOnly <- res_annot[str_detect(res_annot$Symbol, "Olfr"), ]
write_csv(OlfrOnly, file="Output/DESeq10G4_Olfrgenes.csv")

NotOlfr <- res_annot[str_detect(res_annot$Symbol, "Olfr", negate = TRUE), ]
write_csv(NotOlfr, file="Output/DESeq10G4_NOTolfrgenes.csv")

TaarOnly <- res_annot[str_detect(res_annot$Symbol, "Taar"), ]
write_csv(TaarOnly, file="Output/DESeq10G4_Taargenes.csv")

#06 - Import Greek Island info (modified for mm39 positional info) & Dorsal/Ventral Index----
Greek_Islands <- read_table("Inputs/GImm39_Input", col_names = FALSE) %>%
  rename(Chromosome = X1, Start_GI = X2, End_GI = X3, Greek_Island = X4)

Greek_Islands2 <- read_table("Inputs/GImm39_Input2", col_names = FALSE) %>%
  rename(Chromosome = X1, Start_GI = X2, End_GI = X3, Greek_Island = X4) %>%
  mutate(Chromosome = str_remove(Chromosome, "chr"))

Greek_Islands2 <- Greek_Islands2 %>%
  bind_rows(tibble(Chromosome = "7", Start_GI = 102158975, End_GI = 102160824, Greek_Island = "J Element"))

DVI <- read_csv("Inputs/DVI.csv")
  
#07 - Generating combined datasets for plotting----
MeanCounts <- count_table %>%
  dplyr::rename(GeneID = gene_id) %>%
  mutate(WTmeanCount = (`4163-M` + `4166-M` + `4169-M` + `4171-M`)/4,
         MutantmeanCount = (`4164-M` + `4165-M` + `4170-M`)/3) %>%
  dplyr::select(GeneID, WTmeanCount, MutantmeanCount)

a1 <- res_annot %>%
  mutate(GeneType = case_when(str_detect(Symbol, "Olfr") ~ "Olfr", 
                              str_detect(Symbol, "Taar") ~ "Taar",
                              str_detect(Symbol, "Ms4a") ~ "Ms4a",
                              TRUE ~ "Other"),
         Genes = case_when(GeneType == "Olfr" & padj < 0.05 & abs(log2FoldChange) > 1 ~ "Olfr DEG",
                           GeneType == "Taar" & padj < 0.05 & abs(log2FoldChange) > 1 ~ "Taar DEG",
                           GeneType %in% c("Other", "Ms4a") & padj < 0.05 & abs(log2FoldChange) > 1 ~ "DEG", 
                           TRUE ~ "non-DEG"),
         Corrected_padj = padj + subset(res_annot, Symbol == "Olfr711")$padj,
         Genes = factor(Genes, levels = c("DEG", "Olfr DEG", "Taar DEG", "non-DEG"))) %>% #for graphing purposes
  left_join(MeanCounts, by = "GeneID") %>%
  left_join(DVI, by = "Symbol") %>%
  mutate(DVI = factor(DVI, levels = c("1.05", "1.25", "1.4", "1.5", "1.6", "1.7", "1.8", "1.9", 
                                      "2", "2.1", "2.2", "2.3", "2.4", "2.5", "2.6", "2.7", "2.8", "2.9", 
                                      "3.05", "3.2", "3.3", "3.4", "3.5", "3.6", "3.7", "3.8", "3.9", "4",
                                      "4.1", "4.2", "4.3", "4.4", "4.5", "4.6", "4.7", "4.8", "4.9", "5",
                                      "low expression", "unusual")),
         Class = case_when(GeneType == "Olfr" & Chromosome == 7 & Start > 102000000 & End < 106000000 ~ "Class I",
                           GeneType == "Olfr" ~ "Class II",
                           GeneType == "Taar" ~ "Class TAAR")) %>%
  left_join(Greek_Islands2, by = "Chromosome", relationship = "many-to-many") %>% #This line replaced the original Greek_Islands with the updated version and the code was rerun to make a4
  rowwise() %>% 
  mutate(GI_Distance = min(c(abs(End - End_GI), abs(Start - End_GI), abs(End - Start_GI), abs(Start - Start_GI)))) %>%
  ungroup() %>% 
  arrange(Greek_Island, GI_Distance) %>% 
  group_by(Greek_Island) %>% 
  mutate(GI_Prox_Rank = rank(GI_Distance, ties.method = "first")) %>%
  ungroup()

a2 <- a1 %>% 
  filter(GeneType == "Olfr") %>% 
  arrange(Greek_Island, GI_Distance) %>% 
  group_by(Greek_Island) %>%  
  mutate(GI_Prox_RankvsOlfrs = rank(GI_Distance, ties.method = "first")) %>%
  ungroup() %>%
  dplyr::select(Symbol, Greek_Island, GI_Prox_RankvsOlfrs) %>%
  right_join(a1, by = c("Symbol", "Greek_Island"))


a3 <- a2 %>%
  group_by(Symbol) %>%
  mutate(SmallestGIDistance = min(GI_Distance), 
         SmallestGIRankvsOlfrs = min(GI_Prox_RankvsOlfrs)) %>%
  ungroup() %>%
  group_by(Symbol) %>%
  arrange(GI_Distance) %>%
  slice_head() %>%
  ungroup() %>%
  dplyr::select(GeneID, Greek_Island, SmallestGIDistance) %>%
  dplyr::rename(Closest_GI = Greek_Island) %>%
  right_join(a2, by = "GeneID") %>%
  mutate(GI_Distance_Category = case_when(GI_Distance < 100000 ~ "<100k",
                                          GI_Distance < 500000 ~ "100k-499999",
                                          GI_Distance < 1000000 ~ "500k-999999",
                                          GI_Distance < 2000000 ~ "1kk-1999999",
                                          GI_Distance >= 2000000 ~ "2kk+"),
         Closest_GI_Distance_Category = case_when(SmallestGIDistance < 100000 ~ "<100k",
                                                  SmallestGIDistance < 500000 ~ "100k-499999",
                                                  SmallestGIDistance < 1000000 ~ "500k-999999",
                                                  SmallestGIDistance < 2000000 ~ "1kk-1999999",
                                                  SmallestGIDistance >= 2000000 ~ "2kk+")) %>%
  mutate(Closest_GI_Distance_Category = factor(Closest_GI_Distance_Category, levels = c("<100k", "100k-499999", "500k-999999", "1kk-1999999", "2kk+")),
         GI_Distance_Category = factor(GI_Distance_Category, levels = c("<100k", "100k-499999", "500k-999999", "1kk-1999999", "2kk+"))) %>% 
  mutate(WTmeanCount_Category = case_when(WTmeanCount < 10 ~ "<10",
                                          WTmeanCount < 50 ~ "10-49",
                                          WTmeanCount < 100 ~ "50-99",
                                          WTmeanCount < 250 ~ "100-249",
                                          WTmeanCount < 500 ~ "250-499",
                                          WTmeanCount < 1000 ~ "500-999",
                                          WTmeanCount < 5000 ~ "1k-4999",
                                          WTmeanCount < 10000 ~ "5k-9999",
                                          WTmeanCount >= 10000 ~ "10k+", 
                                          TRUE ~ "No_Counts_Detected"),
         MutantmeanCount_Category = case_when(MutantmeanCount < 10 ~ "<10",
                                              MutantmeanCount < 50 ~ "10-49",
                                              MutantmeanCount < 100 ~ "50-99",
                                              MutantmeanCount < 250 ~ "100-249",
                                              MutantmeanCount < 500 ~ "250-499",
                                              MutantmeanCount < 1000 ~ "500-999",
                                              MutantmeanCount < 5000 ~ "1k-4999",
                                              MutantmeanCount < 10000 ~ "5k-9999",
                                              MutantmeanCount >= 10000 ~ "10k+", 
                                              TRUE ~ "No_Counts_Detected")) %>%
  mutate(WTmeanCount_Category = factor(WTmeanCount_Category, levels = c("<10", "10-49", "50-99", "100-249", "250-499", "500-999", "1k-4999", "5k-9999", "10k+", "No_Counts_Detected")),
         MutantmeanCount_Category = factor(MutantmeanCount_Category, levels = c("<10", "10-49", "50-99", "100-249", "250-499", "500-999", "1k-4999", "5k-9999", "10k+", "No_Counts_Detected")))

a4 <- a3 %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                    "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                                    "X", "Y", "GL456211.1", "GL456221.1", "JH584295.1", 
                                                    "GL456212.1", "GL456210.1", "JH584304.1"))) %>%
  arrange(Chromosome, Start) %>%
  mutate(Symbol = factor(Symbol, levels = unique(Symbol))) %>%
  mutate(Result = case_when(Genes == "Olfr DEG" & Class == "Class I" ~ "Class I DEG",
                            Genes == "non-DEG" & Class == "Class I" ~ "Class I non-DEG",
                            Genes == "Olfr DEG" & Class == "Class II" ~ "Class II DEG",
                            Genes == "non-DEG" & Class == "Class II" ~ "Class II non-DEG",
                            Genes == "Taar DEG" & Class == "Class TAAR" ~ "Taar DEG",
                            Genes == "non-DEG" & Class == "Class TAAR" ~ "Taar non-DEG",
                            TRUE ~ "Other")) %>%
  mutate(Result = factor(Result, levels = c("Class I DEG", "Class II DEG", "Taar DEG", "Class I non-DEG", "Class II non-DEG", "Taar non-DEG", "Other"))) 

ExportRTqPCR <- a4 %>%
  filter(Symbol %in% c("Olfr358", "Olfr390", "Olfr510", "Olfr596", "Olfr603", "Olfr690", "Olfr1154")) %>%
  dplyr::select(Symbol, log2FoldChange, Corrected_padj) %>%
  unique()
  
write_csv(ExportRTqPCR, file="Output/DESeq10G4_ExportRTqPCR.csv")

#08 - Volcano plots----
#Volcano plot for all Genes filtering out low counts (< 100)

#Volcano plot similar to above but with WT mean counts under 100 removed
a3 %>%
  filter(WTmeanCount >= 100) %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 4, alpha = 0.75, stat = "unique") + 
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Olfr mRNA is significantly reduced") +
  coord_cartesian(xlim = c(-9, 9)) +
  scale_color_manual(breaks = c("DEG", "Olfr DEG", "Taar DEG", "non-DEG"),
                     values=c("red", "blue", "green3", "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_label_repel(data = subset(a3, log2FoldChange > 1 & Genes == "Olfr DEG" & WTmeanCount >= 100),  
                     aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), size = 2.5, stat = "unique", min.segment.length = 0, seed = 1) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.82, 0.82),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(size = 26, hjust = 0.5, face = "bold"), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 18, face = "bold")) + 
  guides(colour = guide_legend(override.aes = list(size=7)))

ggsave("Output/VP_all_MF2A.png", width = 20, height = 20, units = "cm")

#Volcano for Taar/Ms4a/GcD genes only  

# caption = "DEG limited to a change of at least one log2. 
# To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values"

a3 %>%
  filter(GeneType %in% c("Ms4a") | Symbol %in% c("Gucy1b2", "Gucy2d")) %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 7, alpha = 0.85, stat = "unique") + 
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Ms4a and Gucy genes") +
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  scale_color_manual(breaks = c("DEG", "Olfr DEG", "Taar DEG", "non-DEG"),
                     values=c("red", "blue", "green", "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_label_repel(aes(label = Symbol), size = 4, hjust = 0, nudge_y = 0.15, stat = "unique") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.12, 0.82),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(face="bold"),
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold")) + 
  guides(colour = guide_legend(override.aes = list(size=7)))

ggsave("Output/Fig_VP_Ms4aGucySmaller_SF2B.png", width = 15, height = 15, units = "cm")

#09 - DVI-specific plots-----
GOI2 <- c(subset(a3, log2FoldChange > 0 & Genes == "Olfr DEG" & !is.na(DVI) & DVI %notin% c("low expression", "unusual"))$Symbol) %>%
  unique()

#Same as above, but with WTmean count >= 100
a3 %>%
  filter(WTmeanCount >= 100) %>%
  filter(!is.na(DVI)) %>%
  filter(DVI %notin% c("low expression", "unusual")) %>%
  mutate(DVI = as.numeric(as.character(DVI))) %>%
  dplyr::select(GeneType, Symbol, log2FoldChange, Corrected_padj, Genes, DVI) %>%
  unique() %>%
  ggplot(aes(DVI, log2FoldChange)) + 
  geom_jitter(aes(color = Genes), alpha = 0.75) + 
  geom_smooth(se = FALSE, color = "red") +
  scale_color_manual(breaks = c("DEG", "Olfr DEG", "Taar DEG", "non-DEG"),
                     values=c("red", "blue", "green3", "grey")) +
  labs(x="Dorsal-medial/Ventral-lateral Zone Index (Inferred)", 
       y="Log2 Fold Change", 
       title = "Olfr L2FC values trend lower in higher DVIs") +
  scale_y_continuous(breaks=seq(-8,4,2)) +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dashed", col = 'red', size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2) +
  geom_hline(yintercept = -1, linetype = "dashed", col = 'red', size = 1.2) + 
  theme(legend.title = element_blank(),
        legend.position = "none",
        legend.position.inside = c(0.92, 0.82),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"))

ggsave("Output/DVI_L2FC_MF2C.png", width = 20, height = 16, units = "cm")

#10 - Olfr L2FC by Chromosome Plot----
a4 %>%
  filter(WTmeanCount >= 100) %>%
  filter(GeneType %in% c("Olfr", "Taar")) %>% 
  filter(Symbol != "Olfr151") %>%
  ggplot(aes(Chromosome, log2FoldChange, color = Result, shape = Result)) + 
  geom_jitter(size = 3, stat = "unique", alpha = 0.8) +
  scale_color_manual(values=c("skyblue2", "blue1", "green3", "skyblue2", "blue1", "green3")) +
  scale_shape_manual(values=c(16, 16, 16, 1, 1, 1)) +
  labs(title = "Log2 Fold Change of Olfr genes across Chromosomes", 
       y = "Log2 Fold Change") +
  scale_y_continuous(breaks=seq(-8,4,2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 10, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.5, 0.08), 
        legend.box = "horizontal",
        legend.key = element_blank(),
        legend.box.background = element_rect(color="black", linewidth = 0.8, fill = alpha("white", 0.7)),
        legend.background = element_blank()) +
  geom_vline(xintercept = 1.5:18.5, size=0.25) +
  geom_hline(yintercept = 1, linetype = "dashed", col = 'red', size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2) +
  geom_hline(yintercept = -1, linetype = "dashed", col = 'red', size = 1.2) + 
  guides(colour = guide_legend(ncol = 6, override.aes = list(size = 4)))

ggsave("Output/OlfrL2FC_byChrom_MF2B.png", width = 25, height = 15, units = "cm")


#11 - WTcount, Mutantcount, L2FC Plots-------
Top1 <- a4 %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, GeneType, Class) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  arrange(desc(WTmeanCount)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  filter(Symbol %notin% c("OR10G4", "Olfr151")) %>%
  head(100) %>%
  ggplot() + 
  geom_col(aes(Symbol, WTmeanCount), fill = "red", alpha = 0.5) +
  geom_col(aes(Symbol, MutantmeanCount), fill = "blue", alpha = 0.5) +
  annotate("rect", xmin=c(70), xmax=c(72), ymin=c(25000) , ymax=c(30000), alpha=0.5, color="red", fill="red") +
  annotate("rect", xmin=c(70), xmax=c(72), ymin=c(20000) , ymax=c(25000), alpha=0.5, color="purple", fill="purple") +
  annotate("label", x = 72.5, y = 27500, label = "Mean of Raw WT Counts", hjust = 0, color = "red") +
  annotate("label", x = 72.5, y = 22500, label = "Mean of Raw Mutant Counts", hjust = 0, color = "purple") +
  labs(title = "The Top 100 Olfrs by WT mean count", 
       y = "Raw Mean Counts") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45000)) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 7, face = "bold"),
        axis.ticks.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.92, 0.88), 
        legend.title = element_blank())

Bottom1 <- a4 %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, GeneType, Class, Result, Corrected_padj) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  arrange(desc(WTmeanCount)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  filter(Symbol %notin% c("Olfr10G4", "Olfr151")) %>%
  mutate(p.signif = ifelse(Result %in% c("Class II non-DEG", "Class I non-DEG") & Corrected_padj < 0.05, "*", "")) %>%
  head(100) %>%
  ggplot(aes(Symbol, log2FoldChange, fill = Result)) + 
  geom_col() +
  geom_text(aes(Symbol, y = log2FoldChange - 0.45, label = p.signif), size = 5) +
  annotate("label", x = 15, y = -7, label = "Non-DEG Olfr genes with an * have an adusted p value < 0.05.", size = 2.5) +
  labs(y = "Log2 Fold Change") +
  scale_y_continuous(expand = c(0, 0), limits = c(-8.5, 0),breaks=seq(-8,2,2)) +
  scale_fill_manual(values=c("skyblue2", "blue1", "black", "grey")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 7), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        legend.position = c(0.70, 0.18),
        legend.direction = "horizontal", 
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(color="black", linewidth = 0.8, fill = alpha("white", 0.7)),
        legend.background = element_blank())

Top1/Bottom1

ggsave("Output/Top100_Olfr_MF2D.png", width = 35, height = 20, units = "cm")



ChromColors <- c("blue", "red", "green", "purple", "cyan", "orange", "coral", "magenta", "hotpink", "grey50", 
                 "black", "green4", "skyblue", "gold4", "red4", "khaki", "deepskyblue2", "yellow2", "orangered", "darkviolet",
                 "bisque3", "slateblue1", "firebrick1", "mediumseagreen")

#12 - Exploring positive L2FC ----

#Identifying the DE Olfrs with positive L2FC 
a4 %>% 
  filter(Genes == "Olfr DEG" & log2FoldChange >=1) %>% 
  ggplot(aes(log2FoldChange, -log10(Corrected_padj), fill = WTmeanCount_Category)) + 
  geom_point(stat = "unique") +
  geom_label_repel(aes(label = Symbol), size = 4, hjust = 0.5, stat = "unique", min.segment.length = 0, max.overlaps = 15, nudge_y = 2.5) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Differentially Expressed Olfrs with Positive L2FC mostly have relatively low WT counts") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85, 0.20),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold")) + 
  guides(fill = guide_legend(title = "Raw WT mean Counts")) 

ggsave("Output/Supp_Fig2_VP_PositiveOlfrDE_SF2A.png", width = 30, height = 20, units = "cm")

#13 - Graphs with GI Info-----
morecolor <- colorRampPalette(brewer.pal(8, "Dark2"))(24)

Top2 <- a4 %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, GeneType, Corrected_padj, Class, Result, Closest_GI, Closest_GI_Distance_Category) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  filter(Symbol %notin% c("OR10G4", "Olfr151")) %>%
  arrange(desc(WTmeanCount)) %>%
  head(100) %>%
  arrange(Closest_GI, Closest_GI_Distance_Category) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  mutate(p.signif = case_when(Result %in% c("Class II non-DEG", "Class I non-DEG") & Corrected_padj < 0.05 ~ "*",
                              Result %in% c("Class II DEG", "Class I DEG") ~ "",
                              Result %in% c("Class II non-DEG", "Class I non-DEG") ~ "N")) %>%
  ggplot() + 
  geom_col(aes(Symbol, WTmeanCount, fill = Closest_GI_Distance_Category)) +
  geom_text(aes(Symbol, y = WTmeanCount + 1000, label = p.signif), size = 3.5) +
  labs(title = "Analysis of Top 100 Olfrs by WT mean count: Distance to nearest Greek Island", 
       subtitle = "All Olfr Genes are DE except those marked with * (p < 0.05) or N (Not significant)",
       y = "Raw WT Mean Counts") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45000)) +
  scale_fill_discrete(name="Distance from Nearest\nGreek Island") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 7, face = "bold"),
        axis.ticks.x = element_blank(),
        legend.position = c(0.92, 0.68), 
        legend.key = element_blank(),
        legend.box.background = element_rect(color="black", linewidth = 0.8, fill = alpha("white", 0.7)),
        legend.background = element_blank())

Bottom2 <- a4 %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, GeneType, Class, Result, Corrected_padj, Closest_GI_Distance_Category, Closest_GI) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  filter(Symbol %notin% c("OR10G4", "Olfr151")) %>%
  arrange(desc(WTmeanCount)) %>%
  head(100) %>%
  arrange(Closest_GI, Closest_GI_Distance_Category) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  ggplot(aes(Symbol, log2FoldChange, fill = Closest_GI_Distance_Category)) + 
  geom_col() +
  geom_text(aes(label = Closest_GI), vjust = 0.5, hjust = 1.1, angle = 90, size = 3) +
  labs(y = "Log2 Fold Change") +
  scale_y_continuous(expand = c(0, 0), limits = c(-9, 0),breaks=seq(-8,2,2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 7), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        legend.position = "none",
        legend.direction = "horizontal", 
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(color="black", linewidth = 0.8, fill = alpha("white", 0.7)),
        legend.background = element_blank())

Top2 / Bottom2

ggsave("Output/GIDversusL2FC_Top100v2_MF2E.png", width = 50, height = 20, units = "cm")

#14 - Additional Plots----
library(scales)

# Generate 20 colors using a perceptually uniform colormap
colors <- scales::gradient_n_pal(c("red", "yellow", "cyan", "blue"))(seq(0, 1, length.out = 20))

values <- c(
  seq(0, 0.1, length.out = 8),  # More points in 0-10
  seq(0.12, 0.30, length.out = 6), # Moderate spacing in 10-30
  seq(0.35, 1.00, length.out = 6) # Wider spacing for 30-100
)

library(RColorBrewer)

# Define 6 distinct colors
fill_colors <- brewer.pal(6, "Set2") 

a4 %>%
  filter(Symbol != "OR10G4") %>%
  filter(str_detect(Symbol, "Olfr")) %>% 
  select(Symbol, log2FoldChange, WTmeanCount) %>%
  unique() %>%
  arrange(desc(WTmeanCount)) %>%
  mutate(Group = rep(1:ceiling(n()/5), each = 5, length.out = n())) %>%
  group_by(Group) %>%
  summarise(NonNegativeDE = sum(log2FoldChange > -1)/n() * 100,
            AvgWTcount = mean(WTmeanCount),
            AvgL2FC = mean(log2FoldChange)) %>%
  mutate(NonNegativeDEGroup = factor(NonNegativeDE, levels = seq(0, 100, by = 20))) %>%
  ggplot(aes(Group, AvgL2FC, fill = NonNegativeDEGroup)) +
  geom_col() +
  scale_x_continuous(expand = c(0.005, 0)) +
  labs(y = "Average L2FC for each OR group",
       x = "OR group (5 sequential ORs by WT Count arranged from Max to Min)",
       title = "Grouped by WT count, the top few ORs might resist silencing better than most, but this is not a unique property",
       subtitle = "The groups with the lowest WT counts had large variability in resistance",
       fill = "% Non-negative DEG") +
  scale_fill_manual(values = c("black", "red", "yellow", "green", "blue", "purple")) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        plot.subtitle = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.95,0.2),
        legend.key = element_rect(color = "black", linewidth = 1),
        legend.box.background = element_rect(color = "black", linewidth = 1))

ggsave("Output/Grouped_L2FC_SF2F.png", width = 65, height = 20, units = "cm")  

Raw_Count_Totals <- read_excel("Inputs/Raw_Count_Totals.xlsx")

Raw_Count_Totals %>% 
  select(Sample, Genotype, OR_prop, OR10G4_prop) %>% 
  pivot_longer(cols = c("OR_prop", "OR10G4_prop"), names_to = "Gene", values_to = "Prop") %>% 
  mutate(Gene = str_remove(Gene, "_prop")) %>% 
  mutate(Genotype = ifelse(Genotype == "WT", "WT", "OR10G4"), 
         Gene = ifelse(Gene == "OR10G4", "10G4", "OR")) %>% 
  mutate(Genotype_Gene = paste0(Genotype, "_", Gene)) %>% 
  filter(Genotype_Gene != "WT_10G4") %>% 
  mutate(SampleV2 = paste0(Sample, " (", Genotype, ")")) %>% 
  ggplot(aes(SampleV2, Prop, fill = Gene)) + 
  geom_col() +
  labs(x = "Sample (Genotype)", y = "Proportion of Total Counts",
       title = "10G4 mRNA does not compensate for Endogenous OR mRNA loss") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.position = "top",
        legend.key = element_rect(color = "black", linewidth = 1),
        legend.box.background = element_rect(color = "black", linewidth = 1))

ggsave("Output/OR_Totals_SF2D.png", width = 30, height = 15, units = "cm")  


HKgenes <- a4 %>% 
  filter(Symbol %in% c("Actb", "Gapdh", "Ywhaz", "Ubc", "Ppia")) %>% 
  select(GeneID, Symbol, log2FoldChange, WTmeanCount, MutantmeanCount) %>% 
  unique() %>% 
  left_join(count_table, by = c("GeneID" = "gene_id")) %>%
  select(-c(GeneID, log2FoldChange, WTmeanCount, MutantmeanCount))

Raw_Count_Totals %>%
  select(Sample, All) %>%
  pivot_wider(names_from = Sample, values_from = All) %>%
  bind_cols(tibble(Symbol = "Total")) %>%
  bind_rows(HKgenes) %>%
  pivot_longer(cols = -Symbol, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = Symbol, values_from = "Count") %>%
  pivot_longer(cols = 3:7, names_to = "Symbol", values_to = "Count") %>%
  mutate(Genotype = ifelse(Sample %in% c("4164-M", "4165-M", "4170-M"), "OR10G4", "WT")) %>%
  mutate(Prop = Count/Total) %>%
  ggplot(aes(Genotype, Prop, fill = Sample)) +
  geom_col(position = "dodge") +
  facet_wrap(~Symbol, scales = "free_y", nrow = 1) +
  labs(y = "Proportion of Total Counts",
       title = "No notable difference in HouseKeeping Genes") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 10),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.key = element_rect(color = "black", linewidth = 1),
        legend.box.background = element_rect(color = "black", linewidth = 1))  +
  guides(fill = guide_legend(nrow = 1))
  
ggsave("Output/HouseKeeping_Genotype_SF2E.png", width = 20, height = 10, units = "cm")  

DVI_forCompare <- a4 %>%
  filter(!is.na(DVI)) %>%
  select(Symbol, log2FoldChange, DVI) %>%
  unique()

df <- ggpubr::compare_means(log2FoldChange ~ DVI, DVI_forCompare, method = "t.test") %>% 
  select(group1, group2, p.adj)

# Make symmetric matrix by adding reversed comparisons
df_symmetric <- df %>%
  bind_rows(df %>% rename(group1 = group2, group2 = group1))

# Add diagonal (self-comparisons)
groups <- union(df$group1, df$group2)
df_diag <- tibble(group1 = groups, group2 = groups, p.adj = NA_real_)

# Full dataset
df_full <- bind_rows(df_symmetric, df_diag)

DVI_sums <- DVI_forCompare %>%
  count(DVI) %>%
  mutate(n = as.character(n))

# Add a column to indicate significance
df_full2 <- df_full %>%
  left_join(DVI_sums, by = c("group1" = "DVI")) %>%
  filter(group1 <= group2) %>%
  mutate(significant = case_when(
    is.na(p.adj) ~ "Group Size", 
    p.adj <= 0.05 ~ "Significant",
    TRUE ~ "Not significant"),
    group1 = factor(group1),
    group2 = factor(group2),
    entry = ifelse(is.na(p.adj), n, as.character(sprintf("%.2f", p.adj))))

# Plot binary heatmap
ggplot(df_full2, aes(x = group1, y = group2, fill = significant)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("Significant" = "red", "Not significant" = "grey90", "Group Size" = "skyblue"),
    drop = FALSE) +
  geom_text(aes(label = entry), size = 3) +
  scale_y_discrete(limits=rev) +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = "Pairwise Significance Heatmap (padj ≤ 0.05)", fill = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.75, 0.75))

ggsave("Output/DVI_Comparison_stats_SF2C.png", width = 30, height = 20, units = "cm")  


Enhancer14 <- a4 %>%
  filter(Chromosome == "14") %>%
  select(Greek_Island, Start_GI) %>%
  unique()

EnhancersTesting <- a4 %>% 
  select(Chromosome, Greek_Island, Start_GI, End_GI) %>% 
  unique()

TripleInput <- a4 %>% 
  filter(str_detect(Symbol, "Olfr")) %>% 
  select(Symbol, WTmeanCount, log2FoldChange, Closest_GI, SmallestGIDistance) %>% 
  unique() %>% 
  arrange(desc(SmallestGIDistance)) %>% 
  mutate(Symbol = factor(Symbol, levels = Symbol)) 

LeftPlot <- TripleInput %>% 
  ggplot() + 
  geom_point(aes(log10(SmallestGIDistance), Symbol), color = "green3") +
  labs(x = "Nearest OR Enhancer in bp (Log10 Scale)", y = "OR Genes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 14))

MiddlePlot <- TripleInput %>% 
  ggplot() + 
  geom_point(aes(log2FoldChange, Symbol), color = "blue") +
  labs(x = "Log2 Fold Change") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 14))


RightPlot <- TripleInput %>% 
  ggplot() + 
  geom_point(aes(log10(WTmeanCount), Symbol), color = "red2") +
  labs(x = "Mean Raw WT Counts (Log10 Scale)", y = "OR Genes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 14))

LeftPlot + MiddlePlot + RightPlot + 
  plot_annotation(title = "No Global relationship between OR Enhancer Proximity and L2FC or Counts", 
                  theme =  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

ggsave("Output/GIDversusL2FCCount_AllOR_SF2H.png", width = 30, height = 15, units = "cm")

GI_GI_forplot <- Greek_Islands2 %>% 
  left_join(Greek_Islands2, by = "Chromosome") %>% 
  filter(Greek_Island.x != Greek_Island.y) %>%
  filter(abs(Start_GI.x - End_GI.y) <= 500000 | abs(End_GI.x - Start_GI.y) <= 500000) %>%
  mutate(GI_GI = ifelse(Start_GI.x < Start_GI.y, Start_GI.y - End_GI.x, End_GI.y - Start_GI.x)) %>%
  rename(Greek_Island = Greek_Island.x) %>%
  arrange(Chromosome, Start_GI.x) %>%
  mutate(Greek_Island = factor(Greek_Island, levels = unique(Greek_Island)))

chroms <- c(as.character(c(1:4,6:11, 13:17, 19)), "X")
chrom_colors <- setNames(
  c(
    RColorBrewer::brewer.pal(7, "Set1"),
    RColorBrewer::brewer.pal(7, "Dark2"),
    RColorBrewer::brewer.pal(3, "Set2")
  )[1:length(chroms)],chroms)

# One color per Enhancer strip
Greek_Islands2 <- Greek_Islands2 %>%
  mutate(Chromosome = factor(Chromosome, levels = chroms)) %>%
  arrange(Chromosome)
strip_info <- unique(Greek_Islands2[c("Greek_Island", "Chromosome")])
strip_info <- strip_info[match(unique(Greek_Islands2$Greek_Island), strip_info$Greek_Island), ]
strip_fills <- chrom_colors[strip_info$Chromosome]

# Build ggh4x strip fill and text theme
strip_fill_elements <- lapply(strip_fills, function(col) element_rect(fill = col))
strip_text_elements <- lapply(seq_along(strip_fills), function(i) element_text(color = "white"))

chrom_legend_df <- unique(Greek_Islands2[c("Chromosome")])
chrom_legend_df$x <- 0
chrom_legend_df$y <- 0

GI_GI_forplotHfocus <- Greek_Islands2 %>%
  filter(Greek_Island == "H") %>%
  left_join(Greek_Islands2, by = "Chromosome") %>% 
  filter(Greek_Island.x != Greek_Island.y) %>%
  filter(abs(Start_GI.x - End_GI.y) <= 2000000 | abs(End_GI.x - Start_GI.y) <= 2000000) %>%
  mutate(GI_GI = ifelse(Start_GI.x < Start_GI.y, Start_GI.y - End_GI.x, End_GI.y - Start_GI.x)) %>%
  rename(Greek_Island = Greek_Island.x) %>%
  arrange(Chromosome, Start_GI.x)

NonOlfr_H_Neighbors <- a4 %>%
  filter(!str_detect(Symbol, "Olfr")) %>%
  filter(WTmeanCount >= 1) %>%
  filter(Greek_Island == "H") %>%
  filter(GI_Distance <= 2000000) %>%
  mutate(GeneGroup = ifelse(str_detect(Symbol, "Olfr"), "OR", "NonOR")) %>%
  mutate(R_GI_Distance = ifelse(Start < Start_GI, -GI_Distance, GI_Distance))
    
a4 %>%
  filter(str_detect(Symbol, "Olfr")) %>%
  filter(WTmeanCount >= 72) %>%
  filter(Greek_Island == "H") %>%
  filter(GI_Distance <= 2000000) %>%
  mutate(R_GI_Distance = ifelse(Start < Start_GI, -GI_Distance, GI_Distance)) %>%
  arrange(Chromosome, Start_GI) %>%
  mutate(Greek_Island = factor(Greek_Island, levels = unique(Greek_Island))) %>%
  ggplot(aes(R_GI_Distance, log2FoldChange)) +
  geom_vline(data = GI_GI_forplotHfocus, aes(xintercept = GI_GI), color = "green3") +
  geom_text_repel(data = GI_GI_forplotHfocus, aes(x = GI_GI, y = -5, label = Greek_Island.y), seed = 43, angle = 90) +
  geom_vline(xintercept = 0, color = "red") +
  geom_hline(yintercept = 0, color = "cyan3", linetype = "dashed") +
  geom_hline(yintercept = -2, color = "violet", linetype = "dashed") +
  geom_line() +
  geom_point(aes(color = WTmeanCount), size = 5) +
  geom_point(data = NonOlfr_H_Neighbors, aes(color = WTmeanCount), shape = 17, size = 3) +
  geom_smooth(data = NonOlfr_H_Neighbors, se = FALSE, linewidth = 1) +
  geom_text_repel(aes(label = Symbol),min.segment.length = 0, seed = 42, box.padding = 0.5) +
  scale_color_gradientn(
    colours = c("blue3", "cyan4", "green3", "yellow3", "orange2", "red3"),
    trans = "log",
    breaks = c(1, 100, 1000, 3000, 10000, 40000),
    labels = scales::comma_format()) +
  labs(title = "The H element and surrounding GI Enhancers impact on OR Counts and L2FC",
       x = "Relative position to H element in bp", y = "Log2 Fold Change") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2))

ggsave("Output/GIDversusL2FCCount_Hfocusv2_SF2G.png", width = 40, height = 25, units = "cm")
