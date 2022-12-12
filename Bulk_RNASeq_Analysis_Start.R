#This was started in a clean-ish session. I basically had these packages installed, but not loaded. 
#Also, the location of the count matrix is specific to my machine.
#Heavily based on https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#DESeq2 and basic checks based on above
#Get all the libraries active in your session before loading the data and subsequent steps. Use install and library functions one at a time.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

BiocManager::install("vsn")

library(vsn)

#Preemptive side step in prep for later visualization: Log Fold Change Shrinkage with the apeglm method
BiocManager::install("apeglm")

library(apeglm)


if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

if (!require(pheatmap)) install.packages('pheatmap')
library(pheatmap)

if (!require(biomaRt)) install.packages('biomaRt')
library(biomaRt)

if (!require(ggrepel)) install.packages('ggrepel')
library(ggrepel)

if (!require(cowplot)) install.packages('cowplot')
library(cowplot)

if (!require(ggh4x)) install.packages('ggh4x')
library(ggh4x)

#Import and prepare the data to make a cts matrix and a coldata table. The file directory will be different. Also, import done via "Import Dataset" and "tab-delimited"
count_table <- read_delim("raw_counts_genes.clean.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
View(count_table)
#Extract the gene ids to use as rownames instead
genes <- count_table$gene_id

#Select specific sample columns and remove gene_id column. I'm ignoring the OB samples. Make the counts a matrix and add rownames
cts <- count_table[, c(2:8)]
cts <- as.matrix(cts)
rownames(cts) <- genes

#Making the coldata object with the condition of the samples, samples as rownames, and condition as a factor
coldata <- data.frame(condition = c("untreated", "treated", "treated", "untreated", "untreated", "treated", "untreated"))
rownames(coldata) <- colnames(cts)
coldata$condition <- as.factor(coldata$condition)

#A check on matching the rownames and colnames, aka the samples. Should be true.
all(rownames(coldata) == colnames(cts))

#The rest of the process is kinda simple enough. The (experimental) design here is simple with only one feature to select from the coldata table: condition
dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = coldata, 
                              design = ~ condition)

#We will apply a small pre-filtering step to the data, which in this case will take dds from 55k to about 30k genes.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Here the untreated condition is established as the reference
dds$condition <- relevel(dds$condition, ref = "untreated")

#run the function and see the results.
dds <- DESeq(dds)
res <- results(dds)
res

resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resLFC

#A few quick modified looks at the results data. Note, the results were generated with a default alpha of 0.1.
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
sum(res$padj < 0.05, na.rm = TRUE)
sum(res$padj < 0.01, na.rm = TRUE)
resOrdered <- res[order(res$pvalue),]
resOrdered

#repeating the results with a modified alpha value of 0.05 and 0.01
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
res01 <- results(dds, alpha=0.01)
summary(res01)
sum(res01$padj < 0.01, na.rm=TRUE)

#PLOTS! This will use the res05 data when possible
plotMA(res05, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))
plotCounts(dds, gene=which.min(res05$padj), intgroup="condition") #Looks at the norm counts of the gene with the smallest padj

#plotCounts using ggplot
d <- plotCounts(dds, gene=which.min(res05$padj), intgroup="condition", 
                returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#Data transformations and visualizations: Note- I do not know the value of these transformations. You may be asked to install the hexbin package when you run meanSdPlot. Say yes.
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

ntd <- normTransform(dds)

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#Data QA via sample clustering and visualization. I couldn't figure out how to pheatmap without sizeFactor, but it seems to be a small issue.
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","sizeFactor")])

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)

plotPCA(vsd, intgroup=c("condition")) #This is the simple version of the PCA, below is the ggplot version with more control.
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#Some more QAQC plots
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
plotDispEsts(dds)


##Annotation with biomaRt
listMarts()
## set up connection to ensembl database
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
# list the available datasets (species)
listDatasets(ensembl) %>% 
  filter(str_detect(description, "Mouse"))
# specify a data set to use
ensembl = useDataset("mmusculus_gene_ensembl", mart=ensembl)
listFilters(ensembl) %>% 
  filter(str_detect(name, "ensembl"))
# Set the filter type and values
ourFilterType <- "ensembl_gene_id"
filterValues <- rownames(res05)
# check the available "attributes" - things you can retrieve
listAttributes(ensembl) %>% 
  head(50)
# Set the list of attributes
attributeNames <- c('ensembl_gene_id', 'external_gene_name', 'description', 'chromosome_name', 'start_position', 'end_position')
# run the query
annot2 <- getBM(attributes=attributeNames, 
                filters = ourFilterType, 
                values = filterValues, 
                mart = ensembl)

#With Chromosome Zero, I need to add it to the annot2 table manually.
First <- data.frame('ensembl_gene_id' = "PFMUSG0001", 
                    'external_gene_name' = "Olfr10G4", 
                    'description' = "10G4 Transgene", 
                    'chromosome_name' = "3", 
                    'start_position' = 60436556, 
                    'end_position' = 60487531)


#Weird Genes that don't have a home...
OddOneOut1 <- data.frame('ensembl_gene_id' = "ENSMUSG00000106198", 
                     'external_gene_name' = "Gm42504", 
                     'description' = "predicted gene 42504", 
                     'chromosome_name' = "5", 
                     'start_position' = 143379444,
                     'end_position' = 143383212)

OddOneOut2 <- data.frame('ensembl_gene_id' = "ENSMUSG00000079170", 
                         'external_gene_name' = "Gm13941", 
                         'description' = "predicted gene 13941", 
                         'chromosome_name' = "2", 
                         'start_position' = 108214986,
                         'end_position' = 108234110)

OddOneOut3 <- data.frame('ensembl_gene_id' = "ENSMUSG00000079171", 
                         'external_gene_name' = "Gm13942", 
                         'description' = "predicted gene 13942", 
                         'chromosome_name' = "2", 
                         'start_position' = 110881665,
                         'end_position' = 110887919)

OddOneOut4 <- data.frame('ensembl_gene_id' = "ENSMUSG00000085431", 
                         'external_gene_name' = "4930440I19Rik", 
                         'description' = "RIKEN cDNA 4930440I19 gene", 
                         'chromosome_name' = "2", 
                         'start_position' = 77881521,
                         'end_position' = 78024774)

annot3 <- bind_rows(annot2, OddOneOut1) %>%
  bind_rows(OddOneOut2) %>%
  bind_rows(OddOneOut3) %>%
  bind_rows(OddOneOut4)
annot3 <- annot3 %>% arrange(ensembl_gene_id)
annot3 <- bind_rows(First, annot3)

# %notin% : Convenience function that provides a vectorized != logical for filter(); found online, unknown source ----
`%notin%` <- Negate(`%in%`)
rownames(dds)[which(rownames(dds) %notin% annot3$ensembl_gene_id)]

#Prepare the annotation
all(rownames(dds) == annot3$ensembl_gene_id)
colnames(annot3) <- c("GeneID", "Symbol", "Description", "Chromosome", "Start", "End")

#Combine it with the res05 object
annot05 <- as.data.frame(res05) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(annot3, "GeneID") 

#reordering the table a few ways
annot05Ordered <- annot05[order(annot05$pvalue),]
annot05Ordered
write_csv(annot05Ordered, 
          file="DESeq10G4_30kgenesv2.csv")

OlfrOnly <- annot05Ordered[str_detect(annot05Ordered$Symbol, "Olfr"), ]
OlfrOnly
write_csv(OlfrOnly, 
          file="DESeq10G4_Olfrgenesv2.csv")

NotOlfr <- annot05Ordered[str_detect(annot05Ordered$Symbol, "Olfr", negate = TRUE), ]
NotOlfr
write_csv(NotOlfr, 
          file="DESeq10G4_NOTolfrgenesv2.csv")

TaarOnly <- annot05Ordered[str_detect(annot05Ordered$Symbol, "Taar"), ]
TaarOnly
write_csv(TaarOnly, 
          file="DESeq10G4_Taargenesv2.csv")