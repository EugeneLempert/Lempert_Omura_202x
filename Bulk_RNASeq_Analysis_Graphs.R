#------Plots V1. source of a1 set------
#Plots Final?
cutoff <- annot05Ordered$pvalue[10]
a3 <- annot05Ordered %>% 
  mutate(TopGeneLabel=ifelse(pvalue<=cutoff, Symbol, ""))

ggplot(a3, aes(x = log2(baseMean), y=log2FoldChange)) + 
  geom_point(aes(colour=padj < 0.05), shape=20, size=0.5) +
  geom_text(aes(label=TopGeneLabel), size = 2.5) +
  labs(x="mean of normalised counts", y="log fold change")

a1 <- annot05Ordered %>% mutate(GeneType = case_when(
  str_detect(Symbol, "Olfr") ~ "Olfr", 
  str_detect(Symbol, "Taar") ~ "Taar",
  TRUE ~ "Other")
)

ggplot(a1, aes(x= log2(baseMean), y = log2FoldChange)) + 
  geom_point(aes(color = GeneType), shape = 20, size = 0.5) + 
  geom_text(data = subset(a1, Symbol == "Olfr151"), aes(log2(baseMean), log2FoldChange, label = Symbol)) +
  labs(x="log2(mean of normalized counts)", y="log 2 fold change", title = "10G4 RNA-Seq Analysis") +
  theme_bw()

a1 <- a1 %>% mutate(Genes = case_when(
  GeneType == "Olfr" & padj < 0.05 & log2FoldChange < -1 ~ "Olfr DEG",
  GeneType == "Olfr" & padj < 0.05 & log2FoldChange > 1 ~ "Olfr DEG",
  GeneType == "Taar" & padj < 0.05 & log2FoldChange < -1 ~ "Taar DEG",
  GeneType == "Taar" & padj < 0.05 & log2FoldChange > 1 ~ "Taar DEG",
  GeneType == "Other" & padj < 0.05 & log2FoldChange < -1~ "DEG",
  GeneType == "Other" & padj < 0.05 & log2FoldChange > 1~ "DEG",
  TRUE ~ "non-DEG"
))

a1 <- a1 %>% mutate(Corrected_padj = padj + subset(a1, Symbol == "Olfr711")$padj)

Genelevels <- c("DEG", "Olfr DEG", "Taar DEG", "non-DEG")
a1$Genes <- factor(a1$Genes, levels = Genelevels)

GOI <- c(subset(a1, log2FoldChange > 7 & -log10(Corrected_padj) < 50)$Symbol, 
         subset(a1, log2FoldChange > 2.5 & -log10(Corrected_padj) > 40)$Symbol,
         subset(a1, log2FoldChange > 0 & -log10(Corrected_padj) > 70)$Symbol,
         subset(a1, log2FoldChange > -1 & -log10(Corrected_padj) > 200)$Symbol,
         subset(a1, log2FoldChange < -7.5 & -log10(Corrected_padj) > 0)$Symbol)
GOI <- unique(GOI)

ggplot(a1, aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 0.5) + 
  labs(x="Log2 Fold Change", y="-log10(Adjusted P-Value)", title = "Most Olfr and Taar mRNAs are significantly reduced in 10G4 samples", caption = "To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-12, 12)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_text(data = subset(a1, Symbol %in% GOI), aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), size = 4) +
  theme_bw() +
  theme(legend.title = element_blank()) + 
  guides(colour = guide_legend(override.aes = list(size=7)))

ggplot(a1, aes(x= log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 0.5) + 
  labs(x="Log2 Fold Change", y="-log10(Adjusted P-Value)") +
  coord_cartesian(xlim = c(-10, 10)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash") +
  geom_text(data = subset(a1, Symbol == "Olfr151"), aes(log2FoldChange, -log10(padj), label = Symbol)) +
  theme_bw() 

OlfrTaar <- subset(a1, GeneType %in% c("Olfr", "Taar"))
NonOlfrTaar <- subset(a1, GeneType == "Other")
GOI2 <- GOI[str_detect(GOI, "Olfr")]

ggplot(OlfrTaar, aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 0.5) + 
  geom_text(data = subset(a1, Symbol %in% GOI2), aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), size = 4) +
  labs(x="Log2 Fold Change", y="-log10(Adjusted P-Value)", title = "Most Olfr and Taar mRNAs are significantly reduced in 10G4 samples", caption = "To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-12, 12)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.3) +
  theme_bw() +
  theme(legend.title = element_blank()) + 
  guides(colour = guide_legend(override.aes = list(size=7)))

ggplot(NonOlfrTaar, aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 0.5) + 
  labs(x="Log2 Fold Change", y="-log10(Adjusted P-Value)", title = "10G4 impact on Non-Olfr & Non-Taar mRNAs", caption = "To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-10, 10)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.3) +
  theme_bw() +
  theme(legend.title = element_blank()) + 
  guides(colour = guide_legend(override.aes = list(size=7)))

ggplot(a1, aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 0.5) + 
  labs(x="Log2 Fold Change", y="-log10(Adjusted P-Value)", title = "Most Olfr and Taar mRNAs are significantly reduced in 10G4 samples", caption = "To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-18, 18), ylim = c(0,300)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash") +
  geom_text(data = subset(a1, Symbol %in% GOI2), aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), size = 4) +
  theme_bw() 

#Step to generate some intermediate datasets used in the following plots
Olfrx <- OlfrOnly %>%
  dplyr::select(GeneID, Symbol)

Countx <- count_table %>%
  rename(GeneID = gene_id)

Countx1 <- Countx %>%
  inner_join(Olfrx)

Countx2 <- Countx1 %>%
  dplyr::select(-GeneID, -Symbol)

CountxG <- Countx1 %>%
  filter(Symbol == "Olfr10G4") %>%
  dplyr::select(-GeneID, -Symbol)

Countx3 <- Countx2 %>%
  colSums(na.rm = TRUE)

Countx4 <- CountxG/Countx3

Countx5 <- bind_rows(CountxG, Countx3, Countx4)

Countxcol <- data.frame(Value = c("OR10G4", "All_ORs", "Ratio"))
Countx6 <- bind_cols(Countxcol, Countx5)

write_csv(Countx6, 
          file="DESeq10G4_rawratio.csv")

Countx %>%
  dplyr::select(-GeneID) %>%
  colSums()

annotx1 <- annot05Ordered  %>%
  dplyr::select(GeneID, Symbol)

Countx <- count_table %>%
  rename(GeneID = gene_id)

annotx2 <- Countx %>%
  inner_join(annotx1)

annotx3 <- annotx2 %>%
  dplyr::select(Symbol, GeneID, `4163-M`, `4166-M`, `4169-M`, `4171-M`)

annotx3$Mean <- rowMeans(annotx3[,3:6])

annotx4 <- annotx3 %>%
  arrange(desc(Mean))

annotx4 <- annotx4[1:50,1]
annotx4 <- unlist(annotx4)


annotx5 <- annotx3[str_detect(annotx3$Symbol, "Olfr"),]
annotx5 <- annotx5 %>%
  arrange(desc(Mean))

annotx6 <- annotx5[1:50,1]
annotx6 <- unlist(annotx6)

annotx7 <- unlist(annotx5$Symbol)

annotx3 <- annotx2 %>%
  dplyr::select(Symbol, GeneID, `4164-M`, `4165-M`,`4170-M`)

annotx3$Mean <- rowMeans(annotx3[,3:5])

annotx8 <- annotx3[str_detect(annotx3$Symbol, "Olfr"),]
annotx8 <- annotx8 %>%
  arrange(desc(Mean))

#Final steps will include saving the data and maybe some kind of rich report package work.

ggplot(a1, aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 0.5) + 
  labs(x="Log2 Fold Change", y="-log10(Adjusted P-Value)", title = "Most Olfr and Taar mRNAs are significantly reduced in 10G4 samples", caption = "To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-18, 18), ylim = c(0,300)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash") +
  geom_text(data = subset(a1, Symbol %in% annotx4), aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), size = 4) +
  theme_bw() 

ggplot(a1, aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 0.5) + 
  labs(x="Log2 Fold Change", y="-log10(Adjusted P-Value)", title = "Most Olfr and Taar mRNAs are significantly reduced in 10G4 samples", caption = "To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-18, 18), ylim = c(0,300)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash") +
  geom_text(data = subset(a1, Symbol %in% annotx6), aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), size = 4) +
  theme_bw() 


#------Source of Big Data set-----

DVI2 <- read_csv("../../Supportive Files/DVI2.csv")

All_info <- a1 %>% 
  left_join(annotx3, by = c("Symbol", "GeneID")) %>%
  left_join(DVI2, by = "Symbol") %>%
  mutate(DVI = factor(DVI, levels = c("1.05", "1.25", "1.4", "1.5", "1.6", "1.7", "1.8", "1.9", 
                                      "2", "2.1", "2.2", "2.3", "2.4", "2.5", "2.6", "2.7", "2.8", "2.9", 
                                      "3.05", "3.2", "3.3", "3.4", "3.5", "3.6", "3.7", "3.8", "3.9", "4",
                                      "4.1", "4.2", "4.3", "4.4", "4.5", "4.6", "4.7", "4.8", "4.9", "5",
                                      "low expression", "unusual")))

WTOlfrRanking <- All_info %>% arrange(desc(WTmeanCount)) %>% filter(GeneType == "Olfr") %>% head(50)

annotx1 <- a1  %>%
  dplyr::select(GeneID, Symbol)
Countx <- count_table %>%
  dplyr::rename(GeneID = gene_id)
annotx2 <- Countx %>%
  inner_join(annotx1)
annotx3 <- annotx2 %>%
  dplyr::select(Symbol, GeneID)
annotx3$WTmeanCount <- rowMeans(annotx2[,c(2,5,6,8)])  
annotx3$MutantmeanCount <- rowMeans(annotx2[,c(3,4,7)])

All_info <- All_info %>%
  mutate(Top50WT = case_when(GeneType == "Olfr" & Symbol %in% WTOlfrRanking$Symbol ~ TRUE, 
                             GeneType == "Olfr" & Symbol %notin% WTOlfrRanking$Symbol ~ FALSE))

ClassI <- All_info %>% filter(GeneType =="Olfr") %>% filter(Chromosome == 7, Start > 102000000, End < 106000000)

All_info <- All_info %>%
  mutate(Class = case_when(GeneType == "Olfr" & Symbol %in% ClassI$Symbol ~ "Class I", 
                           GeneType == "Olfr" & Symbol %notin% ClassI$Symbol ~ "Class II",
                           GeneType == "Taar" ~ "Class TAAR"))

#This imports a Greek Island (enhancer) list with positions that come from the same mouse genome: mm39
GKIm39 <- read_delim("../gxAeke3sdzxUnjBB.bed", 
              delim = "\t", escape_double = FALSE, 
              col_names = FALSE, trim_ws = TRUE)

names(GKIm39) <- c("Chromosome", "GI_Start", "GI_End", "Greek_Island")  

All_info <- All_info %>%
  left_join(GKIm39, by = "Chromosome")

Big_Data <- All_info %>% 
  rowwise() %>% 
  mutate(GI_Distance = min(c(abs(End - GI_Start), abs(Start - GI_End)))) %>%
  ungroup()

Big_Data <- Big_Data %>% 
  arrange(Greek_Island, GI_Distance) %>% 
  group_by(Greek_Island) %>% 
  mutate(GI_Prox_Rank = rank(GI_Distance, ties.method = "first")) %>%
  ungroup()

Big_Olfr <- Big_Data %>% 
  filter(GeneType == "Olfr") %>% 
  arrange(Greek_Island, GI_Distance) %>% 
  group_by(Greek_Island) %>%  
  mutate(GI_Prox_RankvsOlfrs = rank(GI_Distance, ties.method = "first")) %>%
  ungroup() %>%
  dplyr::select(Symbol, Greek_Island, GI_Prox_RankvsOlfrs)

Big_Data <- Big_Data %>%
  left_join(Big_Olfr, by = c("Symbol", "Greek_Island"))

Big_Data <- Big_Data %>% 
  mutate(GeneType = case_when(
    str_detect(Symbol, "Olfr") ~ "Olfr", 
    str_detect(Symbol, "Taar") ~ "Taar",
    str_detect(Symbol, "Ms4a") ~ "Ms4a",
    TRUE ~ "Other"))

#------Olfr subset-----
GOI <- c(subset(Big_Data, log2FoldChange > 2.5)$Symbol, 
         subset(Big_Data, log2FoldChange < -7.5)$Symbol,
         subset(Big_Data, -log10(Corrected_padj) > 250)$Symbol)

GOI2 <- unique(GOI)

Olfr_subset <- Big_Data %>%
  filter(GeneType == "Olfr") %>%
  filter(Symbol != "Olfr10G4") %>%
  dplyr::select(GeneType, Symbol, log2FoldChange, Corrected_padj, Genes) %>%
  unique()

Olfr_subset %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 1, alpha = 0.75) + 
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Most Olfr mRNAs are significantly reduced in 10G4 samples", 
       caption = "DEG limited to a change of at least one log2. To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-8, 5)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_text(data = subset(Olfr_subset, Symbol %in% GOI2[str_detect(GOI2, "Olfr")]), 
            aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), 
            size = 4, position = position_jitter(height = 3.7, width = 0)) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.92, 0.13),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))

#------ClassI General----
GOI3 <- c(subset(Big_Data, log2FoldChange > 2)$Symbol, 
          subset(Big_Data, log2FoldChange < -6)$Symbol,
          subset(Big_Data, -log10(Corrected_padj) > 150)$Symbol)

GOI3 <- unique(GOI3)

ClassI_subset <- Big_Data %>%
  filter(Class == "Class I") %>%
  dplyr::select(GeneType, Symbol, log2FoldChange, Corrected_padj, Genes) %>%
  unique()

ClassI_subset %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 1, alpha = 0.75) + 
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Most Class I Olfr mRNAs are significantly reduced in 10G4 samples", 
       caption = "DEG limited to a change of at least one log2. To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-8, 5)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_text(data = subset(ClassI_subset, Symbol %in% GOI3[str_detect(GOI3, "Olfr")]), 
            aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), 
            size = 4, position = position_jitter(height = 3.7, width = 0)) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.92, 0.13),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))

#------ClassII general-----
GOI4 <- c(subset(Big_Data, log2FoldChange > 2.8)$Symbol, 
          subset(Big_Data, log2FoldChange < -6.9)$Symbol,
          subset(Big_Data, -log10(Corrected_padj) > 245)$Symbol)

GOI4 <- unique(GOI4)

ClassII_subset <- Big_Data %>%
  filter(Class == "Class II") %>%
  filter(Symbol != "Olfr10G4") %>%
  dplyr::select(GeneType, Symbol, log2FoldChange, Corrected_padj, Genes) %>%
  unique()

ClassII_subset %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 1, alpha = 0.75) + 
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Most Class II Olfr mRNAs are significantly reduced in 10G4 samples", 
       caption = "DEG limited to a change of at least one log2. To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-8, 5)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_text(data = subset(ClassII_subset, Symbol %in% GOI4[str_detect(GOI4, "Olfr")]), 
            aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), 
            size = 4, position = position_jitter(height = 3.7, width = 0)) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.92, 0.13),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))

#------Taar Subset------
ClassTAAR_subset <- Big_Data %>%
  filter(Class == "Class TAAR") %>%
  filter(Symbol != "Olfr10G4") %>%
  dplyr::select(GeneType, Symbol, log2FoldChange, Corrected_padj, Genes) %>%
  unique()

ClassTAAR_subset %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 2, alpha = 0.75) + 
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Almost all Taar mRNAs are significantly reduced in 10G4 samples", 
       caption = "DEG limited to a change of at least one log2. To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-7, 1)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_text(aes(label = Symbol), 
            size = 4, position = position_jitter(height = 1, width = 0)) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.92, 0.93),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))

#------MS4a subset-------
TypeMs4a_subset <- Big_Data %>%
  filter(GeneType == "Ms4a") %>%
  filter(Symbol != "Olfr10G4") %>%
  dplyr::select(GeneType, Symbol, log2FoldChange, Corrected_padj, Genes) %>%
  unique()

TypeMs4a_subset %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 2, alpha = 0.75) + 
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Some Ms4a-associated mRNAs are slightly, but significantly increased in 10G4 samples", 
       caption = "DEG limited to a change of at least one log2. To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-2, 2)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_text(aes(label = Symbol), 
            size = 4) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.92, 0.13),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))

#------Taar-MS4a-Gcd----
TypeMs4aTaarGcd_subset <- Big_Data %>%
  filter(GeneType %in% c("Ms4a", "Taar") | Symbol %in% c("Gucy1b2", "Gucy2d")) %>%
  filter(Symbol != "Olfr10G4") %>%
  dplyr::select(GeneType, Symbol, log2FoldChange, Corrected_padj, Genes) %>%
  unique()

GOI12 <- c(subset(Big_Data, GeneType == "Taar")$Symbol, 
           subset(Big_Data, Symbol %in% c("Gucy1b2", "Gucy2d"))$Symbol) %>%
  unique()

TypeMs4aTaarGcd_subset %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 5, alpha = 0.75) + 
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Almost all Taar mRNAs are significantly reduced in 10G4 samples",
       subtitle = "Some Ms4a-associated mRNAs, unlabeled below, are slightly, but significantly increased in 10G4 samples",
       caption = "DEG limited to a change of at least one log2. To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-9, 9)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_text(data = subset(TypeMs4aTaarGcd_subset, Symbol %in% GOI12), 
            aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), 
            size = 3) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.92, 0.13),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))

#------ClassI DVI1.05-------
GOI5 <- c(subset(Big_Data, log2FoldChange > 1.25)$Symbol, 
          subset(Big_Data, log2FoldChange < -6.2)$Symbol,
          subset(Big_Data, -log10(Corrected_padj) > 125)$Symbol)

GOI5 <- unique(GOI5)

ClassIDVI1_subset <- Big_Data %>%
  filter(Class == "Class I") %>%
  filter(DVI == "1.05") %>%
  dplyr::select(GeneType, Symbol, log2FoldChange, Corrected_padj, Genes, DVI) %>%
  unique()

ClassIDVI1_subset %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 1, alpha = 0.75) + 
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Most Class I DVI 1.05 Olfr mRNAs are significantly reduced in 10G4 samples", 
       caption = "DEG limited to a change of at least one log2. To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-8, 5)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_text(data = subset(ClassIDVI1_subset, Symbol %in% GOI5[str_detect(GOI5, "Olfr")]), 
            aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), 
            size = 4, position = position_jitter(height = 3.7, width = 0)) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.92, 0.13),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))
#------ClassII DVI------
GOI6 <- subset(Big_Data, log2FoldChange < -1)$Symbol

GOI6 <- unique(GOI6)

ClassIIDVI_subset <- Big_Data %>%
  filter(Class == "Class II") %>%
  filter(Symbol != "Olfr10G4") %>%
  dplyr::select(GeneType, Symbol, log2FoldChange, Corrected_padj, Genes, DVI) %>%
  unique()

ClassIIDVI_subset %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 1.5, alpha = 0.5) + facet_wrap(~DVI) +
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Most Class II Olfr mRNAs across DVIs are significantly reduced in 10G4 samples",
       subtitle = "The vast majority of significant increases come from the unclassified DVI category Low Expression",
       caption = "DEG limited to a change of at least one log2. To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-8, 8)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  #  geom_text(data = subset(ClassII_subset, Symbol %in% GOI4[str_detect(GOI4, "Olfr")]), 
  #            aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), 
  #            size = 4, position = position_jitter(height = 3.7, width = 0)) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.92, 0.08),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))

ClassIIDVIunusual_subset <- ClassIIDVI_subset %>%
  filter(DVI == "unusual") 

ClassIIDVIunusual_subset %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 2, alpha = 1) +
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Only about 50% of Class II Olfr mRNAs with unusual DVIs are significantly reduced in 10G4 samples",
       caption = "DEG limited to a change of at least one log2. To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_text(data = subset(ClassIIDVIunusual_subset, Symbol %in% GOI6[str_detect(GOI6, "Olfr")]), 
            aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), 
            size = 4) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.92, 0.08),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))
#------NonOlfr-NonTaar-Graphs----
GOI7 <- c("Gucy1b2", "Gucy2d")

OtherMS4a_subset <- Big_Data %>%
  filter(GeneType %in% c("Other", "MS4a")) %>%
  dplyr::select(GeneType, Symbol, log2FoldChange, Corrected_padj, Genes) %>%
  unique()

OtherMS4a_subset %>%
  #  filter(!is.na(Corrected_padj)) %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 1.5, alpha = 0.5) +
  geom_point(data = subset(OtherMS4a_subset, Symbol %in% GOI7), 
             aes(log2FoldChange, -log10(Corrected_padj)), color = "green") +
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Most non-Olfr/Taar mRNAs are not significantly affected in 10G4 samples",
       subtitle = "Colored Green are two Guanylyl Cyclase genes: Gucy1b2 and Gucy2d",
       caption = "All rows with missing padj were removed. DEG limited to a change of at least one log2. To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-10, 10)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  #    geom_text(data = subset(OtherMS4a_subset, Symbol %in% GOI7), 
  #              aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), 
  #              size = 4, position = position_jitter(height = 5, width = 2)) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.92, 0.92),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))


GOI8 <- c(subset(Big_Data, log2FoldChange < -5)$Symbol,
          subset(Big_Data, log2FoldChange > 5)$Symbol,
          subset(Big_Data, -log10(Corrected_padj) > 118)$Symbol) %>% 
  unique()

OtherMS4a_subset %>%
  #  filter(!is.na(Corrected_padj)) %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 1.5, alpha = 0.5) +
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Most non-Olfr/Taar mRNAs are not significantly affected in 10G4 samples",
       caption = "All rows with missing padj were removed. DEG limited to a change of at least one log2. To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-10, 10)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_text(data = subset(OtherMS4a_subset, Symbol %in% GOI8), 
            aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), 
            size = 4, position = position_jitter(height = 10, width = 0)) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.92, 0.92),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))

#------Unfiltered Graphs------
UniqueBigData <- Big_Data %>%
  dplyr::select(GeneType, Symbol, log2FoldChange, Corrected_padj, Genes) %>%
  unique()

GOI9 <- c(subset(Big_Data, log2FoldChange < -7.5)$Symbol,
          subset(Big_Data, log2FoldChange > 7.5)$Symbol,
          subset(Big_Data, log2FoldChange > 2.5 & -log10(Corrected_padj) > 45)$Symbol, 
          subset(Big_Data, log2FoldChange > -2 & -log10(Corrected_padj) > 200)$Symbol,
          subset(Big_Data, -log10(Corrected_padj) > 250)$Symbol) %>% 
  unique()

UniqueBigData %>%
  filter(Symbol != "Olfr10G4") %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 1.5, alpha = 0.8) +
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Olfr mRNA are clearly significantly reduced in 10G4 samples",
       caption = "All rows with missing padj were removed. DEG limited to a change of at least one log2. To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-10, 10)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_text(data = subset(UniqueBigData, Symbol %in% GOI9), 
            aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), 
            size = 4, position = position_jitter(height = 10, width = 0)) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.92, 0.5),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))
#------OlfrTaar------
OlfrTaar_subset <- Big_Data %>%
  filter(GeneType %in% c("Olfr", "Taar")) %>%
  dplyr::select(GeneType, Symbol, log2FoldChange, Corrected_padj, Genes) %>%
  unique()

GOI10 <- c(subset(Big_Data, log2FoldChange < -7.5)$Symbol,
           subset(Big_Data, log2FoldChange > 7.5)$Symbol,
           subset(Big_Data, log2FoldChange > 2.5 & -log10(Corrected_padj) > 45)$Symbol, 
           subset(Big_Data, -log10(Corrected_padj) > 250)$Symbol, 
           subset(Big_Data, GeneType == "Taar")$Symbol) %>% 
  unique()

OlfrTaar_subset %>%
  filter(Symbol != "Olfr10G4") %>%
  ggplot(aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Genes), shape = 20, size = 1.5, alpha = 0.5) +
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Olfr and Taar mRNA are significantly reduced in 10G4 samples",
       caption = "All rows with missing padj were removed. DEG limited to a change of at least one log2. To graph p-values of zero,  appx. 1e-281 was added to all Adjusted P-Values") +
  coord_cartesian(xlim = c(-10, 10)) +
  scale_colour_brewer(palette = "Set1") +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash", size = 0.4) +
  geom_text(data = subset(OlfrTaar_subset, Symbol %in% GOI10), 
            aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), 
            size = 3, position = position_jitter(height = 10, width = 0)) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.92, 0.5),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))
#------DVI-specific plots-----
Big_Data %>%
  filter(!is.na(DVI)) %>%
  filter(DVI %notin% c("low expression", "unusual")) %>%
  mutate(DVI = as.numeric(as.character(DVI))) %>%
  dplyr::select(GeneType, Symbol, log2FoldChange, Corrected_padj, Genes, DVI) %>%
  unique() %>%
  ggplot(aes(DVI, log2FoldChange)) + 
  geom_point(alpha = 0.75) + 
  geom_smooth(se = FALSE) +
  labs(x="Dorsal-medial/Ventral-lateral Zone Index (Inferred)", 
       y="Log2 Fold Change", 
       title = "Olfr and Taar mRNA L2FC values appear more significantly reduced in higher DVIs") +
  theme_bw() +
  theme(plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5))

#------WTcount, Mutantcount, L2FC-------
Big_Data %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, Top50WT, GeneType, Class) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  arrange(desc(MutantmeanCount)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  filter(Symbol %notin% c("Olfr10G4", "Olfr151")) %>%
  head(50) %>%
  mutate(Top50WT = ifelse(Top50WT == TRUE, "Top50WT", NA)) %>%
  ggplot(aes(Symbol, log2FoldChange, fill = Top50WT)) + 
  geom_col() +
  labs(y="Log2 Fold Change", 
       title = "The Log2 Fold Change of the Top 50 Olfrs by Mutant Sample count in decreasing order by Mutant Sample count", 
       subtitle = "In red, 21 of the Top 50 Olfrs by WT Sample count remain within the Top 50 Olfrs by Mutant Sample count") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 65, hjust = 1), 
        axis.text = element_text(size = 7), 
        axis.title.x = element_blank(), 
        legend.position = "none")

Top1 <- Big_Data %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, Top50WT, GeneType, Class) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  arrange(desc(WTmeanCount)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  filter(Symbol %notin% c("Olfr10G4", "Olfr151")) %>%
  head(50) %>%
  ggplot(aes(Symbol, WTmeanCount, fill = Class)) + 
  geom_col() +
  labs(title = "The Top 50 Olfrs by mean count from WT samples still experience a negative L2FC", 
       y = "Mean of the Raw Count in WT Samples") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45000)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 7),
        axis.ticks.x = element_blank(),
        legend.position = c(0.92, 0.88), 
        legend.title = element_blank())

Bottom1 <- Big_Data %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, Top50WT, GeneType, Class) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  arrange(desc(WTmeanCount)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  filter(Symbol %notin% c("Olfr10G4", "Olfr151")) %>%
  head(50) %>%
  ggplot(aes(Symbol, log2FoldChange, fill = Class)) + 
  geom_col() +
  labs(y = "Log2 Fold Change") +
  scale_y_continuous(expand = c(0, 0), limits = c(-8.5, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1), 
        axis.text = element_text(size = 7), 
        axis.title.x = element_blank(), 
        legend.position = "none")

plot_grid(Top1,
          Bottom1,
          nrow = 2, 
          align = "v")
##--
Top2 <- Big_Data %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, Top50WT, GeneType, Class) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  arrange(desc(WTmeanCount)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  filter(Symbol %notin% c("Olfr10G4", "Olfr151")) %>%
  head(100) %>%
  ggplot(aes(Symbol, WTmeanCount, fill = Class)) + 
  geom_col() +
  labs(title = "The Top 100 Olfrs by mean count from WT samples still experience a negative L2FC", 
       y = "Mean of the Raw Count in WT Samples") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45000)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 7),
        axis.ticks.x = element_blank(),
        legend.position = c(0.92, 0.88), 
        legend.title = element_blank())

Bottom2 <- Big_Data %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, Top50WT, GeneType, Class) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  arrange(desc(WTmeanCount)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  filter(Symbol %notin% c("Olfr10G4", "Olfr151")) %>%
  head(100) %>%
  ggplot(aes(Symbol, log2FoldChange, fill = Class)) + 
  geom_col() +
  labs(y = "Log2 Fold Change") +
  scale_y_continuous(expand = c(0, 0), limits = c(-8.5, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1), 
        axis.text = element_text(size = 7), 
        axis.title.x = element_blank(), 
        legend.position = "none")

plot_grid(Top2,
          Bottom2,
          nrow = 2, 
          align = "v")

##---
Top4 <- Big_Data %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, Top50WT, GeneType, Class) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  arrange(desc(WTmeanCount)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  filter(Symbol %notin% c("Olfr10G4", "Olfr151")) %>%
  head(100) %>%
  ggplot(aes(Symbol, WTmeanCount)) + 
  geom_col(alpha = 0.4, fill = "red") +
  geom_col(aes(Symbol, MutantmeanCount), alpha = 0.4, fill = "blue") +
  labs(title = "The Majority of Top 100 Olfrs by WT mean count experience major count reduction in Mutant Samples",
       subtitle = "The Top 100 Olfrs still experience a negative L2FC",
       y = "WT mean Count (Red), Mutant mean Count (Purple)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45000)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 7),
        legend.title = element_blank())

Bottom4 <- Big_Data %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, Top50WT, GeneType, Class) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  arrange(desc(WTmeanCount)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  filter(Symbol %notin% c("Olfr10G4", "Olfr151")) %>%
  head(100) %>%
  ggplot(aes(Symbol, log2FoldChange, fill = Class)) + 
  geom_col() +
  labs(y = "Log2 Fold Change") +
  scale_y_continuous(expand = c(0, 0), limits = c(-8.5, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1), 
        axis.text = element_text(size = 7), 
        axis.title.x = element_blank(), 
        legend.position = c(0.92, 0.18))

plot_grid(Top4,
          Bottom4,
          nrow = 2, 
          align = "v")

#---
Top5 <- Big_Data %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, Top50WT, GeneType, Class) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  arrange(desc(WTmeanCount)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  filter(Symbol %notin% c("Olfr10G4", "Olfr151")) %>%
  head(250) %>%
  ggplot(aes(Symbol, WTmeanCount)) + 
  geom_col(alpha = 0.4, fill = "red") +
  geom_col(aes(Symbol, MutantmeanCount), alpha = 0.4, fill = "blue") +
  labs(title = "The Majority of Top 250 Olfrs by WT mean count experience major count reduction in Mutant Samples",
       subtitle = "The Top 250 Olfrs still experience a negative L2FC",
       y = "WT mean Count (Red)/Mutant mean Count (Purple)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45000)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 7),
        legend.title = element_blank())

Bottom5 <- Big_Data %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, Top50WT, GeneType, Class) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  arrange(desc(WTmeanCount)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  filter(Symbol %notin% c("Olfr10G4", "Olfr151")) %>%
  head(250) %>%
  ggplot(aes(Symbol, log2FoldChange, fill = Class)) + 
  geom_col() +
  labs(y = "Log2 Fold Change") +
  scale_y_continuous(expand = c(0, 0), limits = c(-8.5, 0.6)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1), 
        axis.text = element_text(size = 7), 
        axis.title.x = element_blank(), 
        legend.position = c(0.92, 0.18))

plot_grid(Top5,
          Bottom5,
          nrow = 2, 
          align = "v")

#------Bigger Big Data-----
Smaller_Data_GreekIsland <- Big_Data %>%
  group_by(Symbol) %>%
  mutate(SmallestGIDistance = min(GI_Distance), 
         SmallestGIRankvsOlfrs = min(GI_Prox_RankvsOlfrs)) %>%
  ungroup() %>%
  group_by(Symbol) %>%
  arrange(GI_Distance) %>%
  slice_head() %>%
  ungroup() %>%
  dplyr::select(GeneID, Greek_Island, SmallestGIDistance, SmallestGIRankvsOlfrs) %>%
  dplyr::rename(Closest_GI = Greek_Island)

Big_Data2 <- Big_Data %>%
  left_join(Smaller_Data_GreekIsland, by = "GeneID")

Big_Data2 <- Big_Data2 %>%
  mutate(GI_Distance_Category = case_when(GI_Distance < 10000 ~ "<10k",
                                          GI_Distance < 100000 ~ "10k-99999",
                                          GI_Distance < 500000 ~ "100k-499999",
                                          GI_Distance < 1000000 ~ "500k-999999", 
                                          GI_Distance >= 1000000 ~ "1kk+"),
         Closest_GI_Distance_Category = case_when(SmallestGIDistance < 10000 ~ "<10k",
                                                  SmallestGIDistance < 100000 ~ "10k-99999",
                                                  SmallestGIDistance < 500000 ~ "100k-499999",
                                                  SmallestGIDistance < 1000000 ~ "500k-999999", 
                                                  SmallestGIDistance >= 1000000 ~ "1kk+"))

Big_Data2 <- Big_Data2 %>%
  mutate(Closest_GI_Distance_Category = factor(Closest_GI_Distance_Category, levels = c("<10k", "10k-99999", "100k-499999", "500k-999999", "1kk+")),
         GI_Distance_Category = factor(GI_Distance_Category, levels = c("<10k", "10k-99999", "100k-499999", "500k-999999", "1kk+")))

Big_Data3 <- Big_Data2 %>% 
  mutate(WTmeanCount_Category = case_when(WTmeanCount < 10 ~ "<10",
                                          WTmeanCount < 100 ~ "10-99",
                                          WTmeanCount < 250 ~ "100-249",
                                          WTmeanCount < 500 ~ "250-499",
                                          WTmeanCount < 1000 ~ "500-999",
                                          WTmeanCount < 5000 ~ "1k-4999",
                                          WTmeanCount < 10000 ~ "5k-9999",
                                          WTmeanCount >= 10000 ~ "10k+", 
                                          TRUE ~ "No_Counts_Detected"),
         MutantmeanCount_Category = case_when(MutantmeanCount < 10 ~ "<10",
                                              MutantmeanCount < 100 ~ "10-99",
                                              MutantmeanCount < 250 ~ "100-249",
                                              MutantmeanCount < 500 ~ "250-499",
                                              MutantmeanCount < 1000 ~ "500-999",
                                              MutantmeanCount < 5000 ~ "1k-4999",
                                              MutantmeanCount < 10000 ~ "5k-9999",
                                              MutantmeanCount >= 10000 ~ "10k+", 
                                              TRUE ~ "No_Counts_Detected"))         

Big_Data3 <- Big_Data3 %>%
  mutate(WTmeanCount_Category = factor(WTmeanCount_Category, levels = c("<10", "10-99", "100-249", "250-499", "500-999", "1k-4999", "5k-9999", "10k+", "No_Counts_Detected")),
         MutantmeanCount_Category = factor(MutantmeanCount_Category, levels = c("<10", "10-99", "100-249", "250-499", "500-999", "1k-4999", "5k-9999", "10k+", "No_Counts_Detected")))

#------Graphs with GI Info-----
Top3 <- Big_Data2 %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, GeneType, Closest_GI_Distance_Category) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  arrange(desc(WTmeanCount)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  filter(Symbol %notin% c("Olfr10G4", "Olfr151")) %>%
  head(100) %>%
  ggplot(aes(Symbol, WTmeanCount, fill = Closest_GI_Distance_Category)) + 
  geom_col() +
  labs(title = "The Top 100 Olfrs by mean count from WT samples still experience a negative L2FC",
       subtitle = "Greek Island Distance not an obvious indicator of resistance to change outside ~Top 10 Olfrs",
       y = "Mean of the Raw Count in WT Samples") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45000)) +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 7),
        axis.ticks.x = element_blank(),
        legend.position = c(0.85, 0.75))

Bottom3 <- Big_Data2 %>%
  dplyr::select(Symbol, WTmeanCount, MutantmeanCount, log2FoldChange, GeneType, Closest_GI_Distance_Category, Closest_GI) %>%
  unique() %>%
  filter(GeneType == "Olfr") %>%
  arrange(desc(WTmeanCount)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  filter(Symbol %notin% c("Olfr10G4", "Olfr151")) %>%
  head(100) %>%
  ggplot(aes(Symbol, log2FoldChange, fill = Closest_GI_Distance_Category)) + 
  geom_col() +
  geom_text(aes(label = Closest_GI), vjust = 0.5, hjust = 1.1, angle = 90, size = 3) +
  labs(y = "Log2 Fold Change") +
  scale_y_continuous(expand = c(0, 0), limits = c(-8.8, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1), 
        axis.text = element_text(size = 7), 
        axis.title.x = element_blank(), 
        legend.position = "none")

plot_grid(Top3,
          Bottom3,
          nrow = 2, 
          align = "v")

#---
#Proportion Plot: All Genes. Identifies shifts in Count Category across GI Distance for each Greek Island.
Big_Data3 %>% 
  ungroup() %>% 
  count(Chromosome, Greek_Island, GI_Distance_Category, WTmeanCount_Category) %>% 
  ggplot(aes(GI_Distance_Category, n, fill = WTmeanCount_Category)) + 
  geom_bar(position = "fill", stat = "identity") + 
  facet_wrap(~Greek_Island) +
  labs(y = "Proportion of Genes in each category grouping", 
       x = "Distance between a Gene and a Greek Island", 
       title = "An overview of the relationship between Greek Islands, the distance to a gene, and the count in WT samples") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(size = 7)) +
  guides(fill=guide_legend(title="WT Mean Count"))

#Proportion Plot: All Genes. Identifies shifts in Count Category across GI Distance for each Chromosome (Combined Greek Islands values).
Big_Data3 %>% 
  ungroup() %>% 
  count(Chromosome, Greek_Island, GI_Distance_Category, WTmeanCount_Category) %>% 
  filter(Chromosome %notin% c("GL456210.1", "GL456211.1", "GL456212.1", "GL456221.1", "JH584304.1", "JH584295.1")) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"))) %>%
  ggplot(aes(GI_Distance_Category, n, fill = WTmeanCount_Category)) + 
  geom_bar(position = "fill", stat = "identity") + 
  facet_wrap(~Chromosome) +
  labs(y = "Proportion of Genes in each category grouping", 
       x = "Distance between a Gene and a Greek Island", 
       title = "An overview of the relationship between Greek Islands, the distance to a gene, and the count in WT samples", 
       subtitle = "Segregated by Chromosome") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(size = 7)) +
  guides(fill=guide_legend(title="WT Mean Count"))

#Proportion Plot: Olfr Genes. Identifies shifts in Count Category across GI Distance for each Chromosome (Combined Greek Islands values).
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  ungroup() %>% 
  count(Chromosome, Greek_Island, GI_Distance_Category, WTmeanCount_Category) %>% 
  filter(Chromosome %notin% c("GL456210.1", "GL456211.1", "GL456212.1", "GL456221.1", "JH584304.1", "JH584295.1")) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"))) %>%
  ggplot(aes(GI_Distance_Category, n, fill = WTmeanCount_Category)) + 
  geom_bar(position = "fill", stat = "identity") + 
  facet_wrap(~Chromosome) +
  labs(y = "Proportion of Olfr Genes in each category grouping", 
       x = "Distance between an Olfr Gene and a Greek Island", 
       title = "An overview of the relationship between Greek Islands, the distance to an Olfr gene, and the count in WT samples", 
       subtitle = "Segregated by Chromosome") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(size = 7)) +
  guides(fill=guide_legend(title="WT Mean Count"))

#Proportion Plot: Olfr Genes. Identifies shifts in Count Category across GI Distance categories (Combined Greek Islands values).
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  ungroup() %>% 
  count(Chromosome, Greek_Island, GI_Distance_Category, WTmeanCount_Category) %>% 
  filter(Chromosome %notin% c("GL456210.1", "GL456211.1", "GL456212.1", "GL456221.1", "JH584304.1", "JH584295.1")) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"))) %>%
  ggplot(aes(GI_Distance_Category, n, fill = WTmeanCount_Category)) + 
  geom_bar(position = "fill", stat = "identity") + 
  labs(y = "Proportion of Olfr Genes in each category grouping", 
       x = "Distance between an Olfr Gene and a Greek Island", 
       title = "An overview of the relationship between Greek Islands, the distance to an Olfr gene, and the count in WT samples") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 8)) +
  guides(fill=guide_legend(title="WT Mean Count"))

#Proportion Plot: Olfr Genes. Identifies shifts in Count Category across GI Distance for each Greek Island.
Big_Data3 %>% 
  filter(GeneType == "Olfr") %>%
  ungroup() %>% 
  count(Chromosome, Greek_Island, GI_Distance_Category, WTmeanCount_Category) %>% 
  ggplot(aes(GI_Distance_Category, n, fill = WTmeanCount_Category)) + 
  geom_bar(position = "fill", stat = "identity") + 
  facet_wrap(~Greek_Island) +
  labs(y = "Proportion of Olfr genes in each category grouping", 
       x = "Distance between a Gene and a Greek Island", 
       title = "An overview of the relationship between Greek Islands, the distance to an Olfr gene, and the count in WT samples") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(size = 7)) +
  guides(fill=guide_legend(title="WT Mean Count"))

#L2FC Plot: All Genes. Identifies shifts in L2FC across Count Category x GI Distance for each Greek Island.
Big_Data3 %>% 
  ungroup() %>% 
  group_by(Chromosome, Greek_Island, GI_Distance_Category, WTmeanCount_Category) %>% 
  summarize(MeanL2FC = mean(log2FoldChange)) %>% 
  ggplot(aes(GI_Distance_Category, MeanL2FC, fill = WTmeanCount_Category)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  facet_wrap(~Greek_Island) +
  labs(y = "Mean Log2 Fold Change in each category grouping", 
       x = "Distance between a Gene and a Greek Island", 
       title = "An overview of the relationship between Greek Islands, the distance to a gene, WT mean Count, and the Log2 Fold Change in mRNA levels") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(size = 7)) +
  guides(fill=guide_legend(title="WT Mean Count"))

#L2FC Plot: Olfr Genes. Identifies shifts in L2FC across Count Category x GI Distance categories (Combined Greek Islands values).
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  ungroup() %>% 
  group_by(GI_Distance_Category, WTmeanCount_Category) %>% 
  summarize(MeanL2FC = mean(log2FoldChange)) %>% 
  ggplot(aes(WTmeanCount_Category, MeanL2FC, fill = GI_Distance_Category)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  labs(title = "Within a certain distance, Greek Islands might provide cis Olfrs with some resistance to 10G4 impact", 
       y = "Mean Log2 Fold Change in each category grouping", 
       x = "Mean WT count", 
       subtitle = "An overview of the relationship between Greek Islands, the distance to a gene, WT mean Count, and the Log2 Fold Change in Olfr mRNA levels") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(size = 7)) +
  guides(fill=guide_legend(title="Distance btwn Olfr/GI"))

#L2FC Plot: Olfr Genes. Identifies shifts in L2FC across Count Category x GI Distance for each Greek Island.
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  ungroup() %>% 
  group_by(Chromosome, Greek_Island, GI_Distance_Category, WTmeanCount_Category) %>% 
  summarize(MeanL2FC = mean(log2FoldChange)) %>% 
  ggplot(aes(GI_Distance_Category, MeanL2FC, fill = WTmeanCount_Category)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  facet_wrap(~Greek_Island) +
  labs(y = "Mean Log2 Fold Change in each category grouping", 
       x = "Distance between a Gene and a Greek Island", 
       title = "An overview of the relationship between Greek Islands, the distance to a Olfr gene, WT mean Count, and the mean Log2 Fold Change in mRNA levels") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(size = 7)) +
  guides(fill=guide_legend(title="WT Mean Count"))

#L2FC Plot: Olfr Genes. Rank within Top3 for GI Proximity vs Olfrs. Not plotting Count < 10. Identifies shifts in L2FC cross Count Category x GI Distance for each Greek Island.
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  filter(GI_Prox_RankvsOlfrs <= 3) %>%
  filter(WTmeanCount_Category != "<10") %>%
  ungroup() %>% 
  group_by(Greek_Island, GI_Distance_Category, WTmeanCount_Category) %>% 
  summarize(MeanL2FC = mean(log2FoldChange)) %>% 
  ggplot(aes(WTmeanCount_Category, MeanL2FC, fill = GI_Distance_Category)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  facet_wrap(~Greek_Island) +
labs(title = "Within certain categories, Greek Islands might provide Olfr genes with some resistance to 10G4 impact", 
     y = "Mean Log2 Fold Change in each category grouping", 
     x = "Mean WT count", 
     subtitle = "Comparisons limited to Mean WT Count >= 10 and within top 3 Olfrs for proximity to any Greek Island") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 7),
        legend.position = "bottom",
        legend.box = "horizontal") +
  guides(fill=guide_legend(title="Distance btwn Olfr/GI"))

#L2FC Plot: Olfr Genes. Top Rank for GI Proximity vs Olfrs. Not plotting Count < 10. Identifies shifts in L2FC across Count Category x GI Distance categories (Combined GI values).
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  filter(GI_Prox_RankvsOlfrs == 1) %>%
  filter(WTmeanCount_Category != "<10") %>%
  ungroup() %>% 
  group_by(Greek_Island, GI_Distance_Category, WTmeanCount_Category) %>% 
  summarize(MeanL2FC = mean(log2FoldChange)) %>% 
  ggplot(aes(WTmeanCount_Category, MeanL2FC, fill = GI_Distance_Category)) + 
  geom_bar(position = "dodge", stat = "identity") +
  labs(title = "Within certain categories, Greek Islands might provide Olfr genes with some resistance to 10G4 impact", 
       y = "Mean Log2 Fold Change in each category grouping", 
       x = "Mean WT count", 
       subtitle = "Comparisons limited to Mean WT Count >= 10 and Top Olfr Gene for proximity to any Greek Island") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 7),
        legend.position = c(0.92, 0.15)) +
  guides(fill=guide_legend(title="Distance btwn Olfr/GI"))

#Corrected Adjusted P value Plot: All Genes. Identifies shifts in CadjP across Count Category x GI Distance for each Greek Island.
Big_Data3 %>% 
  ungroup() %>% 
  group_by(Chromosome, Greek_Island, GI_Distance_Category, WTmeanCount_Category) %>% 
  summarize(MedianCpadj = median(Corrected_padj, na.rm = TRUE)) %>% 
  ggplot(aes(GI_Distance_Category, -log10(MedianCpadj), fill = WTmeanCount_Category)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  facet_wrap(~Greek_Island) +
  labs(y = "-Log10 of the Median Adjusted P value in each category grouping", 
       x = "Distance between a Gene and a Greek Island", 
       title = "An overview of the relationship between Greek Islands, the distance to a gene, WT mean Count, and the Median P value") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(size = 7)) +
  guides(fill=guide_legend(title="WT Mean Count"))

#Corrected Adjusted P value Plot: Olfr Genes. Identifies shifts in CadjP across Count Category x GI Distance for each Greek Island.
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  ungroup() %>% 
  group_by(Chromosome, Greek_Island, GI_Distance_Category, WTmeanCount_Category) %>% 
  summarize(MedianCpadj = median(Corrected_padj, na.rm = TRUE)) %>% 
  ggplot(aes(GI_Distance_Category, -log10(MedianCpadj), fill = WTmeanCount_Category)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  facet_wrap(~Greek_Island) +
  labs(y = "-Log10 of the Median Adjusted P value in each category grouping", 
       x = "Distance between an Olfr Gene and a Greek Island", 
       title = "An overview of the relationship between Greek Islands, the distance to a Olfr gene, WT mean Count, and the Median P value") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(size = 7)) +
  guides(fill=guide_legend(title="WT Mean Count"))

#WT mean Count: Olfr Genes. Plots WTmeanCount across Proximity Ranks for every Greek Island, with Chromosomes colored. 
GI_Resorted <- Big_Data3 %>% 
  dplyr::select(Greek_Island, Chromosome) %>% 
  unique() %>% 
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"))) %>% 
  arrange(Chromosome)

Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y")),
         Greek_Island = factor(Greek_Island, levels = GI_Resorted$Greek_Island)) %>%
  ggplot(aes(GI_Prox_RankvsOlfrs, log2(WTmeanCount), fill = Chromosome)) + 
  geom_col() + 
  facet_nested_wrap(vars(Chromosome, Greek_Island), dir = "v") +
  labs(title = "Only some Greek Islands seem to be related to higher Mean WT Count",
       y = "Log2 of Mean WT Count", 
       x = "Relative Proximity of an Olfr Gene to a specific Greek Island", 
       subtitle = "An overview of the relationship between Greek Islands, the relative proximity of Olfrs to Greek Islands and Mean WT Count") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 7), 
        legend.position = "none")

#As above, but filtered to Top 100 ranks
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  filter(GI_Prox_RankvsOlfrs <= 100) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y")),
         Greek_Island = factor(Greek_Island, levels = GI_Resorted$Greek_Island)) %>%
  ggplot(aes(GI_Prox_RankvsOlfrs, log2(WTmeanCount), fill = Chromosome)) + 
  geom_col() + 
  facet_nested_wrap(vars(Chromosome, Greek_Island), dir = "v") +
  labs(title = "Only some Greek Islands seem to be related to higher Mean WT Count",
       y = "Log2 of Mean WT Count", 
       x = "Relative Proximity of an Olfr Gene to a specific Greek Island", 
       subtitle = "An overview of the relationship between Greek Islands, the Top 100 relative proximity of Olfrs to Greek Islands and Mean WT Count") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 7), 
        legend.position = "none")

#As above, but filtered to Top 50 ranks
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  filter(GI_Prox_RankvsOlfrs <= 50) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y")),
         Greek_Island = factor(Greek_Island, levels = GI_Resorted$Greek_Island)) %>%
  ggplot(aes(GI_Prox_RankvsOlfrs, log2(WTmeanCount), fill = Chromosome)) + 
  geom_col() + 
  facet_nested_wrap(vars(Chromosome, Greek_Island), dir = "v") +
  labs(title = "Only some Greek Islands seem to be related to higher Mean WT Count",
       y = "Log 2 of Mean WT Count", 
       x = "Relative Proximity of an Olfr Gene to a specific Greek Island", 
       subtitle = "An overview of the relationship between Greek Islands, the Top 50 relative proximity of Olfrs to Greek Islands and Mean WT Count") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 7), 
        legend.position = "none")

#As above, but filtered to Top 50 ranks AND Recolored for actual GI distance category
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  filter(GI_Prox_RankvsOlfrs <= 50) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y")),
         Greek_Island = factor(Greek_Island, levels = GI_Resorted$Greek_Island)) %>%
  ggplot(aes(GI_Prox_RankvsOlfrs, log2(WTmeanCount), fill = GI_Distance_Category)) + 
  geom_col() + 
  facet_nested_wrap(vars(Chromosome, Greek_Island), dir = "v") +
  labs(title = "Only some Greek Islands seem to be related to higher Mean WT Count",
       y = "Log2 of Mean WT Count", 
       x = "Relative Proximity of an Olfr Gene to a specific Greek Island (with Chromosome)", 
       subtitle = "An overview of the relationship between Greek Islands, the Top 50 relative proximity of Olfrs to Greek Islands and Mean WT Count") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 7), 
        legend.position = c(0.95, 0.1)) +
  guides(fill=guide_legend(title="GI Distance"))

#As above, but filtered to Top 50 ranks AND Recolored for actual GI distance category AND filtered for GI Distance < 1kk
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  filter(GI_Prox_RankvsOlfrs <= 50) %>%
  filter(GI_Distance < 1000000) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y")),
         Greek_Island = factor(Greek_Island, levels = GI_Resorted$Greek_Island)) %>%
  ggplot(aes(GI_Prox_RankvsOlfrs, log2(WTmeanCount), fill = GI_Distance_Category)) + 
  geom_col() + 
  facet_nested_wrap(vars(Chromosome, Greek_Island), dir = "v") +
  labs(title = "Only some Greek Islands seem to be related to higher Mean WT Count",
       y = "Log2 of Mean WT Count", 
       x = "Relative Proximity of an Olfr Gene to a specific Greek Island (with Chromosome)", 
       subtitle = "An overview of the relationship between Greek Islands, the Top 50 relative proximity of Olfrs to Greek Islands and Mean WT Count") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 7), 
        legend.position = c(0.95, 0.1)) +
  guides(fill=guide_legend(title="GI Distance"))

#As above, but filtered to Top 50 ranks AND Recolored for actual GI distance category AND filtered for GI Distance < 1kk AND log2 the Mean WT Count
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  filter(GI_Prox_RankvsOlfrs <= 50) %>%
  filter(GI_Distance < 1000000) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y")),
         Greek_Island = factor(Greek_Island, levels = GI_Resorted$Greek_Island)) %>%
  ggplot(aes(GI_Prox_RankvsOlfrs, log2(WTmeanCount), fill = GI_Distance_Category)) + 
  geom_col() + 
  facet_nested_wrap(vars(Chromosome, Greek_Island), dir = "v") +
  labs(title = "Only some Greek Islands seem to be related to higher Mean WT Count",
       y = "Log2 of the Mean WT Count", 
       x = "Relative Proximity of an Olfr Gene to a specific Greek Island (with Chromosome)", 
       subtitle = "An overview of the relationship between Greek Islands, the Top 50 relative proximity of Olfrs to Greek Islands and Mean WT Count") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 7), 
        legend.position = c(0.95, 0.1)) +
  guides(fill=guide_legend(title="GI Distance"))

#Similar to above, but filtered to Top 50 ranks AND Recolored for actual GI distance category AND filtered for GI Distance < 1kk AND Mean WT Count > 10
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  filter(GI_Prox_RankvsOlfrs <= 50) %>%
  filter(GI_Distance < 1000000) %>%
  filter(WTmeanCount > 10) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y")),
         Greek_Island = factor(Greek_Island, levels = GI_Resorted$Greek_Island)) %>%
  ggplot(aes(GI_Prox_RankvsOlfrs, log2(WTmeanCount), fill = GI_Distance_Category)) + 
  geom_col() + 
  facet_nested_wrap(vars(Chromosome, Greek_Island), dir = "v") +
  labs(title = "Only some Greek Islands seem to be related to higher Mean WT Count",
       y = "Log2 of Mean WT Count", 
       x = "Relative Proximity of an Olfr Gene to a specific Greek Island (with Chromosome)", 
       subtitle = "An overview of the relationship between Greek Islands, the Top 50 relative proximity of Olfrs to Greek Islands and Mean WT Count (>10)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 7), 
        legend.position = c(0.95, 0.1)) +
  guides(fill=guide_legend(title="GI Distance"))

#Similar to above, but filtered to Top 30 ranks AND Recolored for actual GI distance category AND filtered for GI Distance < 1kk AND Mean WT Count > 100
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  filter(GI_Prox_RankvsOlfrs <= 30) %>%
  filter(GI_Distance < 1000000) %>%
  filter(WTmeanCount > 100) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y")),
         Greek_Island = factor(Greek_Island, levels = GI_Resorted$Greek_Island)) %>%
  ggplot(aes(GI_Prox_RankvsOlfrs, log2(WTmeanCount), fill = GI_Distance_Category)) + 
  geom_col() + 
  facet_nested_wrap(vars(Chromosome, Greek_Island), dir = "v") +
  labs(title = "Only some Greek Islands seem to be related to higher Mean WT Count",
       y = "Log2 of Mean WT Count", 
       x = "Relative Proximity of an Olfr Gene to a specific Greek Island (with Chromosome)", 
       subtitle = "An overview of the relationship between Greek Islands, the Top 30 relative proximity of Olfrs to Greek Islands and Mean WT Count (>100)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 7), 
        legend.position = c(0.95, 0.1)) +
  guides(fill=guide_legend(title="GI Distance"))

#Similar to above, but filtered to Top 30 ranks AND Recolored for actual GI distance category AND filtered for GI Distance < 1kk AND Mean WT Count > 250 and log2 of MeanWTCount
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  filter(GI_Prox_RankvsOlfrs <= 30) %>%
  filter(GI_Distance < 1000000) %>%
  filter(WTmeanCount > 250) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y")),
         Greek_Island = factor(Greek_Island, levels = GI_Resorted$Greek_Island)) %>%
  ggplot(aes(GI_Prox_RankvsOlfrs, log2(WTmeanCount), fill = GI_Distance_Category)) + 
  geom_col() + 
  facet_nested_wrap(vars(Chromosome, Greek_Island), dir = "v") +
  labs(title = "Only some Greek Islands seem to be related to higher Mean WT Count (< 1 Million bp Away)",
       y = "Log2 of the Mean WT Count", 
       x = "Relative Proximity of an Olfr Gene to a specific Greek Island (with Chromosome)", 
       subtitle = "An overview of the relationship between Greek Islands, the Top 30 relative proximity of Olfrs to Greek Islands and Mean WT Count (>250)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 7), 
        legend.position = c(0.95, 0.1)) +
  guides(fill=guide_legend(title="GI Distance"))

#Similar to above, but filtered to Top 30 ranks AND Recolored for actual GI distance category AND Mean WT Count > 250 and log2 of MeanWTCount
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  filter(GI_Prox_RankvsOlfrs <= 30) %>%
  filter(WTmeanCount > 250) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y")),
         Greek_Island = factor(Greek_Island, levels = GI_Resorted$Greek_Island)) %>%
  ggplot(aes(GI_Prox_RankvsOlfrs, log2(WTmeanCount), fill = GI_Distance_Category)) + 
  geom_col() + 
  facet_nested_wrap(vars(Chromosome, Greek_Island), dir = "v") +
  labs(title = "Only some Greek Islands seem to be related to higher Mean WT Count",
       y = "Log2 of the Mean WT Count", 
       x = "Relative Proximity of an Olfr Gene to a specific Greek Island (with Chromosome)", 
       subtitle = "An overview of the relationship between Greek Islands, the Top 30 relative proximity of Olfrs to Greek Islands and Mean WT Count (>250)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 7), 
        legend.position = c(0.95, 0.1)) +
  guides(fill=guide_legend(title="GI Distance"))

#Similar to above, but filtered to Top 30 ranks AND Recolored for actual GI distance category AND Mean WT Count > 250
Big_Data3 %>%
  filter(GeneType == "Olfr") %>%
  filter(GI_Prox_RankvsOlfrs <= 30) %>%
  filter(WTmeanCount > 250) %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y")),
         Greek_Island = factor(Greek_Island, levels = GI_Resorted$Greek_Island)) %>%
  ggplot(aes(GI_Prox_RankvsOlfrs, log2(WTmeanCount), fill = GI_Distance_Category)) + 
  geom_col() + 
  facet_nested_wrap(vars(Chromosome, Greek_Island), dir = "v") +
  labs(title = "Only some Greek Islands seem to be related to higher Mean WT Count",
       y = "Log2 of Mean WT Count", 
       x = "Relative Proximity of an Olfr Gene to a specific Greek Island (with Chromosome)", 
       subtitle = "An overview of the relationship between Greek Islands, the Top 30 relative proximity of Olfrs to Greek Islands and Mean WT Count (>250)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 7), 
        legend.position = c(0.95, 0.1)) +
  guides(fill=guide_legend(title="GI Distance"))
# Other graphs------
Big_Data %>% 
  filter(Top50WT == TRUE) %>% 
  dplyr::select(WTmeanCount, MutantmeanCount, Symbol, Class) %>% 
  unique() %>% 
  arrange(desc(WTmeanCount)) %>% 
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>% 
  ggplot(aes(Symbol, WTmeanCount)) + 
  geom_col(fill = "blue", alpha = 0.5) + 
  geom_col(aes(Symbol, MutantmeanCount), fill = "red", alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Big_Data %>% 
  filter(GeneType == "Olfr", GI_Prox_RankvsOlfrs <= 3) %>%
  filter(Symbol != "Olfr10G4") %>%
  dplyr::select(Top50WT, log2FoldChange, WTmeanCount, MutantmeanCount, Symbol, Class) %>% 
  unique() %>% 
  arrange(desc(WTmeanCount)) %>% 
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>% 
  ggplot(aes(Symbol, WTmeanCount, fill = Class)) + 
  geom_col() + 
  labs(title = "Olfr mRNA Raw Count, sorted by descending raw count from WT samples, Identified as within 3 Olfrs to an Greek Island") +
  theme(axis.text.x=element_blank(), axis.title.x = element_blank())

Big_Data %>% 
  filter(GeneType == "Olfr", GI_Prox_RankvsOlfrs <= 3) %>%
  filter(Symbol != "Olfr10G4") %>%
  dplyr::select(Top50WT, log2FoldChange, WTmeanCount, MutantmeanCount, Symbol, Class) %>% 
  unique() %>% 
  arrange(desc(WTmeanCount)) %>% 
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>% 
  ggplot(aes(Symbol, log2FoldChange, fill = Class)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#------Finalize Plots----
#Relabeling Results
a5 <- a1 %>% mutate(Result = case_when(
  GeneType == "Olfr" & padj < 0.05 & log2FoldChange < -1 ~ "Olfr DEG",
  GeneType == "Olfr" & padj < 0.05 & log2FoldChange > 1 ~ "Olfr DEG",
  padj < 0.05 & log2FoldChange < -1~ "DEG",
  padj < 0.05 & log2FoldChange > 1~ "DEG",
  TRUE ~ "non-DEG"))

a5 <- a5 %>%
  filter(Symbol != "Olfr10G4") %>%
  mutate(Result = factor(Result, levels = c("Olfr DEG", "DEG", "non-DEG")))

RT_Olfrs <- c("Olfr358", "Olfr390", "Olfr510", "Olfr596", "Olfr603", "Olfr690", "Olfr1154")

#A labeled volcano plot using the RT-qPCR targets
Volcano <- ggplot(a5, aes(x= log2FoldChange, y = -log10(Corrected_padj))) + 
  geom_point(aes(color = Result), shape = 20, size = 0.8, alpha = 0.80) + 
  labs(x="Log2 Fold Change", 
       y="-log10(Adjusted P-Value)", 
       title = "Olfr mRNAs are significantly reduced in the original 10G4 line 7 samples") +
  coord_cartesian(xlim = c(-7.8, 7.8), ylim = c(0,280)) +
  scale_color_manual(values=c("skyblue1", "red3", "yellow2")) +
  geom_vline(xintercept = c(-1, 1), linetype = "longdash") +
  geom_point(data = subset(a5, Symbol %in% RT_Olfrs),
             shape = 5, size = 2) +
  geom_text_repel(data = subset(a5, Symbol %in% RT_Olfrs), 
                  aes(log2FoldChange, -log10(Corrected_padj), label = Symbol), 
                  size = 4) +
  theme_bw() + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.85, 0.15),
        legend.background = element_rect(fill = "white", color = "black"), 
        plot.title = element_text(size = 20, hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size=7)))

a6 <- a5 %>% mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                               "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                                               "X", "Y", "GL456211.1", "GL456221.1", "JH584295.1", 
                                                               "GL456212.1", "GL456210.1", "JH584304.1")))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


a6 %>%
  filter(GeneType == "Olfr") %>%
  arrange(Chromosome, Start) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  ggplot(aes(Symbol, log2FoldChange, color = Chromosome)) + 
  geom_point(size = 3, shape = 21) +
  geom_hline(yintercept = 1, linetype = "dashed", col = 'red') +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -1, linetype = "dashed", col = 'red') +
  scale_colour_manual(values=rep(cbPalette, 3)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=5)))

#Attempting to reproduce plots more similar to paul's
Smaller_Data3 <- Big_Data4 %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                    "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                                    "X", "Y", "GL456211.1", "GL456221.1", "JH584295.1", 
                                                    "GL456212.1", "GL456210.1", "JH584304.1"))) %>%
  arrange(Chromosome, Start) %>%
  dplyr::select(Genes, GeneType, Symbol, log2FoldChange, Chromosome, DVI, Class, Corrected_padj) %>%
  filter(GeneType %in% c("Olfr", "Taar")) %>%
  unique() %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  mutate(Result = case_when(Genes == "Olfr DEG" & Class == "Class I" ~ "Class I DEG",
                            Genes == "non-DEG" & Class == "Class I" ~ "Class I non-DEG",
                            Genes == "Olfr DEG" & Class == "Class II" ~ "Class II DEG",
                            Genes == "non-DEG" & Class == "Class II" ~ "Class II non-DEG",
                            Genes == "Taar DEG" & Class == "Class TAAR" ~ "Taar DEG",
                            Genes == "non-DEG" & Class == "Class TAAR" ~ "Taar non-DEG"
  )) %>%
  mutate(Result = factor(Result, levels = c("Class I DEG", "Class II DEG", "Taar DEG", "Class I non-DEG", "Class II non-DEG", "Taar non-DEG"))) 

#This graph attempts to reproduce a similar plot as in Khan 2011 and Paul's paper, minus the DVI and class detail and p-value (correction needed?)
Smaller_Data3 %>%
  ggplot(aes(Symbol, log2FoldChange, color = Chromosome)) + 
  geom_point(size = 3, shape = 21) +
  geom_hline(yintercept = 1, linetype = "dashed", col = 'red') +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -1, linetype = "dashed", col = 'red') +
  scale_colour_manual(values=rep(cbPalette, 3)) +
  labs(title = "Log2 Fold Change of Olfr genes") +
  geom_point(data = subset(Smaller_Data3, Symbol %in% RT_Olfrs),
             aes(Symbol, log2FoldChange),
             shape = 16, size = 2, color = "black") +
  geom_text_repel(data = subset(Smaller_Data3, Symbol %in% RT_Olfrs), 
                  aes(Symbol, log2FoldChange, label = Symbol), 
                  size = 4, color = 'black') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.line = element_line(colour = "black")) +
  guides(colour = guide_legend(override.aes = list(size=5)))

Smaller_Data4 <- Big_Data4 %>%
  mutate(Chromosome = factor(Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                    "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                                    "X", "Y", "GL456211.1", "GL456221.1", "JH584295.1", 
                                                    "GL456212.1", "GL456210.1", "JH584304.1"))) %>%
  arrange(Chromosome, Start) %>%
  dplyr::select(Genes, GeneType, Symbol, log2FoldChange, Chromosome, Start, DVI, Class, Corrected_padj) %>%
  filter(GeneType %in% c("Olfr", "Taar")) %>%
  unique() %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%
  mutate(Result = case_when(Genes == "Olfr DEG" & Class == "Class I" ~ "Class I DEG",
                            Genes == "Olfr DEG" & Class == "Class II" ~ "Class II DEG",
                            Genes == "Taar DEG" & Class == "Class TAAR" ~ "Taar DEG",
                            Genes == "non-DEG" & Class == "Class I" ~ "Class I non-DEG",
                            Genes == "non-DEG" & Class == "Class II" ~ "Class II non-DEG",
                            Genes == "non-DEG" & Class == "Class TAAR" ~ "Taar non-DEG")) %>%
  mutate(Result = factor(Result, levels = c("Class I DEG", "Class I non-DEG", "Class II DEG", "Class II non-DEG", "Taar DEG", "Taar non-DEG"))) %>%
  arrange(Chromosome, Start)

cbbPalette2 <- c("#56B4E9", "#E69F00", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ChromosomeL2FC <- Smaller_Data4 %>%
  ggplot(aes(Chromosome, log2FoldChange, color = Result, shape = Result)) + 
  geom_jitter(size = 1.3) +
  scale_color_manual(values=c("skyblue1", "skyblue1", "black", "black", "red2", "red2")) +
  scale_shape_manual(values=c(16, 1, 16, 1, 16, 1)) +
  labs(title = "Log2 Fold Change of Olfr genes in 10G4 Line 7 by Chromosome and Class", 
       y = "Log2 Fold Change") +
  scale_y_continuous(breaks=seq(-8,4,2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
        axis.text = element_text(size = 15, face = "bold"), 
        axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.5, 0.04), 
        legend.box = "horizontal",
        legend.key = element_blank(),
        legend.box.background = element_rect(color="black", size=1.2)) +
  geom_vline(xintercept = 1.5:18.5, size=0.25) +
  geom_hline(yintercept = 1, linetype = "dashed", col = 'red', size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2) +
  geom_hline(yintercept = -1, linetype = "dashed", col = 'red', size = 1.2) + 
  guides(colour = guide_legend(ncol = 6, override.aes = list(size = 4)))