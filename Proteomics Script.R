# This is a script to generate Proteomics Figures and Tables for the 10G4 Paper.

# Set your working directory to the "Proteomics 10G4" folder

# 01 - Negated %in% filter function----
`%notin%` <- Negate(`%in%`)

# 02 - Required packages----
if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

if (!require(ggpubr)) install.packages('ggpubr')
library(ggpubr)

if (!require(data.table)) install.packages('data.table')
library(data.table)

if (!require(grid)) install.packages('grid')
library(grid)

if (!require(gridExtra)) install.packages('gridExtra')
library(gridExtra)

if (!require(patchwork)) install.packages('patchwork')
library(patchwork)

if (!require(igraph)) install.packages('igraph')
library(igraph)

if (!require(ggraph)) install.packages('ggraph')
library(ggraph)

if (!require(ggrepel)) install.packages('ggrepel')
library(ggrepel)

if (!require(gtable)) install.packages('gtable')
library(gtable)

if (!require(ggridges)) install.packages('ggridges')
library(ggridges)

library(scales)

library(ggVennDiagram)

# library(ggplotify) #Used to make panels, might be obsolete with current approach
# 
# library(VennDiagram) #Obsolete due to not working as well with 2 groups
# library(RColorBrewer) #Used for color purposes with above


# 03 - Importing the datasets----

#The somewhat processed peptide report
PF1586_peptide_report <- read_csv("Input/PF1586_peptide_report.csv")

#The full processed protein report. Makes unknown decisions about aggregation and specific gene
PF1586_Report <- read_csv("Input/PF1586_Report.csv")

#Information about many different receptors from uniprot
Protein_seqs <- read_delim("uniprotkb_Olfactory_receptor_AND_mus_mu_2024_10_14.tsv", 
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# 04 - Quality Checks----
PF1586_peptide_report %>%
  pivot_longer(cols = 7:18, names_to = "Sample", values_to = "Count") %>%
  filter(!str_detect(Sample, "Homo|2006")) %>% #Removes Hom samples and also Sample 2006
  mutate(Sample = ifelse(str_detect(Sample, "10G4"), str_remove(Sample, "Het"), Sample)) %>% #renames 10G47Het as 10G4
  ggplot() +
  geom_violin(aes(Sample, log10(Count + 1))) +
  geom_boxplot(aes(Sample, log10(Count + 1)), width = 0.2) +
  labs(title = "Combination boxplot/violinplot show that all samples have a similar spread of values") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))

PR <- PF1586_peptide_report %>%
  pivot_longer(cols = 7:18, names_to = "Sample", values_to = "Count") %>%
  filter(!str_detect(Sample, "Homo|2006"))

compare_means(Count ~ Sample, data = PR, method = "t.test")

ggsave("Output2/Overview_QC1_V2_SF10A.png", width = 30, height = 20, units = "cm")

PR2 <- PR %>%
  mutate(Present = Count > 0) %>%
  group_by(Sequence) %>%
  mutate(PresentTotal = sum(Present)) %>%
  ungroup()

PR2 %>%
  select(Sequence, PresentTotal) %>%
  unique() %>%
  count(PresentTotal)

PR2 %>%
  mutate(Sample = ifelse(str_detect(Sample, "10G4"), str_remove(Sample, "Het"), Sample)) %>%
  mutate(Zero_Check = ifelse(Count == 0, "Zero", "Non-zero")) %>%
  ggplot() +
  geom_bar(aes(x = Sample, fill = Zero_Check), stat="count") +
  labs(title = "Evaluating undetected peptide in each sample versus the collective peptide pool") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))

ggsave("Output2/Overview_QC2_V2_SF10B.png", width = 30, height = 20, units = "cm")

PR2 %>%
  mutate(Sample = ifelse(str_detect(Sample, "10G4"), str_remove(Sample, "Het"), Sample)) %>%
  ungroup()%>%
  arrange(PresentTotal, Count) %>%
  mutate(Sequence = factor(Sequence, levels = unique(Sequence))) %>%
  ggplot(aes(x = Sequence, y = Sample, fill = Present)) +
  geom_tile() +
  labs(title = "A sequence-specific comparison of detection for each sample") +
  scale_fill_manual(values = c("FALSE" = "grey90", "TRUE" = "steelblue"),
                    name = "Detected") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", size = 40, hjust = 0.5),
        axis.text = element_text(size = 30),
        axis.title = element_text(size = 35),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30))

ggsave("Output2/Overview_QC3_SF10C.png", width = 100, height = 40, units = "cm")

# 05 - OR10G4 and Olfrs----

#This pulls all OR/OR10G4 proteins detected in the aggregated protein file

#Total protein counts detected... not total peptide counts, so mind the mismatch
#A published figure below
Protein_sums <- PF1586_Report %>%
  select(1:18) %>%
  pivot_longer(cols = 7:18, names_to = "Sample", values_to = "Count") %>%
  group_by(Sample) %>% 
  summarise(Sum_Count = sum(Count, na.rm = TRUE)) %>% 
  separate(Sample, c("Genotype", "Sample")) %>%
  mutate(Genotype = str_remove(Genotype, "7Het"))

PF1586_Report %>%
  select(1:18) %>%
  pivot_longer(cols = 7:18, names_to = "Sample", values_to = "Count") %>% 
  filter(!str_detect(Sample, "Homo|2006")) %>%
  mutate(Sample = ifelse(str_detect(Sample, "10G4"), str_remove(Sample, "7Het"), Sample)) %>%
  filter(Genes %in% ORs) %>%
  mutate(OR_Group = ifelse(Genes == "OR10G4", "OR10G4", "Olfr")) %>%
  filter(!(OR_Group == "OR10G4" & str_detect(Sample, "WT"))) %>%
  group_by(Sample, OR_Group) %>%
  summarise(OR_sums = sum(Count, na.rm = TRUE)) %>% 
  separate(Sample, c("Genotype", "Sample")) %>%
  left_join(Protein_sums, by = c("Genotype", "Sample")) %>%
  mutate(OR_prop = OR_sums/Sum_Count * 100) %>%
  mutate(Sample2 = paste0(Genotype, " (", Sample, ")")) %>%
  ggplot(aes(Sample2, OR_prop, fill = OR_Group)) +
  geom_bar(position="stack", stat="identity") +
  labs(y = "Percent of Total Protein Counts",
       title = "OR10G4 partially compensates for loss of Olfr counts",
       x = "Samples") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 13, face = "bold"),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.8),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) + 
  guides(fill = guide_legend(nrow = 1)) + 
  scale_fill_manual(values = c("red2", "green4")) +
  scale_y_continuous(breaks = seq(0:0.15, by = 0.01))

ggsave("Output2/F07_OR10G4_Olfr_SinglePlot_MF4B.png", width = 20, height = 15, units = "cm")


  
Olfr_prop <- PF1586_Report %>%
  select(1:18) %>%
  pivot_longer(cols = 7:18, names_to = "Sample", values_to = "Count") %>%
  filter(!str_detect(Sample, "Homo|2006")) %>%
  mutate(Sample = ifelse(str_detect(Sample, "10G4"), str_remove(Sample, "7Het"), Sample)) %>%
  filter(str_detect(Genes, "Olfr")) %>%
  group_by(Sample) %>% 
  summarise(Olfr_Count = sum(Count, na.rm = TRUE)) %>% 
  separate(Sample, c("Genotype", "Sample")) %>%
  left_join(Protein_sums, by =  c("Genotype", "Sample")) %>%
  mutate(Olfr_Prop = Olfr_Count/Sum_Count) 

#OR10G4_prop unused due to removal of hom
OR10G4_Prop <- PF1586_Report %>%
  select(1:18) %>%
  pivot_longer(cols = 7:18, names_to = "Sample", values_to = "Count") %>%
  filter(!str_detect(Sample, "Homo|2006")) %>%
  mutate(Sample = ifelse(str_detect(Sample, "10G4"), str_remove(Sample, "7Het"), Sample)) %>%
  filter(str_detect(Genes, "10G4")) %>%
  filter(!str_detect(Sample, "2006|WT")) %>%
  separate(Sample, c("Genotype", "Sample")) %>%
  left_join(Protein_sums, by =  c("Genotype", "Sample")) %>%
  mutate(OR10G4_Prop = Count/Sum_Count) 

Olfr_mean_prop <- Olfr_prop %>%
  group_by(Genotype) %>%
  summarise(Mean_OP = mean(Olfr_Prop))

WT_omp <- Olfr_mean_prop %>%
  filter(Genotype == "WT") %>%
  rename(Mean_WTOP = Mean_OP)

NOP <- Olfr_prop %>%
  bind_cols(WT_omp$Mean_WTOP) %>%
  rename(Mean_WTop = `...6`) %>%
  mutate(N_Olfr_Prop = Olfr_Prop/Mean_WTop) #This is normalized to the mean of WT Olfr_Prop... unnecessary here

NOP2 <- NOP %>%
  filter(Sample != "2006") %>%
  mutate(Percent = Olfr_Prop * 100)

compare_means(Percent ~ Genotype, NOP2, method = "t.test")

NOP2 %>%
  ggplot(aes(Genotype, Percent)) +
  geom_point(size = 5, alpha = 0.85) +
  geom_signif(y_position = c(0.148), xmin = c(1), xmax = c(2),
              annotation = c("*"), tip_length = 0, textsize = 12) +
  coord_cartesian(ylim = c(0.01, 0.17)) +
  labs(y = "Percent of Total Protein Counts",
       title = "OR10G4 samples have less Olfr Protein",
       x = "Samples") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 13, face = "bold"),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.8),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) 

ggsave("Output2/F02_OR10G4_Olfr_v2_MF4A.png", width = 15, height = 15, units = "cm")

# 06 - Looking at sequences----

#Protein_seqs contain sequence information about Receptor proteins from uniprot
#This file is filtered specifically for Olfr protein, see Olfr_protein
#Olfr_protein is further filtered to only keep ORs that have a matching sequence from mass spec aka Detected_Olfrs

Olfr_protein <- Protein_seqs %>%
  filter(str_detect(`Protein names`, "Olfactory receptor"))

#Keeps only Olfr proteins that were possibly detected via mass spec fragments
#What a pain in the... code. 235 groups of ORs, 246 unique OR proteins.
#Here's the catch, 5 of those also had an -ps1 suffix and regular names so got combined, not Olfr1534, only ps1, hence the mismatch
Detected_Olfrs <- Olfr_protein %>%
  rowwise() %>%
  filter(any(str_detect(Sequence, PF1586_peptide_report$Sequence))) %>%
  ungroup() #Also eliminates any sequences from other mouse species

#Taking those 60k peptides and asking which genotype detects it
Genotype_detected <- PF1586_peptide_report %>%
  pivot_longer(cols = 7:18, names_to = "Sample", values_to = "Count") %>%
  filter(!str_detect(Sample, "Homo|2006")) %>%
  mutate(Sample = ifelse(str_detect(Sample, "10G4"), str_remove(Sample, "7Het"), Sample)) %>%
  separate(Sample, c("Genotype", "Sample")) %>%
  group_by(Sequence, Genotype) %>%
  summarise(NonZero = sum(Count) > 0)

#Just the peptide sequences detected. updated, but seems like consistent? 60426 in both cases
peptides <- PF1586_peptide_report %>%
  pivot_longer(cols = 7:18, names_to = "Sample", values_to = "Count") %>%
  filter(!str_detect(Sample, "Homo|2006")) %>%
  pull(Sequence) %>%
  unique()

#I do something weird here and for every row, I add a filtered list of peptides
#based on whether the peptide is detected in the uniprot-derived detected OR protein sequence. 
#I used the peptides to keep only possibly detected Olfrs, and now I am adding what those peptides are to each Olfr.
D_Ov2 <- Detected_Olfrs %>%
  rowwise() %>%
  mutate(MS_peptides = list(peptides[str_detect(Sequence, peptides)])) %>%
  ungroup() %>%
  unnest(MS_peptides) 

#The Venn Diagrams

#This will be a file that checks on the distribution of split ORs (77) from the PF_Report file.
Protein_77_ORs_Nested <- PF1586_Report %>%
  select(1:18) %>%
  pivot_longer(cols = 7:18, names_to = "Sample", values_to = "Count") %>%
  filter(!str_detect(Sample, "Homo|2006")) %>%
  mutate(Sample = ifelse(str_detect(Sample, "10G4"), str_remove(Sample, "7Het"), Sample)) %>%
  separate(Sample, c("Genotype", "Sample")) %>%
  filter(str_detect(Genes, "Olfr")) %>%
  mutate(GeneNamesSeparated = str_split(Genes, ";")) %>%
  unnest(GeneNamesSeparated) %>% 
  mutate(Count = ifelse(is.na(Count), 0, Count)) %>%
  group_by(Genotype, GeneNamesSeparated) %>%
  summarize(AbsentOR = sum(Count) == 0) %>%
  filter(AbsentOR == FALSE) %>%
  select(Genotype, GeneNamesSeparated) %>%
  group_by(Genotype) %>%
  nest()

x <- list(OR10G4 = Protein_77_ORs_Nested[[2]][[1]]$GeneNamesSeparated,
          WT = Protein_77_ORs_Nested[[2]][[2]]$GeneNamesSeparated)

ggVennDiagram(x, label_alpha = 0, label_size = 15, set_size = 15) + 
  scale_fill_gradient(low="grey90",high = "red") +
  labs(title = "OR10G4 has 49% fewer unique OR Proteins than WT") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "none") +
  coord_flip()

ggsave("Output2/F5_Protein_77_ORs_Nested_MF4C.png", width = 25, height = 20, units = "cm")


PeptideSeqs_103_seqs_Nested <- D_Ov2 %>%
  left_join(Genotype_detected, by = c("MS_peptides" = "Sequence")) %>%
  filter(NonZero == TRUE) %>%
  filter(Length > 276 & str_detect(Sequence, "^M")) %>%
  select(Genotype, MS_peptides) %>%
  unique() %>%
  group_by(Genotype) %>%
  nest()

x <- list(OR10G4 = PeptideSeqs_103_seqs_Nested[[2]][[1]]$MS_peptides,
          WT = PeptideSeqs_103_seqs_Nested[[2]][[2]]$MS_peptides)

ggVennDiagram(x, label_alpha = 0, label_size = 15, set_size = 15) + 
  scale_fill_gradient(low="grey90",high = "red") +
  labs(title = "OR10G4 has 56% fewer unique OR Peptide kinds than WT") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "none") +
  coord_flip()

ggsave("Output2/F5_PeptideSeqs_103_seqs_SF10D.png", width = 25, height = 20, units = "cm")

MSpeptides_uniProteins <- D_Ov2 %>%
  mutate(GeneNamesSeparated = str_split(`Gene Names`, " ")) %>%
  unnest(GeneNamesSeparated) %>%
  filter(str_starts(GeneNamesSeparated, "Olfr")) %>%
  filter(Length > 276 & str_detect(Sequence, "^M")) %>%
  select(MS_peptides, GeneNamesSeparated) %>%
  unique() %>%
  group_by(MS_peptides) %>%
  summarize(Combined_Olfrs = str_c(GeneNamesSeparated, collapse = ";")) %>%
  ungroup()

#Seems like I shift around MS and ditch 2006? Info is made
MSuni_Counts <- MSpeptides_uniProteins %>%
  mutate(Number = str_count(Combined_Olfrs, "Olfr"), 
         Info = paste0(MS_peptides, " (Olfrs: ", Number, ")")) %>%
  left_join(PF1586_peptide_report, by = c("MS_peptides" = "Sequence")) %>%
  pivot_longer(cols = 10:21, names_to = "Sample", values_to = "Count") %>%
  filter(!str_detect(Sample, "Homo|2006")) %>%
  mutate(Sample = ifelse(str_detect(Sample, "10G4"), str_remove(Sample, "7Het"), Sample)) %>%
  separate(Sample, c("Genotype", "Sample")) %>%
  left_join(Protein_sums, by = c("Sample", "Genotype")) %>%
  mutate(Percent_Total = Count/Sum_Count * 100)

CompareMeans_Info <- compare_means(Percent_Total ~ Genotype, MSuni_Counts, method = "t.test", group = "Info") %>%
  arrange(p) %>%
  select(Info, p, p.signif)

MSuni_Counts %>%
  left_join(CompareMeans_Info, by = "Info") %>%
  arrange(p) %>%
  mutate(Info = factor(Info, levels = unique(Info))) %>%
  ggplot() +
  geom_boxplot(aes(Info, Percent_Total, color = Genotype), alpha = 0.7, position = position_dodge(width = 0.8)) +  
  geom_point(aes(Info, Percent_Total, color = Genotype), alpha = 0.7, size = 1.8, position = position_dodge(width = 0.8)) +
  geom_text(aes(Info, label = p.signif),y = -0.0003, size = 3) +
  labs(y = "Percent of Total Peptide Counts", x = "OR-associated peptide sequences and the number of associated OR proteins",
       title = "Approximately 40% of detected Olfr-associated peptides are differentially affected in 10G4-expressing mice") +
  theme_bw() +
  theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.95),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 1))

ggsave("Output2/F04_OR10G4_Olfr_peptide_DEP_final_MF4D.png", width = 45, height = 25, units = "cm")

WTxG4_peptide <- intersect(PeptideSeqs_103_seqs_Nested[[2]][[1]]$MS_peptides,
                           PeptideSeqs_103_seqs_Nested[[2]][[2]]$MS_peptides)

MSuni_Counts %>%
  left_join(CompareMeans_Info, by = "Info") %>%
  filter(MS_peptides %in% WTxG4_peptide) %>%
  arrange(p) %>%
  mutate(Info = factor(Info, levels = unique(Info))) %>%
  ggplot() +
  geom_boxplot(aes(Info, Percent_Total, color = Genotype), alpha = 0.7, position = position_dodge(width = 0.8)) +  
  geom_point(aes(Info, Percent_Total, color = Genotype), alpha = 0.7, size = 1.8, position = position_dodge(width = 0.8)) +
  geom_text(aes(Info, label = p.signif),y = -0.0003, size = 3) +
  labs(y = "Percent of Total Peptide Counts", x = "OR-associated peptide sequences and the number of associated OR proteins",
       title = "Approximately 30% of Olfr-associated peptides detected in both Genotypes are differentially expressed") +
  theme_bw() +
  theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.95),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 1))

ggsave("Output2/F04_OR10G4_Olfr_peptide_OVERLAP_DEP_final_SF10E.png", width = 40, height = 20, units = "cm")


# - End----
sessionInfo()