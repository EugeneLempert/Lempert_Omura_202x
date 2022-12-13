# Required packages ----
if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

if (!require(ggpubr)) install.packages('ggpubr')
library(ggpubr)

if (!require(grid)) install.packages('grid')
library(grid)

if (!require(gridExtra)) install.packages('gridExtra')
library(gridExtra)

# Import all the raw data; most data not used in this analysis as all qPCR analyses use the same import code----
Input_list <-
  list.files(path = "../../Olfr-RT-qPCR/Input",
             pattern = "*.txt", 
             full.names = TRUE) %>% 
  lapply(read_plus)

Input_list_names <-   
  list.files(path = "../../Olfr-RT-qPCR/Input",
             pattern = "*.txt") %>%
  str_replace(pattern = ".txt", 
              replacement = "")

names(Input_list) <- Input_list_names

# Import the primer efficiency file
AA10 <- read_csv("../../Olfr-RT-qPCR/Code/AA10.csv")

# Labeling and cleaning the data ----
# Label the data that compares the OR10G4 lines and generate some simple QC values and plots (7, 7J, and 10) This is my first run
Raw1 <- Label_qPCR_data(CqFile = Input_list$`220502_Exp_Ct`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c("MSS-4164", "MSS-4165", "MSS-4170", "PFH-610", "PFH-613", "PFH-617", "PFH-733", 
                                                       "MSS-4361", "MSS-4362", "MSS-4364", "MSS-4365", "OR10G4L7_gDNA"), 15), 
                                                 rep("UPWater", 15)),
                                 GeneNames = c(rep(c("Acsm4", "Acss2", "Slc25a35", "OR10G4_CDS", "OR10G4_mRNA"), each = 36),
                                               rep(c("Acsm4", "Acss2", "Slc25a35", "OR10G4_CDS", "OR10G4_mRNA"), each = 3)),
                                 Mutation = c(rep(c("Line 7", "Line 7", "Line 7", "Line 7J", "Line 7J", "Line 7J", "Line 7J", 
                                                    "Line 10", "Line 10", "Line 10", "Line 10", "OR10G4L7_gDNA"), 15), 
                                              rep("UPWater", 15)),
                                 Plate = rep("Run1", 195),   
                                 ReplicateNames = c(rep(c("One", "Two", "Three"), times = 5, each = 12),
                                                    rep(c("One", "Two", "Three"), times = 5)))
Raw1_QE <- Quality_Evaluation(Raw1, 
                                       MeltFile = Input_list$`220502_Exp_Melt`,
                                       RepDiffAllow = 0.6, CqLower = 15)
Raw1_Plot <- Raw1_QE$QE_table %>% Plot_QE(MutationColumn = "Mutation")

grid.arrange(Raw1_Plot$SdMean, Raw1_Plot$TmHeight, Raw1_Plot$RawRQ, Raw1_Plot$Combo)

# This is the second attempt at comparing 10G4 line mRNA expression. Mostly a repeat, but not completely due the GCaMP primer and L10_gDNA sample
Raw2 <- Label_qPCR_data(CqFile = Input_list$`220818_Exp_Ct`, 
                                          Style = "Custom", 
                                          SampleNames = c(rep(c("MSS-4164", "MSS-4165", "MSS-4170", "PFH-610", "PFH-613", "PFH-617", "PFH-733", 
                                                                "MSS-4361", "MSS-4362", "MSS-4364", "MSS-4365", "OR10G4L7_gDNA"), 18), 
                                                          rep(c("OR10G4L10_gDNA", "UPWater"), each = 18)),
                                          GeneNames = c(rep(c("Acsm4", "Acss2", "Slc25a35", "OR10G4_CDS", "OR10G4_mRNA", "GCaMP"), each = 36),
                                                        rep(c("Acsm4", "Acss2", "Slc25a35", "OR10G4_CDS", "OR10G4_mRNA", "GCaMP"), each = 3, times = 2)),
                                          Mutation = c(rep(c("Line 7", "Line 7", "Line 7", "Line 7J", "Line 7J", "Line 7J", "Line 7J", 
                                                             "Line 10", "Line 10", "Line 10", "Line 10", "OR10G4L7_gDNA"), 18), 
                                                       rep(c("OR10G4L10_gDNA", "UPWater"), 18)),
                                          Plate = rep("Run2", 252),
                                          ReplicateNames = c(rep(c("One", "Two", "Three"), times = 6, each = 12),
                                                             rep(c("One", "Two", "Three"), times = 12)))
Raw2_QE <- Quality_Evaluation(Raw2, 
                                                MeltFile = Input_list$`220818_Exp_Melt`,
                                                RepDiffAllow = 0.6, CqLower = 15)
Raw2_Plot <- Raw2_QE$QE_table %>% Plot_QE(MutationColumn = "Mutation")

grid.arrange(Raw2_Plot$SdMean, Raw2_Plot$TmHeight, Raw2_Plot$RawRQ, Raw2_Plot$Combo)

# After reviewing the raw data, removing unwanted columns and samples
Clean1 <- Raw1_QE$QE_table %>% 
  filter(Sample %notin% c("OR10G4L7_gDNA", "UPWater")) %>%
  filter(Peak == 1)

Clean2 <- Raw2_QE$QE_table %>% 
  filter(Sample %notin% c("OR10G4L7_gDNA","OR10G4L10_gDNA", "UPWater")) %>%
  filter(Peak == 1)

# Merging the two data sets for analysis
Clean_Input <- Clean1 %>%
  bind_rows(Clean2)

# Analysis start, converting raw Cq into relative values, though I switch to using Ct instead of Cq----
Output1 <- Clean_Input %>%
  ungroup() %>%
  group_by(Plate, Gene, Sample, Mutation) %>%
  summarize(MeanSampleCt = mean(Cq, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(AA10, by = "Gene") %>%
  select(Plate, Gene, Sample, Mutation, MeanSampleCt, Primer_Efficiency) %>%
  rename(Genotype = Mutation) %>%
  mutate(MeanSCtxlog2Eff = MeanSampleCt*log2(Primer_Efficiency))

Calibrators <- Output1 %>%
  filter(Gene %in% c("Acsm4", "Acss2", "Slc25a35")) %>%
  select(Sample, Gene, Plate, MeanSCtxlog2Eff) %>%
  pivot_wider(names_from = Gene, values_from = MeanSCtxlog2Eff) %>%
  rename(MeanSCtxlog2Eff_Acsm4 = Acsm4, 
         MeanSCtxlog2Eff_Acss2 = Acss2,
         MeanSCtxlog2Eff_Slc25a35 = Slc25a35)

Output2 <- Output1 %>% 
  left_join(Calibrators, by = c("Sample", "Plate")) %>%
  mutate(DCtSample = (MeanSCtxlog2Eff_Acsm4 + MeanSCtxlog2Eff_Acss2 + MeanSCtxlog2Eff_Slc25a35)/3 - MeanSCtxlog2Eff)

#Taking the mean of the DCtSample values from each plate to unify across plates
Output_Sample <- Output2 %>%
  ungroup() %>% 
  group_by(Sample, Genotype, Gene) %>% 
  summarize(MeanDCtSample = mean(DCtSample))

#Calculating a DCt for each genotype/gene pair
Output_Genotype <- Output_Sample %>%
  ungroup() %>% 
  group_by(Genotype, Gene) %>% 
  summarize(MeanDCtGenotype = mean(MeanDCtSample),
            DCtGenotypeSize = length(MeanDCtSample),
            DCtGenotypeSEM = sd(MeanDCtSample)/sqrt(DCtGenotypeSize))

#Calculating the DDCt for each genotype/gene pair using Line 10 as the reference point
Output_GenotypeL10 <- Output_Genotype %>%
  filter(Genotype == "Line 10") %>%
  rename(DCtL10Size = DCtGenotypeSize, 
         MeanDCtL10 = MeanDCtGenotype, 
         DCtL10SEM = DCtGenotypeSEM) %>%
  ungroup() %>%
  select(-Genotype)

Output_Genotype_DDCt <- Output_Genotype %>%
  left_join(Output_GenotypeL10, by = "Gene") %>%
  mutate(MeanDDCt = MeanDCtGenotype - MeanDCtL10, 
         DDCtSEM = sqrt(DCtGenotypeSEM^2 + DCtL10SEM^2),
         DDCtCI = qt(0.975, df = DCtGenotypeSize + DCtL10Size - 2) * DDCtSEM)

#DDCt for each sample
Output_Sample_DDCt <- Output_Sample %>%
  ungroup() %>%
  left_join(Output_GenotypeL10, by = "Gene") %>%
  select(Sample, Gene, Genotype, MeanDCtSample, MeanDCtL10) %>%
  mutate(SampleDDCt = MeanDCtSample - MeanDCtL10)

DDCtSampleGenotype_Filtered <- Output_Sample_DDCt %>%
  left_join(Output_Genotype_DDCt, by = c("Gene", "Genotype")) %>%
  filter(Gene %notin% c("Acsm4", "Acss2", "Slc25a35")) %>%
  mutate(Gene = case_when(Gene == "OR10G4_CDS" ~ "Coding exon: OR region",
                          Gene == "GCaMP" ~ "Coding exon: GCaMP region",
                          TRUE ~ "5'UTR/Coding exon splice junction")) %>%
  mutate(Genotype = factor(Genotype, levels = c("Line 10", "Line 7", "Line 7J")))

# Graphs and statistics----
DDCtSampleGenotype_Filtered %>%
  filter(Gene != "5'UTR/Coding exon splice junction") %>%
  ggplot(aes(Gene, SampleDDCt, color = Genotype)) +
  geom_point(position = position_dodge(width = 0.6), size = 4) + 
  geom_errorbar(aes(ymin = MeanDDCt - DDCtCI, ymax = MeanDDCt + DDCtCI), width = 0.4, position = position_dodge(0.6)) +
  labs(title = "Log2 Fold Change in OR10G4 mRNA for three 10G4 lines relative to Line 10", 
       subtitle = "Line 7/7J produce twice as much total 10G4 mRNA as Line 10.", 
       caption = "DDCt calculated using corrected primer efficiencies and reference primers Acsm4, Acss2, and Slc25a35. Each Sample and Genotype DDCt calibrated using the mean of Line 10 spliced mRNA for each target. \nEach sample/target pair loaded in triplicate. At least 3 biological replicates per genotype. Plotted L2FC of each Sample with the 95% confidence interval based on the mean and single-propagated error of each Genotype", 
       y = expression("Log2 Fold Change (-" ~ Delta*Delta ~ "Ct)")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = c(0.87, 0.21),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank(),
        plot.caption = element_text(size = 7, face = "italic"), 
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5))

#Quick stats
compare_means(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered, method = "anova", group.by = "Gene")
compare_means(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered, method = "t.test", group.by = "Gene")

#Intra-Assay CV: 0.442
Clean_Input %>% 
  ungroup() %>% 
  select(Plate, Sample, Gene, Cq) %>%
  distinct() %>%
  group_by(Plate, Sample, Gene) %>%
  summarize(Plate = Plate, 
            Sample = Sample, 
            Gene = Gene, 
            RepMean = mean(Cq, na.rm = TRUE),
            RepSd = sd(Cq, na.rm = TRUE)) %>% 
  ungroup() %>% 
  distinct() %>%
  select(Plate, Sample, Gene, RepMean, RepSd) %>% 
  distinct() %>%
  mutate(CV = RepSd/RepMean * 100) %>%
  group_by(Plate) %>%
  summarize(MeanCV = mean(CV, na.rm = TRUE)) %>%
  summarize(mean = mean(MeanCV))


#Inter-Assay CV: 8.77
Clean_Input %>% 
  ungroup() %>% 
  select(Plate, Sample, Gene, Cq) %>%
  distinct() %>%
  group_by(Plate, Sample, Gene) %>%
  summarize(Plate = Plate, 
            Sample = Sample, 
            Gene = Gene, 
            RepMean = mean(Cq, na.rm = TRUE),
            RepSd = sd(Cq, na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(Plate, Sample, Gene, RepMean, RepSd) %>% 
  distinct() %>%
  group_by(Sample, Gene) %>%
  summarize(PlateMeans = mean(RepMean, na.rm = TRUE), 
            PlateSd = sd(RepMean, na.rm = TRUE)) %>% 
  ungroup() %>% 
  summarise(CV = PlateSd/PlateMeans * 100) %>%
  summarize(meanCV = mean(CV, na.rm = TRUE))
