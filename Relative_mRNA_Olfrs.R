# Required packages ----
if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

if (!require(ggpubr)) install.packages('ggpubr')
library(ggpubr)

if (!require(grid)) install.packages('grid')
library(grid)

if (!require(gridExtra)) install.packages('gridExtra')
library(gridExtra)

if (!require(car)) install.packages('car')
library(car)

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

# Labeling and set-specific cleaning the data that compares the OR10G4 lines and generate some simple QC values and plots (7, 7J, and 10) ----
    #The QE Plots are turned into comments, so run the two lines of code fora dataset to see the QE plots
Plate1 <- Label_qPCR_data(CqFile = Input_list$`220221_Exp_Cq`, 
                                 Style = "Simple", 
                                 SampleNames = c("PFH-609", "PFH-610", "PFH-611", "PFH-613", "PFH-614", "PFH-617", "PFH-618", "PFH-732", "PFH-733", "PFH-734",
                                                 "MSS-4163", "MSS-4164", "MSS-4165", "MSS-4166", "MSS-4169", "MSS-4170", "MSS-4171", "WT_gDNA", "UPwater"), 
                                 GeneNames = c("Olfr358", "Olfr609", "Olfr690", "Olfr1154", "Acsm4", "Olfr596_603", "Olfr390", "Olfr510"),
                                 Mutation = c("WT", "L7J", "WT", "L7J", "WT", "L7J", "WT", "WT", "L7J", "WT", 
                                              "WT", "L7", "L7", "WT", "WT", "L7", "WT", "WT_gDNA", "Blank"),
                                 ReplicateNames = c("R1", "R2"),
                                 Plate = "P1")
  Plate1_QE <- Quality_Evaluation(Plate1, 
                                MeltFile = Input_list$`220221_Exp_Melt`,
                                RepDiffAllow = 0.6)
#  Plate1_Plot <- Plate1_QE$QE_table %>% Plot_QE(MutationColumn = "Mutation")
#  grid.arrange(Plate1_Plot$SdMean, Plate1_Plot$TmHeight, Plate1_Plot$RawRQ, Plate1_Plot$Combo)
  Plate1_Clean <- Plate1_QE$QE_table %>%
    filter(Pos != "N17") %>%
    filter(Gene != "Olfr609")

Plate2 <- Label_qPCR_data(CqFile = Input_list$`220228_Exp_Cq_1`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c(rep(c("MSS-4359", "MSS-4360", "MSS-4361", "MSS-4362", "MSS-4363", "MSS-4364", "MSS-4365", "MSS-4368"), 3),
                                                       rep(c("PFH-609", "PFH-610", "PFH-611", "PFH-613", "PFH-614", "PFH-617", "PFH-618", "PFH-733"), 3)), 6), 
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35"), each = 48),
                                               rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35"), each = 6)),
                                 Mutation =    c(rep(c(rep(c("WT", "WT", "L10", "L10", "WT", "L10", "L10", "WT"), 3),
                                                       rep(c("WT", "L7J", "WT", "L7J", "WT", "L7J", "WT", "L7J"), 3)), 6),
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 Plate = rep("P2", times = 324),
                                 ReplicateNames = c(rep(c("R1", "R2", "R3"), times = 12, each = 8),
                                                    rep(c("R1", "R2", "R3"), times = 12)))
  Plate2_QE <- Quality_Evaluation(Plate2, 
                                      MeltFile = Input_list$`220228_Exp_Melt_1`,
                                      RepDiffAllow = 0.6)
#  Plate2_Plot <- Plate2_QE$QE_table %>% Plot_QE(MutationColumn = "Mutation")
#  grid.arrange(Plate2_Plot$SdMean, Plate2_Plot$TmHeight, Plate2_Plot$RawRQ, Plate2_Plot$Combo)

Plate3 <- Label_qPCR_data(CqFile = Input_list$`220228_Exp_Cq_2`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c(rep(c("MSS-4359", "MSS-4360", "MSS-4361", "MSS-4362", "MSS-4363", "MSS-4364", "MSS-4365", "MSS-4368"), 3),
                                                       rep(c("PFH-609", "PFH-610", "PFH-611", "PFH-613", "PFH-614", "PFH-617", "PFH-618", "PFH-733"), 3)), 6), 
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35"), each = 48),
                                               rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35"), each = 6)),
                                 Mutation =    c(rep(c(rep(c("WT", "WT", "L10", "L10", "WT", "L10", "L10", "WT"), 3),
                                                       rep(c("WT", "L7J", "WT", "L7J", "WT", "L7J", "WT", "L7J"), 3)), 6),
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 Plate = rep("P3", times = 324),
                                 ReplicateNames = c(rep(c("R1", "R2", "R3"), times = 12, each = 8),
                                                    rep(c("R1", "R2", "R3"), times = 12)))
  Plate3_QE <- Quality_Evaluation(Plate3, 
                                      MeltFile = Input_list$`220228_Exp_Melt_2`,
                                      RepDiffAllow = 0.6)
#  Plate3_Plot <- Plate3_QE$QE_table %>% Plot_QE(MutationColumn = "Mutation")
#  grid.arrange(Plate3_Plot$SdMean, Plate3_Plot$TmHeight, Plate3_Plot$RawRQ, Plate3_Plot$Combo)
  
Plate4 <- Label_qPCR_data(CqFile = Input_list$`220228_Exp_Cq_3`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c(rep(c("MSS-4359", "MSS-4360", "MSS-4361", "MSS-4362", "MSS-4363", "MSS-4364", "MSS-4365", "MSS-4368"), 3),
                                                       rep(c("PFH-609", "PFH-610", "PFH-611", "PFH-613", "PFH-614", "PFH-617", "PFH-618", "PFH-733"), 3)), 6), 
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35"), each = 48),
                                               rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35"), each = 6)),
                                 Mutation =    c(rep(c(rep(c("WT", "WT", "L10", "L10", "WT", "L10", "L10", "WT"), 3),
                                                       rep(c("WT", "L7J", "WT", "L7J", "WT", "L7J", "WT", "L7J"), 3)), 6),
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 Plate = rep("P4", times = 324),
                                 ReplicateNames = c(rep(c("R1", "R2", "R3"), times = 12, each = 8),
                                                    rep(c("R1", "R2", "R3"), times = 12)))
  Plate4_QE <- Quality_Evaluation(Plate4, 
                                      MeltFile = Input_list$`220228_Exp_Melt_3`,
                                      RepDiffAllow = 0.6)
#  Plate4_Plot <- Plate4_QE$QE_table %>% Plot_QE(MutationColumn = "Mutation")
#  grid.arrange(Plate4_Plot$SdMean, Plate4_Plot$TmHeight, Plate4_Plot$RawRQ, Plate4_Plot$Combo)
  Plate4_Clean <- Plate4_QE$QE_table %>% 
    filter(Pos %notin% c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "E17", "E18", "E19", "E20", "E21", "E22", "E23", "E24", "L24")) 

Plate5 <- Label_qPCR_data(CqFile = Input_list$`220301_Exp_Cq_1`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c(rep(c("MSS-4359", "MSS-4360", "MSS-4361", "MSS-4362", "MSS-4363", "MSS-4364", "MSS-4365", "MSS-4368"), 3),
                                                       rep(c("PFH-609", "PFH-610", "PFH-611", "PFH-613", "PFH-614", "PFH-617", "PFH-618", "PFH-733"), 3)), 6), 
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 GeneNames = c(rep(c("Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 48),
                                               rep(c("Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 6)),
                                 Mutation =    c(rep(c(rep(c("WT", "WT", "L10", "L10", "WT", "L10", "L10", "WT"), 3),
                                                       rep(c("WT", "L7J", "WT", "L7J", "WT", "L7J", "WT", "L7J"), 3)), 6),
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 Plate = rep("P5", times = 324),
                                 ReplicateNames = c(rep(c("R1", "R2", "R3"), times = 12, each = 8),
                                                    rep(c("R1", "R2", "R3"), times = 12)))
  Plate5_QE <- Quality_Evaluation(Plate5, 
                                      MeltFile = Input_list$`220301_Exp_Melt_1`,
                                      RepDiffAllow = 0.6)
#  Plate5_Plot <- Plate5_QE$QE_table %>% Plot_QE(MutationColumn = "Mutation")
#  grid.arrange(Plate5_Plot$SdMean, Plate5_Plot$TmHeight, Plate5_Plot$RawRQ, Plate5_Plot$Combo)

Plate6 <- Label_qPCR_data(CqFile = Input_list$`220301_Exp_Cq_2`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c(rep(c("MSS-4359", "MSS-4360", "MSS-4361", "MSS-4362", "MSS-4363", "MSS-4364", "MSS-4365", "MSS-4368"), 3),
                                                       rep(c("PFH-609", "PFH-610", "PFH-611", "PFH-613", "PFH-614", "PFH-617", "PFH-618", "PFH-733"), 3)), 6), 
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 GeneNames = c(rep(c("Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 48),
                                               rep(c("Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 6)),
                                 Mutation =    c(rep(c(rep(c("WT", "WT", "L10", "L10", "WT", "L10", "L10", "WT"), 3),
                                                       rep(c("WT", "L7J", "WT", "L7J", "WT", "L7J", "WT", "L7J"), 3)), 6),
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 Plate = rep("P6", times = 324),
                                 ReplicateNames = c(rep(c("R1", "R2", "R3"), times = 12, each = 8),
                                                    rep(c("R1", "R2", "R3"), times = 12)))
  Plate6_QE <- Quality_Evaluation(Plate6, 
                                      MeltFile = Input_list$`220301_Exp_Melt_2`,
                                      RepDiffAllow = 0.6)
#  Plate6_Plot <- Plate6_QE$QE_table %>% Plot_QE(MutationColumn = "Mutation")
#  grid.arrange(Plate6_Plot$SdMean, Plate6_Plot$TmHeight, Plate6_Plot$RawRQ, Plate6_Plot$Combo)

Plate7 <- Label_qPCR_data(CqFile = Input_list$`220301_Exp_Cq_3`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c(rep(c("MSS-4359", "MSS-4360", "MSS-4361", "MSS-4362", "MSS-4363", "MSS-4364", "MSS-4365", "MSS-4368"), 3),
                                                       rep(c("PFH-609", "PFH-610", "PFH-611", "PFH-613", "PFH-614", "PFH-617", "PFH-618", "PFH-733"), 3)), 6), 
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 GeneNames = c(rep(c("Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 48),
                                               rep(c("Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 6)),
                                 Mutation =    c(rep(c(rep(c("WT", "WT", "L10", "L10", "WT", "L10", "L10", "WT"), 3),
                                                       rep(c("WT", "L7J", "WT", "L7J", "WT", "L7J", "WT", "L7J"), 3)), 6),
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 Plate = rep("P7", times = 324),
                                 ReplicateNames = c(rep(c("R1", "R2", "R3"), times = 12, each = 8),
                                                    rep(c("R1", "R2", "R3"), times = 12)))
  Plate7_QE <- Quality_Evaluation(Plate7, 
                                      MeltFile = Input_list$`220301_Exp_Melt_3`,
                                      RepDiffAllow = 0.6)
#  Plate7_Plot <- Plate7_QE$QE_table %>% Plot_QE(MutationColumn = "Mutation")
#  grid.arrange(Plate7_Plot$SdMean, Plate7_Plot$TmHeight, Plate7_Plot$RawRQ, Plate7_Plot$Combo)

Plate8 <- Label_qPCR_data(CqFile = Input_list$`220412_Exp_Ct`, 
                         Style = "Custom", 
                         SampleNames = c(rep(c("MSS-4163", "MSS-4164", "MSS-4165", "MSS-4166", "MSS-4169", "MSS-4170", "MSS-4171", "WT_gDNA"), 27), 
                                         rep("UPWater", 27)),
                         GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 24),
                                       rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 3)),
                         Mutation = c(rep(c("WT", "L7", "L7", "WT", "WT", "L7", "WT", "WT_gDNA"), 27), 
                                      rep("UPWater", 27)),
                         Plate = rep("P8", 243),   
                         ReplicateNames = c(rep(c("R1", "R2", "R3"), times = 9, each = 8),
                                            rep(c("R1", "R2", "R3"), times = 9)))
  Plate8_QE <- Quality_Evaluation(Plate8, 
                                MeltFile = Input_list$`220412_Exp_Melt`,
                                RepDiffAllow = 0.6, CqLower = 17)
#  Plate8_Plot <- Plate8_QE$QE_table %>% Plot_QE(MutationColumn = "Mutation")
#  grid.arrange(Plate8_Plot$SdMean, Plate8_Plot$TmHeight, Plate8_Plot$RawRQ, Plate8_Plot$Combo)

Plate9 <- Label_qPCR_data(CqFile = Input_list$`220413_Exp_Ct`, 
                         Style = "Custom", 
                         SampleNames = c(rep(c("MSS-4163", "MSS-4164", "MSS-4165", "MSS-4166", "MSS-4169", "MSS-4170", "MSS-4171", "WT_gDNA"), 27), 
                                         rep("UPWater", 27)),
                         GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 24),
                                       rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 3)),
                         Mutation = c(rep(c("WT", "L7", "L7", "WT", "WT", "L7", "WT", "WT_gDNA"), 27), 
                                      rep("UPWater", 27)),
                         Plate = rep("P9", 243),   
                         ReplicateNames = c(rep(c("R1", "R2", "R3"), times = 9, each = 8),
                                            rep(c("R1", "R2", "R3"), times = 9)))
  Plate9_QE <- Quality_Evaluation(Plate9, 
                                MeltFile = Input_list$`220413_Exp_Melt`,
                                RepDiffAllow = 0.6, CqLower = 17)
#  Plate9_Plot <- Plate9_QE$QE_table %>% Plot_QE(MutationColumn = "Mutation")
#  grid.arrange(Plate9_Plot$SdMean, Plate9_Plot$TmHeight, Plate9_Plot$RawRQ, Plate9_Plot$Combo)

Plate10 <- Label_qPCR_data(CqFile = Input_list$`220414_Exp_Ct`, 
                         Style = "Custom", 
                         SampleNames = c(rep(c("MSS-4163", "MSS-4164", "MSS-4165", "MSS-4166", "MSS-4169", "MSS-4170", "MSS-4171", "WT_gDNA"), 27), 
                                         rep("UPWater", 27)),
                         GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 24),
                                       rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 3)),
                         Mutation = c(rep(c("WT", "L7", "L7", "WT", "WT", "L7", "WT", "WT_gDNA"), 27), 
                                      rep("UPWater", 27)),
                         Plate = rep("P10", 243),   
                         ReplicateNames = c(rep(c("R1", "R2", "R3"), times = 9, each = 8),
                                            rep(c("R1", "R2", "R3"), times = 9)))
  Plate10_QE <- Quality_Evaluation(Plate10, 
                                MeltFile = Input_list$`220414_Exp_Melt`,
                                RepDiffAllow = 0.6, CqLower = 17)
#  Plate10_Plot <- Plate10_QE$QE_table %>% Plot_QE(MutationColumn = "Mutation")
#  grid.arrange(Plate10_Plot$SdMean, Plate10_Plot$TmHeight, Plate10_Plot$RawRQ, Plate10_Plot$Combo)

# Cleaning and Combining plates----
OR10G4_Input <- Plate1_Clean %>%
    bind_rows(Plate2_QE$QE_table) %>%
    bind_rows(Plate3_QE$QE_table) %>%
    bind_rows(Plate4_Clean) %>%
    bind_rows(Plate5_QE$QE_table) %>%
    bind_rows(Plate6_QE$QE_table) %>%
    bind_rows(Plate7_QE$QE_table) %>%
    bind_rows(Plate8_QE$QE_table) %>%
    bind_rows(Plate9_QE$QE_table) %>%
    bind_rows(Plate10_QE$QE_table) %>%
    filter(FalseCq == FALSE) %>%
    filter(NoPeak == FALSE) %>%
    filter(MaxCq == FALSE) %>%
    filter(Sample %notin% c("WT_gDNA", "Blank", "UPwater", "UPWater")) %>%
    filter(Peak == 1) %>%
    filter(NoCq == FALSE) %>%
    filter(NotWithinCqRange == FALSE) %>% 
    filter(Cq < 33) %>% 
    filter(PoorReplicate == FALSE)

# Analysis start, converting raw Cq into relative values, though I switch to using Ct instead of Cq----
Output1 <- OR10G4_Input %>%
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
  mutate(Refs = 3 - (is.na(MeanSCtxlog2Eff_Acsm4) + is.na(MeanSCtxlog2Eff_Acss2) + is.na(MeanSCtxlog2Eff_Slc25a35)),
         MeanSCtxlog2Eff_Acsm4 = ifelse(is.na(MeanSCtxlog2Eff_Acsm4), 0, MeanSCtxlog2Eff_Acsm4),
         MeanSCtxlog2Eff_Acss2 = ifelse(is.na(MeanSCtxlog2Eff_Acss2), 0, MeanSCtxlog2Eff_Acss2),
         MeanSCtxlog2Eff_Slc25a35 = ifelse(is.na(MeanSCtxlog2Eff_Slc25a35), 0, MeanSCtxlog2Eff_Slc25a35),
         DCtSample = (MeanSCtxlog2Eff_Acsm4 + MeanSCtxlog2Eff_Acss2 + MeanSCtxlog2Eff_Slc25a35)/Refs - MeanSCtxlog2Eff)

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

#Calculating the DDCt for each genotype/gene pair using WT as the reference point
Output_GenotypeWT <- Output_Genotype %>%
  filter(Genotype == "WT") %>%
  rename(DCtWTSize = DCtGenotypeSize, 
         MeanDCtWT = MeanDCtGenotype, 
         DCtWTSEM = DCtGenotypeSEM) %>%
  ungroup() %>%
  select(-Genotype)

Output_Genotype_DDCt <- Output_Genotype %>%
  left_join(Output_GenotypeWT, by = "Gene") %>%
  mutate(MeanDDCt = MeanDCtGenotype - MeanDCtWT, 
         DDCtSEM = sqrt(DCtGenotypeSEM^2 + DCtWTSEM^2),
         DDCtCI = qt(0.975, df = DCtGenotypeSize + DCtWTSize - 2) * DDCtSEM)

#DDCt for each sample
Output_Sample_DDCt <- Output_Sample %>%
  ungroup() %>%
  left_join(Output_GenotypeWT, by = "Gene") %>%
  select(Sample, Gene, Genotype, MeanDCtSample, MeanDCtWT) %>%
  mutate(SampleDDCt = MeanDCtSample - MeanDCtWT)

DDCtSampleGenotype_Filtered <- Output_Sample_DDCt %>%
  left_join(Output_Genotype_DDCt, by = c("Gene", "Genotype")) %>%
  filter(Gene %notin% c("Acsm4", "Acss2", "Slc25a35")) %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "L10", "L7", "L7J")))


# Graphs and statistics----
DDCtSampleGenotype_Filtered %>%
  ggplot() +
  geom_point(aes(Gene, SampleDDCt, color = Genotype, group = Genotype), position = position_dodge(width = 0.6), size = 3) +
  scale_color_manual(values=c("WT" = "yellow2", "L7" = "skyblue1", "L7J" = "red3", "L10" = "darkorchid1")) +
  geom_errorbar(aes(x = Gene, group = Genotype, ymin = MeanDDCt - DDCtCI, ymax = MeanDDCt + DDCtCI), 
                width = 0.4, position = position_dodge(0.6), color = "black") +
  labs(title = "Log2 Fold Change in Ct values for Olfr mRNA in Three 10G4 lines", 
       y = expression("Log2 Fold Change (-" ~ Delta*Delta ~ "Ct)") ) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(face="bold", size = 15),
        axis.title.y = element_text(size = 17),
        legend.position = c(0.11, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(face="bold"),
        legend.title = element_blank()) +
  geom_vline(xintercept = 2.5, linetype="dotted", color = "grey40", linewidth=1.5) +
  scale_y_continuous(breaks = c(-7, -6, -5, -4, -3, -2, -1, 0, 1)) +
  annotate("text", label = "Class I", x = 1.5, y = 1.35, size = 8, colour = "blue", fontface = 2) +
  annotate("text", label = "Class II", x = 4.5, y = 1.35, size = 8, colour = "red", fontface = 2) +
  annotate("text", label = "DVI: 1.05", x = 1, y = -7.20, size = 4, colour = "black") +
  annotate("text", label = "DVI: 1.05", x = 2, y = -7.20, size = 4, colour = "black") +
  annotate("text", label = "DVI: 1.05", x = 3, y = -7.20, size = 4, colour = "black") +
  annotate("text", label = "DVI: 1.05", x = 4, y = -7.20, size = 4, colour = "black") +
  annotate("text", label = "DVI: 1.5", x = 5, y = -7.20, size = 4, colour = "black") +
  annotate("text", label = "DVI: 3.6", x = 6, y = -7.20, size = 4, colour = "black") 

compare_means(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered, method = "anova", group.by = "Gene")

compare_means(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered, method = "t.test", group.by = "Gene")
#----More stats----
#Let's explore this data set with the idea that I have three Genotypes and a column of DDCt values for each of them
#Let's plot the data points as a histogram for each Genotype
ggplot(DDCtSampleGenotype_Filtered, aes(SampleDDCt, fill = Gene)) + geom_histogram(binwidth = 0.5) +  facet_wrap(~Genotype)
ggplot(DDCtSampleGenotype_Filtered, aes(SampleDDCt, Genotype)) + geom_boxplot()

#Let's run a one-way ANOVA and then check the residuals
AOV1 <- aov(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered)
#A histogram for the residuals
ggplot(DDCtSampleGenotype_Filtered, aes(AOV1$residuals)) + geom_histogram(binwidth = 0.25) + facet_wrap(~Genotype)
#A multiplot of the residuals and other metrics
par(mfrow=c(2,2))
plot(AOV1, lwd = 2)
par(mfrow=c(1,1))
#leveneTest and barlete test for homoscedasticity
leveneTest(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered)

bartlett.test(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered)

fligner.test(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered)

#tests for normality
shapiro.test(AOV1$residuals)

#non-parametric one-way anova 
kruskal.test(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered)

#multiple comparisons using wilcox and t.test
compare_means(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered, method = "t.test")
compare_means(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered, method = "wilcox.test")

compare_means(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered, method = "t.test") %>%
  arrange(p) %>%
  mutate(row = row_number(),
         HolmCorrection = p * (6 - row + 1),
         Pass_DunnSidakCorrection = p < 1 - (1 - 0.05)^(1/6))

#Now let's run Gene-split comparisons
ggplot(DDCtSampleGenotype_Filtered, aes(SampleDDCt)) + geom_histogram(binwidth = 0.25) +  facet_grid(rows = vars(Gene), cols = vars(Genotype))
ggplot(DDCtSampleGenotype_Filtered, aes(SampleDDCt, Genotype)) + geom_boxplot() + facet_wrap(~Gene)

#A quick ANOVA and krusal test for each gene group
aov(SampleDDCt ~ Genotype * Gene, data = DDCtSampleGenotype_Filtered) %>% summary()

compare_means(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered, method = "anova", group.by = c("Gene"))
compare_means(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered, method = "kruskal.test", group.by = c("Gene"))

#Making subgroups and running the toughest homoscedasticity test 
ForStats_358 <- DDCtSampleGenotype_Filtered %>% filter(Gene == "Olfr358")
ForStats_390 <- DDCtSampleGenotype_Filtered %>% filter(Gene == "Olfr390")
ForStats_510 <- DDCtSampleGenotype_Filtered %>% filter(Gene == "Olfr510")
ForStats_596_603 <- DDCtSampleGenotype_Filtered %>% filter(Gene == "Olfr596_603")
ForStats_690 <- DDCtSampleGenotype_Filtered %>% filter(Gene == "Olfr690")
ForStats_1154 <- DDCtSampleGenotype_Filtered %>% filter(Gene == "Olfr1154")

fligner.test(SampleDDCt ~ Genotype, data = ForStats_358)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_390)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_510)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_596_603)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_690)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_1154)

bartlett.test(SampleDDCt ~ Genotype, data = ForStats_358)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_390)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_510)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_596_603)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_690)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_1154)

leveneTest(SampleDDCt ~ Genotype, data = ForStats_358)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_390)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_510)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_596_603)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_690)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_1154)

#ANOVA for residuals for normality
AOV358 <- aov(SampleDDCt ~ Genotype, data = ForStats_358)
shapiro.test(AOV358$residuals)
AOV390 <- aov(SampleDDCt ~ Genotype, data = ForStats_390)
shapiro.test(AOV390$residuals)
AOV510 <- aov(SampleDDCt ~ Genotype, data = ForStats_510)
shapiro.test(AOV510$residuals)
AOV596_603 <- aov(SampleDDCt ~ Genotype, data = ForStats_596_603)
shapiro.test(AOV596_603$residuals)
AOV690 <- aov(SampleDDCt ~ Genotype, data = ForStats_690)
shapiro.test(AOV690$residuals)
AOV1154 <- aov(SampleDDCt ~ Genotype, data = ForStats_1154)
shapiro.test(AOV1154$residuals)

#Performing parametric and non-parametric pairwise comparisons
compare_means(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered, method = "wilcox.test", group.by = c("Gene"))

compare_means(SampleDDCt ~ Genotype, data = DDCtSampleGenotype_Filtered, method = "t.test", group.by = c("Gene")) %>%
  select(Gene, group1, group2, p) %>%
  arrange(p) %>%
  mutate(row = row_number(),
         HolmAdjustment = 0.05/(36 - row + 1), 
         PassHA = p < HolmAdjustment, 
         HolmCorrection = p * (36 - row + 1),
         PassHC = HolmCorrection < 0.05, 
         Pass_DunnSidakCorrection = p < 1 - (1 - 0.05)^(1/36)) %>%
  mutate(NEWp = p.adjust(p, method = "holm"))

#ExTra plots
par(mfrow=c(2,2))
plot(AOV358, lwd = 2)
plot(AOV390, lwd = 2)
plot(AOV510, lwd = 2)
plot(AOV596_603, lwd = 2)
plot(AOV690, lwd = 2)
plot(AOV1154, lwd = 2)
par(mfrow=c(1,1))

#----CV values for this entire dataset----
RawForCV <- Plate1_QE$QE_table %>%
  bind_rows(Plate2_QE$QE_table) %>%
  bind_rows(Plate3_QE$QE_table) %>%
  bind_rows(Plate4_QE$QE_table) %>%
  bind_rows(Plate5_QE$QE_table) %>%
  bind_rows(Plate6_QE$QE_table) %>%
  bind_rows(Plate7_QE$QE_table) %>%
  bind_rows(Plate8_QE$QE_table) %>%
  bind_rows(Plate9_QE$QE_table) %>%
  bind_rows(Plate10_QE$QE_table)
  
#Raw value Intra-Assay CV: 0.785
RawForCV %>% 
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
  summarize(mean = mean(MeanCV, na.rm = TRUE))

#Filtered value Intra-Assay CV: 0.497
OR10G4_Input %>% 
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
  summarize(mean = mean(MeanCV, na.rm = TRUE))

#Raw value Inter-Assay CV: 2.10
RawForCV %>% 
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

#Filtered value Inter-Assay CV: 1.32
OR10G4_Input %>% 
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