# This is a script to generate Bioassay Figures and Tables for the 10G4 Paper.

# Set your working directory to the "Bioassay 10G4" folder

# 01 - Negated %in% filter function----
`%notin%` <- Negate(`%in%`)

# 02 - Prereq for vectorized read_csv function----
read_plus <- function(flnm) {
  read_csv(flnm, col_names = FALSE)
}

# 03 - Required packages----
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

# 04 - MyAssays.com 4PL standard curve interpretation----
myAssayOutput <- read_csv("MyAssay_outputs_simple.csv")
# 05 - Liquid One-dose 50uM Bioassays----
Label_list_L1 <- list(tibble(Exp = "Exp42b",
                         Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 42),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 42),
                         Ligand = rep(c("DMSO", "FSK", "Guaiacol", "Vanillin", "Acetovanillone", "Vanillyl butyl ether", "Ethyl vanillin", "Dehydrodivanillin",
                                        "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde", "o-Cresol", "2-Ethylphenol"), each = 24)),
                  tibble(Exp = "Exp46a",
                         Sample = rep(c("PFH1721", "PFH1741", "PFH1742", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 42),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 42),
                         Ligand = rep(c("DMSO", "FSK", "Guaiacol", "Vanillin", "Acetovanillone", "Vanillyl butyl ether", "Ethyl vanillin", "Dehydrodivanillin",
                                        "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde", "o-Cresol", "2-Ethylphenol"), each = 24)),
                  tibble(Exp = "Exp51Newa",
                         Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 36),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 36),
                         Ligand = rep(c("DMSO", "FSK", "4-Fluoroguaiacol", "5-Fluoroguaiacol", "4-Bromoguaiacol", "5-Bromoguaiacol", "4-Chloroguaiacol", "4-Iodoguaiacol",
                                        "4-Methoxyguaiacol", "4-Methylguaiacol", "2-Ethoxyphenol", "Guaiacol"), each = 24)),
                  tibble(Exp = "Exp54NewNew",
                         Sample = rep(c("PFH0887", "PFH0888", "PFH1790","PFH1885", "PFH1789", "PFH1792", "PFH1793","PFH1884"), 36),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 36),
                         Ligand = rep(c("DMSO", "FSK", "4-Fluoroguaiacol", "5-Fluoroguaiacol", "4-Bromoguaiacol", "5-Bromoguaiacol", "4-Chloroguaiacol", "4-Iodoguaiacol",
                                        "4-Methoxyguaiacol", "5-Methoxyguaiacol", "2-Ethoxyphenol", "4-Methylguaiacol"), each = 24)),
                  tibble(Exp = "Exp54OldNew",
                         Sample = rep(c("PFH0887", "PFH0888", "PFH1790","PFH1885", "PFH1789", "PFH1792", "PFH1793","PFH1884"), 42),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 42),
                         Ligand = rep(c("DMSO", "FSK", "Guaiacol", "Vanillin", "Acetovanillone", "Vanillyl butyl ether", "Ethyl vanillin", "Dehydrodivanillin",
                                        "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde", "o-Cresol", "2-Ethylphenol"), each = 24)),
                  tibble(Exp = "Exp54NewOld",
                         Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 36),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 36),
                         Ligand = rep(c("DMSO", "FSK", "4-Fluoroguaiacol", "5-Fluoroguaiacol", "4-Bromoguaiacol", "5-Bromoguaiacol", "4-Chloroguaiacol", "4-Iodoguaiacol",
                                        "4-Methoxyguaiacol", "5-Methoxyguaiacol", "2-Ethoxyphenol", "4-Methylguaiacol"), each = 24)),
                  tibble(Exp = "Exp54OldOld",
                         Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 42),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 42),
                         Ligand = rep(c("DMSO", "FSK", "Guaiacol", "Vanillin", "Acetovanillone", "Vanillyl butyl ether", "Ethyl vanillin", "Dehydrodivanillin",
                                        "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde", "o-Cresol", "2-Ethylphenol"), each = 24)))

Input_list_L1 <-
  list.files(path = "Input/Liquid One Dose", 
             pattern = "*.csv", 
             full.names = TRUE) %>% 
  lapply(read_plus)

Input_list_L1_names <-   
  list.files(path = "Input/Liquid One Dose", 
             pattern = "*.csv") %>%
  str_replace(pattern = ".csv", 
              replacement = "")

names(Input_list_L1) <- Input_list_L1_names

# A rather convoluted piece of code that produces a first table of data for each experiment via an unlabeled intermediate
Input_tibble_L1 <- tibble(Assay = Input_list_L1, Exp = Input_list_L1_names) %>%
  left_join(myAssayOutput, by = "Exp") %>%
  mutate(Background = (Blank1 + Blank2)/2, 
         Measure = map2(Assay, Background, ~ c(t(as.matrix(.x))) - .y),
         Concentration = pmap(list(a, b, c, d, Measure), .f = function(a, b, c, d, Measure) {c * ((a - d)/(Measure - d) - 1)^(1/b)} ),
         Position = map(Assay, ~ 1:length(unlist(.x))), 
         Unbound = map(Concentration, ~ ifelse(is.na(.x), TRUE, FALSE)), 
         Concentration2 = pmap(list(Concentration, Measure, a, d), .f = function(Concentration, Measure, a, d) 
         {case_when(is.na(Concentration) & Measure > d ~ min(Concentration, na.rm= TRUE),
                    is.na(Concentration) & Measure < a ~ max(Concentration, na.rm = TRUE),
                    TRUE ~ Concentration)}),
         Combined = pmap(list(Measure, Position, Concentration, Unbound, Concentration2), 
                         .f = function (Measure, Position, Concentration, Unbound, Concentration2) {
           tibble(Measure = Measure, Position = Position, Concentration = Concentration, 
                  Unbound = Unbound, Concentration2 = Concentration2)}))

Exp_L1 <- tibble(Labels = Label_list_L1) %>%
  mutate(Exp = map_chr(Labels, ~unique(.x$Exp))) %>%
  inner_join(Input_tibble_L1, by = "Exp") %>% 
  mutate(Full_Table = map2(Combined, Labels, bind_cols)) %>%
  select(Exp, Full_Table, a, b, c, d) %>%
  mutate(QC = map(Full_Table, function(.x) {.x %>% group_by(Sample, Ligand) %>%
      nest() %>%
      mutate(MeanC = map(data, function(x) mean(x$Concentration2)),
             SD = map(data, function(x) sd(x$Concentration2))) %>%
      unnest(c(data, MeanC, SD))  %>%
      ungroup() %>%
      mutate(CV = SD/MeanC * 100,
             SEM = SD/sqrt(3))}))

Exp_L1$QC[["Exp54OldNew"]] <- Exp_L1$QC[["Exp54OldNew"]] %>%
  filter(Position != 41)
Exp_L1$QC[["Exp54OldOld"]] <- Exp_L1$QC[["Exp54OldOld"]] %>%
  filter(Position != 22)

# Stupid large file, possibly smarter to split it with multiple separate map calls? Works fine. 
Exp_L1_Processed <- Exp_L1 %>%
  mutate(Concentration_Triplicate_mean = pmap(list(QC, a, b, c, d, Exp), function(QC, a, b, c, d, Exp) {QC %>% group_by(Sample, Mutation, Ligand) %>%
      summarize(S_Measure = mean(Measure)) %>%
      mutate(Concentration_S_Measure = c * ((a - d)/(S_Measure - d) - 1)^(1/b)) %>%
      ungroup() %>%
      bind_cols(tibble(Experiment = Exp))}),
      All_Analysis_Values = map(Concentration_Triplicate_mean, function(.x){
        .x <- .x %>% rename(Concentration = Concentration_S_Measure)
        
        ExpDMSO <- .x %>%
          filter(Ligand == "DMSO") %>%
          select(Sample, Concentration) %>%
          rename(DMSOconcentration = Concentration)
        
        WithDMSO <- .x %>%
          left_join(ExpDMSO, by = "Sample") %>%
          mutate(CCviaDMSO = Concentration - DMSOconcentration) 
        
        ExpFSK <- WithDMSO %>%
          filter(Ligand == "FSK") %>%
          select(Sample, Concentration, CCviaDMSO) %>%
          rename(FSKconcentration = Concentration,
                 FSKCCviaDMSO = CCviaDMSO)
        
        WithDMSOandFSK <- WithDMSO %>%
          left_join(ExpFSK, by = "Sample") %>%
          mutate(RAviaFSK = Concentration/FSKconcentration,
                 RAviaDMSO = Concentration/DMSOconcentration,
                 CCRAviaFSKCC = CCviaDMSO/FSKCCviaDMSO, 
                 CCRAviaDMSO = CCviaDMSO/DMSOconcentration)
        
        WTmean <- WithDMSOandFSK %>%
          filter(Mutation == "WT") %>%
          group_by(Ligand) %>%
          summarize(mCCRAviaDMSO = mean(CCRAviaDMSO))
        
        WithDMSOandFSK2 <- WithDMSOandFSK %>%
          left_join(WTmean, by = "Ligand") %>%
          mutate(CCRAviaDMSO2 = CCRAviaDMSO - mCCRAviaDMSO)
        
        Output <- WithDMSOandFSK2 %>%
          select(Sample, Mutation, Ligand, Concentration, DMSOconcentration, FSKconcentration, 
                 RAviaDMSO, CCviaDMSO, RAviaFSK, CCRAviaFSKCC, CCRAviaDMSO, CCRAviaDMSO2) %>%
          group_by(Mutation, Ligand) %>%
          nest() %>%
          mutate(MeanRAviaDMSO = map(data, function(x) mean(x$RAviaDMSO)),
                 SDRAviaDMSO = map(data, function(x) sd(x$RAviaDMSO)),
                 MeanRAviaFSK = map(data, function(x) mean(x$RAviaFSK)),
                 SDRAviaFSK = map(data, function(x) sd(x$RAviaFSK)),
                 MeanCCRAviaFSKCC = map(data, function(x) mean(x$CCRAviaFSKCC)),
                 SDCCRAviaFSKCC = map(data, function(x) sd(x$CCRAviaFSKCC)),
                 MeanCCRAviaDMSO = map(data, function(x) mean(x$CCRAviaDMSO)),
                 SDCCRAviaDMSO = map(data, function(x) sd(x$CCRAviaDMSO)),
                 SEMCCRAviaDMSO = map(data, function(x) sd(x$CCRAviaDMSO)/sqrt(nrow(x))),
                 CICCRAviaDMSO = map(data, function(x) qt(0.975, df = nrow(x) - 2) * sd(x$CCRAviaDMSO)/sqrt(nrow(x))),
                 MeanCCRAviaDMSO2 = map(data, function(x) mean(x$CCRAviaDMSO2)),
                 SDCCRAviaDMSO2 = map(data, function(x) sd(x$CCRAviaDMSO2)),
                 SEMCCRAviaDMSO2 = map(data, function(x) sd(x$CCRAviaDMSO2)/sqrt(nrow(x))),
                 CICCRAviaDMSO2 = map(data, function(x) qt(0.975, df = nrow(x) - 2) * sd(x$CCRAviaDMSO2)/sqrt(nrow(x)))) %>%
          unnest(c(MeanRAviaDMSO, SDRAviaDMSO, MeanRAviaFSK, SDRAviaFSK, MeanCCRAviaFSKCC, SDCCRAviaFSKCC, 
                   MeanCCRAviaDMSO, SDCCRAviaDMSO, SEMCCRAviaDMSO, CICCRAviaDMSO, 
                   MeanCCRAviaDMSO2, SDCCRAviaDMSO2, SEMCCRAviaDMSO2, CICCRAviaDMSO2, data)) %>%
          ungroup()
      }),
      Anova_Results_CCRAviaDMSO = map2(All_Analysis_Values, Exp, function(.x, .y) {
        Trim <- .x %>% filter(Ligand %notin% c("FSK", "DMSO"))
        compare_means(CCRAviaDMSO2 ~ Mutation, Trim, method = "anova") %>%
          bind_cols(tibble(Exp = .y))
      }),
      Ttest_Results_CCRAviaDMSO = map2(All_Analysis_Values, Exp, function(.x, .y){
        Trim <- .x %>% filter(Ligand %notin% c("FSK", "DMSO"))
        compare_means(CCRAviaDMSO2 ~ Mutation, Trim, method = "t.test", group.by = "Ligand") %>%
          bind_cols(tibble(Exp = .y))
      }),
      First_plots = map2(All_Analysis_Values, Exp, function(.x, .y){
        Plotoutput <- .x %>%
          filter(Ligand %notin% c("DMSO", "FSK")) %>%
          ggplot(aes(Ligand, CCRAviaDMSO2, color = Mutation)) + 
          geom_point(position = position_dodge(width = 0.3), size = 2.5) +
          geom_errorbar(aes(x = Ligand, group = Mutation, ymin = MeanCCRAviaDMSO2 - SDCCRAviaDMSO2, ymax = MeanCCRAviaDMSO2 + SDCCRAviaDMSO2), 
                        width = 0.2, color = "black", position = position_dodge(width = 0.3)) +
          labs(title = .y, y = bquote((cAMP[Ligand] - cAMP[DMSO])/cAMP[DMSO]), 
               caption = "Mean +/- SD.") +
          theme_classic() +
          theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
                axis.title.x = element_blank(),
                axis.text = element_text(face="bold", size = 10),
                axis.text.x = element_text(angle = 15, hjust = 1),
                axis.title.y = element_text(size = 14),
                legend.box = "horizontal",
                legend.background = element_rect(fill = "white", color = "black"),
                legend.text = element_text(face="bold"),
                legend.title = element_blank())
      }))

# I require an intermediate to combine separate experiments properly...
Exp_L1_wError <- rbindlist(Exp_L1_Processed$All_Analysis_Values, use.names=TRUE) %>%
  group_by(Sample, Ligand, Mutation) %>%
  summarize(CCRAviaDMSO = mean(CCRAviaDMSO)) %>% #Combining all technical replicates so I have a single value for each Bio Rep/Odor
  ungroup() %>%
  group_by(Ligand, Mutation) %>%
  nest() %>%
  mutate(MeanCCRAviaDMSO = map(data, function(x) mean(x$CCRAviaDMSO)),
         SDCCRAviaDMSO = map(data, function(x) sd(x$CCRAviaDMSO)),
         SEMCCRAviaDMSO = map(data, function(x) sd(x$CCRAviaDMSO)/sqrt(nrow(x))),
         CICCRAviaDMSO = map(data, function(x) qt(0.975, df = nrow(x) - 2) * sd(x$CCRAviaDMSO)/sqrt(nrow(x)))) %>%
  unnest(c(MeanCCRAviaDMSO, SDCCRAviaDMSO, SEMCCRAviaDMSO, CICCRAviaDMSO, data)) %>%
  ungroup()

Stats_Exp_L1 <- compare_means(CCRAviaDMSO ~ Mutation, Exp_L1_wError, method = "t.test", group.by = "Ligand") %>%
  filter(Ligand != "FSK") %>%
  select(Ligand, p) %>%
  arrange(p) %>%
  mutate(row = row_number(),
         HolmAdjustment = 0.05/(22 - row + 1), 
         PassHA = p < HolmAdjustment, 
         HolmCorrection = p * (22 - row + 1),
         PassHC = HolmCorrection < 0.05, 
         Pass_DunnSidakCorrection = p < 1 - (1 - 0.05)^(1/22)) %>%
  mutate(padj = case_when(HolmCorrection <= 0.0001 ~ "****", HolmCorrection <= 0.001~ "***", 
                          HolmCorrection <= 0.01 ~ "**", HolmCorrection <= 0.05 ~"*", TRUE ~ "ns")) %>%
  select(Ligand, padj) #selecting for left join

OdorOrder <- c("5-Methoxyguaiacol", "4-Bromoguaiacol", "4-Methylguaiacol", "4-Chloroguaiacol", "4-Methoxyguaiacol", "Guaiacol", "Ethyl vanillin", 
               "Vanillin", "4-Iodoguaiacol", "Dehydrodivanillin", "2-Hydroxy-4-Methoxybenzaldehyde", "Acetovanillone", "Vanillyl butyl ether",
               "5-Bromoguaiacol", "5-Fluoroguaiacol", "2-Ethoxyphenol", "4-Fluoroguaiacol", "2-Ethylphenol", "o-Cresol","Eugenol",
               "Menthoxypropanediol", "Salicylaldehyde")

Exp_L1_Final <- Exp_L1_wError %>%
  filter(Ligand %notin% c("FSK", "DMSO")) %>%
  group_by(Ligand) %>%
  nest() %>%
  mutate(Max = map(data, function(x) max(x$CCRAviaDMSO))) %>%
  unnest(c(Max, data)) %>%
  left_join(Stats_Exp_L1, by = "Ligand") %>%
  mutate(padj = ifelse(Mutation == "WT", "", padj)) %>%
  mutate(Ligand = factor(Ligand, levels = c("DMSO", "FSK", OdorOrder))) %>%
  ungroup()

# Figure 5A
#"Seven BioAssay experiments total, each odor tested at least three times. 7-8 biological replicates for each Genotype, tested in triplicate per experiment.
#p values Holm-corrected for multiple comparisons. * 0.05, ** 0.01, *** 0.001, **** 0.0001, ns not significant. Mean +/- SD."

Fig5A <- ggplot() + 
  geom_point(data = Exp_L1_Final, aes(Ligand, CCRAviaDMSO, color = Mutation), position = position_dodge(width = 0.3), size = 2, alpha = 0.8, stat = "unique") +  
  geom_errorbar(data = Exp_L1_Final, aes(Ligand, CCRAviaDMSO, color = Mutation, ymin = MeanCCRAviaDMSO - SDCCRAviaDMSO, ymax = MeanCCRAviaDMSO + SDCCRAviaDMSO), 
                position = position_dodge(width = 0.3), width = 0.3, stat = "unique") + 
  geom_text(data = Exp_L1_Final, aes(Ligand, y = Max + 0.4, label = padj), stat = "unique", size = 4) +
  labs(subtitle = "Liquid-Phase Responses for Guaiacol and Vanillin-related Ligands", y = "Activation Index (DMSO-Corrected)") +
  geom_vline(xintercept = 13.5, linetype = "dashed") +
  annotate("text", x = 7, y = 7, label = "Vanilla-like Odor quality", size = 5, fontface = "bold") +
  annotate("text", x = 18, y = 7, label = "Other Odor qualities", size = 5, fontface = "bold") +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(face="bold", size = 10),
        axis.text.x = element_text(angle = 15, hjust = 0.7, vjust = 0.85),
        axis.title.y = element_text(size = 12),
        legend.position = "top",
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(face="bold"),
        legend.title = element_blank())

# 06 - Liquid Dose Response Bioassays----
Label_list_LDR <- list(tibble(Exp = "Exp43b",
                              Sample = rep(c("PFH1721", "PFH1741", "PFH1742", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 45),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Guaiacol",168), rep("o-Cresol", 144)),
                              Dilution = rep(c(0, 5e-7, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9), each = 24)),
                       tibble(Exp = "Exp45a",
                              Sample = rep(c("PFH1721", "PFH1741", "PFH1742", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 45),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Guaiacol",168), rep("o-Cresol", 144)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), each = 24)),
                       tibble(Exp = "Exp47a",
                              Sample = rep(c("PFH1721", "PFH1741", "PFH1742", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 45),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Guaiacol",168), rep("o-Cresol", 144)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), each = 24)),
                       tibble(Exp = "Exp48_1a",
                              Sample = rep(c("PFH1721", "PFH1741", "PFH1742", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 45),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Vanillin",168), rep("2-Ethylphenol", 144)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), each = 24)),
                       tibble(Exp = "Exp48_2a",
                              Sample = rep(c("PFH1721", "PFH1741", "PFH1742", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 45),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Ethyl vanillin",168), rep("Dehydrodivanillin", 144)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), each = 24)),
                       tibble(Exp = "Exp51DSa",
                              Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 42),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 42),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Vanillin",144), rep("Ethyl vanillin", 144)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), each = 24)),
                       tibble(Exp = "Exp52DS1a",
                              Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 45),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("5-Bromoguaiacol",168), rep("4-Bromoguaiacol", 144)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), each = 24)),
                       tibble(Exp = "Exp52DS2a",
                              Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 45),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("2-Ethoxyphenol",168), rep("5-Methoxyguaiacol", 144)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), each = 24)),
                       tibble(Exp = "Exp53DRa",
                              Sample = rep(c("PFH0887", "PFH0888", "PFH1790", "PFH1885", "PFH1789", "PFH1792", "PFH1793", "PFH1884"), 42),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 42),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Vanillin",144), rep("Ethyl vanillin", 144)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), each = 24)),
                       tibble(Exp = "Exp53DRb",
                              Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 42),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 42),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Vanillin",144), rep("Ethyl vanillin", 144)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), each = 24)),
                       tibble(Exp = "Exp55a",
                              Sample = rep(c("PFH0887", "PFH0888", "PFH1790", "PFH1885", "PFH1789", "PFH1792", "PFH1793", "PFH1884"), 45),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("2-Ethylphenol",144), rep("2-Ethoxyphenol", 168)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9), each = 24)),
                       tibble(Exp = "Exp55b",
                              Sample = rep(c("PFH0887", "PFH0888", "PFH1790", "PFH1885", "PFH1789", "PFH1792", "PFH1793", "PFH1884"), 45),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("4-Bromoguaiacol",144), rep("5-Bromoguaiacol", 168)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9), each = 24)),
                       tibble(Exp = "Exp55c",
                              Sample = rep(c("PFH0887", "PFH0888", "PFH1790", "PFH1885", "PFH1789", "PFH1792", "PFH1793", "PFH1884"), 45),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Dehydrodivanillin",144), rep("5-Methoxyguaiacol", 168)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9), each = 24)),
                       tibble(Exp = "Exp56a",
                              Sample = rep(c("PFH0887", "PFH0888", "PFH1790", "PFH1885", "PFH1789", "PFH1792", "PFH1793", "PFH1884"), 45),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("4-Bromoguaiacol",144), rep("5-Bromoguaiacol", 168)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9), each = 24)),
                       tibble(Exp = "Exp56b",
                              Sample = rep(c("PFH0887", "PFH0888", "PFH1790", "PFH1885", "PFH1789", "PFH1792", "PFH1793", "PFH1884"), 45),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("2-Ethylphenol",144), rep("2-Ethoxyphenol", 168)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9), each = 24)),
                       tibble(Exp = "Exp56c",
                              Sample = rep(c("PFH0887", "PFH0888", "PFH1790", "PFH1885", "PFH1789", "PFH1792", "PFH1793", "PFH1884"), 45),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Dehydrodivanillin",144), rep("5-Methoxyguaiacol", 168)),
                              Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9), each = 24)))

Input_list_LDR <-
  list.files(path = "Input/Liquid Dose Response", 
             pattern = "*.csv", 
             full.names = TRUE) %>% 
  lapply(read_plus)

Input_list_LDR_names <-   
  list.files(path = "Input/Liquid Dose Response", 
             pattern = "*.csv") %>%
  str_replace(pattern = ".csv", 
              replacement = "")

names(Input_list_LDR) <- Input_list_LDR_names

# A rather convoluted piece of code that produces a first table of data for each experiment via an unlabeled intermediate
Input_tibble_LDR <- tibble(Assay = Input_list_LDR, Exp = Input_list_LDR_names) %>%
  left_join(myAssayOutput, by = "Exp") %>%
  mutate(Background = (Blank1 + Blank2)/2, 
         Measure = map2(Assay, Background, ~ c(t(as.matrix(.x))) - .y),
         Concentration = pmap(list(a, b, c, d, Measure), .f = function(a, b, c, d, Measure) {c * ((a - d)/(Measure - d) - 1)^(1/b)} ),
         Position = map(Assay, ~ 1:length(unlist(.x))), 
         Unbound = map(Concentration, ~ ifelse(is.na(.x), TRUE, FALSE)), 
         Concentration2 = pmap(list(Concentration, Measure, a, d), .f = function(Concentration, Measure, a, d) 
         {case_when(is.na(Concentration) & Measure > d ~ min(Concentration, na.rm= TRUE),
                    is.na(Concentration) & Measure < a ~ max(Concentration, na.rm = TRUE),
                    TRUE ~ Concentration)}),
         Combined = pmap(list(Measure, Position, Concentration, Unbound, Concentration2), 
                         .f = function (Measure, Position, Concentration, Unbound, Concentration2) {
                           tibble(Measure = Measure, Position = Position, Concentration = Concentration, 
                                  Unbound = Unbound, Concentration2 = Concentration2)}))

Exp_LDR <- tibble(Labels = Label_list_LDR) %>%
  mutate(Exp = map_chr(Labels, ~unique(.x$Exp))) %>%
  inner_join(Input_tibble_LDR, by = "Exp") %>% 
  mutate(Full_Table = map2(Combined, Labels, bind_cols)) %>%
  select(Exp, Full_Table, a, b, c, d) %>%
  mutate(QC = map(Full_Table, function(.x) {.x %>% group_by(Sample, Ligand, Dilution) %>%
      nest() %>%
      mutate(MeanC = map(data, function(x) mean(x$Concentration2)),
             SD = map(data, function(x) sd(x$Concentration2))) %>%
      unnest(c(data, MeanC, SD))  %>%
      ungroup() %>%
      mutate(CV = SD/MeanC * 100,
             SEM = SD/sqrt(3))}))

#The cleanup requires manual investigation and selection; there is no "correct" choice in all cases. 
Exp_LDR$QC[["Exp43b"]] <- Exp_LDR$QC[["Exp43b"]] %>%
  filter(Position != 170)
Exp_LDR$QC[["Exp47a"]] <- Exp_LDR$QC[["Exp47a"]] %>%
  filter(Position != 25)
Exp_LDR$QC[["Exp48_1a"]] <- Exp_LDR$QC[["Exp48_1a"]] %>%
  filter(Position != 27)
Exp_LDR$QC[["Exp51DSa"]] <- Exp_LDR$QC[["Exp51DSa"]] %>%
  filter(Position %notin% c(74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120))
Exp_LDR$QC[["Exp53DRb"]] <- Exp_LDR$QC[["Exp53DRb"]] %>%
  filter(Position != 22)
Exp_LDR$QC[["Exp55a"]] <- Exp_LDR$QC[["Exp55a"]] %>%
  filter(Sample != "PFH1793")

#This next step takes the DR files to a final form. The DR files will be exported for use with prism to extrapolate the DR curves
Exp_LDR_Export <- Exp_LDR %>%
  mutate(Concentration_Triplicate_mean = pmap(list(QC, a, b, c, d, Exp), function(QC, a, b, c, d, Exp) {QC %>% group_by(Sample, Mutation, Ligand, Dilution) %>%
      summarize(S_Measure = mean(Measure)) %>% #I take the mean of MEASURE instead of the mean of CONCENTRATION. 
      mutate(Concentration_S_Measure = c * ((a - d)/(S_Measure - d) - 1)^(1/b)) %>%
      ungroup() %>%
      bind_cols(tibble(Experiment = Exp))}),
      Corrected_CTm_Relative_DMSO = map(Concentration_Triplicate_mean, function(.x){
        ExpDMSO <- .x %>%
          filter(Ligand == "DMSO") %>%
          select(Sample, Concentration_S_Measure) %>%
          rename(DMSOconcentration = Concentration_S_Measure)
        
        CCRAviaDMSO_Table <- .x %>%
          left_join(ExpDMSO, by = "Sample") %>%
          mutate(CCRAviaDMSO = (Concentration_S_Measure - DMSOconcentration)/DMSOconcentration) %>%
          filter(Ligand %notin% c("DMSO", "FSK")) %>%
          select(Experiment, Sample, Mutation, Ligand, Dilution, CCRAviaDMSO)
      }))

Exp_LDR_Export <- rbindlist(Exp_LDR_Export[[9]], use.names=TRUE)

#Keeps the relative distance between WT and 10G4 while dropping to the lowest value to zero.
Zero_Baseline_LDR <- Exp_LDR_Export %>%
  group_by(Sample, Mutation, Ligand, Dilution) %>%
  summarize(meanCCRAviaDMSO = mean(CCRAviaDMSO)) %>% #Taking the mean of Samples that fit into the same group_by grouping so as to avoid over-representing a cilia sample
  ungroup() %>% 
  group_by(Mutation, Ligand, Dilution) %>% 
  summarise(meanCC = mean(meanCCRAviaDMSO)) %>% 
  ungroup() %>% 
  group_by(Ligand) %>% 
  summarise(minmeanCC = min(meanCC)) %>%
  ungroup()

#Export file to be entered into prism to generate 4PL curve results including EC50.
Exp_LDR_Prism_Export_Final <- Exp_LDR_Export %>%
  group_by(Sample, Mutation, Ligand, Dilution) %>%
  summarize(SmeanCCRAviaDMSO = mean(CCRAviaDMSO)) %>% #Taking the mean of Samples that fit into the same group_by grouping so as to avoid over-representing a cilia sample
  ungroup() %>%
  left_join(Zero_Baseline_LDR, by = "Ligand") %>%
  mutate(meanBaselineCCRAviaDMSO = SmeanCCRAviaDMSO - minmeanCC) %>%
  rename(SMBL_CCRAviaDMSO = meanBaselineCCRAviaDMSO)

write.csv(Exp_LDR_Prism_Export_Final, file = "Output/LDR_Prism_Input.csv")

# Import Prism results 
g <- read_csv("DR_10G4_4Pcurve_prism.csv") %>%
  mutate(Ligand = factor(Ligand, levels = c("DMSO", "FSK", "Guaiacol", "Vanillin",  "Ethyl vanillin", "4-Bromoguaiacol", "5-Bromoguaiacol", "5-Methoxyguaiacol", 
                                            "2-Ethoxyphenol", "2-Ethylphenol", "o-Cresol", "Dehydrodivanillin", "Acetovanillone", "Vanillyl butyl ether", 
                                            "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde",
                                            "5-Fluoroguaiacol", "4-Fluoroguaiacol", "4-Chloroguaiacol", "4-Iodoguaiacol",
                                            "4-Methoxyguaiacol", "4-Methylguaiacol"))) %>%
  filter(Ligand != "Dehydrodivanillin") %>%
  select(Ligand, Top, EC50) %>%
  rename(`Max Signal` = Top,
         `EC50 (Liquid)` = EC50) %>%
  arrange(Ligand) %>%
  tableGrob(rows=NULL, theme = ttheme_minimal()) 

g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
Fig5D <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(g))

#Error calculations for the means
Exp_LDR_wError <- Exp_LDR_Prism_Export_Final %>%
  group_by(Mutation, Ligand, Dilution) %>%
  nest() %>%
  mutate(MeanCCRAviaDMSO = map(data, function(x) mean(x$SMBL_CCRAviaDMSO)),
         SDCCRAviaDMSO = map(data, function(x) sd(x$SMBL_CCRAviaDMSO)),
         SEMCCRAviaDMSO = map(data, function(x) sd(x$SMBL_CCRAviaDMSO)/sqrt(nrow(x))),
         CICCRAviaDMSO = map(data, function(x) qt(0.975, df = nrow(x) - 2) * sd(x$SMBL_CCRAviaDMSO)/sqrt(nrow(x)))) %>%
  unnest(c(MeanCCRAviaDMSO, SDCCRAviaDMSO, SEMCCRAviaDMSO, CICCRAviaDMSO)) %>%
  ungroup()

#Making stat tables to combine with graphs.
Stat_Test_LDR <- Exp_LDR_Prism_Export_Final %>%
  group_by(Ligand) %>%
  nest() %>%
  mutate(Ttest_Dilution_pairs = map(data, ~ compare_means(SMBL_CCRAviaDMSO ~ Mutation, .x, method = "t.test", group.by = "Dilution"))) %>%
  ungroup() %>%
  unnest(Ttest_Dilution_pairs) %>%
  select(Ligand, Dilution, p.signif) %>%
  mutate(p.signif = ifelse(p.signif == "ns", "", p.signif))

Fig5B <- Exp_LDR_wError %>%
  left_join(Stat_Test_LDR, by = c("Ligand", "Dilution")) %>%
  mutate(p.signif = ifelse(Mutation == "WT", "", p.signif)) %>%
  mutate(Ligand = factor(Ligand, levels = c("DMSO", "FSK", "Guaiacol",  "Ethyl vanillin", "Vanillin", "o-Cresol", "2-Ethylphenol", "Dehydrodivanillin", 
                                            "Acetovanillone", "Vanillyl butyl ether", "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde",
                                            "5-Fluoroguaiacol", "4-Fluoroguaiacol", "5-Bromoguaiacol", "4-Bromoguaiacol", "4-Chloroguaiacol", "4-Iodoguaiacol",
                                            "5-Methoxyguaiacol", "4-Methoxyguaiacol", "4-Methylguaiacol", "2-Ethoxyphenol"))) %>%
  filter(Ligand %in% c("Guaiacol", "5-Bromoguaiacol")) %>% 
  ggplot() + 
  geom_point(aes(log10(Dilution), MeanCCRAviaDMSO, color = Mutation, group = Mutation), size = 3) + 
  geom_errorbar(aes(x = log10(Dilution), color = Mutation, group = Mutation, ymin = MeanCCRAviaDMSO - SDCCRAviaDMSO, ymax = MeanCCRAviaDMSO + SDCCRAviaDMSO), width = 0.2) + 
  geom_text(aes(x=log10(Dilution), y = MeanCCRAviaDMSO + SDCCRAviaDMSO + 0.08, label = p.signif), size = 3) +
  scale_x_continuous(breaks = c(-10, -8, -6, -4, -2, -0)) +
  facet_wrap(~Ligand, nrow = 1) +
  labs(subtitle = "Liquid-Phase Dose Response",
       y = "Activation Index (DMSO-Corrected)",
       x = "Ligand Concentration log10(M)") +
  theme_bw() +
  theme(plot.subtitle = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(face="bold", size = 10),
        axis.text.x = element_text(hjust = 1),
        axis.title = element_text(size = 12),
        legend.position="none",
        strip.text = element_text(size = 12))

Supp_toFig5B_01 <- Exp_LDR_wError %>%
  left_join(Stat_Test_LDR, by = c("Ligand", "Dilution")) %>%
  mutate(p.signif = ifelse(Mutation == "WT", "", p.signif)) %>%
  mutate(Ligand = factor(Ligand, levels = c("DMSO", "FSK", "Guaiacol", "Vanillin",  "Ethyl vanillin", "4-Bromoguaiacol", "5-Bromoguaiacol", "5-Methoxyguaiacol", 
                                            "2-Ethoxyphenol", "2-Ethylphenol", "o-Cresol", "Dehydrodivanillin", "Acetovanillone", "Vanillyl butyl ether", 
                                            "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde",
                                            "5-Fluoroguaiacol", "4-Fluoroguaiacol", "4-Chloroguaiacol", "4-Iodoguaiacol",
                                            "4-Methoxyguaiacol", "4-Methylguaiacol"))) %>%
  ggplot() + 
  geom_point(aes(log10(Dilution), MeanCCRAviaDMSO, color = Mutation, group = Mutation), size = 3) + 
  geom_errorbar(aes(x = log10(Dilution), color = Mutation, group = Mutation, ymin = MeanCCRAviaDMSO - SDCCRAviaDMSO, ymax = MeanCCRAviaDMSO + SDCCRAviaDMSO), width = 0.2) + 
  geom_text(aes(x=log10(Dilution), y = MeanCCRAviaDMSO + SDCCRAviaDMSO + 0.08, label = p.signif), size =3) +
  scale_x_continuous(breaks = c(-10, -8, -6, -4, -2, -0)) +
  facet_wrap(~Ligand, nrow = 1) +
  labs(subtitle = "Liquid-Phase Dose Response",
       y = "Activation Index (DMSO-Corrected)",
       x = "Ligand Concentration log10(M)") +
  theme_bw() +
  theme(plot.subtitle = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(face="bold", size = 10),
        axis.text.x = element_text(hjust = 1),
        axis.title = element_text(size = 12),
        legend.position="none",
        strip.text = element_text(size = 12))

# 07 - Vapor Dose Response Bioassays----
Label_list_VDR <- list(tibble(Exp = "Exp61",
                              Sample = rep(c("PFH887", "PFH888", "PFH1722", "PFH1790", "PFH1885", "MSS4541", 
                                             "PFH1547", "PFH1789", "PFH1792", "PFH1793", "PFH1884", "PFH1886"), 24),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "WT", "WT",
                                               "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 24),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 36), rep("GUA", 180), rep("Mineral Oil", 36)),
                              Delivery = c(rep("Liq", 108), rep("Vap", 180)),
                              Dilution = rep(c(7.040*1e-4, 1e-7, 5*1e-5, 1e-1, 1e-2, 1e-3, 1e-4, 0), each = 36)),
                       tibble(Exp = "Exp62",
                              Sample = rep(c("PFH1722", "PFH1722", "PFH1790", "PFH1790", "PFH1885", "PFH1885", 
                                             "PFH1547", "PFH1547", "PFH1884", "PFH1884", "PFH1886", "PFH1886",
                                             "PFH1722", "PFH1790", "PFH1885", "PFH1547", "PFH1884", "PFH1886"), 8),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "WT", "WT",
                                               "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het",
                                               "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 8),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 18), rep("5BG", 90), rep("Mineral Oil", 18)),
                              Delivery = c(rep("Liq", 54), rep("Vap", 90)),
                              Dilution = rep(c(7.040*1e-4, 1e-7, 5*1e-5, 1e-1, 1e-2, 1e-3, 1e-4, 0), each = 18)),
                       tibble(Exp = "Exp63",
                              Sample = rep(c("PFH1722", "PFH1722", "PFH1790", "PFH1790", "PFH1885", "PFH1885", 
                                             "PFH1547", "PFH1547", "PFH1884", "PFH1884", "PFH1886", "PFH1886",
                                             "PFH1722", "PFH1790", "PFH1885", "PFH1547", "PFH1884", "PFH1886"), 8),
                              Mutation = rep(c("WT", "WT", "WT", "WT", "WT", "WT",
                                               "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het",
                                               "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 8),
                              Ligand = c(rep(c("DMSO", "FSK"), each = 18), rep("5FG", 90), rep("Mineral Oil", 18)),
                              Delivery = c(rep("Liq", 54), rep("Vap", 90)),
                              Dilution = rep(c(7.040*1e-4, 1e-7, 5*1e-5, 1e-1, 1e-2, 1e-3, 1e-4, 0), each = 18)))

Input_list_VDR <-
  list.files(path = "Input/Vapor Dose Response", 
             pattern = "*.csv", 
             full.names = TRUE) %>% 
  lapply(read_plus)

Input_list_VDR_names <-   
  list.files(path = "Input/Vapor Dose Response", 
             pattern = "*.csv") %>%
  str_replace(pattern = ".csv", 
              replacement = "")

names(Input_list_VDR) <- Input_list_VDR_names

# A rather convoluted piece of code that produces a first table of data for each experiment via an unlabeled intermediate
Input_tibble_VDR <- tibble(Assay = Input_list_VDR, Exp = Input_list_VDR_names) %>%
  left_join(myAssayOutput, by = "Exp") %>%
  mutate(Background = (Blank1 + Blank2)/2, 
         Measure = map2(Assay, Background, ~ c(t(as.matrix(.x))) - .y),
         Concentration = pmap(list(a, b, c, d, Measure), .f = function(a, b, c, d, Measure) {c * ((a - d)/(Measure - d) - 1)^(1/b)} ),
         Position = map(Assay, ~ 1:length(unlist(.x))), 
         Unbound = map(Concentration, ~ ifelse(is.na(.x), TRUE, FALSE)), 
         Concentration2 = pmap(list(Concentration, Measure, a, d), .f = function(Concentration, Measure, a, d) 
         {case_when(is.na(Concentration) & Measure > d ~ min(Concentration, na.rm= TRUE),
                    is.na(Concentration) & Measure < a ~ max(Concentration, na.rm = TRUE),
                    TRUE ~ Concentration)}),
         Combined = pmap(list(Measure, Position, Concentration, Unbound, Concentration2), 
                         .f = function (Measure, Position, Concentration, Unbound, Concentration2) {
                           tibble(Measure = Measure, Position = Position, Concentration = Concentration, 
                                  Unbound = Unbound, Concentration2 = Concentration2)}))

Exp_VDR <- tibble(Labels = Label_list_VDR) %>%
  mutate(Exp = map_chr(Labels, ~unique(.x$Exp))) %>%
  inner_join(Input_tibble_VDR, by = "Exp") %>% 
  mutate(Full_Table = map2(Combined, Labels, bind_cols)) %>%
  select(Exp, Full_Table, a, b, c, d) %>%
  mutate(QC = map(Full_Table, function(.x) {.x %>% group_by(Sample, Ligand, Dilution, Delivery) %>%
      nest() %>%
      mutate(MeanC = map(data, function(x) mean(x$Concentration2)),
             SD = map(data, function(x) sd(x$Concentration2))) %>%
      unnest(c(data, MeanC, SD))  %>%
      ungroup() %>%
      mutate(CV = SD/MeanC * 100,
             SEM = SD/sqrt(3))}))

#This next step takes the DR files to a final form. The DR files will be exported for use with prism to extrapolate the DR curves
Exp_VDR_Export <- Exp_VDR %>%
  mutate(Concentration_Triplicate_mean = pmap(list(QC, a, b, c, d, Exp), function(QC, a, b, c, d, Exp) {QC %>% group_by(Sample, Mutation, Ligand, Dilution, Delivery) %>%
      summarize(S_Measure = mean(Measure)) %>% #I take the mean of MEASURE instead of the mean of CONCENTRATION. 
      mutate(Concentration_S_Measure = c * ((a - d)/(S_Measure - d) - 1)^(1/b)) %>%
      ungroup() %>%
      bind_cols(tibble(Experiment = Exp))}),
      Corrected_CTm_Relative_DMSO = map(Concentration_Triplicate_mean, function(.x){
        ExpMO <- .x %>%
          filter(Ligand == "Mineral Oil") %>%
          select(Sample, Concentration_S_Measure) %>%
          rename(MOconcentration = Concentration_S_Measure)
        ExpFSK <- .x %>%
          filter(Ligand == "FSK") %>%
          select(Sample, Concentration_S_Measure) %>%
          rename(FSKconcentration = Concentration_S_Measure)
        ExpDMSO <- .x %>%
          filter(Ligand == "DMSO") %>%
          select(Sample, Concentration_S_Measure) %>%
          rename(DMSOconcentration = Concentration_S_Measure)
        CCviaMO_RAviaFSK_Table <- .x %>%
          left_join(ExpMO, by = "Sample") %>%
          left_join(ExpFSK, by = "Sample") %>%
          left_join(ExpDMSO, by = "Sample") %>%
          mutate(CCviaMO_RAviaFSK = (Concentration_S_Measure - MOconcentration)/(FSKconcentration - DMSOconcentration)) %>%
          filter(Ligand %notin% c("DMSO", "FSK", "Mineral Oil")) %>%
          select(Experiment, Sample, Mutation, Ligand, Dilution, Delivery, CCviaMO_RAviaFSK)
      })) 

Exp_VDR_Export <- rbindlist(Exp_VDR_Export[[9]], use.names=TRUE)%>%
  filter(Ligand %notin% c("DMSO", "FSK", "Mineral Oil")) %>%
  filter(Delivery == "Vap")

#Keeps the relative distance between WT and 10G4 while dropping to the lowest value to zero.
Zero_Baseline_VDR <- Exp_VDR_Export %>%
  group_by(Sample, Mutation, Ligand, Dilution) %>%
  summarize(meanCCviaMO_RAviaFSK = mean(CCviaMO_RAviaFSK)) %>% #Taking the mean of Samples that fit into the same group_by grouping so as to avoid over-representing a cilia sample
  ungroup() %>% 
  group_by(Mutation, Ligand, Dilution) %>% 
  summarise(meanCC = mean(meanCCviaMO_RAviaFSK)) %>% 
  ungroup() %>% 
  group_by(Ligand) %>% 
  summarise(minmeanCC = min(meanCC)) %>%
  ungroup()

#Export file to be entered into prism to generate 4PL curve results including EC50.
Exp_VDR_Prism_Export_Final <- Exp_VDR_Export %>%
  group_by(Sample, Mutation, Ligand, Dilution) %>%
  summarize(SmeanCCviaMO_RAviaFSK = mean(CCviaMO_RAviaFSK)) %>% #Taking the mean of Samples that fit into the same group_by grouping so as to avoid over-representing a cilia sample
  ungroup() %>%
  left_join(Zero_Baseline_VDR, by = "Ligand") %>%
  mutate(meanBaselineCCviaMO_RAviaFSK = SmeanCCviaMO_RAviaFSK - minmeanCC) %>%
  rename(SMBL_CCviaMO_RAviaFSK = meanBaselineCCviaMO_RAviaFSK)

write.csv(Exp_VDR_Prism_Export_Final, file = "Output/VDR_Prism_Input.csv")

#Error calculations for the means
Exp_VDR_wError <- Exp_VDR_Prism_Export_Final %>%
  group_by(Mutation, Ligand, Dilution) %>%
  nest() %>%
  mutate(MeanCCviaMO_RAviaFSK = map(data, function(x) mean(x$SMBL_CCviaMO_RAviaFSK)),
         SDCCviaMO_RAviaFSK = map(data, function(x) sd(x$SMBL_CCviaMO_RAviaFSK)),
         SEMCCviaMO_RAviaFSK = map(data, function(x) sd(x$SMBL_CCviaMO_RAviaFSK)/sqrt(nrow(x))),
         CICCviaMO_RAviaFSK = map(data, function(x) qt(0.975, df = nrow(x) - 2) * sd(x$SMBL_CCviaMO_RAviaFSK)/sqrt(nrow(x)))) %>%
  unnest(c(MeanCCviaMO_RAviaFSK, SDCCviaMO_RAviaFSK, SEMCCviaMO_RAviaFSK, CICCviaMO_RAviaFSK)) %>%
  ungroup()

#Making stat tables to combine with graphs.
Stat_Test_VDR <- Exp_VDR_Prism_Export_Final %>%
  group_by(Ligand) %>%
  nest() %>%
  mutate(Ttest_Dilution_pairs = map(data, ~ compare_means(SMBL_CCviaMO_RAviaFSK ~ Mutation, .x, method = "t.test", group.by = "Dilution"))) %>%
  ungroup() %>%
  unnest(Ttest_Dilution_pairs) %>%
  select(Ligand, Dilution, p.signif) %>%
  mutate(p.signif = ifelse(p.signif == "ns", "", p.signif))

Fig5C <- Exp_VDR_wError %>%
  left_join(Stat_Test_VDR, by = c("Ligand", "Dilution")) %>%
  mutate(p.signif = ifelse(Mutation == "WT", "", p.signif)) %>%
  mutate(Ligand = case_when(Ligand == "GUA" ~ "Guaiacol", 
                   Ligand == "5BG" ~ "5-Bromoguaiacol",
                   Ligand == "5FG" ~ "5-Fluoroguaiacol")) %>%
  mutate(Ligand = factor(Ligand, levels = c("DMSO", "FSK", "Guaiacol",  "Ethyl vanillin", "Vanillin", "o-Cresol", "2-Ethylphenol", "Dehydrodivanillin", 
                                            "Acetovanillone", "Vanillyl butyl ether", "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde",
                                            "5-Fluoroguaiacol", "4-Fluoroguaiacol", "5-Bromoguaiacol", "4-Bromoguaiacol", "4-Chloroguaiacol", "4-Iodoguaiacol",
                                            "5-Methoxyguaiacol", "4-Methoxyguaiacol", "4-Methylguaiacol", "2-Ethoxyphenol"))) %>%
  filter(Ligand %in% c("Guaiacol", "5-Bromoguaiacol")) %>% 
  ggplot() + 
  geom_point(aes(log10(Dilution), MeanCCviaMO_RAviaFSK, color = Mutation, group = Mutation), size = 3) + 
  geom_errorbar(aes(x = log10(Dilution), color = Mutation, group = Mutation, ymin = MeanCCviaMO_RAviaFSK - SDCCviaMO_RAviaFSK, ymax = MeanCCviaMO_RAviaFSK + SDCCviaMO_RAviaFSK), width = 0.2) + 
  geom_text(aes(x=log10(Dilution), y = MeanCCviaMO_RAviaFSK + SDCCviaMO_RAviaFSK + 0.05, label = p.signif), size = 3) +
  scale_x_continuous(breaks = c(-10, -8, -6, -4, -2, -0)) +
  facet_wrap(~Ligand, nrow = 1) +
  labs(subtitle = "Vapor-Phase Dose Response",
       y = "Activation Index (MO- & FSK-Corrected)",
       x = "Ligand Concentration (Pre-Gauze) log10(M)") +
  theme_bw() +
  theme(plot.subtitle = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(face="bold", size = 10),
        axis.text.x = element_text(hjust = 1),
        axis.title = element_text(size = 12),
        legend.position="none", 
        strip.text = element_text(size = 12))

# 08 - Figure 5: Bioassay Assembled & Supplemental of Liquid Dose Response----

# Fig5A / (Fig5B | Fig5C | Fig5D) +
#   plot_annotation(tag_levels = 'A',
#                   title = "Our Bioassay reveals ligands with higher Activation Indeces for OR10G4 than previously validated Ligands",
#                   theme = theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))) & 
#   theme(plot.tag = element_text(size = 16, face = "bold"))
# 
# ggsave("Output/Composite_MF5ABCD.png", width = 35, height = 25, units = "cm")

Fig5A
ggsave("Output/BioAssay_MF5A.png", width = 35, height = 12, units = "cm")
Fig5B
ggsave("Output/BioAssay_MF5B.png", width = 13, height = 12, units = "cm")
Fig5C
ggsave("Output/BioAssay_MF5C.png", width = 13, height = 12, units = "cm")
Fig5D
ggsave("Output/BioAssay_MF5D.png", plot = Fig5D, width = 9, height = 8, units = "cm")

Supp_toFig5B_01

ggsave("Output/BioAssay_SF11.png", width = 45, height = 15, units = "cm")

# 09 - Comparing Odor Properties----
Odor_Table <- read_csv("Odor_Table.csv") 

Odor_Table_long <- Odor_Table %>%
  separate(OdorQuality, into = c("OQ1", "OQ2", "OQ3", "OQ4", "OQ5")) %>%
  pivot_longer(cols = 7:11, names_to = "OQ", values_to = "OdorQuality", values_drop_na = TRUE) %>%
  select(-OQ)

Nodes <- Odor_Table_long %>% 
  group_by(OdorQuality) %>% 
  mutate(medmedAI = median(medianAI), 
         instances = n()) %>%
  select(OdorQuality, medmedAI, instances) %>%
  unique() %>%
  filter(!is.na(OdorQuality))

Qualities <- Nodes %>%
  select(OdorQuality)

Pairs <- Qualities %>%
  rename(OQ2 = OdorQuality) %>%
  cross_join(Qualities) %>%
  rename(OQ1 = OdorQuality) %>%
  filter(OQ1 != OQ2) %>%
  filter(OQ1 < OQ2) %>%
  select(OQ1, OQ2)

Edges <- Odor_Table %>%
  cross_join(Pairs) %>%
  filter(str_detect(OdorQuality, OQ1) & str_detect(OdorQuality, OQ2)) %>%
  select(OQ1, OQ2, everything())

net <- graph_from_data_frame(d=Edges, vertices=Nodes, directed=T) 

ggraph(net, layout = "kk") +
  geom_edge_fan(aes(color = medianAI), linewidth = 0.75) +            
  geom_node_point(aes(color = medmedAI), size = 10) +
  geom_node_label(aes(label = name), size = 8, repel = TRUE) +
  scale_edge_color_gradientn(name = "Median Ligand A.I. (Edge)", colors = c("black", "purple", "blue", "green", "yellow", "red")) +
  scale_color_gradientn("Median Odor Quality A.I. (Node)", colors = c("black", "purple", "blue", "green", "yellow", "red")) +
  labs(title = "The Smoky odor quality has the highest median Activation Index for OR10G4 cilia") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggsave("Output/BioAssay_Network_MF5F.png", width = 30, height = 30, units = "cm")

lm(log10(OdorThreshold) ~ medianAI, data = Odor_Table) %>% summary()

Odor_Table %>%
  ggplot(aes(medianAI, log10(OdorThreshold))) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5, color = "red") +
  geom_point(size = 3) +
  geom_text(aes(label = Ligand), nudge_y = 0.1, size = 7) + 
  scale_x_continuous(limits = c(0,3)) +
  labs(title = "Activation Index is negatively correlated with Vapor Odor Thresholds", 
       subtitle = "Log10(OT) = -1.481 * AI + 0.8177; Adjusted R squared = 0.8756; p = 0.0001293",
       x = "Median Activation Index", 
       y = "log10(Odor Threshold (ng/L air))") +
  theme_classic() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(face="bold", size = 18),
        axis.title = element_text(size = 20),
        legend.position = "top",
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(face="bold"),
        legend.title = element_blank())

ggsave("Output/BioAssay_MF5E.png", width = 30, height = 30, units = "cm")

# 10 - End----
sessionInfo()