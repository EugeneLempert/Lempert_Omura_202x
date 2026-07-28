# This is a script to generate RTqPCR Figures and Tables for the 10G4 Paper.

# Set your working directory to the "RTqPCR 10G4" folder

#01 - Required Packages----
if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

if (!require(broom)) install.packages('broom')
library(broom)

if (!require(magrittr)) install.packages('magrittr')
library(magrittr)

if (!require(ggpubr)) install.packages('ggpubr')
library(ggpubr)

if (!require(scales)) install.packages('scales')
library(scales)

if (!require(grid)) install.packages('grid')
library(grid)

if (!require(gridExtra)) install.packages('gridExtra')
library(gridExtra)

if (!require(RColorBrewer)) install.packages('RColorBrewer')
library(RColorBrewer)

#02 - Functions----
`%notin%` <- Negate(`%in%`)

read_plus <- function(flnm) {
  read_delim(flnm, , delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE, skip = 1)
}



Quality_Evaluation <- function(Labeled_qPCR_data, 
                               MeltFile = NULL, 
                               MaxCq = 45, 
                               CqLower = 18, 
                               CqUpper = 37,
                               RepDiffAllow = 0.3,
                               gdnaTMallow = 0.5){
  
  message("The MaxCq is to set to ", MaxCq, " and Cq range is set between ", 
          CqLower, " and ", CqUpper, " and the allowed difference between replicates is set to ", RepDiffAllow)
  
  if(!is.null(MeltFile)){
    
    if(!identical(Labeled_qPCR_data$Pos, MeltFile$Pos)) stop ("Joining by Pos will fail. Review Files and/or File names.")
    
    Labeled_qPCR_data <- Labeled_qPCR_data %>%
      left_join(MeltFile, by = "Pos") %>%
      select(-Include, 
             -Color, 
             -Name, 
             -Status)
  }
  
  QC_table1 <- Labeled_qPCR_data %>%
    rename(Cq = Cp) %>%
    mutate(MaxCq = Cq == MaxCq, 
           NoCq = is.na(Cq),
           NotWithinCqRange = !between(Cq, CqLower, CqUpper))
  
  QC_table2 <- QC_table1 %>% 
    select(Sample, Gene, Cq) %>%
    group_by(Sample, Gene) %>%
    summarize(Sample = Sample, 
              Gene = Gene, 
              RepMean = mean(Cq, na.rm = TRUE),
              RepSd = sd(Cq, na.rm = TRUE),
              RepNA = sum(is.na(Cq))/length(Cq),
              RepMax = max(Cq, na.rm = TRUE), 
              RepMin = min(Cq, na.rm = TRUE),
              RepMed = median(Cq, na.rm = TRUE))
  
  QC_table1 <- QC_table1 %>% 
    full_join(QC_table2, by = c("Sample", "Gene"))
  
  QC_table1 <- QC_table1 %>%
    mutate(DiffRepMax = Cq - RepMax, 
           DiffRepMed = Cq - RepMed,
           DiffRepMin = Cq - RepMin, 
           PoorReplicate = ((abs(DiffRepMin) > RepDiffAllow) + (abs(DiffRepMed) > RepDiffAllow) + (abs(DiffRepMax) > RepDiffAllow)) > 1)
  
  if(!is.null(MeltFile)){
    
    QC_table1 <- QC_table1 %>%
      mutate(NoPeak = is.na(Tm1),
             FalseCq = !is.na(Cq) & is.na(Tm1))
    
    QC_table1 <- QC_table1 %>% 
      group_by(Gene) %>% 
      nest() %>% 
      mutate(gDNATM1diff = map(data, function(x) mean(x$Tm1[which(x$Sample == "WT_gDNA")], na.rm = TRUE) - x$Tm1)) %>% 
      unnest(c(data, gDNATM1diff)) %>% 
      mutate(BadTM1 = abs(gDNATM1diff) > gdnaTMallow)
    
    if("Tm2" %in% colnames(QC_table1)){
      QC_table1 <-QC_table1 %>%
        mutate(TwoPeaks = !is.na(Tm2))
    }
    
    QC_table1 <- QC_table1 %>%
      pivot_longer(cols = c(starts_with("Tm"), starts_with("Area"), starts_with("Height"), starts_with("Width")), 
                   names_to = c(".value", "Peak"), 
                   names_pattern = "(.*)(.)")
  }
  QC_table1 <- QC_table1 %>% distinct() %>% ungroup()
  
  QE_summary <- QC_table1 %>% 
    select(Pos, Sample, Gene, Cq, where(is.logical)) %>%
    distinct() %>% 
    summary()
  
  QE_sum2 <- QC_table1 %>% 
    select(Sample, Gene, where(is.logical)) %>% 
    distinct() %>% 
    group_by(Sample, Gene) %>% 
    summarize(across(everything(), ~sum(.x, na.rm = TRUE)))
  
  QE_list <- list(QE_table = QC_table1,
                  QE_summary = QE_summary, 
                  QE_detail = QE_sum2)
  
  return(QE_list)
}


Label_qPCR_data <- function(CqFile, 
                            Style,
                            SampleNames, 
                            GeneNames,
                            Replicates = NULL,
                            ReplicateNames = 1:Replicates, 
                            ...) {
  
  #The ... is a variable length input for all "Conditions" to incorporate
  
  cDNAObject <- CqFile %>%
    select(-Include, 
           -Color, 
           -Name, 
           -Status, 
           -Concentration, 
           -Standard)
  
  #The style options are Simple for continuous lines, Standard for blocks, and Custom for customized continuous line labeling
  
  R <- length(ReplicateNames)
  S <- length(SampleNames)
  G <- length(GeneNames)
  Z <- floor(24/S)
  groups <- trunc(G/Z) 
  extra <- G%%Z
  
  if(Style == "Custom"){
    cDNAObject$Sample <- SampleNames
    cDNAObject$Gene <- GeneNames
    cDNAObject$Replicate <- ReplicateNames
    cDNAObject <- cDNAObject %>% bind_cols(data.frame(...))
    
  } else if(Style == "Simple"){
    #The ... related code should convert vectors of conditions into the appropriate length columns to bind
    cDNAObject$Sample <- rep(SampleNames, times = R * G)
    cDNAObject$Gene <- rep(GeneNames, each = S * R)
    cDNAObject$Replicate <- rep(rep(ReplicateNames, each = S), times = G)
    cDNAObject <- cDNAObject %>% bind_cols(as.data.frame(map(list(...), rep, times = R * G)))
    
  } else if(Style == "Standard"){
    
    #The code in standard was taken and modified from stackoverflow by Tony Breyal.
    #It allows the columns/labels to work when using the "standard" style, which does not pattern nicely when the plate is loaded this way.
    
    f1 <- as.character(sort(rep(1:groups, Z)))
    f <- as.character(c(f1, rep("Last", extra)))
    g <- split(GeneNames, f)
    
    if(extra > 0){
      g.names <- names(g[-length(g)])
      g.names.ordered <- as.character(sort(as.numeric(g.names)))
      g.names.ordered <- c(g.names.ordered, "Last")
    } else {
      g.names.ordered <- as.character(sort(as.numeric(names(g))))
    }
    
    #The g.names.ordered object allows the proper repetition of the Gene names. Using map was my inspired thought.
    
    cDNAObject$Sample <- rep(SampleNames, times = R * G)
    cDNAObject$Gene <- (g[g.names.ordered]) %>%
      map(~rep(.x, each = S, times = R)) %>% 
      unlist() %>%
      unname()
    cDNAObject$Replicate <- c(rep(ReplicateNames, each = S* Z, times = groups), rep(ReplicateNames, each = S * extra))
    cDNAObject <- cDNAObject %>% bind_cols(as.data.frame(map(list(...), rep, times = R * G)))
  }
  
  return(cDNAObject)
}

#03 - Import data, DESeq results, and DVI----
Input_list_G4 <- list.files(path = "Input/OR10G4", pattern = "*.txt", full.names = TRUE) %>%
  lapply(read_plus)
Input_list_G4_names <- list.files(path = "Input/OR10G4", pattern = "*.txt") %>%
  str_replace(pattern = ".txt", replacement = "")
names(Input_list_G4) <- Input_list_G4_names

Input_list_P <- list.files(path = "Input/Primers", pattern = "*.txt", full.names = TRUE) %>%
  lapply(read_plus)
Input_list_P_names <- list.files(path = "Input/Primers", pattern = "*.txt") %>%
  str_replace(pattern = ".txt", replacement = "")
names(Input_list_P) <- Input_list_P_names

Input_list_C <- list.files(path = "Input/Compare", pattern = "*.txt", full.names = TRUE) %>%
  lapply(read_plus)
Input_list_C_names <- list.files(path = "Input/Compare", pattern = "*.txt") %>%
  str_replace(pattern = ".txt", replacement = "")
names(Input_list_C) <- Input_list_C_names

DVI <- read_csv("Input/DVI.csv")

DE <- read_csv("Input/DESeq10G4_ExportRTqPCR.csv") %>%
  select(-Corrected_padj)

#04 - Primer Efficiency & Plots----
#Resulting plots and tables not in figures, but exist in supplemental excel file.
PE1 <- Input_list_P$`Primer-Cq-220201` %>%
  select(Cp) %>%
  bind_cols(Input_list_P$`Primer-Melt-220201`) %>%
  select(!c(Include, Pos, Color, Name, Status))

PE1$Target <- c(rep("Acsm4", 12), rep("Acss2", 12), rep("Acsm4", 12), rep("Acss2", 12),
                rep("Olfr358", 30), rep("Olfr609", 30), rep("Olfr690", 30), rep("Olfr1154", 30), 
                rep("Olfr596_603", 30), rep("Olfr390", 30), rep("Olfr510", 30))

PE1$Sample <- c(rep(c(rep("MSS-4163", 9), "Blank", "MSS-3184", "MSS-3256") ,4), 
                rep(c(rep("MSS-4166", 12), "Blank", "MSS-3254", "MSS-3257") ,8), 
                rep(c(rep("MSS-4169", 12), "Blank", "MSS-3255", "MSS-3258") ,6))

PE1$Dilution <- c(rep(c(1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 0, 1e-1, 1e-1), 4), 
                  rep(c(2^(0), 2^(-2), 2^(-4), 2^(-6), 2^(-7), 2^(-8), 2^(-10), 2^(-9), 2^(-11), 2^(-12), 2^(-13), 2^(-14), 0, 1e-1, 1e-1), 8), 
                  rep(c(2^(0), 2^(-2), 2^(-4), 2^(-6), 2^(-7), 2^(-8), 2^(-9), 2^(-10), 2^(-11), 2^(-12), 2^(-13), 2^(-14), 0, 1e-1, 1e-1), 6))

PE1$Genotype <- c(rep(c(rep("WT", 9), "Blank", "WT_1.1.2", "OR_1.1.2") ,4), 
                  rep(c(rep("WT", 12), "Blank", "WT_1.1.2", "OR_1.1.2") ,8), 
                  rep(c(rep("WT", 12), "Blank", "WT_1.1.2", "OR_1.1.2") ,6))

PE2 <- tibble(Dilution = c(1, 1, 0.1, 0.01, 0.001, 0.001, 0.0001, 0.0001, 0.00005, 0.00005, 0.000025, 0.000025, 0.0000125, 0.0000125),
              Target = rep("Slc25a35", 14),
              Cp = c(18.98, 18.80, 21.90, 25.41, 28.99, 28.73, 31.96, 31.01, 32.89, 32.48, 34.28, 33.61, 35.23, 34.28))

PE3 <- Input_list_P$`220317_Primer_Cq` %>%
  select(Cp) %>%
  bind_cols(Input_list_P$`220317_Primer_Melt`) %>%
  select(!c(Include, Pos, Color, Name, Status))

PE3$Target <- c(rep("OR1A1_m", 24), rep("OR5AN1", 24), rep("GCaMP", 24))

PE3$Sample <- c(rep(c(rep("MSS-3257", 11), "Blank"), 2), 
                          rep(c(rep("MSS-3278", 11), "Blank"), 2), 
                          rep(c(rep("MSS-4361", 11), "Blank"), 2))

PE3$Dilution <- c(rep(c(1, 0.1, 0.025, 0.00625, 0.0015625, 0.0015625/2, 0.0015625/4, 0.0015625/8, 0.0015625/16, 0.0015625/32, 0.0015625/64, 0), 6))

PE3$Genotype <- c(rep(c(rep("OR1A1_V1.1.2", 11), "Blank"), 2), 
                            rep(c(rep("OR5AN1 L15", 11), "Blank"), 2), 
                            rep(c(rep("Or10G4 L10", 11), "Blank"), 2))

PE3 <- PE3 %>%
  filter(Sample != "Blank") %>%
  filter(Cp <= 34) %>%
  select(Dilution, Target, Cp)

PE4 <- Input_list_P$`220429_Primer_Cq` %>%
  select(Cp) %>%
  bind_cols(Input_list_P$`220429_Primer_Melt`) %>%
  select(!c(Include, Pos, Color, Name, Status))

PE4$Target <- c(rep("OR10G4_CDS", 54), rep("OR10G4_mRNA", 54))

PE4$Sample <- c(rep(c(rep("MSS-4164", 8), "Blank", rep("MSS-4165", 8), "Blank"), 6))

PE4$Dilution <- c(rep(c(1, 2^(-2), 2^(-4), 2^(-6), 2^(-8), 2^(-10), 2^(-12), 2^(-14), 0), 12))

PE4$Genotype <- c(rep(c(rep("OR10G4 Line 7", 8), "Blank"), 12))

PE4 <- PE4 %>%
  select(Sample, Genotype, Dilution, Target, Cp) %>%
  filter(Sample != "Blank") %>%
  filter(Cp <= 34)

MeanDifferenceCorrection <- PE4 %>% 
  group_by(Sample, Genotype, Dilution, Target) %>% 
  summarize(MeanCp = mean(Cp)) %>% 
  pivot_wider(names_from = Sample, values_from = MeanCp) %>% 
  mutate(SampleDifference = `MSS-4164` - `MSS-4165`) %>% 
  arrange(Target, desc(Dilution)) %>% 
  ungroup() %>% 
  group_by(Target) %>% 
  summarize(MeanDifference = mean(SampleDifference))

PE4 <- PE4 %>%
  left_join(MeanDifferenceCorrection, by = "Target") %>%
  mutate(NewCp = ifelse(Sample == "MSS-4164", Cp - MeanDifference, Cp)) %>%
  select(-c(Cp, MeanDifference)) %>%
  rename(Cp = NewCp) %>%
  select(Dilution, Target, Cp)

PE5 <- PE1 %>%
  select(Sample, Genotype, Dilution, Target, Cp, Tm1) %>%
  filter(Genotype == "WT" & !is.na(Tm1) & !is.na(Cp)) %>%
  select(Dilution, Target, Cp) %>%
  bind_rows(PE2) %>%
  bind_rows(PE3) %>%
  bind_rows(PE4) %>%
  filter(Cp < 33) %>%
  mutate(Target = factor(Target, levels = c("Acsm4", "Acss2", "Slc25a35", 
                                            "Olfr358", "Olfr390", "Olfr510", "Olfr596_603", "Olfr609", "Olfr690", "Olfr1154", 
                                            "OR5AN1", "OR1A1_m", "GCaMP", "OR10G4_CDS", "OR10G4_mRNA"))) %>%
  filter(Target %notin% c("OR5AN1", "OR1A1_m", "OR10G4_mRNA"))

PE_Plot <- PE5 %>%
  rename(Gene = Target) %>%
  ggplot(aes(log10(Dilution), Cp, color = Gene)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_manual(values=c(brewer.pal(8, "Dark2"), brewer.pal(7, "Set2"))) +
  theme_bw() +
  labs(title = "Primer Efficiency calculated from RT-qPCR dilution curves", 
       x = "Log10 of Relative Template Concentration", 
       y = "RT-qPCR Crossing Point value") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.15, 0.15),
        legend.background = element_rect(fill = "white", color = "black")) +
  guides(colour = guide_legend(ncol = 2))

PE6 <- PE5 %>%
  group_by(Target) %>%
  nest() %>%
  mutate(fit = map(data, ~lm(Cp ~ log10(Dilution), data = .x))) %>%
  mutate(rsq = map_dbl(fit, ~ summary(.x)[["r.squared"]])) %>%
  mutate(coef =  map(fit, tidy, conf.int = TRUE)) %>%
  unnest(coef) %>%
  filter(term == "log10(Dilution)") %>%
  mutate(conflow_eff = 10^(-1/conf.low),
         estimate_efficiency = 10^(-1/estimate),
         confhigh_eff = 10^(-1/conf.high)) %>%
  select(Target, rsq, estimate, estimate_efficiency, std.error) %>%
  rename(Gene = Target, R_Squared = rsq, Slope_Estimate = estimate, Primer_Efficiency = estimate_efficiency, Std_error = std.error)

PE7 <- PE5 %>%
  group_by(Target) %>%
  nest() %>%
  mutate(fit = map(data, ~lm(Cp ~ log10(Dilution), data = .x))) %>%
  mutate(rsq = map_dbl(fit, ~ summary(.x)[["r.squared"]])) %>%
  mutate(coef =  map(fit, tidy, conf.int = TRUE)) %>%
  unnest(coef) %>%
  filter(term == "(Intercept)") %>%
  select(Target, estimate) %>%
  rename(Gene = Target, Intercept_Estimate = estimate) %>%
  full_join(PE6) %>%
  select(Gene, R_Squared, Intercept_Estimate, Slope_Estimate, Primer_Efficiency, Std_error)

PE8 <- PE5 %>% 
  group_by(Target) %>% 
  summarize(Max_Cp = max(Cp)) %>% 
  rename(Gene = Target) %>% 
  full_join(PE7)

PE_Table <- PE8 %>%
  mutate(across(3:7, ~round(.x, 3))) %>%
  mutate(Gene = factor(Gene, levels = c("Acsm4", "Acss2", "Slc25a35", 
                                        "Olfr358", "Olfr390", "Olfr510", "Olfr596_603", "Olfr609", "Olfr690", "Olfr1154",
                                        "OR5AN1", "OR1A1_m", "GCaMP", "OR10G4_CDS", "OR10G4_mRNA"))) %>%
  arrange(Gene) %>%
  tableGrob(rows=NULL)

grid.arrange(PE_Plot, PE_Table, heights = c(1.0, 0.7), as.table = TRUE)

g <- arrangeGrob(PE_Plot, PE_Table, nrow=2, heights = c(1.2, 0.7))

ggsave("Output/Supp_Fig1_Primer_Efficiency.pdf",plot = g, width = 30, height = 30, units = "cm")

Primer_Efficiency <- PE8 %>% 
  select(Gene, Primer_Efficiency)

Primer_Efficiency_Table <- PE8 %>%
  mutate(across(3:7, ~round(.x, 3))) %>%
  mutate(Gene = factor(Gene, levels = c("Acsm4", "Acss2", "Slc25a35", 
                                        "Olfr358", "Olfr390", "Olfr510", "Olfr596_603", "Olfr609", "Olfr690", "Olfr1154",
                                        "OR5AN1", "OR1A1_m", "GCaMP", "OR10G4_CDS", "OR10G4_mRNA"))) %>%
  arrange(Gene)

openxlsx::write.xlsx(Primer_Efficiency_Table, file = "Output/Primer_Efficiency.xlsx")

#05 - Generating dataset----
#220228 sets
OR10G4_7J10_1 <- Label_qPCR_data(CqFile = Input_list_G4$`220228_Exp_Cq_1`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c(rep(c("MSS-4359", "MSS-4360", "MSS-4361", "MSS-4362", "MSS-4363", "MSS-4364", "MSS-4365", "MSS-4368"), 3),
                                                       rep(c("PFH-609", "PFH-610", "PFH-611", "PFH-613", "PFH-614", "PFH-617", "PFH-618", "PFH-733"), 3)), 6), 
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35"), each = 48),
                                               rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35"), each = 6)),
                                 Mutation =    c(rep(c(rep(c("L10_WT", "L10_WT", "L10_Het", "L10_Het", "L10_WT", "L10_Het", "L10_Het", "L10_WT"), 3),
                                                       rep(c("L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het"), 3)), 6),
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 Mouse_line = c(rep(c(rep("L10", 24),  rep("L7J", 24)), 6), rep("Ctrl", 36)),
                                 Plate = rep("One", times = 324),
                                 ReplicateNames = c(rep(c("One", "Two", "Three"), times = 12, each = 8), rep(c("One", "Two", "Three"), times = 12)))

OR10G4_7J10_QE1 <- Quality_Evaluation(OR10G4_7J10_1, MeltFile = Input_list_G4$`220228_Exp_Melt_1`, RepDiffAllow = 0.6)

OR10G4_7J10_2 <- Label_qPCR_data(CqFile = Input_list_G4$`220228_Exp_Cq_2`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c(rep(c("MSS-4359", "MSS-4360", "MSS-4361", "MSS-4362", "MSS-4363", "MSS-4364", "MSS-4365", "MSS-4368"), 3),
                                                       rep(c("PFH-609", "PFH-610", "PFH-611", "PFH-613", "PFH-614", "PFH-617", "PFH-618", "PFH-733"), 3)), 6), 
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35"), each = 48),
                                               rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35"), each = 6)),
                                 Mutation =    c(rep(c(rep(c("L10_WT", "L10_WT", "L10_Het", "L10_Het", "L10_WT", "L10_Het", "L10_Het", "L10_WT"), 3),
                                                       rep(c("L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het"), 3)), 6),
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 Mouse_line = c(rep(c(rep("L10", 24), rep("L7J", 24)), 6), rep("Ctrl", 36)),
                                 Plate = rep("Two", times = 324),
                                 ReplicateNames = c(rep(c("One", "Two", "Three"), times = 12, each = 8), rep(c("One", "Two", "Three"), times = 12)))

OR10G4_7J10_QE2 <- Quality_Evaluation(OR10G4_7J10_2,  MeltFile = Input_list_G4$`220228_Exp_Melt_2`, RepDiffAllow = 0.6)

OR10G4_7J10_3 <- Label_qPCR_data(CqFile = Input_list_G4$`220228_Exp_Cq_3`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c(rep(c("MSS-4359", "MSS-4360", "MSS-4361", "MSS-4362", "MSS-4363", "MSS-4364", "MSS-4365", "MSS-4368"), 3),
                                                       rep(c("PFH-609", "PFH-610", "PFH-611", "PFH-613", "PFH-614", "PFH-617", "PFH-618", "PFH-733"), 3)), 6), 
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35"), each = 48),
                                               rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35"), each = 6)),
                                 Mutation =    c(rep(c(rep(c("L10_WT", "L10_WT", "L10_Het", "L10_Het", "L10_WT", "L10_Het", "L10_Het", "L10_WT"), 3),
                                                       rep(c("L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het"), 3)), 6),
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 Mouse_line = c(rep(c(rep("L10", 24),  rep("L7J", 24)), 6), rep("Ctrl", 36)),
                                 Plate = rep("Three", times = 324),
                                 ReplicateNames = c(rep(c("One", "Two", "Three"), times = 12, each = 8), rep(c("One", "Two", "Three"), times = 12)))

OR10G4_7J10_QE3 <- Quality_Evaluation(OR10G4_7J10_3,   MeltFile = Input_list_G4$`220228_Exp_Melt_3`, RepDiffAllow = 0.6)

OR10G4_7J10_QE3x <- OR10G4_7J10_QE3$QE_table %>% filter(Pos %notin% c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "E17", "E18", "E19", "E20", "E21", "E22", "E23", "E24", "L24")) 

#220301 sets
OR10G4_7J10_4 <- Label_qPCR_data(CqFile = Input_list_G4$`220301_Exp_Cq_1`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c(rep(c("MSS-4359", "MSS-4360", "MSS-4361", "MSS-4362", "MSS-4363", "MSS-4364", "MSS-4365", "MSS-4368"), 3),
                                                       rep(c("PFH-609", "PFH-610", "PFH-611", "PFH-613", "PFH-614", "PFH-617", "PFH-618", "PFH-733"), 3)), 6), 
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 GeneNames = c(rep(c("Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 48),
                                               rep(c("Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 6)),
                                 Mutation =    c(rep(c(rep(c("L10_WT", "L10_WT", "L10_Het", "L10_Het", "L10_WT", "L10_Het", "L10_Het", "L10_WT"), 3),
                                                       rep(c("L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het"), 3)), 6),
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 Mouse_line = c(rep(c(rep("L10", 24), rep("L7J", 24)), 6), rep("Ctrl", 36)),
                                 Plate = rep("Four", times = 324),
                                 ReplicateNames = c(rep(c("One", "Two", "Three"), times = 12, each = 8), rep(c("One", "Two", "Three"), times = 12)))

OR10G4_7J10_QE4 <- Quality_Evaluation(OR10G4_7J10_4,  MeltFile = Input_list_G4$`220301_Exp_Melt_1`, RepDiffAllow = 0.6)

OR10G4_7J10_5 <- Label_qPCR_data(CqFile = Input_list_G4$`220301_Exp_Cq_2`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c(rep(c("MSS-4359", "MSS-4360", "MSS-4361", "MSS-4362", "MSS-4363", "MSS-4364", "MSS-4365", "MSS-4368"), 3),
                                                       rep(c("PFH-609", "PFH-610", "PFH-611", "PFH-613", "PFH-614", "PFH-617", "PFH-618", "PFH-733"), 3)), 6), 
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 GeneNames = c(rep(c("Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 48),
                                               rep(c("Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 6)),
                                 Mutation =    c(rep(c(rep(c("L10_WT", "L10_WT", "L10_Het", "L10_Het", "L10_WT", "L10_Het", "L10_Het", "L10_WT"), 3),
                                                       rep(c("L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het"), 3)), 6),
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 Mouse_line = c(rep(c(rep("L10", 24),  rep("L7J", 24)), 6), rep("Ctrl", 36)),
                                 Plate = rep("Five", times = 324),
                                 ReplicateNames = c(rep(c("One", "Two", "Three"), times = 12, each = 8), rep(c("One", "Two", "Three"), times = 12)))

OR10G4_7J10_QE5 <- Quality_Evaluation(OR10G4_7J10_5, MeltFile = Input_list_G4$`220301_Exp_Melt_2`, RepDiffAllow = 0.6)

OR10G4_7J10_6 <- Label_qPCR_data(CqFile = Input_list_G4$`220301_Exp_Cq_3`, 
                                 Style = "Custom", 
                                 SampleNames = c(rep(c(rep(c("MSS-4359", "MSS-4360", "MSS-4361", "MSS-4362", "MSS-4363", "MSS-4364", "MSS-4365", "MSS-4368"), 3),
                                                       rep(c("PFH-609", "PFH-610", "PFH-611", "PFH-613", "PFH-614", "PFH-617", "PFH-618", "PFH-733"), 3)), 6), 
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 GeneNames = c(rep(c("Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 48),
                                               rep(c("Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 6)),
                                 Mutation =    c(rep(c(rep(c("L10_WT", "L10_WT", "L10_Het", "L10_Het", "L10_WT", "L10_Het", "L10_Het", "L10_WT"), 3),
                                                       rep(c("L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het", "L7J_WT", "L7J_Het"), 3)), 6),
                                                 rep(c(rep("WT_gDNA", 3), rep("Blank", 3)), 6)),
                                 Mouse_line = c(rep(c(rep("L10", 24),  rep("L7J", 24)), 6), rep("Ctrl", 36)),
                                 Plate = rep("Six", times = 324),
                                 ReplicateNames = c(rep(c("One", "Two", "Three"), times = 12, each = 8), rep(c("One", "Two", "Three"), times = 12)))

OR10G4_7J10_QE6 <- Quality_Evaluation(OR10G4_7J10_6, MeltFile = Input_list_G4$`220301_Exp_Melt_3`, RepDiffAllow = 0.6)

#220412 set
Data1 <- Label_qPCR_data(CqFile = Input_list_G4$`220412_Exp_Ct`, 
                         Style = "Custom", 
                         SampleNames = c(rep(c("MSS-4163", "MSS-4164", "MSS-4165", "MSS-4166", "MSS-4169", "MSS-4170", "MSS-4171", "WT_gDNA"), 27), 
                                         rep("UPWater", 27)),
                         GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 24),
                                       rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 3)),
                         Mutation = c(rep(c("WT", "Line 7", "Line 7", "WT", "WT", "Line 7", "WT", "WT_gDNA"), 27), rep("UPWater", 27)),
                         Plate = rep("Seven", 243),   
                         ReplicateNames = c(rep(c("One", "Two", "Three"), times = 9, each = 8), rep(c("One", "Two", "Three"), times = 9)))

Data_QE_1 <- Quality_Evaluation(Data1, MeltFile = Input_list_G4$`220412_Exp_Melt`,  RepDiffAllow = 0.6, CqLower = 17)

#220413 set
Data2 <- Label_qPCR_data(CqFile = Input_list_G4$`220413_Exp_Ct`, 
                         Style = "Custom", 
                         SampleNames = c(rep(c("MSS-4163", "MSS-4164", "MSS-4165", "MSS-4166", "MSS-4169", "MSS-4170", "MSS-4171", "WT_gDNA"), 27), 
                                         rep("UPWater", 27)),
                         GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 24),
                                       rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 3)),
                         Mutation = c(rep(c("WT", "Line 7", "Line 7", "WT", "WT", "Line 7", "WT", "WT_gDNA"), 27),  rep("UPWater", 27)),
                         Plate = rep("Eight", 243),   
                         ReplicateNames = c(rep(c("One", "Two", "Three"), times = 9, each = 8), rep(c("One", "Two", "Three"), times = 9)))

Data_QE_2 <- Quality_Evaluation(Data2,  MeltFile = Input_list_G4$`220413_Exp_Melt`, RepDiffAllow = 0.6, CqLower = 17)

#220414
Data3 <- Label_qPCR_data(CqFile = Input_list_G4$`220414_Exp_Ct`, 
                         Style = "Custom", 
                         SampleNames = c(rep(c("MSS-4163", "MSS-4164", "MSS-4165", "MSS-4166", "MSS-4169", "MSS-4170", "MSS-4171", "WT_gDNA"), 27), 
                                         rep("UPWater", 27)),
                         GeneNames = c(rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 24),
                                       rep(c("Olfr358", "Olfr690", "Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"), each = 3)),
                         Mutation = c(rep(c("WT", "Line 7", "Line 7", "WT", "WT", "Line 7", "WT", "WT_gDNA"), 27), rep("UPWater", 27)),
                         Plate = rep("Nine", 243),   
                         ReplicateNames = c(rep(c("One", "Two", "Three"), times = 9, each = 8), rep(c("One", "Two", "Three"), times = 9)))

Data_QE_3 <- Quality_Evaluation(Data3, MeltFile = Input_list_G4$`220414_Exp_Melt`, RepDiffAllow = 0.6, CqLower = 17)

#----This section to clean the Data
#The ...Fullx version is for CV calculations only?
OR10G4_7J10_Fullx <- OR10G4_7J10_QE1$QE_table %>%
  bind_rows(OR10G4_7J10_QE2$QE_table) %>%
  bind_rows(OR10G4_7J10_QE3$QE_table) %>%
  bind_rows(OR10G4_7J10_QE4$QE_table) %>%
  bind_rows(OR10G4_7J10_QE5$QE_table) %>%
  bind_rows(OR10G4_7J10_QE6$QE_table)

#This version removes some bad rows from Run 3, whereby I loaded twice and that replicate set is inappropriate. This is the correct one for the rest of the analysis.
OR10G4_7J10_Full <- OR10G4_7J10_QE1$QE_table %>%
  bind_rows(OR10G4_7J10_QE2$QE_table) %>%
  bind_rows(OR10G4_7J10_QE3x) %>%
  bind_rows(OR10G4_7J10_QE4$QE_table) %>%
  bind_rows(OR10G4_7J10_QE5$QE_table) %>%
  bind_rows(OR10G4_7J10_QE6$QE_table)

OR10G4_7J10_Full <- OR10G4_7J10_Full %>%
  filter(FalseCq == FALSE) %>%
  filter(NoPeak == FALSE) %>%
  filter(MaxCq == FALSE) %>%
  filter(Sample %notin% c("WT_gDNA", "Blank")) %>%
  filter(Peak == 1)

Data_PreInput1 <- Data_QE_1$QE_table %>%
  filter(NoCq == FALSE) %>%
  filter(MaxCq == FALSE) %>%
  filter(NotWithinCqRange == FALSE) %>% 
  filter(FalseCq == FALSE) %>% 
  filter(Cq < 33) %>% 
  filter(Sample != "WT_gDNA") %>%
  filter(PoorReplicate == FALSE) %>%
  filter(Peak == 1)

Data_PreInput2 <- Data_QE_2$QE_table %>%
  filter(NoCq != TRUE) %>%
  filter(MaxCq != TRUE) %>%
  filter(NoPeak != TRUE) %>%
  filter(NotWithinCqRange != TRUE) %>%
  filter(PoorReplicate != TRUE) %>%
  filter(Cq < 33) %>%
  filter(Sample != "WT_gDNA") %>%
  filter(Peak == 1) 

Data_PreInput3 <- Data_QE_3$QE_table %>%
  filter(NoCq != TRUE) %>%
  filter(MaxCq != TRUE) %>% 
  filter(FalseCq != TRUE) %>%
  filter(Cq < 33) %>%
  filter(PoorReplicate != TRUE) %>%
  filter(Sample != "WT_gDNA") %>%
  filter(Peak == 1)

Data_PreInput <- bind_rows(Data_PreInput1, Data_PreInput2, Data_PreInput3)

# Generating DDCt
Input_Start <- bind_rows(OR10G4_7J10_Full, Data_PreInput) %>%
  ungroup() %>%
  group_by(Plate, Gene, Sample, Mutation) %>%
  summarize(MeanSampleCt = mean(Cq, na.rm = TRUE)) %>%
  mutate(Gene = factor(Gene, levels = c("Acsm4", "Acss2", "Slc25a35", 
                                        "Olfr358", "Olfr390", "Olfr510", "Olfr596_603", "Olfr609", "Olfr690", "Olfr1154")),
         Plate = factor(Plate, levels = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")),
         Genotype = case_when(Mutation == "Line 7" ~ "Line 7",
                              Mutation == "L7J_Het" ~ "Line 7Jax",
                              Mutation == "L10_Het" ~ "Line 10",
                              TRUE ~ "WT")) %>%
  ungroup() %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "Line 7", "Line 7Jax", "Line 10"))) %>%
  left_join(Primer_Efficiency, by = "Gene") %>%
  select(Plate, Gene, Sample, Genotype, MeanSampleCt, Primer_Efficiency) %>%
  mutate(MeanSCtxlog2Eff = MeanSampleCt*log2(Primer_Efficiency))

Output1 <- Input_Start %>%
  filter(Gene %in% c("Acsm4", "Acss2", "Slc25a35")) %>%
  select(Sample, Gene, Plate, MeanSCtxlog2Eff) %>%
  pivot_wider(names_from = Gene, values_from = MeanSCtxlog2Eff) %>%
  rename(MeanSCtxlog2Eff_Acsm4 = Acsm4, 
         MeanSCtxlog2Eff_Acss2 = Acss2,
         MeanSCtxlog2Eff_Slc25a35 = Slc25a35)

Output2 <- Input_Start %>% 
  left_join(Output1, by = c("Sample", "Plate")) %>%
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

#Calculating the DDCt for each genotype/gene pair
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

#06 - Comparing Lines----
OR10G4_Analysis_081822 <- Label_qPCR_data(CqFile = Input_list_C$`220818_Exp_Ct`, 
                                          Style = "Custom", 
                                          SampleNames = c(rep(c("MSS-4164", "MSS-4165", "MSS-4170", "PFH-610", "PFH-613", "PFH-617", "PFH-733", "MSS-4361", 
                                                                "MSS-4362", "MSS-4364", "MSS-4365", "OR10G4L7_gDNA"), 18), 
                                                          rep(c("OR10G4L10_gDNA", "UPWater"), each = 18)),
                                          GeneNames = c(rep(c("Acsm4", "Acss2", "Slc25a35", "OR10G4_CDS", "OR10G4_mRNA", "GCaMP"), each = 36),
                                                        rep(c("Acsm4", "Acss2", "Slc25a35", "OR10G4_CDS", "OR10G4_mRNA", "GCaMP"), each = 3, times = 2)),
                                          Mutation = c(rep(c("Line 7", "Line 7", "Line 7", "Line 7J", "Line 7J", "Line 7J", "Line 7J", "Line 10", "Line 10", 
                                                             "Line 10", "Line 10", "OR10G4L7_gDNA"), 18), 
                                                       rep(c("OR10G4L10_gDNA", "UPWater"), 18)),
                                          ReplicateNames = c(rep(c("One", "Two", "Three"), times = 6, each = 12), rep(c("One", "Two", "Three"), times = 12)))

OR10G4_Analysis_081822_QE <- Quality_Evaluation(OR10G4_Analysis_081822, MeltFile = Input_list_C$`220818_Exp_Melt`, RepDiffAllow = 0.6, CqLower = 15)

OR10G4_Analysis <- OR10G4_Analysis_081822_QE$QE_table %>% 
  filter(Sample %notin% c("OR10G4L7_gDNA","OR10G4L10_gDNA", "UPWater")) %>%
  filter(Gene != "OR10G4_mRNA") %>%
  filter(Peak == 1) %>%
  select(Gene, Sample, Mutation, Replicate, Cq) %>%
  group_by(Gene, Sample, Mutation) %>%
  summarize(MeanCq = mean(Cq, na.rm = TRUE)) %>%
  left_join(Primer_Efficiency, by = "Gene") %>%
  mutate(MeanSCtxlog2Eff = MeanCq*log2(Primer_Efficiency))

Compare1 <- OR10G4_Analysis %>%
  filter(Gene %in% c("Acsm4", "Acss2", "Slc25a35")) %>%
  select(Sample, Gene, MeanSCtxlog2Eff) %>%
  pivot_wider(names_from = Gene, values_from = MeanSCtxlog2Eff) %>%
  rename(MeanSCtxlog2Eff_Acsm4 = Acsm4, 
         MeanSCtxlog2Eff_Acss2 = Acss2,
         MeanSCtxlog2Eff_Slc25a35 = Slc25a35)

Compare2 <- OR10G4_Analysis %>% 
  left_join(Compare1, by = c("Sample")) %>%
  mutate(DCtSample = (MeanSCtxlog2Eff_Acsm4 + MeanSCtxlog2Eff_Acss2 + MeanSCtxlog2Eff_Slc25a35)/3 - MeanSCtxlog2Eff)

#Taking the mean of the DCtSample values from each plate to unify across plates
Output_Compare <- Compare2 %>%
  ungroup() %>% 
  group_by(Sample, Mutation, Gene) %>% 
  summarize(MeanDCtSample = mean(DCtSample))

#Calculating a DCt for each genotype/gene pair
Output_Compare_Geno <- Output_Compare %>%
  ungroup() %>% 
  group_by(Mutation, Gene) %>% 
  summarize(MeanDCtGenotype = mean(MeanDCtSample),
            DCtGenotypeSize = length(MeanDCtSample),
            DCtGenotypeSEM = sd(MeanDCtSample)/sqrt(DCtGenotypeSize))

#Calculating the DDCt for each genotype/gene pair
Output_Compare_Geno10 <- Output_Compare_Geno %>%
  filter(Mutation == "Line 10") %>%
  rename(DCt10Size = DCtGenotypeSize, 
         MeanDCt10 = MeanDCtGenotype, 
         DCt10SEM = DCtGenotypeSEM) %>%
  ungroup() %>%
  select(-Mutation)

Output_Compare_Geno_DDCt <- Output_Compare_Geno %>%
  left_join(Output_Compare_Geno10, by = "Gene") %>%
  mutate(MeanDDCt = MeanDCtGenotype - MeanDCt10, 
         DDCtSEM = sqrt(DCtGenotypeSEM^2 + DCt10SEM^2),
         DDCtCI = qt(0.975, df = DCtGenotypeSize + DCt10Size - 2) * DDCtSEM)

#DDCt for each sample
Output_Compare_DDCt <- Output_Compare %>%
  ungroup() %>%
  left_join(Output_Compare_Geno10, by = "Gene") %>%
  select(Sample, Gene, Mutation, MeanDCtSample, MeanDCt10) %>%
  mutate(SampleDDCt = MeanDCtSample - MeanDCt10) %>%
  left_join(Output_Compare_Geno_DDCt, by = c("Gene", "Mutation", "MeanDCt10"))

#06 - Plots----
#Olfr genes
DDCtSampleGenotype_Filtered1 <- Output_Sample_DDCt %>%
  left_join(Output_Genotype_DDCt, by = c("Gene", "Genotype")) %>%
  mutate(Gene = factor(Gene, levels = c("Acsm4", "Acss2", "Slc25a35", "Olfr358", "Olfr390", "Olfr510", "Olfr596_603", "Olfr609", "Olfr690", "Olfr1154"))) %>%
  filter(Gene %notin% c("Acsm4", "Acss2", "Slc25a35")) %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "Line 7", "Line 7Jax", "Line 10"))) %>%
  mutate(Class = case_when(Gene %in% c("Olfr596_603", "Olfr690") ~ "Class I", 
                           TRUE ~ "Class II")) %>%
  left_join(DVI, by = c("Gene" = "Symbol")) %>%
  mutate(DVI = ifelse(is.na(DVI), 1.05, DVI)) %>%
  mutate(Gene = factor(Gene, levels = c("Olfr596_603", "Olfr690", "Olfr510", "Olfr1154", "Olfr358", "Olfr390")))

ggplot(data = DDCtSampleGenotype_Filtered1) +
  geom_point(aes(Gene, SampleDDCt, color = Genotype, group = Genotype), position = position_dodge(width = 0.6), size = 5) +
  scale_color_manual(values=c("WT" = "black", "Line 7" = "skyblue1", "Line 7Jax" = "red2", "Line 10" = "green3")) +
  geom_errorbar(aes(x = Gene, group = Genotype, ymin = MeanDDCt - DDCtCI, ymax = MeanDDCt + DDCtCI), 
                width = 0.6, position = position_dodge(0.6), color = "black") +
  scale_y_continuous(breaks=seq(-6,1,1)) +
  labs(title = "Olfr mRNA L2FC in Three 10G4 lines", 
       y = expression("Log2 Fold Change (-" ~ Delta*Delta ~ "Ct)") ) +
  theme_bw() +
  theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(face="bold", size = 18),
        axis.title.y = element_text(size = 20),
        legend.position = c(0.50, 0.25),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(face="bold", size = 12),
        legend.title = element_blank()) +
  geom_vline(xintercept = 2.5, linetype="dotted", color = "grey40", size=1.5) +
  annotate("label", label = "Class I", x = 1.5, y = 1.25, size = 8, colour = "blue", fontface = 2) +
  annotate("label", label = "Class II", x = 4.5, y = 1.25, size = 8, colour = "red", fontface = 2) +
  annotate("label", label = "DVI: 1.05", x = 1, y = -6.20, size = 5, colour = "black") +
  annotate("label", label = "DVI: 1.05", x = 2, y = -6.20, size = 5, colour = "black") +
  annotate("label", label = "DVI: 1.05", x = 3, y = -6.20, size = 5, colour = "black") +
  annotate("label", label = "DVI: 1.05", x = 4, y = -6.20, size = 5, colour = "black") +
  annotate("label", label = "DVI: 1.5", x = 5, y = -6.20, size = 5, colour = "black") +
  annotate("label", label = "DVI: 3.6", x = 6, y = -6.20, size = 5, colour = "black") 

ggsave("Output/RTqPCR_SF1G.png", width = 30, height = 15, units = "cm")

#10G4 transgene
compare_means(SampleDDCt ~ Mutation, data = subset(Output_Compare_DDCt, subset = Gene %in% c("GCaMP", "OR10G4_CDS")), group.by = c("Gene"), method = "t.test")

# caption = "DDCt calculated using corrected primer efficiencies and reference primers Acsm4, Acss2, and Slc25a35. 
#        Each Sample and Genotype DDCt calibrated using the mean Line 10 DCt. Ct and DCt in triplicate. 
#        Plotted L2FC of each Sample with the 95% confidence interval based on the mean and single-propagated error of each Genotype
#        t.test p values not corrected for multiple comparisons.",

Output_Compare_DDCt %>%
  filter(Gene %in% c("GCaMP", "OR10G4_CDS")) %>%
  mutate(Mutation = ifelse(Mutation == "Line 7J", "Line 7Jax", Mutation)) %>%
  mutate(Mutation = factor(Mutation, levels = c("WT", "Line 7", "Line 7Jax", "Line 10"))) %>%
  ggplot(aes(Gene, SampleDDCt, color = Mutation, group = Mutation)) +
  geom_errorbar(aes(ymin = MeanDDCt - DDCtCI, ymax = MeanDDCt + DDCtCI), width = 0.4, 
                position = position_dodge(0.6), size = 1.5, stat = "unique", color = "black") +
  geom_point(position = position_dodge(width = 0.6), size = 7, alpha = 0.95) + 
  labs(title = "10G4-IRES-GCaMP mRNA L2FC in Three 10G4 lines", 
       y = expression("Log2 Fold Change (-" ~ Delta*Delta ~ "Ct)")) +
  scale_color_manual(values=c("WT" = "black", "Line 7" = "skyblue1", "Line 7Jax" = "red2", "Line 10" = "green3")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 22),
        legend.position = c(0.50, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(face="bold", size = 12),
        plot.caption = element_text(size = 8, face = "italic"), 
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), 
        axis.text = element_text(face="bold", size = 20),
        legend.title = element_blank()) +
  geom_signif(
    annotation = 0.033, textsize = 8,
    y_position = 2.2, xmin = 0.8, xmax = 1.2,
    tip_length = c(0.02, 0.02), color = "black") +
  geom_signif(
    annotation = 0.012, textsize = 8,
    y_position = 2, xmin = 1, xmax = 1.2,
    tip_length = c(0.02, 0.02), color = "black") +
  geom_signif(
    annotation = 0.020, textsize = 8,
    y_position = 1.45, xmin = 1.8, xmax = 2.2,
    tip_length = c(0.02, 0.02), color = "black") +
  geom_signif(
    annotation = 0.0060, textsize = 8,
    y_position = 1.25, xmin = 2, xmax = 2.2,
    tip_length = c(0.02, 0.02), color = "black")

ggsave("Output/RTqPCR_SF1F.png", width = 25, height = 25, units = "cm")

#07 - Tables----
#Significance Comparison
Table10G4 <- compare_means(SampleDDCt ~ Genotype, Output_Sample_DDCt, method = "t.test", group.by = c("Gene")) %>% 
  filter(Gene %notin% c("Acsm4", "Acss2", "Slc25a35")) %>%
  select(Gene, group1, group2, p) %>% 
  mutate(row = row_number(),
         HolmAdjustment = 0.05/(36 - row + 1), 
         PassHA = p < HolmAdjustment, 
         HolmCorrection = p * (36 - row + 1),
         PassHC = HolmCorrection < 0.05, 
         Pass_DunnSidakCorrection = p < 1 - (1 - 0.05)^(1/36)) %>% 
  select(Gene, group1, group2, p, HolmCorrection, PassHC, Pass_DunnSidakCorrection) %>% 
  rename(p_value = p, Pass_HolmCorrection = PassHC, Group1 = group1, Group2 = group2) %>%
  arrange(Gene) %>%
  mutate(p.signif = case_when(HolmCorrection <= 0.0001 ~ "****",
                              HolmCorrection <= 0.001 ~ "***",
                              HolmCorrection <= 0.01 ~ "**",
                              HolmCorrection <= 0.05 ~ "*",
                              TRUE ~ "NS"),
         Comparison = paste(Group1, "vs", Group2)) %>%
  select(Gene, Comparison, p.signif) %>%
  pivot_wider(names_from = Gene, values_from = p.signif) %>%
  select(Comparison, Olfr596_603, Olfr690, Olfr510, Olfr1154, Olfr358, Olfr390) %>%
  tableGrob(rows = NULL)

Text10G4 <- text_grob(paste("Thirty-six independent t.tests were performed.", 
                            "With alpha = 0.05, Holm-Bonferroni corrections were calculated. ", 
                            "**** = p <= 0.0001, *** = p <= 0.001, ** = p <= 0.01, * = p <= 0.05, NS = Not Significant", sep = "\n"), face = "italic")

grid.arrange(Table10G4, Text10G4, heights = c(0.5, 0.1), as.table = TRUE)

g <- arrangeGrob(Table10G4, Text10G4, nrow=2, heights = c(1, 0.6))

g <- arrangeGrob(Table10G4)

ggsave("Output/Supp_Table1_Significant_SF1H.png",plot = g, width = 18, height = 6, units = "cm")

#Comparing RTqPCR vs DESeq
Comparison <- Output_Genotype_DDCt %>%
  select(Genotype, Gene, MeanDDCt) %>%
  filter(Genotype != "WT" & Gene %notin% c("Acsm4", "Acss2", "Slc25a35")) %>%
  pivot_wider(names_from = Genotype, values_from = MeanDDCt) %>%
  separate_rows(Gene, sep = "_") %>%
  mutate(Gene = ifelse(Gene == "603", "Olfr603", Gene)) %>%
  left_join(DE, by = c("Gene" = "Symbol")) %>%
  rename(`Line 7 RNASeq` = log2FoldChange) %>%
  select(Gene, `Line 7 RNASeq`, everything()) %>% 
  mutate(across(where(is.numeric), round, 2)) %>%
  tableGrob(rows = NULL)

g <- arrangeGrob(Comparison)

ggsave("Output/Supp_Table1_Comparison_ST1.png", plot = g, width = 13, height = 6.5, units = "cm")
