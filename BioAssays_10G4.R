# This is a script designed to import all relevant BioAssay data and produce analyses and graphs
# Read through this and make the appropriate changes in file paths AND if you add new experiments, edit the relevant code for them.----

# Required packages
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

# Import and name all assay files using a function from: https://stackoverflow.com/questions/11433432/how-to-import-multiple-csv-files-at-once
# The Input files are in the Input folder located relative to the set directory. Adjust all file paths as necessary.
read_plus <- function(flnm) {
  read_csv(flnm, col_names = FALSE)
  }

Input_list <-
  list.files(path = "../Input/", #Adjust here
             pattern = "*.csv", 
             full.names = TRUE) %>% 
  lapply(read_plus)

Input_list_names <-   
  list.files(path = "../Input/", #Adjust here
             pattern = "*.csv") %>%
  str_replace(pattern = ".csv", 
              replacement = "")

names(Input_list) <- Input_list_names

# Import MyAssays curve data. In this file path, the csv is located at the same level as the parent folder. Adjust as needed. 
myAssayOutput <- read_csv("../MyAssay_outputs_simple.csv")

# A rather convoluted piece of code that produces a first table of data for each experiment, unlabeled
Input_tibble <- tibble(Assay = Input_list, Exp = Input_list_names) %>%
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
         Combined = pmap(list(Measure, Position, Concentration, Unbound, Concentration2), .f = function (Measure, Position, Concentration, Unbound, Concentration2) {
                                                                                                        tibble(Measure = Measure, Position = Position, Concentration = Concentration, 
                                                                                                               Unbound = Unbound, Concentration2 = Concentration2)}))
# LABELS: A list containing all label information about each experiment; EDIT the list when adding new results----
Label_list <- list(tibble(Exp = "Exp38a", 
                          Sample = rep(c("PFH878", "PFH879","PFH880","PFH886","PFH887","PFH888", "MSS4512", "MSS4577", "MSS4591", "MSS4592"), 18),
                          Mutation = rep(c("OR10G4L7J", "OR10G4L7J", "OR10G4L7J", "OR10G4L7J", "WT", "WT", "OR1A1V1.1", "WT", "OR1A1V1.1", "OR1A1V1.1"), 18),
                          Ligand = c(rep(c("DMSO", "FSK"), each = 10, times = 3), rep(c("Guaiacol_TCI", "Vanillin_TCI"), each = 10, times = 3), rep(c("MND_SC_50uM", "MND_SC_50nM"), each = 10, times = 3))),
                   tibble(Exp = "Exp38b", 
                          Sample = rep(c("PFH878", "PFH879","PFH880","PFH886","PFH887","PFH888", "MSS4512", "MSS4577", "MSS4591", "MSS4592"), 18),
                          Mutation = rep(c("OR10G4L7J", "OR10G4L7J", "OR10G4L7J", "OR10G4L7J", "WT", "WT", "OR1A1V1.1", "WT", "OR1A1V1.1", "OR1A1V1.1"), 18),
                          Ligand = c(rep(c("DMSO", "FSK"), each = 10, times = 3), rep(c("Guaiacol_TCI", "Vanillin_TCI"), each = 10, times = 3), rep(c("MND_SC_50uM", "MND_SC_50nM"), each = 10, times = 3))),
                   tibble(Exp = "Exp39a", 
                          Sample = rep(c("PFH887", "PFH888","MSS4763(1:100)","MSS4763(1:200)","PFH878","PFH879", "PFH880", "PFH886","MSS4764(1:100)", "MSS4764(1:200)", "MSS4765(1:100)", "MSS4765(1:200)", "MSS4766(1:100)", "MSS4766(1:200)"), 9),
                          Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J", "OR10G4L7J", "OR10G4L7J", "OR10G4L7J", "OR10G4V2.1", "OR10G4V2.1", "OR10G4V2.1", "OR10G4V2.1", "OR10G4V2.1", "OR10G4V2.1"), 9),
                          Ligand = c(rep(c("DMSO", "FSK", "Guaiacol_TCI"), each = 42))),
                   tibble(Exp = "Exp39b", 
                          Sample = rep(c("PFH887", "PFH888","MSS4763(1:100)","MSS4763(1:200)","PFH878","PFH879", "PFH880", "PFH886","MSS4764(1:100)", "MSS4764(1:200)", "MSS4765(1:100)", "MSS4765(1:200)", "MSS4766(1:100)", "MSS4766(1:200)"), 9),
                          Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J", "OR10G4L7J", "OR10G4L7J", "OR10G4L7J", "OR10G4V2.1", "OR10G4V2.1", "OR10G4V2.1", "OR10G4V2.1", "OR10G4V2.1", "OR10G4V2.1"), 9),
                          Ligand = c(rep(c("DMSO", "FSK", "Guaiacol_TCI"), each = 42))),
                   tibble(Exp = "Exp39c", 
                          Sample = rep(c("PFH887", "PFH888","MSS4763(1:100)","MSS4763(1:200)","PFH878","PFH879", "PFH880", "PFH886","MSS4764(1:100)", "MSS4764(1:200)", "MSS4765(1:100)", "MSS4765(1:200)", "MSS4766(1:100)", "MSS4766(1:200)"), 9),
                          Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J", "OR10G4L7J", "OR10G4L7J", "OR10G4L7J", "OR10G4V2.1", "OR10G4V2.1", "OR10G4V2.1", "OR10G4V2.1", "OR10G4V2.1", "OR10G4V2.1"), 9),
                          Ligand = c(rep(c("DMSO", "FSK", "Guaiacol_TCI"), each = 42))),
                   tibble(Exp = "Exp41a_a", 
                          Sample = rep(c("PFH1721", "PFH1722", "PFH1741","PFH1742","PFH1743", "PFH1547", "PFH1548", "PFH1550", "PFH1740", "PFH1744", "PFH1745", "PFH1549", "PFH1551", "PFH1553"), 12),
                          Mutation = rep(c("WT", "WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het","OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Hom", "OR10G4L7J_Hom", "OR10G4L7J_Hom"), 12),
                          Ligand = c(rep(c("DMSO", "FSK", "Guaiacol_TCI", "Vanillin_TCI"), each = 42))),
                   tibble(Exp = "Exp41a_b", 
                          Sample = rep(c("PFH1721", "PFH1722", "PFH1741","PFH1742","PFH1743", "PFH1547", "PFH1548", "PFH1550", "PFH1740", "PFH1744", "PFH1745", "PFH1549", "PFH1551", "PFH1553"), 12),
                          Mutation = rep(c("WT", "WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het","OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Hom", "OR10G4L7J_Hom", "OR10G4L7J_Hom"), 12),
                          Ligand = c(rep(c("DMSO", "FSK", "Guaiacol_TCI", "Vanillin_TCI"), each = 42))),
                   tibble(Exp = "Exp41a_c", 
                          Sample = rep(c("PFH1721", "PFH1722", "PFH1741","PFH1742","PFH1743", "PFH1547", "PFH1548", "PFH1550", "PFH1740", "PFH1744", "PFH1745", "PFH1549", "PFH1551", "PFH1553"), 12),
                          Mutation = rep(c("WT", "WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het","OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Hom", "OR10G4L7J_Hom", "OR10G4L7J_Hom"), 12),
                          Ligand = c(rep(c("DMSO", "FSK", "Guaiacol_TCI", "Vanillin_TCI"), each = 42))),
                   tibble(Exp = "Exp41a_d", 
                          Sample = rep(c("PFH1721", "PFH1722", "PFH1741","PFH1742","PFH1743", "PFH1547", "PFH1548", "PFH1550", "PFH1740", "PFH1744", "PFH1745", "PFH1549", "PFH1551", "PFH1553"), 12),
                          Mutation = rep(c("WT", "WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het","OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Hom", "OR10G4L7J_Hom", "OR10G4L7J_Hom"), 12),
                          Ligand = c(rep(c("DMSO", "FSK", "Guaiacol_TCI", "Vanillin_TCI"), each = 42))),
                   tibble(Exp = "Exp41b_a", 
                          Sample = rep(c("PFH1721", "PFH1722", "PFH1741","PFH1742","PFH1743", "PFH1547", "PFH1548", "PFH1550", "PFH1740", "PFH1744", "PFH1745", "PFH1549", "PFH1551", "PFH1553"), 12),
                          Mutation = rep(c("WT", "WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het","OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Hom", "OR10G4L7J_Hom", "OR10G4L7J_Hom"), 12),
                          Ligand = c(rep(c("DMSO", "FSK", "Salicylicaldehyde_TCI", "O-Cresol_TCI"), each = 42))),
                   tibble(Exp = "Exp41b_b", 
                          Sample = rep(c("PFH1721", "PFH1722", "PFH1741","PFH1742","PFH1743", "PFH1547", "PFH1548", "PFH1550", "PFH1740", "PFH1744", "PFH1745", "PFH1549", "PFH1551", "PFH1553"), 12),
                          Mutation = rep(c("WT", "WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het","OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Hom", "OR10G4L7J_Hom", "OR10G4L7J_Hom"), 12),
                          Ligand = c(rep(c("DMSO", "FSK", "Salicylicaldehyde_TCI", "O-Cresol_TCI"), each = 42))),
                   tibble(Exp = "Exp41b_c", 
                          Sample = rep(c("PFH1721", "PFH1722", "PFH1741","PFH1742","PFH1743", "PFH1547", "PFH1548", "PFH1550", "PFH1740", "PFH1744", "PFH1745", "PFH1549", "PFH1551", "PFH1553"), 12),
                          Mutation = rep(c("WT", "WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het","OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Hom", "OR10G4L7J_Hom", "OR10G4L7J_Hom"), 12),
                          Ligand = c(rep(c("DMSO", "FSK", "Salicylicaldehyde_TCI", "O-Cresol_TCI"), each = 42))),
                   tibble(Exp = "Exp41b_d", 
                          Sample = rep(c("PFH1721", "PFH1722", "PFH1741","PFH1742","PFH1743", "PFH1547", "PFH1548", "PFH1550", "PFH1740", "PFH1744", "PFH1745", "PFH1549", "PFH1551", "PFH1553"), 12),
                          Mutation = rep(c("WT", "WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het","OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Hom", "OR10G4L7J_Hom", "OR10G4L7J_Hom"), 12),
                          Ligand = c(rep(c("DMSO", "FSK", "Salicylicaldehyde_TCI", "O-Cresol_TCI"), each = 42))),
                   tibble(Exp = "Exp42a",
                          Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 42),
                          Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 42),
                          Ligand = rep(c("DMSO", "FSK", "Guaiacol", "Vanillin", "Acetovanillone", "Vanillyl butyl ether", "Ethyl vanillin", "Dehydrodivanillin",
                                         "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde", "o-Cresol", "2-Ethylphenol"), each = 24)),
                  tibble(Exp = "Exp42b",
                         Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 42),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 42),
                         Ligand = rep(c("DMSO", "FSK", "Guaiacol", "Vanillin", "Acetovanillone", "Vanillyl butyl ether", "Ethyl vanillin", "Dehydrodivanillin",
                                        "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde", "o-Cresol", "2-Ethylphenol"), each = 24)),
                  tibble(Exp = "Exp42c",
                         Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 42),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 42),
                         Ligand = rep(c("DMSO", "FSK", "Guaiacol", "Vanillin", "Acetovanillone", "Vanillyl butyl ether", "Ethyl vanillin", "Dehydrodivanillin",
                                        "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde", "o-Cresol", "2-Ethylphenol"), each = 24)),
                  tibble(Exp = "Exp43a",
                         Sample = rep(c("PFH1721", "PFH1741", "PFH1742", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 45),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                         Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Guaiacol",168), rep("o-Cresol", 144)),
                         Dilution = rep(c(0, 5e-7, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9), each = 24)),
                  tibble(Exp = "Exp43b",
                         Sample = rep(c("PFH1721", "PFH1741", "PFH1742", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 45),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                         Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Guaiacol",168), rep("o-Cresol", 144)),
                         Dilution = rep(c(0, 5e-7, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9), each = 24)),
                  tibble(Exp = "Exp43c",
                         Sample = rep(c("PFH1721", "PFH1741", "PFH1742", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 45),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                         Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Guaiacol",168), rep("o-Cresol", 144)),
                         Dilution = rep(c(0, 5e-7, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9), each = 24)),
                  tibble(Exp = "Exp45a",
                         Sample = rep(c("PFH1721", "PFH1741", "PFH1742", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 45),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                         Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Guaiacol",168), rep("o-Cresol", 144)),
                         Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), each = 24)),
                  tibble(Exp = "Exp46a",
                         Sample = rep(c("PFH1721", "PFH1741", "PFH1742", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 42),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 42),
                         Ligand = rep(c("DMSO", "FSK", "Guaiacol", "Vanillin", "Acetovanillone", "Vanillyl butyl ether", "Ethyl vanillin", "Dehydrodivanillin",
                                        "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde", "o-Cresol", "2-Ethylphenol"), each = 24)),
                  tibble(Exp = "Exp47a",
                         Sample = rep(c("PFH1721", "PFH1741", "PFH1742", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 45),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 45),
                         Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Guaiacol",168), rep("o-Cresol", 144)),
                         Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), each = 24)),
                  tibble(Exp = "Exp47b",
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
                  tibble(Exp = "Exp49a",
                         Sample = rep(c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20"), 9),
                         WT_cilia = rep(c(1/50, 1/100, 1/200, 1/400, 1/800, 1/50, 1/100, 1/200, 1/400, 1/800, 0, 0, 0, 0, 0, 1/200, 1/200, 1/200, 1/200, 1/200), 9),
                         OR10G4_cilia = rep(c(0, 0, 0, 0, 0, 1/200, 1/200, 1/200, 1/200, 1/200, 1/50, 1/100, 1/200, 1/400, 1/800, 1/50, 1/100, 1/200, 1/400, 1/800), 9),
                         Ligand = c(rep(c("DMSO", "FSK", "Guaiacol"), each = 60))),
                  tibble(Exp = "Exp49b",
                         Sample = rep(c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20"), 9),
                         WT_cilia = rep(c(1/50, 1/100, 1/200, 1/400, 1/800, 1/50, 1/100, 1/200, 1/400, 1/800, 0, 0, 0, 0, 0, 1/200, 1/200, 1/200, 1/200, 1/200), 9),
                         OR10G4_cilia = rep(c(0, 0, 0, 0, 0, 1/200, 1/200, 1/200, 1/200, 1/200, 1/50, 1/100, 1/200, 1/400, 1/800, 1/50, 1/100, 1/200, 1/400, 1/800), 9),
                         Ligand = c(rep(c("DMSO", "FSK", "Guaiacol"), each = 60))),
                  tibble(Exp = "Exp50a",
                         Sample = rep(c("PFH1791", "PFH1792", "PFH1793", "PFH1789", "PFH1790"), 9),
                         Mutation = rep(c("OR10G4L7J Hom", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "WT"), 9),
                         Ligand = c(rep(c("DMSO", "FSK", "Guaiacol"), each = 15))),
                  tibble(Exp = "Exp51DSa",
                         Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 42),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 42),
                         Ligand = c(rep(c("DMSO", "FSK"), each = 24), rep("Vanillin",144), rep("Ethyl vanillin", 144)),
                         Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8), each = 24)),
                  tibble(Exp = "Exp51Newa",
                         Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 36),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 36),
                         Ligand = rep(c("DMSO", "FSK", "4-Fluoroguaiacol", "5-Fluoroguaiacol", "4-Bromoguaiacol", "5-Bromoguaiacol", "4-Chloroguaiacol", "4-Iodoguaiacol",
                                        "4-Methoxyguaiacol", "4-Methylguaiacol", "2-Ethoxyphenol", "Guaiacol"), each = 24)),
                  tibble(Exp = "Exp51Olda",
                         Sample = rep(c("PFH1721", "PFH1742", "PFH1741", "PFH1743", "PFH1548", "PFH1740", "PFH1744", "PFH1745"), 42),
                         Mutation = rep(c("WT", "WT", "WT", "WT", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het", "OR10G4L7J_Het"), 42),
                         Ligand = rep(c("DMSO", "FSK", "Guaiacol", "Vanillin", "Acetovanillone", "Vanillyl butyl ether", "Ethyl vanillin", "Dehydrodivanillin",
                                        "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde", "o-Cresol", "2-Ethylphenol"), each = 24)),
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
                                        "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde", "o-Cresol", "2-Ethylphenol"), each = 24)),
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
                         Dilution = rep(c(0, 5e-7, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9), each = 24)),
                  tibble(Exp = "Exp57a",
                         Sample = rep(c("PFH0308", "PFH0312", "PFH0313", "PFH0197", "PFH0201", "PFH0206", "PFH0216", "PFH0224", "PFH0544", "PFH1321", "PFH0215", "PFH0169"), 30),
                         Mutation = rep(c("OR2S2B", "OR2S2B", "OR2S2B", "OR2S2", "OR2S2", "OR2S2", "Olfr155", "Olfr155", "WT", "WT", "WT", "WT"), 30),
                         Ligand = rep(c("DMSO", "FSK", "Pentenone", "Hexyl Acetate", "Octanal", "Hippuric Acid", "5-Ethyl-2,3-dimethylpyrazine", "Octamethyleneimine",
                                        "2-Isopropyl-3-methoxypyrazine", "2-Isobutyl-3-methoxypyrazine"), each = 36)),
                  tibble(Exp = "Exp57b",
                         Sample = rep(c("PFH1322", "PFH1821", "PFH1823", "PFH0543", "PFH0545", "PFH0546", "PFH0167", "PFH0168", "PFH0172", "PFH0544", "PFH1321", "PFH0215"), 30),
                         Mutation = rep(c("Olfr157", "Olfr157", "Olfr157", "Olr292", "Olr292", "Olr292", "Olr841", "Olr841", "Olr841", "WT", "WT", "WT"), 30),
                         Ligand = rep(c("DMSO", "FSK", "Pentenone", "Hexyl Acetate", "Octanal", "Hippuric Acid", "5-Ethyl-2,3-dimethylpyrazine", "Octamethyleneimine",
                                        "2-Isopropyl-3-methoxypyrazine", "2-Isobutyl-3-methoxypyrazine"), each = 36)),
                  tibble(Exp = "Exp58a",
                         Sample = rep(c("PFH0308", "PFH0312", "PFH0313", "PFH0197", "PFH0201", "PFH0206", "PFH0216", "PFH0224", "PFH0544", "PFH1321", "PFH0215", "PFH0169"), 30),
                         Mutation = rep(c("OR2S2B", "OR2S2B", "OR2S2B", "OR2S2", "OR2S2", "OR2S2", "Olfr155", "Olfr155", "WT", "WT", "WT", "WT"), 30),
                         Ligand = rep(c("DMSO", "FSK", "Pentenone", "Hexyl Acetate", "Octanal", "Hippuric Acid", "5-Ethyl-2,3-dimethylpyrazine", "Octamethyleneimine",
                                        "2-Isopropyl-3-methoxypyrazine", "2-Isobutyl-3-methoxypyrazine"), each = 36)),
                  tibble(Exp = "Exp58b",
                         Sample = rep(c("PFH1322", "PFH1821", "PFH1823", "PFH0543", "PFH0545", "PFH0546", "PFH0167", "PFH0168", "PFH0172", "PFH0544", "PFH1321", "PFH0215"), 30),
                         Mutation = rep(c("Olfr157", "Olfr157", "Olfr157", "Olr292", "Olr292", "Olr292", "Olr841", "Olr841", "Olr841", "WT", "WT", "WT"), 30),
                         Ligand = rep(c("DMSO", "FSK", "Pentenone", "Hexyl Acetate", "Octanal", "Hippuric Acid", "5-Ethyl-2,3-dimethylpyrazine", "Octamethyleneimine",
                                        "2-Isopropyl-3-methoxypyrazine", "2-Isobutyl-3-methoxypyrazine"), each = 36)))

# Continued processing of the input tables and splitting based on the type of experiment ----
# The Label_tibble object is based on the list found in the LABELS section
Label_tibble <- tibble(Labels = Label_list) %>%
  mutate(Exp = map_chr(Labels, ~unique(.x$Exp)))

# Completely labeled tibbles, prior to QC and further cleaning/processing
Label_Input_Tibble <- Input_tibble %>%
  inner_join(Label_tibble, by = "Exp") %>% #Only labeled assays remain due to using inner_join; some Experiments were dropped from Input_tibble and Label_tibble
  mutate(Full_Table = map2(Combined, Labels, bind_cols))

Full_Tibbles <- Label_Input_Tibble %>%
  select(Exp, Full_Table, a, b, c, d)

# Filter "boxes" to work on specific types of BioAssays; add or remove experiments as needed
Dose_Response_10G4 <- c("Exp43b", "Exp45a", "Exp47a", "Exp48_1a", "Exp48_2a", "Exp51DSa", "Exp52DS1a", "Exp52DS2a", "Exp53DRa", "Exp53DRb", 
                        "Exp55a", "Exp55b", "Exp55c", "Exp56a", "Exp56b", "Exp56c")
Static_50uM_10G4 <- c("Exp42b", "Exp46a", "Exp51Newa", "Exp54NewNew", "Exp54NewOld", "Exp54OldNew", "Exp54OldOld")
Static_50uM_PD <- c("Exp57a", "Exp57b", "Exp58a", "Exp58b")

# Splitting the full tibble into specific groups and generating some QC data to analyze
Exp_DR_10G4 <- Full_Tibbles %>%
  filter(Exp %in% Dose_Response_10G4) %>%
  mutate(QC = map(Full_Table, function(.x) {.x %>% group_by(Sample, Ligand, Dilution) %>%
                                                   nest() %>%
                                                   mutate(MeanC = map(data, function(x) mean(x$Concentration2)),
                                                          SD = map(data, function(x) sd(x$Concentration2))) %>%
                                                   unnest(c(data, MeanC, SD))  %>%
                                                   ungroup() %>%
                                                   mutate(CV = SD/MeanC * 100,
                                                          SEM = SD/sqrt(3))} ) )

Exp_Static50uM_10G4 <- Full_Tibbles %>%
  filter(Exp %in% Static_50uM_10G4) %>%
  mutate(QC = map(Full_Table, function(.x) {.x %>% group_by(Sample, Ligand) %>%
                                                    nest() %>%
                                                    mutate(MeanC = map(data, function(x) mean(x$Concentration2)),
                                                           SD = map(data, function(x) sd(x$Concentration2))) %>%
                                                    unnest(c(data, MeanC, SD))  %>%
                                                    ungroup() %>%
                                                    mutate(CV = SD/MeanC * 100,
                                                           SEM = SD/sqrt(3))} ) )
  
Exp_Static50uM_PD <- Full_Tibbles %>%
  filter(Exp %in% Static_50uM_PD) %>%
  mutate(QC = map(Full_Table, function(.x) {.x %>% group_by(Sample, Ligand) %>%
                                                    nest() %>%
                                                    mutate(MeanC = map(data, function(x) mean(x$Concentration2)),
                                                           SD = map(data, function(x) sd(x$Concentration2))) %>%
                                                    unnest(c(data, MeanC, SD))  %>%
                                                    ungroup() %>%
                                                    mutate(CV = SD/MeanC * 100,
                                                           SEM = SD/sqrt(3))} ) )

# These experiments in the DR_10G4 set require some review and clean up prior to analysis and plots.----
map2(Exp_DR_10G4$Exp, Exp_DR_10G4$QC, ~ if(TRUE %in% .y$Unbound | max(.y$CV) > 30) {.x} ) %>% unlist()

# This negate function was taken from somewhere online, useful for filtering away multiple specific rows 
`%notin%` <- Negate(`%in%`)

# The cleanup requires manual investigation and selection; there is no "correct" choice in all cases. 
Exp_DR_10G4$QC[["Exp43b"]] <- Exp_DR_10G4$QC[["Exp43b"]] %>%
  filter(Position != 170)
Exp_DR_10G4$QC[["Exp47a"]] <- Exp_DR_10G4$QC[["Exp47a"]] %>%
  filter(Position != 25)
Exp_DR_10G4$QC[["Exp48_1a"]] <- Exp_DR_10G4$QC[["Exp48_1a"]] %>%
  filter(Position != 27)
Exp_DR_10G4$QC[["Exp51DSa"]] <- Exp_DR_10G4$QC[["Exp51DSa"]] %>%
  filter(Position %notin% c(74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120))
Exp_DR_10G4$QC[["Exp53DRb"]] <- Exp_DR_10G4$QC[["Exp53DRb"]] %>%
  filter(Position != 22)
Exp_DR_10G4$QC[["Exp55a"]] <- Exp_DR_10G4$QC[["Exp55a"]] %>%
  filter(Sample != "PFH1793")

# This next step takes the DR files to a final form. The DR files will be exported for use with prism to extrapolate the DR curves
Exp_DR_10G4_ForExport <- Exp_DR_10G4 %>%
  mutate(Concentration_Triplicate_mean = pmap(list(QC, a, b, c, d, Exp), function(QC, a, b, c, d, Exp) {QC %>% group_by(Sample, Mutation, Ligand, Dilution) %>%
                                                                                           summarize(S_Measure = mean(Measure)) %>%
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

Exp_DR_10G4_ForExport_Final <- rbindlist(Exp_DR_10G4_ForExport[[9]], use.names=TRUE) %>%
  mutate(Mutation = ifelse(Mutation == "OR10G4L7J Het", "OR10G4L7J_Het", Mutation)) #Fixing a slight error in the Mutation names, might already be corrected

# A new step to lower the baseline signal to approximately zero
Zero_Baseline <- Exp_DR_10G4_ForExport_Final %>%
  group_by(Sample, Mutation, Ligand, Dilution) %>%
  summarize(meanCCRAviaDMSO = mean(CCRAviaDMSO)) %>% #Taking the mean of Samples that fit into the same group_by grouping so as to avoid over-representating a cilia sample
  ungroup() %>% group_by(Mutation, Ligand, Dilution) %>% 
  summarise(meanCC = mean(meanCCRAviaDMSO)) %>% 
  ungroup() %>% 
  group_by(Mutation, Ligand) %>% 
  summarise(minmeanCC = min(meanCC)) %>%
  ungroup()

# The table exported and used in prism to estimate a 4PL/DR curve for each ligand
Exp_DR_10G4_ForExport_Final2 <- Exp_DR_10G4_ForExport_Final %>%
  group_by(Sample, Mutation, Ligand, Dilution) %>%
  summarize(meanCCRAviaDMSO = mean(CCRAviaDMSO)) %>% #Taking the mean of Samples that fit into the same group_by grouping so as to avoid over-representating a cilia sample
  ungroup() %>%
  left_join(Zero_Baseline, by = c("Mutation", "Ligand")) %>%
  mutate(meanBaselineCCRAviaDMSO = meanCCRAviaDMSO - minmeanCC) %>%
  select(Ligand, Dilution, Mutation, Sample, meanBaselineCCRAviaDMSO) %>%
  arrange(Ligand, Mutation) %>%
  pivot_wider(names_from = c(Mutation, Sample), values_from = meanBaselineCCRAviaDMSO) %>%
  mutate(Dilution = log10(Dilution)) %>% #Transform the dilutions into Log10 scale
  arrange(Ligand, desc(Dilution))

# A quick visual of the data
Exp_DR_10G4_ForExport_Final %>% 
  ggplot(aes(log10(Dilution), CCRAviaDMSO, color = Mutation)) + 
  geom_point() + 
  facet_wrap(~Ligand)

Exp_DR_10G4_ForExport_Final %>%
  group_by(Sample, Mutation, Ligand, Dilution) %>%
  summarize(meanCCRAviaDMSO = mean(CCRAviaDMSO)) %>% #Taking the mean of Samples that fit into the same group_by grouping so as to avoid over-representating a cilia sample
  ungroup() %>%
  left_join(Zero_Baseline, by = c("Mutation", "Ligand")) %>%
  mutate(meanBaselineCCRAviaDMSO = meanCCRAviaDMSO - minmeanCC) %>%
  ggplot(aes(log10(Dilution), meanBaselineCCRAviaDMSO, color = Mutation)) + 
  geom_point() + 
  facet_wrap(~Ligand)

# Export the Dose Response data to be manually analyzed using the non-linear regression functions of Prism software
write_csv(Exp_DR_10G4_ForExport_Final, "Dose_Respose_10G4.csv")
write_csv(Exp_DR_10G4_ForExport_Final2, "Dose_Respose_10G4_V3.csv") #V2 ignored

# After finding out the various values for the Dose-Response curves, import that data

DR_10G4_curves <- read_csv("DR_10G4_4Pcurve_prism.csv")

# A quick table of those prism results
DR_Grob <- DR_10G4_curves %>%
  mutate(Ligand = factor(Ligand, levels = c("DMSO", "FSK", "Guaiacol",  "Ethyl vanillin", "Vanillin", "o-Cresol", "2-Ethylphenol", "Dehydrodivanillin", 
                                            "Acetovanillone", "Vanillyl butyl ether", "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde",
                                            "5-Fluoroguaiacol", "4-Fluoroguaiacol", "5-Bromoguaiacol", "4-Bromoguaiacol", "4-Chloroguaiacol", "4-Iodoguaiacol",
                                            "5-Methoxyguaiacol", "4-Methoxyguaiacol", "4-Methylguaiacol", "2-Ethoxyphenol"))) %>%
  arrange(Ligand) %>%
  tableGrob(rows=NULL)

grid.arrange(DR_Grob)

# Revert to a different plan of analysis because I don't want to deal with prism for graphs and statistical tests

# Generate the prism-based DR curves to plot on top of experimental data
DR_curves_output <- DR_10G4_curves %>% 
  mutate(Data = pmap(list(Top, Bottom, HillSlope, EC50), function(Top, Bottom, HillSlope, EC50) {
           Input <- 10^seq(-10, 0, length.out = 100)
           Output <- Bottom + (Top - Bottom)/(1 + (Input/EC50)^(-HillSlope))
           tibble(Input = Input, Output = Output)
         }),
         DR_plots = map(Data, function(.x){
           Plotoutput <- .x %>%
             ggplot(aes(log10(Input), Output)) + 
             geom_line() +
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

DR_curves_output_small <- DR_curves_output %>%
  select(Ligand, Data) %>%
  unnest(Data) %>%
  filter(Output < 8)

# Calculating the required means and error bars (SD is used as of the first attempt)
ForPlots <- Exp_DR_10G4_ForExport_Final %>%
  group_by(Sample, Mutation, Ligand, Dilution) %>%
  summarize(SmeanCCRAviaDMSO = mean(CCRAviaDMSO)) %>% #Taking the mean of Samples that fit into the same group_by grouping so as to avoid over-representating a cilia sample
  ungroup() %>%
  left_join(Zero_Baseline, by = c("Mutation", "Ligand")) %>%
  mutate(meanBaselineCCRAviaDMSO = SmeanCCRAviaDMSO - minmeanCC) %>%
  rename(SMBL_CCRAviaDMSO = meanBaselineCCRAviaDMSO) %>%
  group_by(Mutation, Ligand, Dilution) %>%
  nest() %>%
  mutate(MeanCCRAviaDMSO = map(data, function(x) mean(x$SMBL_CCRAviaDMSO)),
         SDCCRAviaDMSO = map(data, function(x) sd(x$SMBL_CCRAviaDMSO)),
         SEMCCRAviaDMSO = map(data, function(x) sd(x$SMBL_CCRAviaDMSO)/sqrt(nrow(x))),
         CICCRAviaDMSO = map(data, function(x) qt(0.975, df = nrow(x) - 2) * sd(x$SMBL_CCRAviaDMSO)/sqrt(nrow(x)))) %>%
  unnest(c(MeanCCRAviaDMSO, SDCCRAviaDMSO, SEMCCRAviaDMSO, CICCRAviaDMSO)) %>%
  ungroup()
  
ForPlots %>%
  left_join(DR_curves_output_small, by = "Ligand") %>%
  mutate(Ligand = factor(Ligand, levels = c("DMSO", "FSK", "Guaiacol",  "Ethyl vanillin", "Vanillin", "o-Cresol", "2-Ethylphenol", "Dehydrodivanillin", 
                                            "Acetovanillone", "Vanillyl butyl ether", "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde",
                                            "5-Fluoroguaiacol", "4-Fluoroguaiacol", "5-Bromoguaiacol", "4-Bromoguaiacol", "4-Chloroguaiacol", "4-Iodoguaiacol",
                                            "5-Methoxyguaiacol", "4-Methoxyguaiacol", "4-Methylguaiacol", "2-Ethoxyphenol"))) %>%
  ggplot() + 
  geom_point(aes(log10(Dilution), MeanCCRAviaDMSO, color = Mutation, group = Mutation)) + 
  geom_errorbar(aes(x = log10(Dilution), color = Mutation, group = Mutation, ymin = MeanCCRAviaDMSO - SDCCRAviaDMSO, ymax = MeanCCRAviaDMSO + SDCCRAviaDMSO), width = 0.2) + 
  geom_line(aes(log10(Input), Output)) +
  scale_x_continuous(breaks = c(-10, -8, -6, -4, -2, -0)) +
  facet_wrap(~Ligand, nrow = 2, scales = "free_y") +
  labs(y = bquote((cAMP[Ligand] - cAMP[DMSO])/cAMP[DMSO]), 
       caption = "Mean +/- SD.") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(face="bold", size = 10),
        axis.text.x = element_text(angle = 15, hjust = 1),
        axis.title.y = element_text(size = 14),
        legend.box = "horizontal",
        legend.position="bottom",
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(face="bold"),
        legend.title = element_blank())

# Making stat tables, leaving the results alone for now. Do not know how to best present them
StatTest <- Exp_DR_10G4_ForExport_Final %>%
  group_by(Sample, Mutation, Ligand, Dilution) %>%
  summarize(SmeanCCRAviaDMSO = mean(CCRAviaDMSO)) %>% #Taking the mean of Samples that fit into the same group_by grouping so as to avoid over-representating a cilia sample
  ungroup() %>%
  left_join(Zero_Baseline, by = c("Mutation", "Ligand")) %>%
  mutate(meanBaselineCCRAviaDMSO = SmeanCCRAviaDMSO - minmeanCC) %>%
  rename(SMBL_CCRAviaDMSO = meanBaselineCCRAviaDMSO) %>%
  group_by(Ligand) %>%
  nest() %>%
  mutate(Ttest_Dilution_pairs = map(data, ~ compare_means(SMBL_CCRAviaDMSO ~ Mutation, .x, method = "t.test", group.by = "Dilution")),
         Anova_Mutation_sets = map(data, ~ compare_means(SMBL_CCRAviaDMSO ~ Dilution, .x, method = "anova", group.by = "Mutation")),
         Ttest_Mutation_sets = map(data, ~ compare_means(SMBL_CCRAviaDMSO ~ Dilution, .x, method = "t.test", group.by = "Mutation"))) %>%
  ungroup()

# CV for DR work ----
# Raw value Intra-Assay CV: 8.58
Exp_DR_10G4 %>%
  select(Full_Table) %>%
  unnest(Full_Table) %>%
  group_by(Exp, Sample, Ligand, Dilution) %>%
  summarize(repMean = mean(Concentration, na.rm = TRUE),
            repSD = sd(Concentration, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(CV = repSD/repMean * 100) %>%
  group_by(Exp) %>%
  summarize(MeanCV = mean(CV, na.rm = TRUE)) %>%
  summarize(mean = mean(MeanCV))

# Filtered value Intra-Assay CV: 8.32
Exp_DR_10G4 %>%
  select(QC) %>%
  unnest(QC) %>%
  group_by(Exp, Sample, Ligand, Dilution) %>%
  summarize(repMean = mean(Concentration, na.rm = TRUE),
            repSD = sd(Concentration, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(CV = repSD/repMean * 100) %>%
  group_by(Exp) %>%
  summarize(MeanCV = mean(CV, na.rm = TRUE)) %>%
  summarize(mean = mean(MeanCV))

##Raw value Inter-Assay CV: 32.2
Exp_DR_10G4 %>%
  select(Full_Table) %>%
  unnest(Full_Table) %>%
  group_by(Exp, Sample, Ligand, Dilution) %>%
  summarize(repMean = mean(Concentration, na.rm = TRUE),
            repSD = sd(Concentration, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Sample, Ligand, Dilution) %>%
  summarize(PlatesMeans = mean(repMean, na.rm = TRUE),
            PlateSd = sd(repMean, na.rm = TRUE)) %>%
  ungroup() %>%
  summarize(CV = PlateSd/PlatesMeans * 100) %>%
  summarize(mean = mean(CV, na.rm = TRUE))

##Filtered value Inter-Assay CV: 31.0
Exp_DR_10G4 %>%
  select(QC) %>%
  unnest(QC) %>%
  group_by(Exp, Sample, Ligand, Dilution) %>%
  summarize(repMean = mean(Concentration, na.rm = TRUE),
            repSD = sd(Concentration, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Sample, Ligand, Dilution) %>%
  summarize(PlatesMeans = mean(repMean, na.rm = TRUE),
            PlateSd = sd(repMean, na.rm = TRUE)) %>%
  ungroup() %>%
  summarize(CV = PlateSd/PlatesMeans * 100) %>%
  summarize(mean = mean(CV, na.rm = TRUE))

# Process a single DR experiment starting from the Post-QC, post-mean triplicate steps: just use the correct Exp value. ----
FilteredExp <- "Exp45a" # Replace with SINGLE experiment.
SingleDR <- Exp_DR_10G4_ForExport %>%
  filter(Exp == FilteredExp) %>%
  select(Concentration_Triplicate_mean) %>%
  unnest(Concentration_Triplicate_mean) %>%
  rename(Concentration = Concentration_S_Measure) %>%
  select(Sample, Mutation, Ligand, Dilution, Concentration)

ExpDMSO <- SingleDR %>%
  filter(Ligand == "DMSO") %>%
  select(Sample, Concentration) %>%
  rename(DMSOconcentration = Concentration)

SingleDR2 <- SingleDR %>%
  left_join(ExpDMSO, by = "Sample") %>%
  mutate(CCviaDMSO = Concentration - DMSOconcentration) 

ExpFSK <- SingleDR2 %>%
  filter(Ligand == "FSK") %>%
  select(Sample, Concentration, CCviaDMSO) %>%
  rename(FSKconcentration = Concentration,
         FSKCCviaDMSO = CCviaDMSO)

SingleDR3 <- SingleDR2 %>%
  left_join(ExpFSK, by = "Sample") %>%
  mutate(RAviaFSK = Concentration/FSKconcentration,
         RAviaDMSO = Concentration/DMSOconcentration,
         CCRAviaFSKCC = CCviaDMSO/FSKCCviaDMSO, 
         CCRAviaDMSO = CCviaDMSO/DMSOconcentration)

SingleDR4 <- SingleDR3 %>%
  select(Sample, Mutation, Ligand, Concentration, Dilution, DMSOconcentration, FSKconcentration, RAviaDMSO, CCviaDMSO, RAviaFSK, CCRAviaFSKCC, CCRAviaDMSO) %>%
  group_by(Mutation, Ligand, Dilution) %>%
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
         CICCRAviaDMSO = map(data, function(x) qt(0.975, df = nrow(x) - 2) * sd(x$CCRAviaDMSO)/sqrt(nrow(x)))) %>%
  unnest(c(MeanRAviaDMSO, SDRAviaDMSO, MeanRAviaFSK, SDRAviaFSK, MeanCCRAviaFSKCC, SDCCRAviaFSKCC, MeanCCRAviaDMSO, SDCCRAviaDMSO, SEMCCRAviaDMSO, CICCRAviaDMSO, data)) %>%
  ungroup()

#Statistics
SingleDR4_trim <- SingleDR4 %>%
  filter(Ligand != "DMSO") %>%
  filter(Ligand != "FSK")

# Specific comparisons were made, most critically the output value CCRAviaDMSO. Change to test other Relative Activation values
compare_means(CCRAviaDMSO ~ Mutation, SingleDR4_trim, method = "anova")
compare_means(CCRAviaDMSO ~ Mutation, SingleDR4_trim, method = "t.test", group.by = c("Ligand", "Dilution")) 

# Graph without DMSO and FSK, graphing the DMSO-corrected concentration relative to DMSO values. CCRAviaDMSO/MeanCCRAviaDMSO/SDCCRAviaDMSO need to be changed together 
SingleDR4_trim %>%
  ggplot() + 
  geom_point(aes(log10(Dilution), CCRAviaDMSO, color = Sample)) +
  geom_line(aes(log10(Dilution), MeanCCRAviaDMSO, color = Mutation)) +
  geom_errorbar(aes(x = log10(Dilution), group = Mutation, ymin = MeanCCRAviaDMSO - SDCCRAviaDMSO, ymax = MeanCCRAviaDMSO + SDCCRAviaDMSO), 
                width = 0.2, color = "black") +
  facet_wrap(~Ligand) +
  labs(title = FilteredExp)

# Continuing on to collectively analyze all the Static 50uM 10G4 experiments, essentially taking them through the SingleDR path, but without the intermediates visible----
# These experiments in the Exp_Static50uM_10G4 set require some review and clean up prior to analysis and plots.
map2(Exp_Static50uM_10G4$Exp, Exp_Static50uM_10G4$QC, ~ if(TRUE %in% .y$Unbound | max(.y$CV) > 30) {.x} ) %>% unlist()

# Removing poor performing wells
Exp_Static50uM_10G4$QC[["Exp54OldNew"]] <- Exp_Static50uM_10G4$QC[["Exp54OldNew"]] %>%
  filter(Position != 41)
Exp_Static50uM_10G4$QC[["Exp54OldOld"]] <- Exp_Static50uM_10G4$QC[["Exp54OldOld"]] %>%
  filter(Position != 22)

# Stupid large file, possibly smarter to split it with multiple separate map calls? Works fine. 
Exp_Static50uM_10G4_Processed <- Exp_Static50uM_10G4 %>%
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

# Run this line to see the anova results
Exp_Static50uM_10G4_Processed$Anova_Results_CCRAviaDMSO

# Run the below line to see some basic Ligand vs CCRAviaDMSO plots
Exp_Static50uM_10G4_Processed$First_plots

# Making a plot for all tested odors, requires manual modification
Aggregateplot <- rbindlist(Exp_Static50uM_10G4_Processed$All_Analysis_Values, use.names=TRUE) %>%
  filter(Ligand %notin% c("FSK", "DMSO")) %>%
  mutate(Sample = factor(Sample, levels = c("PFH1721", "PFH1722", "PFH1741","PFH1742","PFH1743","PFH887", "PFH888", 
                                            "PFH1790","PFH1885", "PFH1789", "PFH1792", "PFH1793","PFH1884", 
                                            "PFH1547", "PFH1548", "PFH1550", "PFH1740", "PFH1744", "PFH1745", "PFH1549", "PFH1551", "PFH1553"))) %>%
  mutate(Ligand = factor(Ligand, levels = c("DMSO", "FSK", "Guaiacol",  "Ethyl vanillin", "Vanillin", "o-Cresol", "2-Ethylphenol", "Dehydrodivanillin", 
                                            "Acetovanillone", "Vanillyl butyl ether", "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde",
                                            "5-Fluoroguaiacol", "4-Fluoroguaiacol", "5-Bromoguaiacol", "4-Bromoguaiacol", "4-Chloroguaiacol", "4-Iodoguaiacol",
                                            "5-Methoxyguaiacol", "4-Methoxyguaiacol", "4-Methylguaiacol", "2-Ethoxyphenol"))) %>%
  ungroup() %>%
  ggplot(aes(Ligand, CCRAviaDMSO2)) + 
  geom_violin(aes(color = Mutation), position = position_dodge(width = 0.3)) +
  geom_point(aes(color = Mutation), position = position_dodge(width = 0.3), size = 0.5, alpha = 0.8) +
  labs(title = "Violin Plot of Aggregated Results from all 50uM-concentration BioAssays", y = bquote(Signal:(cAMP[Ligand] - cAMP[DMSO])/cAMP[DMSO]), 
       caption = "Seven BioAssay experiments total, each odor tested at least three times. Seven to Eight biological replicates for each Genotype, tested in triplicate per experiment.
       Individual Signal values from all experiments were plotted. P values are log10-scale means from all experiments. P values NOT corrected for multiple comparisons; asterisks mark p values that pass with a Holm correction.") +
  scale_y_continuous(limits = c(-0.7, 8.5)) +
  geom_signif(y_position = c(2.9, 2.4, 1.3, 1.9, 2.5, 1.3, 1.2, 1.0, 1.6, 2.4, 2.2, 2.2, 6.3, 3.8, 8.3, 6.8, 3.8, 2.8, 7.0, 3.8, 3.8, 6.4), 
              xmin = c(0.8, 1.8, 2.8, 3.8, 4.8, 5.8, 6.8, 7.8, 8.8, 9.8, 10.8, 11.8, 12.8, 13.8, 14.8, 15.8, 16.8, 17.8, 18.8, 19.8, 20.8, 21.8), 
              xmax = c(1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 7.2, 8.2, 9.2, 10.2, 11.2, 12.2, 13.2, 14.2, 15.2, 16.2, 17.2, 18.2, 19.2, 20.2, 21.2, 22.2), 
              annotation = c("0.0012*", 0.0138, 0.0169, 0.0516, 0.0140, 0.136, 0.293, 0.478, 0.290, 0.279, 0.394, 0.281,
                             "0.0025*", "0.0012*", "0.0021*", "0.0014*", 0.0040, 0.0921, 0.0060, 0.0078, "0.0009*", "0.0025*"), 
              tip_length = 0, color = "black", textsize = 4) +
  theme_classic() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(face="bold", size = 10),
        axis.text.x = element_text(angle = 15, hjust = 1),
        axis.title.y = element_text(size = 14),
        legend.position = c(0.50, 0.90),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(face="bold"),
        legend.title = element_blank())

# See the plot
Aggregateplot

# Generate some text and the t test results
TextGrob <- text_grob(paste("Twenty-two independent t.tests were performed. With alpha = 0.05, Holm and Dunn-Šidák corrections were calculated. ", 
                             "Dunn-Šidák corrected alpha evaluated to 0.00285.", sep = "\n"), 
                       face = "italic")

TableGrob <- rbindlist(Exp_Static50uM_10G4_Processed$Ttest_Results_CCRAviaDMSO, use.names=TRUE) %>%
  mutate(Log10p = log10(p)) %>%
  group_by(Ligand) %>%
  summarize(MeanLog10p = mean(Log10p)) %>%
  mutate(Aggregate_p = 10^MeanLog10p) %>% 
  arrange(Aggregate_p) %>%
  mutate(row = row_number(),
         HolmAdjustment = 0.05/(22 - row + 1), 
         PassHA = Aggregate_p < HolmAdjustment, 
         HolmCorrection = Aggregate_p * (22 - row + 1),
         PassHC = HolmCorrection < 0.05, 
         Pass_DunnSidakCorrection = Aggregate_p < 1 - (1 - 0.05)^(1/22)) %>%
  mutate(Ligand = factor(Ligand, levels = c("DMSO", "FSK", "Guaiacol",  "Ethyl vanillin", "Vanillin", "o-Cresol", "2-Ethylphenol", "Dehydrodivanillin", 
                                            "Acetovanillone", "Vanillyl butyl ether", "2-Hydroxy-4-Methoxybenzaldehyde", "Menthoxypropanediol", "Eugenol", "Salicylaldehyde",
                                            "5-Fluoroguaiacol", "4-Fluoroguaiacol", "5-Bromoguaiacol", "4-Bromoguaiacol", "4-Chloroguaiacol", "4-Iodoguaiacol",
                                            "5-Methoxyguaiacol", "4-Methoxyguaiacol", "4-Methylguaiacol", "2-Ethoxyphenol"))) %>%
  arrange(Ligand) %>%
  tableGrob(rows=NULL)

# Overview of the p value work
grid.arrange(TableGrob, TextGrob, heights = c(0.6, 0.1))