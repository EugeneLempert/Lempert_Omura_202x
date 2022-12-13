# Functions written and used in this publication:
# %notin% : Convenience function that provides a vectorized != logical for filter(); found online, unknown source ----
`%notin%` <- Negate(`%in%`)

# read_plus : Import all Cq and Melt files (https://stackoverflow.com/questions/11433432/how-to-import-multiple-csv-files-at-once), Run before analysis code ----
read_plus <- function(flnm) {
  read_delim(flnm, 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE, skip = 1)
}

# Label_qPCR_data : Label the Cq data, some code from stackoverflow by Tony Breyal ----

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
# Quality_Evaluation : Optionally combine with a corresponding curve Melt file and provide some (default) QC metrics ----
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
# Plot_QE : Some plots to evaluate the qPCR data ----
Plot_QE <- function(QE_table, MutationColumn, ...){
  
  P1 <- QE_table %>% ggplot(aes(RepSd, RepMean, color = Sample)) + geom_point() + facet_wrap(~Gene)
  
  P2 <- QE_table %>% ggplot(aes(Tm, Height, color = Sample)) + geom_point() + facet_wrap(~Gene) + geom_errorbarh(aes(xmax = Tm + Width/2, xmin = Tm - Width/2))
  
  P3 <- QE_table %>% ggplot(aes(.data[[MutationColumn]], RepMean, fill = Sample)) + geom_bar(position = "dodge", stat = "identity") + facet_wrap(vars(Gene))
  
  P4 <- QE_table %>% ggplot(aes(RepMean, Sample, color = Gene, size = RepSd)) + geom_point()
  
  Plot_list <- list(SdMean = P1, 
                    TmHeight = P2, 
                    RawRQ = P3, 
                    Combo = P4)
  
  return(Plot_list)
  
}