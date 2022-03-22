
# PACKAGES ----------
# Data Manipulation
library(tidyverse)
library(reshape2)
library(biomaRt)
library(limma)
library(edgeR)
library(variancePartition)

# Plots
library(ggpubr)
library(pals)
library(RColorBrewer)
library(viridis)
library(svglite)
library(wesanderson)

# Statistics
library(rstatix)
library(emmeans)

# circRNA Sequences
library(FcircSEC)


# BEFORE USING THIS SCRIPT ----------
# User should alter these variables as necessary

# General script
main_dir <- "G:/Noctsol/GitHub Repositories/temp/"
pheno_dir <- paste0(main_dir, "powell/")
circ_dir <- paste0(main_dir, "powell/circ/")

# Input and output directories for FcircSEC
in_dir <- "G:/Where the Big Files Are/Rat Reference Files/in"
out_dir <- "G:/Where the Big Files Are/Rat Reference Files/out"

# Update and set ggplot theme globally
new_theme <- theme_bw(base_size = 12) %+replace%
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))
theme_set(new_theme)


# LOAD IN DATA & RENAME COLUMNS ----------
# NOTE: Sample order must be maintained for accuracy.
# Sample order is always numerical from low (SRS6085261) to high (SRS6085272).

# Load in phenotype file
pheno <- read_csv(paste0(pheno_dir, "powell_pheno.csv", sep = ""))
sample_ids <- pheno$Sample_Treatment # Vector of sample IDs

# Load in circRNA files
# Get list of files, then load in all files
ciri2_files <- list.files(path = circ_dir, pattern = "ciri2", full.names = TRUE)
ce2_files <- list.files(path = circ_dir, pattern = "known", full.names = TRUE)

# Read in above files
list2env(envir = .GlobalEnv,
         lapply(setNames(ciri2_files, make.names(paste0("ciri2_", sample_ids))), 
                read.delim, sep = "\t", header = TRUE,
                col.names = c("CircRNA_ID", "Chromosome", "Start", "End",
                              "Junction_Reads", "SM_MS_SMS", "NonJunction_Reads",
                              "Junction_Reads_Ratio", "CircRNA_Type", "Gene_ID",
                              "Strand", "Junction_Reads_ID")))

list2env(envir = .GlobalEnv,
         lapply(setNames(ce2_files, make.names(paste0("ce2_", sample_ids))), 
                read.delim, sep = "\t", header = FALSE, 
                col.names = c("Chromosome", "Start", "End", "Name", "Score",
                              "Strand", "ThickStart", "ThickEnd", "ItemRGB",
                              "ExonCount", "ExonSizes", "ExonOffsets",
                              "ReadNumber", "CircType", "GeneName",
                              "IsoformName", "Index", "FlankIntron")))


# Create list of lists for easier transformation into data frame
ciri2_sample_list <- list("SRS6085261_CI1cue" = ciri2_SRS6085261_CI1cue,
                          "SRS6085262_CI21cue" = ciri2_SRS6085262_CI21cue,
                          "SRS6085263_CE1cue" = ciri2_SRS6085263_CE1cue,
                          "SRS6085264_CE21cue" = ciri2_SRS6085264_CE21cue,
                          "SRS6085265_CE21cue" = ciri2_SRS6085265_CE21cue,
                          "SRS6085266_CE1cue" = ciri2_SRS6085266_CE1cue,
                          "SRS6085267_CI1cue" = ciri2_SRS6085267_CI1cue,
                          "SRS6085268_CI1cue" = ciri2_SRS6085268_CI1cue,
                          "SRS6085269_CI21cue" = ciri2_SRS6085269_CI21cue,
                          "SRS6085270_CE1cue" = ciri2_SRS6085270_CE1cue,
                          "SRS6085271_CE21cue" = ciri2_SRS6085271_CE21cue,
                          "SRS6085272_CI21cue" = ciri2_SRS6085272_CI21cue) %>%
  map2(., names(.), ~cbind(.x, Sample_ID = .y))

ce2_sample_list <- list("SRS6085261_CI1cue" = ce2_SRS6085261_CI1cue,
                        "SRS6085262_CI21cue" = ce2_SRS6085262_CI21cue,
                        "SRS6085263_CE1cue" = ce2_SRS6085263_CE1cue,
                        "SRS6085264_CE21cue" = ce2_SRS6085264_CE21cue,
                        "SRS6085265_CE21cue" = ce2_SRS6085265_CE21cue,
                        "SRS6085266_CE1cue" = ce2_SRS6085266_CE1cue,
                        "SRS6085267_CI1cue" = ce2_SRS6085267_CI1cue,
                        "SRS6085268_CI1cue" = ce2_SRS6085268_CI1cue,
                        "SRS6085269_CI21cue" = ce2_SRS6085269_CI21cue,
                        "SRS6085270_CE1cue" = ce2_SRS6085270_CE1cue,
                        "SRS6085271_CE21cue" = ce2_SRS6085271_CE21cue,
                        "SRS6085272_CI21cue" = ce2_SRS6085272_CI21cue) %>%
  map2(., names(.), ~cbind(.x, Sample_ID = .y)) %>%
  lapply(mutate, Chromosome = as.character(Chromosome))

# Sanity check: Sample IDs in the sample lists & pheno files are in the same order
identical(names(ciri2_sample_list), sample_ids) # Should be TRUE
identical(names(ce2_sample_list), sample_ids) # Should be TRUE

# Convert list of lists to single data frame
ciri2_original <- Reduce(function(x, y) 
  full_join(x, y), ciri2_sample_list)

ce2_original <- Reduce(function(x, y) 
  full_join(x, y), ce2_sample_list)

# Tidy dataframes for both tools, getting rid of extraneous columns. 
# For CircExplorer2, create a CircRNA_ID column to match CIRI2.
# Note that CIRI2 start points are off by 1 base and must be revised.
# Rename similar columns to have similar names between dataframes.
ciri2 <- ciri2_original %>%
  dplyr::select(Sample_ID, CircRNA_ID, Gene_ID, CircRNA_Type, Junction_Reads, Strand) %>%
  dplyr::rename(CIRI2_Gene_ID = "Gene_ID", CIRI2_CircRNA_Type = "CircRNA_Type",
         CIRI2_Junction_Reads = "Junction_Reads") %>%
  mutate(across(everything(), function(x) na_if(x, "n/a"))) %>%
  mutate(CIRI2_Gene_ID = gsub(",", "", CIRI2_Gene_ID))

ce2 <- ce2_original %>%
  mutate(NewStart = Start + 1) %>%
  unite(CircRNA_ID, c(Chromosome, NewStart), sep = ":", remove = FALSE) %>%
  unite(CircRNA_ID, c(CircRNA_ID, End), sep = "|", remove = FALSE) %>%
  dplyr::select(Sample_ID, CircRNA_ID, GeneName, CircType, ReadNumber, Strand) %>%
  dplyr::rename(CE2_Transcript_ID = "GeneName", CE2_CircRNA_Type = "CircType",
         CE2_Junction_Reads = "ReadNumber")

# Join together into a single dataframe
complete_circs <- inner_join(ciri2, ce2)

# Create dataframe with mapped reads from BWA for each sample
# Mapped reads quantified from BAMs using bamtools stats
bwa_maps <- data.frame(Sample_ID = sample_ids,
                       Mapped_Reads = c(37836418, 41818651, 28744392, 63373025, 
                                        45085767, 41954076, 44177399, 39936764, 
                                        33314139, 49986109, 52199617, 54460142))
# Check mean mapped reads
mean_maps <- mean(bwa_maps$Mapped_Reads)
mean_maps # 44407233
mean(bwa_maps[-9,]$Mapped_Reads) # 45415669

# Merge dataframes and calculate RPM for each circRNA
# RPM = JunctionReads/(Mapped_Reads/1000000)
all_circs <- Reduce(function(x, y)
  full_join(x, y), list(ciri2, ce2, bwa_maps)) %>%
  mutate(CIRI2_RPM = CIRI2_Junction_Reads/(Mapped_Reads/1000000),
         CE2_RPM = CE2_Junction_Reads/(Mapped_Reads/1000000)) %>%
  full_join(pheno, by = c("Sample_ID" = "Sample_Treatment")) %>%
  dplyr::select(c(1:12, 20:26)) 


## PLOT CIRCRNA COUNTS BEFORE FILTERING ---------
# Plot for each sample and the number of circRNAs called by each tool
# First modify dataframe
counts_4plot <- all_circs %>%
  mutate(Sample_ID = strtrim(Sample_ID, 10)) %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Sample_ID = 
           ordered(Sample_ID, 
                  levels = c("SRS6085261", "SRS6085267", "SRS6085268", 
                             "SRS6085263", "SRS6085266", "SRS6085270",
                             "SRS6085262", "SRS6085269", "SRS6085272", 
                             "SRS6085264", "SRS6085265", "SRS6085271"))) %>%
  mutate(Treatment_Group = fct_relevel(case_when(Treatment_Group == "CI1cue" ~ "1d ABS, Isolation",
                                                 Treatment_Group == "CE1cue" ~ "1d ABS, Enrichment",
                                                 Treatment_Group == "CI21cue" ~ "21d ABS, Isolation",
                                                 Treatment_Group == "CE21cue" ~ "21d ABS, Enrichment"),
                                       levels = c("1d ABS, Isolation", "1d ABS, Enrichment",
                                                  "21d ABS, Isolation", "21d ABS, Enrichment"))) %>%
  # mutate(Sample_ID = case_when(Sample_ID == "SRS6085261_CI1cue" ~ "Iso 1d - A",
  #                              Sample_ID == "SRS6085267_CI1cue" ~ "Iso 1d - B",
  #                              Sample_ID == "SRS6085268_CI1cue" ~ "Iso 1d - C",
  #                              Sample_ID == "SRS6085263_CE1cue" ~ "EE 1d - A",
  #                              Sample_ID == "SRS6085266_CE1cue" ~ "EE 1d - B",
  #                              Sample_ID == "SRS6085270_CE1cue" ~ "EE 1d - C",
  #                              Sample_ID == "SRS6085262_CI21cue" ~ "Iso 21d - A",
  #                              Sample_ID == "SRS6085269_CI21cue" ~ "Iso 21d - B",
  #                              Sample_ID == "SRS6085272_CI21cue" ~ "Iso 21d - C",
  #                              Sample_ID == "SRS6085264_CE21cue" ~ "EE 21d - A",
  #                              Sample_ID == "SRS6085265_CE21cue" ~ "EE 21d - B",
  #                              Sample_ID == "SRS6085271_CE21cue" ~ "EE 21d - C")) %>%
  mutate(Identified_By = 
           factor(case_when(!is.na(CIRI2_Junction_Reads) & !is.na(CE2_Junction_Reads) ~ "Both",
                            !is.na(CIRI2_Junction_Reads) & is.na(CE2_Junction_Reads) ~ "CIRI2",
                            is.na(CIRI2_Junction_Reads) & !is.na(CE2_Junction_Reads) ~ "CIRCexplorer2"),
                  levels = c("CIRI2", "CIRCexplorer2", "Both"))) %>%
  pivot_longer(cols = c(CIRI2_Junction_Reads, CE2_Junction_Reads),
               names_to = "Tool", values_to = "Count") %>%
#  mutate(Count = ifelse(Count >= 30, 30, Count)) %>%
  filter(!is.na(Count)) %>%
  ungroup() %>%
  mutate(Count = ifelse(Count >= 30, 30, Count))

# This is zoomed in to just read count >= 60
count_plot <- counts_4plot %>%
  ggplot(aes(x = Count, fill = Identified_By)) +
  geom_histogram(binwidth = 1, color = "black") +
  labs(x = "Read Count", y = "Number of circRNAs Identified") +
  # X axis limits
  scale_x_continuous(labels = c("0", "10", "20", "30+")) +
  scale_fill_manual(name = "Tool", values = c("white", "grey80", "grey30")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color = "grey80")) +
  geom_vline(xintercept = 1.56, lty = "longdash", color = "red") +
  facet_wrap(~ Sample_ID, ncol = 3, scales = "free_y")

count_plot

count_gtable <- ggplot_gtable(ggplot_build(count_plot))
striprt <- which(grepl('strip-r', count_gtable$layout$name) | 
                   grepl('strip-t', count_gtable$layout$name))
fills <- rep(c(rep("#3F826D", 3), rep("#F2B272", 3), 
               rep("#B2CDC5", 3), rep("#FBE8D5", 3)))
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', count_gtable$grobs[[i]]$grobs[[1]]$childrenOrder))
  count_gtable$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
  
grid::grid.draw(count_gtable)

svglite("C:/Users/Annika/Documents/Figures for circRNA Project/Read_counts.svg", width = 8, height = 7)
grid::grid.draw(count_gtable)
dev.off()

svglite("C:/Users/Annika/Documents/Figures for circRNA Project/Read_counts_legend2.svg", width = 8, height = 7)
counts_4plot %>%
  ggplot(aes(x = Sample_ID, y = Count, fill = Treatment_Group)) +
  geom_col(color = "black") +
  scale_fill_manual(name = "Treatment Group", values = c("#FBE8D5", "#B2CDC5",
                                                         "#F2B272", "#3F826D"))
dev.off()

# Get number of unique circRNAs identified in any sample for each tool
all_circs %>%
  dplyr::select(CircRNA_ID, CIRI2_CircRNA_Type, CE2_CircRNA_Type) %>%
  unique() %>%
  mutate(CE2_CircRNA_Type = case_when(CE2_CircRNA_Type == "circRNA" ~ "exon",
                                      CE2_CircRNA_Type == "ciRNA" ~ "intron"),
         Both_ID = case_when(CE2_CircRNA_Type == CIRI2_CircRNA_Type ~ CE2_CircRNA_Type,
                             TRUE ~ NA_character_)) %>%
  pivot_longer(cols = 2:4, names_to = "Tool", values_to = "Type") %>%
  mutate(Tool = fct_relevel(case_when(str_detect(Tool, "Both") ~ "Both",
                                      str_detect(Tool, "CE2") ~ "CIRCexplorer2",
                                      str_detect(Tool, "CIRI2") ~ "CIRI2"),
                            levels = c("CIRI2", "CIRCexplorer2", "Both")),
         Type = ordered(Type, levels = c("exon", "intron", "intergenic_region"))) %>%
  filter(!is.na(Type)) %>%
  group_by(Tool, Type) %>%
  mutate(count = length(Type)) %>%
  ungroup() %>%
  View()
  
# Get number of unique circRNAs identified in any sample for each tool
type_plot <- all_circs %>%
  dplyr::select(CircRNA_ID, CIRI2_CircRNA_Type, CE2_CircRNA_Type) %>%
  unique() %>%
  mutate(CE2_CircRNA_Type = case_when(CE2_CircRNA_Type == "circRNA" ~ "exon",
                                      CE2_CircRNA_Type == "ciRNA" ~ "intron"),
         Both_ID = case_when(CE2_CircRNA_Type == CIRI2_CircRNA_Type ~ CE2_CircRNA_Type,
                             TRUE ~ NA_character_)) %>%
  pivot_longer(cols = 2:4, names_to = "Tool", values_to = "Type") %>%
  mutate(Tool = fct_relevel(case_when(str_detect(Tool, "Both") ~ "Both",
                                      str_detect(Tool, "CE2") ~ "CIRCexplorer2",
                                      str_detect(Tool, "CIRI2") ~ "CIRI2"),
                            levels = c("CIRI2", "CIRCexplorer2", "Both")),
         Type = ordered(Type, levels = c("exon", "intron", "intergenic_region"))) %>%
  filter(!is.na(Type)) %>%
  group_by(Tool, Type) %>%
  mutate(count = length(Type)) %>%
  ungroup() %>%
  ggplot(aes(x = Tool, fill = Type)) +
  geom_bar(position = position_dodge(preserve = "single"), color = "black") +
  scale_y_continuous(limits = c(0, 3550), breaks = seq(0, 3500, 1000), 
                     expand = c(0, 0)) +
   scale_fill_manual(values = c("#A2A475", "#02401B", "#D8B70A"), name = "Span",
                     labels = c("Exon(s)", "Intron(s)", "Intergenic Region")) +
  labs(y = "Number Identified", x = "") +
  theme(plot.title = element_text(size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey70"),
        panel.grid.minor.y = element_line(color = "grey85"),
        axis.ticks.x = element_blank())

type_plot

svglite("C:/Users/Annika/Documents/Figures for circRNA Project/Type_counts.svg", width = 8.2, height = 2)
type_plot
dev.off()



 ## FILTER CIRCRNAS ----------
# How many circRNAs per sample, per tool before filtering?
ciri2_before <- all_circs %>%
  filter(!is.na(CIRI2_Junction_Reads)) %>%
  group_by(Sample_ID) %>%
  summarize(CIRI2_Before = n())

ce2_before <- all_circs %>%
  filter(!is.na(CE2_Junction_Reads)) %>%
  group_by(Sample_ID) %>%
  summarize(CE2_Before = n())

circs_before <- full_join(ciri2_before, ce2_before)
circs_before

# Filtering by just 1 tool
readcount2_1tool_ciri2 <- all_circs %>%
  filter(CIRI2_Junction_Reads >= 2) %>%
  group_by(Sample_ID) %>%
  summarize(CIRI2_Read2 = n())

readcount2_1tool_ce2 <- all_circs %>%
  filter(CE2_Junction_Reads >= 2) %>%
  group_by(Sample_ID) %>%
  summarize(CE2_Read2 = n())


# Get mean RPM per treatment group per circRNA, and filter out circRNAs w/
# a read count of less than 2 across samples for both tools, and calculate the 
# number of samples within a treatment group that have a circRNA
keep_circs <- all_circs %>%
  group_by(Treatment_Group, CircRNA_ID) %>%
  mutate(CIRI2_MeanRPM = mean(CIRI2_RPM),
         CE2_MeanRPM = mean(CE2_RPM)) %>%
  filter(CIRI2_Junction_Reads >= 2, CE2_Junction_Reads >= 2) %>%
  mutate(SampleTreatment_Count = n()) %>%
  group_by(CircRNA_ID) %>%
  mutate(Sample_Count = n()) %>%
  ungroup()

# How many circular RNAs are left after filtering?
all_after <- keep_circs %>%
  group_by(Sample_ID) %>%
  summarize(After_Circs = n())

# Put before/after summary data together
filter_summary <- Reduce(function(x, y) 
  full_join(x, y), list(ciri2_before, ce2_before, readcount2_1tool_ciri2,
                        readcount2_1tool_ce2, all_after)) %>%
  arrange(After_Circs)
filter_summary


# Get circRNAs that are shared by the tools with no required read count threshold
circs_nofilter <- all_circs %>%
  group_by(Treatment_Group, CircRNA_ID) %>%
  mutate(CIRI2_MeanRPM = mean(CIRI2_RPM),
         CE2_MeanRPM = mean(CE2_RPM)) %>%
  filter(!is.na(CIRI2_Junction_Reads), !is.na(CE2_Junction_Reads)) %>%
  mutate(SampleTreatment_Count = n()) %>%
  group_by(CircRNA_ID) %>%
  mutate(Sample_Count = n()) %>%
  ungroup()

circs_nofilter_summary <- circs_nofilter %>%
  group_by(Sample_ID) %>%
  summarize(Both_AnyRead = n())



## FIND CIRCRNAS EXPRESSED IN ALL GROUPS ----------
# Which circRNAs are still present in all 4 groups?
circs4 <- keep_circs %>%
  group_by(CircRNA_ID) %>%
  mutate(Sample_Count = n()) %>%
  filter(Sample_Count == 12) %>%
  pull(CircRNA_ID) %>%
  unique()
circs4 # "2:24583888|24598599"   "3:17741317|17753437"  "3:57821407|57847776" 
# "5:145948264|145951212" "9:105629568|105660854"

# "2:24583888|24598599" = ENSRNOG00000047014 = Homer1
# "3:17741317|17753437" = ENSRNOG00000017583 = Mapkap1
# "3:57821407|57847776" = ENSRNOG00000060479 = Sp3
# "5:145948264|145951212" = ENSRNOG00000006137 = Arid1a
# "9:105629568|105660854" = ENSRNOG00000012733 = Ankrd12

# Check with no filter requirement
circs4_nofilter <- circs_nofilter %>%
  group_by(CircRNA_ID) %>%
  mutate(Sample_Count = n()) %>%
  filter(Sample_Count == 12) %>%
  pull(CircRNA_ID) %>%
  unique()
circs4_nofilter # "2:24583888|24598599"   "3:17741317|17753437"  
# "3:57821407|57847776" "5:102326743|102327450" "5:145948264|145951212" 
# "9:105629568|105660854"

# 3 circRNAs present in all 4 groups
circs4_df <- data.frame(CircRNA_ID = circs4,
                        Gene_ID = c("ENSRNOG00000047014", "ENSRNOG00000017583",
                                    "ENSRNOG00000006137", "ENSRNOG00000060479",
                                    "ENSRNOG00000012733"),
                        Gene = c("Homer1", "Mapkap1", "Arid1a", "Sp3", 
                                 "Ankrd12"))

# Check if there are circRNAs present in all samples except the outlier.
# It's still the same 5 circRNAs.
circs4_v2 <- keep_circs %>%
  filter(Sample_ID != "SRS6085269_CI21cue") %>%
  group_by(CircRNA_ID) %>%
  mutate(Sample_Count = n()) %>%
  filter(Sample_Count == 11) %>%
  pull(CircRNA_ID) %>%
  unique()
circs4_v2 # Same circRNAs

circs4_v2_nofilter <- circs_nofilter %>%
  filter(Sample_ID != "SRS6085269_CI21cue") %>%
  group_by(CircRNA_ID) %>%
  mutate(Sample_Count = n()) %>%
  filter(Sample_Count == 11) %>%
  pull(CircRNA_ID) %>%
  unique()
circs4_v2_nofilter
# [1] "10:89324561|89325932"  "2:24583888|24598599"   "3:17741317|17753437"  
# [4] "3:57821407|57847776"   "5:102326743|102327450" "5:145948264|145951212"
# [7] "9:105629568|105660854"
# New circ: (10)


## STATISTICS FOR CIRCRNA EXPRESSION (ANOVAs) ----------
# Remove outlier and reformat data so that the RPM from each tool can be entered 
# into the ANOVA models separately, along with sequencing lane
anova_format <- keep_circs %>%
  filter(Sample_ID != "SRS6085269_CI21cue") %>%
  filter(CircRNA_ID %in% circs4_v2) %>%
  ungroup() %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Mean_CircExpr = mean(c(CIRI2_RPM, CE2_RPM))) %>%
  pivot_longer(cols = c(CIRI2_RPM, CE2_RPM), 
               names_to = "Tool", values_to = "Tool_RPM") %>%
  ungroup() %>%
  mutate(across(c(Tool, Lane, Housing_Condition, Abstinence_Length), as.factor))

# Shapiro-Wilk tests for normality
circ_shapiro <- function(circ) {
  anova_format %>%
    filter(CircRNA_ID == circ) %>%
    lm(Tool_RPM ~ Tool + Lane + Housing_Condition * Abstinence_Length, .) %>%
    residuals(.) %>%
    shapiro_test(.)
}

circ_shapiro_noToolLane <- function(circ) {
  anova_format %>%
    filter(CircRNA_ID == circ) %>%
    lm(Mean_CircExpr ~ Housing_Condition * Abstinence_Length, .) %>%
    residuals(.) %>%
    shapiro_test(.)
}

circ_shapiro_noLane <- function(circ) {
  anova_format %>%
    filter(CircRNA_ID == circ) %>%
    lm(Tool_RPM ~ Lane + Housing_Condition * Abstinence_Length, .) %>%
    residuals(.) %>%
    shapiro_test(.)
}

# Shapiro-Wilk tests - all not significant
ch1_shapiro <- circ_shapiro("2:24583888|24598599")
cm1_shapiro <- circ_shapiro("3:17741317|17753437")
car1_shapiro <- circ_shapiro("5:145948264|145951212")
cs3_shapiro <- circ_shapiro("3:57821407|57847776")
cah12_shapiro <- circ_shapiro("9:105629568|105660854")


# Perform the ANOVAs using Tool and Lane in the model
circ_anova <- function(circ) {
  options(contrasts=c('contr.sum', 'contr.poly'))
  anova_format %>%
    filter(CircRNA_ID == circ) %>%
    lm(Tool_RPM ~ Tool + Lane + Housing_Condition * Abstinence_Length, .) %>%
    car::Anova(., type = "III")
}

circ_anova_noToolLane <- function(circ) {
  options(contrasts=c('contr.sum', 'contr.poly'))
  anova_format %>%
    filter(CircRNA_ID == circ) %>%
    lm(Mean_CircExpr ~ Housing_Condition * Abstinence_Length, .) %>%
    car::Anova(., type = "III")
}

circ_anova_noLane <- function(circ) {
  options(contrasts=c('contr.sum', 'contr.poly'))
  anova_format %>%
    filter(CircRNA_ID == circ) %>%
    lm(Tool_RPM ~ Tool + Housing_Condition * Abstinence_Length, .) %>%
    car::Anova(., type = "III")
}

# ANOVAs for circ-Homer1, circ-Mapkap1, and circ-Arid1a (in that order)
ch1_anova <- circ_anova("2:24583888|24598599") # Housing, 0.097 Hous*Abs
cm1_anova <- circ_anova("3:17741317|17753437") # Tool, Lane, Housing, Hous*Abs
car1_anova <- circ_anova("5:145948264|145951212") # Tool, 0.06 Lane, Hous, Hous*Abs
cs3_anova <- circ_anova("3:57821407|57847776") # 0.065 Tool
cah12_anova <- circ_anova("9:105629568|105660854") # Tool, Hous*Abs


# Post-hoc Tukey tests (need to check this in SPSS)
circ_tukey <- function(circ) {
  anova_format %>%
    filter(CircRNA_ID == circ) %>%
    aov(Tool_RPM ~ Tool + Lane + Housing_Condition * Abstinence_Length, .) %>%
    emmeans(., list(pairwise ~ Housing_Condition * Abstinence_Length), 
            adjust = "tukey")
}

circ_tukey_noToolLane <- function(circ) {
  anova_format %>%
    filter(CircRNA_ID == circ) %>%
    aov(Mean_CircExpr ~ Housing_Condition * Abstinence_Length, .) %>%
    emmeans(., list(pairwise ~ Housing_Condition * Abstinence_Length), 
            adjust = "tukey")
}

ch1_tukey <- circ_tukey("2:24583888|24598599") # E21 vs. I21, 0.06 E1 vs. I21
cm1_tukey <- circ_tukey("3:17741317|17753437") # I1 vs. I21, E21 vs. I21, 0.07 E1 vs. E21
car1_tukey <- circ_tukey("5:145948264|145951212") # E1 vs. E21, E21 vs. I21, 0.08 I1 vs E21
cs3_tukey <- circ_tukey("3:57821407|57847776") # None
cah12_tukey <- circ_tukey("9:105629568|105660854") #  I1 vs. I21, E21 vs. I21


# Interpretations:
# circ-Homer1: Not quite significant. Seems like Enrichment effect; possible incubation effect with larger sample size
# circ-Mapkap1: Housing and Interaction; Incubation, Enrichment
# circ-Arid1a: Housing and Interaction; Enrichment
# circ-Sp3: Nothing
# circ-Ankrd12: Interaction; Incubation, Enrichment



## PLOT CIRCRNA EXPRESSION FOR 3 CIRCRNAS IN ALL (11) SAMPLES ----------
# First rearrange dataframe for plotting
# Create two versions, with and without outlier sample

# Outlier stays
keep_circs_4plot <- keep_circs %>%
  filter(CircRNA_ID %in% circs4) %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Mean_CircExpr = mean(c(CIRI2_RPM, CE2_RPM))) %>%
  mutate(Housing_Condition = 
           fct_relevel(Housing_Condition,
                       levels = c("Isolation", "Enrichment.")),
         Abstinence_Length = 
           fct_relevel(case_when(Abstinence_Length == 1 ~ "1d",
                                 Abstinence_Length == 21 ~ "21d"),
                       levels = c("1d", "21d"))) %>%
  mutate(CircRNA_ID = 
           ordered(
             case_when(CircRNA_ID == "2:24583888|24598599" ~ "circHomer1",
                       CircRNA_ID == "3:17741317|17753437" ~ "circMapkap1",
                       CircRNA_ID == "5:145948264|145951212" ~ "circArid1a",
                       CircRNA_ID == "3:57821407|57847776" ~ "circSp3",
                       CircRNA_ID == "9:105629568|105660854" ~ "circAnkrd12"),
             levels = c("circAnkrd12", "circArid1a", "circHomer1", 
                        "circMapkap1", "circSp3"))) %>%
  unique()

# Outlier removed
keep_circs_4plot_no_outlier <- keep_circs_4plot %>%
  filter(Sample_ID != "SRS6085269_CI21cue")

# Only outlier remains; All others are NA
keep_circs_4plot_outlier <- keep_circs_4plot %>%
  mutate(Mean_CircExpr = ifelse(Sample_ID == "SRS6085269_CI21cue",
                                Mean_CircExpr, NA_integer_),
         Cue_Reactivity_ALP = ifelse(Sample_ID == "SRS6085269_CI21cue",
                                     Mean_CircExpr, NA_integer_))


# Create plot
expr_plot <- keep_circs_4plot_no_outlier %>%
  filter(!is.na(CircRNA_ID)) %>%
  ggplot(aes(x = as.factor(Abstinence_Length), y = Mean_CircExpr, 
             fill = Housing_Condition)) +
  labs(x = "Abstinence Length", y = "Mean Expression (CPM)") +
  # Bars and errorbars for all points except outliers
  stat_summary(fun.data = "mean_se", geom = "bar", color = "black",
               position = position_dodge(width = 0.8), width = 0.8) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", size = 0.75, width = 0.3, 
               position = position_dodge(width = 0.8)) +
  # Show outlier points
  # geom_point(data = keep_circs_4plot_outlier, position = position_dodge(width = 0.9), 
  #            shape = 4, size = 3, stroke = 1.5, show.legend = FALSE) +
  # All other points
  geom_point(position = position_dodge(width = 0.8), size = 2.5, shape = 19, 
             show.legend = FALSE) +
  scale_fill_manual(name = "Housing", values = c("#F2B272", "#3F826D")) +
  scale_y_continuous(limits = c(0, 1.25), breaks = seq(0, 1.4, 0.2), expand = c(0.035, 0)) +
  theme(legend.position = c(0.905, 0.73),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey80"),
        panel.grid.minor.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9)) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  facet_wrap(~ CircRNA_ID, ncol = 5)

expr_plot

svglite("C:/Users/Annika/Documents/Figures for circRNA Project/Circ_expr.svg", width = 8.1, height = 3)
expr_plot
dev.off()




## CHECK FOR CORRELATIONS WITH CUE REACTIVITY ----------
# Must be present in at least 9 samples (none are significant)
keep_circs %>%
  filter(Sample_ID != "SRS6085269_CI21cue") %>%
  filter(Sample_Count >= 9) %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Mean_CircExpr = mean(c(CIRI2_RPM, CE2_RPM))) %>%
  group_by(CircRNA_ID) %>%
  cor_test(vars = c(Mean_CircExpr, Cue_Reactivity_ALP)) %>%
  arrange(p)

# What does the behavior look like for these animals specifically?
keep_circs_4plot_no_outlier %>%
  ggplot(aes(x = as.factor(Abstinence_Length), y = Cue_Reactivity_ALP, 
             fill = Housing_Condition)) +
  labs(x = "Abstinence Length", y = "Active Lever Presses") +
  # Bars and errorbars for all points except outliers
  stat_summary(fun.data = "mean_se", geom = "bar", 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", size = 1, width = 0.3, 
               position = position_dodge(width = 0.9)) +
  # Show outlier points
  geom_point(data = keep_circs_4plot_outlier, position = position_dodge(width = 0.9), 
             shape = 4, size = 3, stroke = 1.5, show.legend = FALSE) +
  # All other points
  geom_point(position = position_dodge(width = 0.9), size = 3, shape = 16) +
  scale_fill_manual(name = "Housing", values = c("#E64B35", "#4DBBD6")) +
  scale_y_continuous(expand = c(0.035, 0)) +
  theme(
    #legend.position = c(0.07, 0.82),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey80"),
        panel.grid.minor.y = element_blank(),
        axis.ticks.x = element_blank())
  


## CHECK IF ANY CIRCRNAS ARE EXPRESSED EXCLUSIVELY GROUPS/SAMPLES ----------
# CircRNAs present for both tools in 12, 11, 10, etc. samples
# For 9, 6, and 3 samples, check if samples are exclusive to 3, 2, or 1 group(s)
# Slightly different filtering - 0.05 RPM for both tools, otherwise eliminate
# individual samples

# Get number of samples per circRNA
circs_bysample <- keep_circs %>%
  group_by(CircRNA_ID) %>%
  mutate(Sample_Count = n()) %>%
  unique() %>%
  ungroup()

circs_bysample_nofilter <- circs_nofilter %>%
  group_by(CircRNA_ID) %>%
  mutate(Sample_Count = n()) %>%
  unique() %>%
  ungroup()


# Function for general searching by sample number
sample_search <- function(sample_num) {
  circs_bysample %>%
    filter(Sample_Count == sample_num) %>%
    pull(CircRNA_ID) %>%
    unique()
}

# Function for searching by both sample number and group number
group_search <- function(sample_num, group_num) {
  circs_bysample %>%
    filter(Sample_Count == sample_num) %>%
    dplyr::select(CircRNA_ID, Treatment_Group, Sample_Count) %>%
    unique() %>%
    arrange(CircRNA_ID) %>%
    group_by(CircRNA_ID, Treatment_Group) %>%
    mutate(Sample_TreatmentCount = n()) %>%
    group_by(CircRNA_ID) %>%
    mutate(Number_TreatmentGroups = n()) %>%
    filter(Number_TreatmentGroups == group_num) %>%
    pull(CircRNA_ID) %>%
    unique()
}

# Function for searching by both sample number and group number
group_search_nofilter <- function(sample_num, group_num) {
  circs_bysample_nofilter %>%
    filter(Sample_Count == sample_num) %>%
    dplyr::select(CircRNA_ID, Treatment_Group, Sample_Count) %>%
    unique() %>%
    arrange(CircRNA_ID) %>%
    group_by(CircRNA_ID, Treatment_Group) %>%
    mutate(Sample_TreatmentCount = n()) %>%
    group_by(CircRNA_ID) %>%
    mutate(Number_TreatmentGroups = n()) %>%
    filter(Number_TreatmentGroups == group_num) %>%
    pull(CircRNA_ID) %>%
    unique()
}

# Identified in both tools and expressed in 12, 11, 10, etc. samples,
# as a list of character vectors
exclusive_list <- list()

for (i in 1:12) {
  exclusive_list[[i]] <- sample_search(i)
}

seq_names <- seq(1, 12, 1)
names(exclusive_list) <- paste0("hits", seq_names)

# Identified by both tools and expressed exclusively in 3, 2, or 1 group(s)
# Account for the 1 outlier by looking for 3*groups and 3*groups - 1
# Verify those that are in 3*groups - 1 samples

# 3 groups
hits9_groups3 <- group_search(9, 3)
hits9_groups3 # 0
hits8_groups3 <- group_search(8, 3)
hits8_groups3 # 0

hits9_groups3_nofilter <- group_search_nofilter(9, 3)
hits9_groups3_nofilter # 0
hits8_groups3_nofilter <- group_search_nofilter(8, 3)
hits8_groups3_nofilter # 0

# 2 groups
hits6_groups2 <- group_search(6, 2)
hits6_groups2 # 0
hits5_groups2 <- group_search(5, 2)
hits5_groups2 
# [1] "10:30326084|30338951" "4:11695909|11711196"  "5:2849035|2899494" "9:48793273|48794864" 

# Check these for validity - Of those 4, none are valid (none are present in 
# the I21 group, which has only 2 samples left after removing the outlier)
keep_circs %>% 
  filter(CircRNA_ID %in% hits5_groups2) %>%
  select(CircRNA_ID, Treatment_Group) %>%
  unique()


# 1 group
hits3_groups1 <- group_search(3, 1)
hits3_groups1 # "11:50592225|50597061"  "16:6256270|6268891"    "5:97929048|97944520"
hits2_groups1 <- group_search(2, 1)
hits2_groups1
# [1] "1:2277718|2282190"**     "13:39756638|39758247" 
# [3] "17:5177078|5192799"    "2:115743198|115774378"
# [5] "2:194018644|194027413" "2:24414920|24426075"  
# [7] "3:33691214|33699525"   "3:57815992|57847776"  
# [9] "3:93186876|93272406"   "4:67192493|67204879"  
# [11] "5:150376677|150388680" "6:51793529|51801412"  
# [13] "8:62559509|62579508"**   "8:72517667|72534248"  
# [15] "9:61856348|61901553"   "9:61863815|61879316"  
# [17] "X:148094778|148106617" "X:9630483|9632250"  

# Check these for validity - only the ones marked with ** are exclusive to I21
keep_circs %>% 
  filter(CircRNA_ID %in% hits2_groups1) %>%
  select(CircRNA_ID, Treatment_Group) %>%
  unique()


# USE FCIRCSEC TO GET CIRCRNA SEQUENCES ----------
# Get rat genome from Ensembl (mRatBN7.2)
mRatBN7.2 <- useEnsembl("ensembl", dataset = "rnorvegicus_gene_ensembl")

# Create bed file for FcircSEC of all circRNAs detected w/ at least 2 reads in 
# any sample for both tools (i.e., keep_circs dataframe into file). Only need 
# the chromosome, start, and end, but needs to be 0-indexed like CE2.
# Also need to create a txt file with more information, like splices sites.

# Start position needs to be adjusted for FcircSEC
circs_4bed <- keep_circs %>%
  select(CircRNA_ID, CIRI2_CircRNA_Type) %>%
  unique() %>%
  separate(CircRNA_ID, into = c("Chromosome", "Start_End"), sep = ":") %>%
  separate(Start_End, into = c("Start", "End"), sep = "\\|") %>%
  mutate(Start = as.numeric(Start) - 1) %>%
  unite(Start_End, c("Start", "End"), sep = "|") %>%
  unite(CircRNA_ID, c("Chromosome", "Start_End"), sep = ":")

# Get exonic and intronic circRNAs separately
exonic_circs_4bed <- circs_4bed %>%
  filter(CIRI2_CircRNA_Type == "exon") %>%
  pull(CircRNA_ID)
intronic_circs_4bed <- circs_4bed %>%
  filter(CIRI2_CircRNA_Type == "intron") %>%
  pull(CircRNA_ID)

# Use CIRCexplorer2 to get the other data
# Temporary dataframe with most of the necessary information  
temp_circ_df <- ce2_original %>%
  unite(Start_End_Dash, c("Start", "End"), sep = "-", remove = FALSE) %>%
  unite(Start_End_Pipe, c("Start", "End"), sep = "|", remove = FALSE) %>%
  unite(CircRNA_ID_Dash, c("Chromosome", "Start_End_Dash"), 
        sep = ":", remove = FALSE) %>%
  unite(CircRNA_ID_Pipe, c("Chromosome", "Start_End_Pipe"), 
        sep = ":", remove = FALSE) %>%
  filter(CircRNA_ID_Pipe %in% exonic_circs_4bed | 
           CircRNA_ID_Pipe %in% intronic_circs_4bed) %>%
  mutate(CircRNA_Type = case_when(
    CircRNA_ID_Pipe %in% exonic_circs_4bed ~ "exonic",
    CircRNA_ID_Pipe %in% intronic_circs_4bed ~ "other")) %>%
  mutate(ExonSizes2 = str_replace_all(ExonSizes, ",", "+"))

# Calculate splice length based on exon sizes for exonic circRNAs
exonsizes <- temp_circ_df %>% pull(ExonSizes2)
exp <- lapply(exonsizes, FUN = function(x) parse(text = x)[[1]])
eval <- lapply(exp, eval)
eval_unlist <- unlist(eval)

# Add splice length
temp_circ_df2 <- temp_circ_df %>%
  mutate(Splice_Length = eval_unlist) %>%
  select(CircRNA_ID_Dash, Chromosome, Start, End, Strand, Splice_Length,
         CircRNA_Type, ExonCount, ExonSizes, ExonOffsets, IsoformName) %>%
  setNames(c("ID", "chr", "circ_start", "circ_end", "circ_strand", "splice_L",
             "circ_type", "e_count", "e_sizes", "e_offsets", "b_transcript")) %>%
  unique()

# Pull transcripts
b_transcript <- temp_circ_df2 %>% pull(b_transcript)

# Get dataframe with relevant transcript info from biomaRt
# Then join to create full dataframe
circ_df_4txt <- getBM(attributes = c("ensembl_transcript_id",
                                     "strand", "transcript_start",
                                     "transcript_end", "external_gene_name"),
                      filter = "ensembl_transcript_id", 
                      values = b_transcript,
                      mart = mRatBN7.2) %>%
  mutate(across(everything(), function(x) na_if(x, "")),
         strand = ifelse(str_detect(strand, "-1"), "-", "+")) %>%
  setNames(c("b_transcript", "b_strand", "b_trans_start", 
             "b_trans_end", "b_gene")) %>%
  full_join(temp_circ_df2, .)

# Reformat for BED file
circ_df_4bed <- circ_df_4txt %>%
  select(chr, circ_start, circ_end) %>%
  mutate(across(c(2:3), as.numeric))

# Write files for FcircSEC to use
write.table(circ_df_4txt, "G:/Where the Big Files Are/Rat Reference Files/circ_class_1_23_2022.txt", 
            row.names = FALSE, sep = "\t", quote = FALSE)
write.table(circ_df_4bed, "G:/Where the Big Files Are/Rat Reference Files/circ_class_1_23_2022.bed", 
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# START FCIRCSEC
# General input files
ref_genome <- file.path(in_dir, "Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa")
annotation <- file.path(in_dir, "Rattus_norvegicus.mRatBN7.2.105.gtf")
circ_class_bed <- file.path(in_dir, "final_circ_class_1_23_2022.bed")
circ_class_txt <- file.path(in_dir, "final_circ_class_1_23_2022.txt")

# Output file from get.fasta()
circ_fasta <- file.path(out_dir, "circs_from_getfasta_1_23_2022.fa")

# Output file from circSeqExt()
circ_seq <- file.path(out_dir, "circ_seq_1_23_2022.txt")

# Extract transcript information
transcriptExtract(annotation, "other", t_data) # Takes several minutes

# Get FASTAs of transcripts, then full circRNA sequences
get.fasta(ref_genome, circ_class_bed, circ_fasta)
circSeqExt(circ_fasta, circ_class_txt, circ_seq)

# The 2 necessary files (fasta and sequences) should now be written in the 
# correct directory.



## OTHER
more_miranda_info <- read.delim("C:/Users/Annika/Documents/circrna_full_miranda_info_interest.txt")

cah12_conf <- circ_mir_conf %>% filter(circRNA_ID == "9:105629567|105660854") %>% pull(miRNA)
car1_conf <- circ_mir_conf %>% filter(circRNA_ID == "5:145948263|145951212") %>% pull(miRNA)
ch1_conf <- circ_mir_conf %>% filter(circRNA_ID == "2:24583887|24598599") %>% pull(miRNA)
cm1_conf <- circ_mir_conf %>% filter(circRNA_ID == "3:17741316|17753437") %>% pull(miRNA)
cs3_conf <- circ_mir_conf %>% filter(circRNA_ID == "3:57821406|57847776") %>% pull(miRNA)

more_miranda_info %>%
  mutate(miRNA = str_replace(miRNA, ">", "")) %>%
  filter((circRNA == "9:105629567-105660854" & miRNA %in% cah12_conf) |
         (circRNA == "5:145948263-145951212" & miRNA %in% car1_conf) |
         (circRNA == "2:24583887-24598599" & miRNA %in% ch1_conf) |
         (circRNA == "3:17741316-17753437" & miRNA %in% cm1_conf) |
         (circRNA == "3:57821406-57847776" & miRNA %in% cs3_conf)) %>%
  write_clip()










