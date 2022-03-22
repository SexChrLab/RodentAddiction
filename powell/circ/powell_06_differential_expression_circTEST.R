
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



# BEFORE USING THIS SCRIPT ----------
# User should alter these variables as necessary

# Home PCd
main_dir <- "D:/GitHub Repositories/temp/"
pheno_dir <- paste0(main_dir, "powell/")
gene_dir <- paste0(main_dir, "powell/genes/")
circ_dir <- paste0(main_dir, "powell/circ")

# Update and set ggplot theme globally
new_theme <- theme_bw(base_size = 14) %+replace%
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
# Mapped reads quantified with samtools as follows:
# "samtools view -c -F 4 {input.SAM}"
bwa_maps <- data.frame(Sample_ID = sample_ids,
                       Mapped_Reads = c(36503658, 40498013, 27725648, 61025798, 
                                        43444856, 40516853, 42587934, 38625093, 
                                        32042541, 48159370, 50153018, 52294417))
# Check mean mapped reads
mean_maps <- mean(bwa_maps$Mapped_Reads)
mean_maps

Mapped_Reads = c(x, x, x, x, 
                 x, x, x, x, 
                 x, x, x, x))

# Merge dataframes and calculate RPM for each circRNA
# RPM = JunctionReads/(Mapped_Reads/1000000)
all_circs <- Reduce(function(x, y)
  full_join(x, y), list(ciri2, ce2, bwa_maps)) %>%
  mutate(CIRI2_RPM = CIRI2_Junction_Reads/(Mapped_Reads/1000000),
         CE2_RPM = CE2_Junction_Reads/(Mapped_Reads/1000000)) %>%
  full_join(pheno, by = c("Sample_ID" = "Sample_Treatment")) %>%
  dplyr::select(c(1:12, 20:26)) 



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

count_plot <- all_circs %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Sample_ID = factor(Sample_ID,
                            levels = c("SRS6085261_CI1cue", "SRS6085267_CI1cue", "SRS6085268_CI1cue",
                                       "SRS6085263_CE1cue", "SRS6085266_CE1cue", "SRS6085270_CE1cue",
                                       "SRS6085262_CI21cue", "SRS6085269_CI21cue", "SRS6085272_CI21cue",
                                       "SRS6085264_CE21cue", "SRS6085265_CE21cue", "SRS6085271_CE21cue"))) %>%
  mutate(Sample_ID = case_when(Sample_ID == "SRS6085261_CI1cue" ~ "Iso, 1d - A",
                               Sample_ID == "SRS6085267_CI1cue" ~ "Iso, 1d - B",
                               Sample_ID == "SRS6085268_CI1cue" ~ "Iso, 1d - C",
                               Sample_ID == "SRS6085263_CE1cue" ~ "EE, 1d - A",
                               Sample_ID == "SRS6085266_CE1cue" ~ "EE, 1d - B",
                               Sample_ID == "SRS6085270_CE1cue" ~ "EE, 1d - C",
                               Sample_ID == "SRS6085262_CI21cue" ~ "Iso, 21d - A",
                               Sample_ID == "SRS6085269_CI21cue" ~ "Iso, 21d - B",
                               Sample_ID == "SRS6085272_CI21cue" ~ "Iso, 21d - C",
                               Sample_ID == "SRS6085264_CE21cue" ~ "EE, 21d - A",
                               Sample_ID == "SRS6085265_CE21cue" ~ "EE, 21d - B",
                               Sample_ID == "SRS6085271_CE21cue" ~ "EE, 21d - C")) %>%
  mutate(Identified_By = factor(case_when(!is.na(CIRI2_Junction_Reads) & !is.na(CE2_Junction_Reads) ~ "Both",
                                   !is.na(CIRI2_Junction_Reads) & is.na(CE2_Junction_Reads) ~ "CIRI2",
                                   is.na(CIRI2_Junction_Reads) & !is.na(CE2_Junction_Reads) ~ "CE2"),
                                levels = c("CIRI2", "CE2", "Both"))) %>%
  pivot_longer(cols = c(CIRI2_Junction_Reads, CE2_Junction_Reads),
               names_to = "Tool", values_to = "Count") %>%
  filter(!is.na(Count)) %>%
  ggplot(aes(x = Count, fill = Identified_By)) +
  geom_histogram(binwidth = 1, color = "black") +
  labs(x = "Read Count", y = "Number of circRNAs Identified") +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  scale_fill_manual(name = "Tool", values = c("white", "grey80", "grey30")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color = "grey80")) +
  geom_vline(xintercept = 1.56, lty = "dashed", color = "red") +
  facet_wrap(~ Sample_ID, ncol = 3, scales = "free_y")

count_plot
pal <- c(wes_palette("Chevalier1", 4, type = "discrete"))
  
count_gtable <- ggplot_gtable(ggplot_build(count_plot))
striprt <- which(grepl('strip-r', count_gtable$layout$name) | grepl('strip-t', count_gtable$layout$name))
fills <- rep(c(rep("#446455", 3), rep("#FDD262", 3), rep("#D3DDDC", 3), rep("#C7B19C", 3)))
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', count_gtable$grobs[[i]]$grobs[[1]]$childrenOrder))
  count_gtable$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
  
grid::grid.draw(count_gtable)

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
  full_join(x, y), list(ciri2_before, ce2_before, all_after)) %>%
  arrange(After_Circs)
filter_summary






# Rat (Rn_6.0)
rn6 <- useEnsembl("ensembl", dataset = "rnorvegicus_gene_ensembl", 
                  host = "https://may2021.archive.ensembl.org")



# Create bed file for FcircSEC using any and all identified circs
# Start position needs to be adjusted for FcircSEC
circs_4bed <- keep_circs %>%
  select(CircRNA_ID, CIRI2_CircRNA_Type) %>%
  unique() %>%
  separate(CircRNA_ID, into = c("Chromosome", "Start_End"), sep = ":") %>%
  separate(Start_End, into = c("Start", "End"), sep = "\\|") %>%
  mutate(Start = as.numeric(Start) - 1) %>%
  unite(Start_End, c("Start", "End"), sep = "|") %>%
  unite(CircRNA_ID, c("Chromosome", "Start_End"), sep = ":")

exonic_circs_4bed <- circs_4bed %>%
  filter(CIRI2_CircRNA_Type == "exon") %>%
  pull(CircRNA_ID)

intronic_circs_4bed <- circs_4bed %>%
  filter(CIRI2_CircRNA_Type == "intron") %>%
  pull(CircRNA_ID)

# Temp dataframe with most of the info  
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
exonsizes <- test %>% pull(ExonSizes2)
exp <- lapply(exonsizes, FUN = function(x) parse(text = x)[[1]])
eval <- lapply(exp, eval)
eval_unlist <- unlist(eval)

# Add splice length
temp_circ_df <- temp_circ_df %>%
  mutate(Splice_Length = eval_unlist) %>%
  select(CircRNA_ID_Dash, Chromosome, Start, End, Strand, Splice_Length,
         CircRNA_Type, ExonCount, ExonSizes, ExonOffsets, IsoformName) %>%
  setNames(c("ID", "chr", "circ_start", "circ_end", "circ_strand", "splice_L",
             "circ_type", "e_count", "e_sizes", "e_offsets", "b_transcript")) %>%
  unique()
  
# Get rat mart (Rn_6.0)
rn6 <- useEnsembl("ensembl", dataset = "rnorvegicus_gene_ensembl", 
                  host = "https://may2021.archive.ensembl.org")

# Pull transcripts
b_transcript <- temp_circ_df %>% pull(b_transcript)

# Get dataframe with relevant transcript info from biomaRt
# Then join to create full dataframe
circ_df_4txt <- getBM(attributes = c("ensembl_transcript_id",
                                       "strand", "transcript_start",
                                       "transcript_end", "external_gene_name"),
                        filter = "ensembl_transcript_id", 
                        values = b_transcript,
                        mart = rn6) %>%
  mutate(across(everything(), function(x) na_if(x, "")),
         strand = ifelse(str_detect(strand, "-1"), "-", "+")) %>%
  setNames(c("b_transcript", "b_strand", "b_trans_start", 
             "b_trans_end", "b_gene")) %>%
  full_join(temp_circ_df, .)

circ_df_4bed <- circ_df_4txt %>%
  select(chr, circ_start, circ_end) %>%
  mutate(across(c(2:3), as.numeric))

write.table(circ_df_4txt, "G:/Where the Big Files Are/Rat Reference Files/final_circ_class.txt", row.names = FALSE, sep = "\t")
write.table(circ_df_4bed, "G:/Where the Big Files Are/Rat Reference Files/final_circ_class.bed", row.names = FALSE, col.names = FALSE, sep = "\t")














# Which circRNAs are still present in all 4 groups?
circs4 <- keep_circs %>%
  group_by(CircRNA_ID) %>%
  mutate(Sample_Count = n()) %>%
  filter(Sample_Count == 12) %>%
  pull(CircRNA_ID) %>%
  unique()
circs4 # "2:22950023|22964734", "3:13510410|13522529", "5:151944769|151947717"

# 3 circRNAs present in all 4 groups
circs4_df <- data.frame(CircRNA_ID = circs4,
                        Gene_ID = c("ENSRNOG00000047014", "ENSRNOG00000017583",
                                    "ENSRNOG00000006137"),
                        Gene = c("Homer1", "Mapkap1", "Arid1a"))

# One sample is weird
# Check if there are circRNAs present in all samples except that one -
# It's still the same 3 circRNAs.
circs4_v2 <- keep_circs %>%
  filter(Sample_ID != "SRS6085269_CI21cue") %>%
  group_by(CircRNA_ID) %>%
  mutate(Sample_Count = n()) %>%
  filter(Sample_Count == 11) %>%
  pull(CircRNA_ID) %>%
  unique()
circs4_v2 # "2:22950023|22964734", "3:13510410|13522529", "5:151944769|151947717"

# Reformat data so that the RPM from each tool can be entered into the model
# separately; Remove outlier
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

is.ordered(anova_format$Tool)
is.ordered(anova_format$Lane)
is.ordered(anova_format$Housing_Condition)
is.ordered(anova_format$Abstinence_Length)

# Primary ANOVA
options(contrasts=c('contr.sum', 'contr.poly'))
anova_format %>%
  filter(CircRNA_ID == circs4_v2[1]) %>%
  # aov(Tool_RPM ~ Tool + Lane + Housing_Condition * Abstinence_Length, data = .) %>%
  # emmeans(., list(pairwise ~ Housing_Condition * Abstinence_Length), adjust = "tukey")
  aov(Tool_RPM ~ Tool + Lane + Housing_Condition + Abstinence_Length, .) %>%
  anova(.) %>%
  as.data.frame(.)


# Primary ANOVA
# This is correct and matches SPSS results
anova_format %>%
  filter(CircRNA_ID == circs4_v2[1]) %>%
  lm(Tool_RPM ~ Tool + Lane + Housing_Condition * Abstinence_Length, .) %>%
  car::Anova(., type = "III")



anova_format %>%
  filter(CircRNA_ID == circs4_v2[1]) %>%
  lm(Tool_RPM ~ Tool + Lane + Housing_Condition * Abstinence_Length, .) %>%
  aov(.)
  summary(.)
  as.data.frame(.)
  

anova_format %>%
  filter(CircRNA_ID == circs4_v2[1]) %>%
  aov(Tool_RPM ~ Housing_Condition * Abstinence_Length, .) %>%
  anova(.) %>%
  as.data.frame(.)


anova_format %>%
  filter(CircRNA_ID == circs4_v2[1]) %>%
  anova_test(Tool_RPM ~ Housing_Condition * Abstinence_Length, type = "III") %>%
  as.data.frame(.)


# "2:22950023|22964734" (Mapkap1)
# "3:13510410|13522529" (Homer1)
# "5:151944769|151947717" (Arid1a)

# Post-hoc Tukey tests
anova_format %>%
  filter(CircRNA_ID == circs4_v2[1]) %>%
  aov(Tool_RPM ~ Tool + Lane + Housing_Condition * Abstinence_Length, data = .) %>%
  summary()
  # emmeans(., list(pairwise ~ Housing_Condition * Abstinence_Length), adjust = "tukey")

# Shapiro-Wilk tests for each circRNA
# Test is passed for all 3 circRNAs
anova_format %>%
  filter(CircRNA_ID == circs4_v2[1]) %>%
  lm(Tool_RPM ~ Tool + Lane + Housing_Condition * Abstinence_Length, .) %>%
  residuals(.) %>%
  shapiro_test(.)
anova_format %>%
  filter(CircRNA_ID == circs4_v2[2]) %>%
  lm(Tool_RPM ~ Tool + Lane + Housing_Condition * Abstinence_Length, .) %>%
  residuals(.) %>%
  shapiro_test(.)
anova_format %>%
  filter(CircRNA_ID == circs4_v2[3]) %>%
  lm(Tool_RPM ~ Tool + Lane + Housing_Condition * Abstinence_Length, .) %>%
  residuals(.) %>%
  shapiro_test(.)





# 2:22950023|22964734 (Homer1)
# 3:13510410|1352252 (Mapkap1)
# 5:151944769|151947717 (Arid1a)

# Correlations with cue reactivity
# Must be present in at least 9 samples (none are significant)
keep_circs %>%
#  filter(Sample_ID != "SRS6085269_CI21cue") %>%
  filter(Sample_Count >= 9) %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Mean_CircExpr = mean(c(CIRI2_RPM, CE2_RPM))) %>%
  group_by(CircRNA_ID) %>%
  cor_test(vars = c(Mean_CircExpr, Cue_Reactivity_ALP)) %>%
  arrange(p)

keep_circs %>%
  filter(Sample_Count >= 9) %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Mean_CircExpr = mean(c(CIRI2_RPM, CE2_RPM))) %>%
  group_by(CircRNA_ID) %>%
  cor_test(vars = c(Mean_CircExpr, Cue_Reactivity_ALP), method = "spearman") %>%
  arrange(p)

# Plot of 1 circRNA that is trending
keep_circs %>%
  filter(CircRNA_ID == "4:17032161|17033623") %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Mean_CircExpr = mean(c(CIRI2_RPM, CE2_RPM))) %>%
  ggplot(aes(x = Mean_CircExpr, y = Cue_Reactivity_ALP)) +
  geom_point()


keep_circs %>%
  filter(Sample_ID != "SRS6085269_CI21cue") %>%
  filter(CircRNA_ID %in% circs4) %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Mean_CircExpr = mean(c(CIRI2_RPM, CE2_RPM))) %>%
  mutate(Housing_Condition = fct_relevel(Housing_Condition,
                                         levels = c("Isolation", "Enrichment")),
         Abstinence_Length = fct_relevel(case_when(Abstinence_Length == 1 ~ "1d",
                                                   Abstinence_Length == 21 ~ "21d"),
                                         levels = c("1d", "21d"))) %>%
  mutate(CircRNA_ID = case_when(CircRNA_ID == "2:22950023|22964734" ~ "circ-Homer1",
                                CircRNA_ID == "3:13510410|13522529" ~ "circ-Mapkap1",
                                CircRNA_ID == "5:106171074|106171781" ~ "circ-Mltt3",
                                CircRNA_ID == "5:151944769|151947717" ~ "circ-Arid1a")) %>%
  unique() %>%
  ggplot(aes(x = as.factor(Abstinence_Length), y = Mean_CircExpr, fill = Housing_Condition)) +
  labs(x = "Abstinence Length", y = "Mean Counts per Million") +
  stat_summary(fun.data = "mean_se", geom = "bar", 
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", size = 1, width = 0.3, 
               position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9), size = 3, shape = 16) +
  scale_fill_manual(name = "Housing", values = c("#E64B35", "#4DBBD6")) +
  theme(legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey80"),
        panel.grid.minor.y = element_blank()) +
  facet_wrap(~ CircRNA_ID, scales = "free_y")





## CHECK IF ANY CIRCRNAS ARE EXPRESSED EXCLUSIVELY GROUPS/SAMPLES ----------
# CircRNAs present for both tools in 12, 11, 10, etc. samples
# For 9, 6, and 3 samples, check if samples are exclusive to 3, 2, or 1 group(s)
# Slightly different filtering - 0.05 RPM for both tools, otherwise eliminate
# individual samples

# Get number of samples per circRNA
circs_bysample <- keep_circs %>%
  filter(CIRI2_RPM > 0.05, CE2_RPM > 0.05) %>%
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
    select(CircRNA_ID, Treatment_Group, Sample_Count) %>%
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
hits9_groups3 <- group_search(9, 3)
# [1] "10:16724062|16730140"  "8:130762665|130768360"
# "10:16724062|16730140" (rno-Crebrf_0009; Creb3 regulatory factor) - NOT found in I21
# "8:130762665|130768360" (rno-Snrk_0003) - NOT found in I21
# These are still good with no expression filter

hits6_groups2 <- group_search(6, 2)
# [1] "10:29318809|29319752" "8:508994|509712" 
# "10:29318809|29319752" (rno-Pwwp2a_0002) - I1 and E21 only
# "8:508994|509712" (rno-Gucy1a2_0008) - I1 and E1 only 
# These are still good with no expression filter

hits3_groups1 <- group_search(3, 1)
# [1] "1:192952650|192955267" "2:189723524|189729005" "4:8223840|8239403" "5:101648419|101663967"
# "1:192952650|192955267" (rno-Rbbp6_0006) - Exclusive to E21
# "2:189723524|189729005" (rno-Gatad2b_0003) - Exclusive to I1
# "4:8223840|8239403" (rno-Kmt2e_0001) - Exclusive to I1
# "5:101648419|101663967" (rno-AABR07049060_0016) - Exclusive to I1
# These are still good with no expression filter

test <- append(append(hits9_groups3, hits6_groups2), hits3_groups1)

keep_circs %>%
  filter(Sample_ID != "SRS6085269_CI21cue") %>%
  filter(CircRNA_ID %in% append(hits9_groups3, hits6_groups2) | Sample_ID == "SRS6085262_CI21cue", 
         CircRNA_ID != "2:22950023|22964734") %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Mean_CircExpr = mean(c(CIRI2_RPM, CE2_RPM)),
         logMean = log1p(Mean_CircExpr)) %>%
  mutate(Housing_Condition = fct_relevel(Housing_Condition,
                                         levels = c("Isolation", "Enrichment"))) %>%
  mutate(Mean_CircExpr = case_when(Sample_ID == "" ~ 0, TRUE ~ Mean_CircExpr),
         CircRNA_ID = case_when(CircRNA_ID == "5:106171074|106171781" ~ "10:16724062|16730140",
                                TRUE ~ CircRNA_ID)) %>%
  # mutate(CircRNA_ID = case_when(CircRNA_ID == "2:22950023|22964734" ~ "circ-Homer1",
  #                               CircRNA_ID == "3:13510410|13522529" ~ "circ-Mapkap1",
  #                               CircRNA_ID == "5:106171074|106171781" ~ "circ-Mltt3",
  #                               CircRNA_ID == "5:151944769|151947717" ~ "circ-Arid1a")) %>%
  unique() %>%
  ggplot(aes(x = as.factor(Abstinence_Length), y = Mean_CircExpr, fill = Housing_Condition, shape = Treatment_Group)) +
  labs(x = "Abstinence Length", y = "Mean Expression in Counts per Million (CPM)\nacross Both Tools") +
  geom_boxplot(width = 0.7, position = position_dodge(preserve = "single")) +
#  geom_boxplot(width = 0.7) +
  geom_point(position = position_dodge(width = 0.7)) +
  theme(legend.position = "bottom") +
  facet_wrap(~ CircRNA_ID, scales = "free_y")



# ANOVAs and t-tests for circRNAs expressed in only 2 or 3 groups
# ANOVAs for 3 - not significant
keep_circs %>%
  filter(CircRNA_ID %in% hits9_groups3) %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Mean_CircExpr = mean(c(CIRI2_RPM, CE2_RPM))) %>%
  group_by(CircRNA_ID) %>%
  anova_test(Mean_CircExpr ~ Treatment_Group)
# K-Ws for 3 - not significant
keep_circs %>%
  filter(CircRNA_ID %in% hits9_groups3) %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Mean_CircExpr = mean(c(CIRI2_RPM, CE2_RPM))) %>%
  group_by(CircRNA_ID) %>%
  kruskal_test(Mean_CircExpr ~ Treatment_Group)
# T-tests for 2 - not significant
keep_circs %>%
  filter(CircRNA_ID %in% hits6_groups2) %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Mean_CircExpr = mean(c(CIRI2_RPM, CE2_RPM))) %>%
  group_by(CircRNA_ID) %>%
  t_test(Mean_CircExpr ~ Treatment_Group)
# Mann-Whitney U-tests for 2 - not significant
keep_circs %>%
  filter(CircRNA_ID %in% hits6_groups2) %>%
  group_by(Sample_ID, CircRNA_ID) %>%
  mutate(Mean_CircExpr = mean(c(CIRI2_RPM, CE2_RPM))) %>%
  group_by(CircRNA_ID) %>%
  wilcox_test(Mean_CircExpr ~ Treatment_Group)


all_circs %>%
  mutate(CIRI2_CircRNA_Type = as.factor(CIRI2_CircRNA_Type),
         CE2_CircRNA_Type = as.factor(CE2_CircRNA_Type)) %>%
  select(CIRI2_CircRNA_Type, CE2_CircRNA_Type) %>%
  summary()







>rno-Homer1_0004
CUAUGUGAAAAUGGCAAUGCAUUCUGAGCUGGCUCAGCCCUUGGCUCUGAGUUCUGUGUCACAUCGGGUGUUCUCUCAUCAUCUGUCCCAUUGAUACUUUCUGGUGUUAAAGGAGACUGAAGAUCUCCUCCUGCUGAUUCCUGUGAAGGGGUACUGGUCAGUUCCAUCUUCUCCUGCGACUUCUCCUUUGCCAGCCGAGCAGCUUCUUUAAAUUCCUGAAACUUUUCUGCAAAUUUUGAGAGAUGAUGCUCAGAGGAGAAUCCCAGUCCAUAAACAGUGUUUGCCCGGCUAUCAGCCCAUUGGCCAAACUUUUGAGAUGUUUUAGUAAAUGUCAUGUUUGGAGUGAUGGUGCUAUUUAUUAUUGCCUUUGAGCCGUCUAGACUGAUUAUCCUAUACACAUUCCUUGUGCUGUCAUAGAAAUAAGACACAGUAACUGCAUGCUUGCUGGUGGGUACCCAGUUCUUCUUUGUGUUUGGGUCGAUCUGGAAGACAUGAGCUCGAGUGCUGAAGAUAGGUUGUUCC








###########
library(FcircSEC)
out_dir <- "G:/Where the Big Files Are/Rat Reference Files/"

# General input files
ref_genome <- file.path(out_dir, "Rn6_genome.fa")

# Output files from circClassification()
circ_class_bed <- file.path(out_dir, "final_circ_class.bed")
circ_class_txt <- file.path(out_dir, "final_circ_class.txt")

# Output file from get.fasta()
circ_fasta <- file.path(out_dir, "circs_from_getfasta.fa")

# Output file from circSeqExt()
circ_seq <- file.path(out_dir, "circ_seq.txt")

# Get FASTAs of transcripts, then full circRNA sequences
get.fasta(ref_genome, circ_class_bed, circ_fasta)
circSeqExt(circ_fasta, circ_class_txt, circ_seq)

############
Should find accurate sequences for the following, since they have Ns:
>20:27409749-27415954 - single
>4:17032160-17053486 - ~18
>4:17032160-17033623
>6:29032001-29074704
>1:166231692-166259335



# Extract transcript information
transcriptExtract(annotation, "other", t_data)

# Transcript data has some errors in gene column; fix by removing "gene_id " string
t_data <- read.delim(file.path(out_dir, "extracted_transcripts.txt"))
t_data2 <- t_data %>% mutate(gene = str_replace(gene, "gene_id ", ""))
# Save modified version of file
write.table(t_data2, file.path(out_dir, "extracted_transcripts2.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
# Set new variable
t_data2 <- file.path(out_dir, "extracted_transcripts2.txt")

# Get circRNA classification information
circClassification(t_data2, circ_bed, circ_class_txt, circ_class_bed)

# Rewrite text file to get rid of erroneous columns
circ_class_txt <- read.delim(file.path(out_dir, "circ_class.txt"))
circ_class_txt2 <- circ_class_txt[-1, ]
write.table(circ_class_txt2, file.path(out_dir, "circ_class2.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
circ_class_txt2 <- file.path(out_dir, "circ_class2.txt")

# Rewrite bed file to get rid of headers
circ_class_bed <- read.delim(file.path(out_dir, "circ_class.bed"), header = FALSE)
circ_class_bed2 <- circ_class_bed[-1, ] %>% mutate(across(2:3, as.numeric))
write.table(circ_class_bed2, file.path(out_dir, "circ_class2.bed"), col.names = FALSE, row.names = FALSE, sep = "\t")
circ_class_bed2 <- file.path(out_dir, "circ_class2.bed")

# Get FASTAs of transcripts, then full circRNA sequences
get.fasta(ref_genome, circ_class_bed2, circ_fasta)
circSeqExt(circ_fasta, circ_class_txt2, circ_seq)

test <- readDNAStringSet(file.path(out_dir, "circ_seq.txt"))
circRNA_ids <- names(test) %>% head()
circRNA_seqs <- paste(head(test))
test2 <- data.frame(circRNA_ids, circRNA_seqs)



>2:22950023-22964734
GAACAACCTATCTTCAGCACTCGAGCTCATGTCTTCCAGATCGACCCAAACACAAAGAAGAACTGGGTACCCACCAGCAA
GCATGCAGTTACTGTGTCTTATTTCTATGACAGCACAAGGAATGTGTATAGGATAATCAGTCTAGACGGCTCAAAGGCAA
TAATAAATAGCACCATCACTCCAAACATGACATTTACTAAAACATCTCAAAAGTTTGGCCAATGGGCTGATAGCCGGGCA
AACACTGTTTATGGACTGGGATTCTCCTCTGAGCATCATCTCTCAAAATTTGCAGAAAAGTTTCAGGAATTTAAAGAAGC
TGCTCGGCTGGCAAAGGAGAAGTCGCAGGAGAAGATGGAACTGACCAGTACCCCTTCACAGGAATCAGCAGGAGGAGATC
TTCAGTCTCCTTTAACACCAGAAAGTATCAATGGGACAGATGATGAGAGAACACCCGATGTGACACAGAACTCAGAGCCA
AGGGCTGAGCCAGCTCAGAATGCATTGCCATTTTCACATAG


test <- readDNAStringSet(file.path(out_dir, "circ_seq.txt"))
as.data.frame(test)

seq(1:22)
setwd("G:/Where the Big Files Are/Rat Reference Files/")
annotation_file <- data(refGenchr1)
annotation_file <- refGenchr1
write.table(annotation_file, "annotation_file.gtf", row.names = FALSE,
            sep = "/t", quote = FALSE, col.names = FALSE)


circClassification(transcript_data, file.path(out_dir, "test_circ_bed.bed"), file.path(out_dir, "test_circ_txt.txt"), file.path(out_dir, "test_circ_class_bed.bed"))



# Extract transcript information
transcriptExtract(refGenchr1, "other", t_data)

data("refGenchr1")
write.table(refGenchr1, file.path(out_dir, "refGenchr1.gtf"), row.names = FALSE, sep = "\t", quote = FALSE)
refGenchr1 <- file.path(out_dir, "refGenchr1.gtf")

# Extract transcript information (FAILS)
transcriptExtract(refGenchr1, "other", file.path(out_dir, "test_transcript_data.txt"))

data("transcript_data")
write.table(transcript_data, file.path(out_dir, "test_transcript_data.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
transcript_data <- file.path(out_dir, "test_transcript_data.txt")
data("output_CIRI")
write.table(output_CIRI, file.path(out_dir, "output_CIRI.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
output_CIRI <- file.path(out_dir, "output_CIRI.txt")

# circRNA classification (SUCCESS W/ ERROR MESSAGE - NAs INTRODUCED BY COERCION)
circClassification(transcript_data, output_CIRI, 
                   file.path(out_dir, "test_circ_txt.txt"), 
                   file.path(out_dir, "test_circ_class_bed.bed"))



class_bed <- file.path(out_dir, "test_circ_class_bed.bed")
class_txt <- file.path(out_dir, "test_circ_txt.txt")



data("chr1")
dataframe2fas(chr1, file.path(out_dir, "chr1.fa"))
chr1 <- file.path(out_dir, "chr1.fa")

data("circRNA_classb")
class_bed <- circRNA_classb
# write.table(class_bed, file.path(out_dir, "test_circ_class_bed.bed"), col.names = FALSE, row.names = FALSE)
class_bed <- file.path(out_dir, "test_circ_class_bed.bed")

data("circRNA_classt")
class_txt <- circRNA_classt
write.table(class_txt, file.path(out_dir, "test_circ_txt.txt"), row.names = FALSE)
class_txt <- file.path(out_dir, "test_circ_txt.txt")


# Get FASTAs (SUCCESS)
get.fasta(chr1, class_bed, file.path(out_dir, "circ_fasta.fa"))
circ_fasta <- file.path(out_dir, "circ_fasta.fa")

# Get sequences (SUCCESS)
circSeqExt(circ_fasta, class_txt, file.path(out_dir, "test_circ_seq.txt"))


data("circRNA_genomic_sequence")

