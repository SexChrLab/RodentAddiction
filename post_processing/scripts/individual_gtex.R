## LOAD PACKAGES ----------
# Data manipulation
library(tidyverse)
library(rlist)
library(readxl)

# Plots and visualizations
library(grid)   
library(viridis)
library(ggrepel)
library(ggpubr)
library(svglite)

# Statistics
library(rstatix)
# library(emmeans)
# library(multcomp)
# library(multcompView)

# Matching samples
library(MatchIt)

# Other
library(clipr)



## BEFORE USING THIS SCRIPT ----------
# User should alter these variables as necessary
main_dir <- "G:/Annika/Github Repositories/RodentAddiction/"
de_dir <- paste0(main_dir, "post_processing/results/gene_info/")
gtex_dir <- paste0(main_dir, "post_processing/downloaded_data/gtex/")

# Reduce chances of scientific notation
options(scipen = 9999)

# Update and set ggplot theme globally
new_theme <- theme_bw(base_size = 14) %+replace%
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))
theme_set(new_theme)

# Set dplyr functions
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename



## LOAD IN GENE LISTS AND CONSERVATION ANALYSIS ----------
# All in one list
all_genes_list <- readRDS(file = paste0(main_dir, "post_processing/results/other/all_genes_list.RDS"))

# Conservation analysis
# This is only necessary to load in for creating a combined dataframe with expression data
full_cons_df <- readRDS(file = paste0(main_dir, "post_processing/results/other/full_cons_df.RDS"))

# Separate list into multiple objects
list2env(all_genes_list, envir = .GlobalEnv)



## LOAD IN GTEX DATA (OVERALL MEANS/MEDIANS) -----------
# Get list of files
gtex_all_files <- list.files(path = paste0(gtex_dir, "counts_mean_median"),
                             pattern = "*_overall_counts_mean_median.txt", 
                             full.names = TRUE)

# Get dataframe names
gtex_all_names <- gtex_all_files %>%
  gsub(paste0(gtex_dir, "counts_mean_median/"), "", .) %>%
  gsub("_overall_counts_mean_median.txt", "", .) %>%
  gsub("-", ".", .)

# Read in all GTEx files
list2env(envir = .GlobalEnv,
         lapply(setNames(gtex_all_files, make.names(gtex_all_names)), 
                read.table, sep = "\t", header = TRUE, fill = TRUE,
                check.names = TRUE, row.names = 1))

# Create list of lists for easier transformation into data frame
all_tissue_list <- list(
  "ADP-SBC" = Adipose.Subcutaneous, 
  "ADP-VSO" = Adipose.Visceral_Omentum,
  "ADRNGL" = AdrenalGland,
  "ART-AOR" = Artery.Aorta,
  "ART-CRN" = Artery.Coronary,
  "ART-TB" = Artery.Tibial,
  "BLADDER" = Bladder,
  "BRN-AMY" = Brain.Amygdala,
  "BRN-ACC" = Brain.Anteriorcingulatecortex_BA24,
  "BRN-CAU" = Brain.Caudate_basalganglia,
  "BRN-CB-a"= Brain.CerebellarHemisphere,
  "BRN-CB-b" = Brain.Cerebellum,
  "BRN-CTX-b" = Brain.Cortex,
  "BRN-CTX-a" = Brain.FrontalCortex_BA9,
  "BRN-HIPP" = Brain.Hippocampus,
  "BRN-HYP" = Brain.Hypothalamus,
  "BRN-NAC" = Brain.Nucleusaccumbens_basalganglia,
  "BRN-PUT" = Brain.Putamen_basalganglia,
  "BRN-SPC" = Brain.Spinalcord_cervicalc.1,
  "BRN-SN" = Brain.Substantianigra,
  "BREAST" = Breast.MammaryTissue,
  "CELL-FB" = Cells.Culturedfibroblasts,
  "CELL-LYM" = Cells.EBV.transformedlymphocytes,
  "CVX-ECT" = Cervix.Ectocervix,
  "CVX-END" = Cervix.Endocervix,
  "CLN-SIG" = Colon.Sigmoid,
  "CLN-TRN" = Colon.Transverse,
  "ESP-GEJ" = Esophagus.GastroesophagealJunction,
  "ESP-MCS" = Esophagus.Mucosa,
  "ESP-MSL" = Esophagus.Muscularis,
  "FLPTB" = FallopianTube,
  "HRT-AA" = Heart.AtrialAppendage,
  "HRT-LV" = Heart.LeftVentricle,
  "KDY-CTX" = Kidney.Cortex,
  "KDY-MDL" = Kidney.Medulla,
  "LIVER" = Liver,
  "LUNG" = Lung,
  "SALGL" = MinorSalivaryGland,
  "MSC-SK" = Muscle.Skeletal,
  "NRV-TB" = Nerve.Tibial,
  "OVARY" = Ovary,
  "PANCREAS" = Pancreas,
  "BRN-PTRY" = Pituitary,
  "PROSTATE" = Prostate,
  "SKN-NSP" = Skin.NotSunExposed_Suprapubic,
  "SKN-SLL" = Skin.SunExposed_Lowerleg,
  "SIN-TIL" = SmallIntestine.TerminalIleum,
  "SPLEEN" = Spleen,
  "STOMACH" = Stomach,
  "TESTIS" = Testis,
  "THYROID" = Thyroid,
  "UTERUS" = Uterus,
  "VAGINA" = Vagina,
  "WHLBLD" = WholeBlood)

# Convert list of lists to data frames; rename columns to match tissues; pivot
gtex_prep <- function(list_of_lists, cols_to_convert, new_col) {
  # Create dataframe
  df <- Reduce(function(x, y) 
    full_join(x, y, by = c("Name", "Description")), list_of_lists) %>%
    select("Name", "Description", contains(cols_to_convert))
  # Rename columns to match tissues
  names(df) <- c("Human_ID", "Human_Symbol", names(list_of_lists))
  # Pivot longer
  pivot_longer(df, cols = 3:56, names_to = "Tissue", values_to = new_col)
}

# Prepare overall mean/median dataframes
gtex_all_means <- gtex_prep(all_tissue_list, "mean", "Mean_TPM")
gtex_all_medians <- gtex_prep(all_tissue_list, "median", "Median_TPM")

# Merge dataframes & get rid of version numbers on ENSEMBL Gene IDs
gtex_all <- full_join(gtex_all_means, gtex_all_medians) %>%
  mutate(Human_ID = sub("\\.[0-9]+", "", Human_ID))

# Get temporary dataframe without Gene_Group labeled
gtex_temp <- gtex_all %>%
  mutate(System = factor(case_when(Tissue == "BRN-CB-b" ~ NA_character_,
                                   Tissue == "BRN-CTX-b" ~ NA_character_,
                                   # Remove brain replicates from this &
                                   # Put all other brain tissues in CNS
                                   grepl("BRN", Tissue) == TRUE ~ "CNS",
                                   # Remove cell lines from this calculation
                                   Tissue == "CELL-FB" ~ NA_character_,
                                   Tissue == "CELL-LYM" ~ NA_character_,
                                   TRUE ~ "Body")),
         .after = "Tissue") %>%
  # Calculate sum of mean TPMs for each gene for brain vs. non-brain
  group_by(Human_ID, System) %>%
  mutate(System_Total = sum(Mean_TPM), .after = "System") %>%
  ungroup()

# Get separate dataframes for each dataset and then join
carp_gtex <- gtex_temp %>%
  filter(Human_ID %in% sig_carp_hs_ids) %>%
  mutate(Gene_Group = "Carpenter")
walk_gtex <- gtex_temp %>%
  filter(Human_ID %in% sig_walk_hs_ids) %>%
  mutate(Gene_Group = "Walker")

# Join for supplementary material
gtex <- Reduce(function(x, y) 
  full_join(x, y), list(carp_gtex, walk_gtex, gtex_temp)) %>%
  mutate(Gene_Group = fct_relevel(ifelse(is.na(Gene_Group), 
                                         "All Other Genes", Gene_Group),
                                  levels = c("Carpenter", "Walker", 
                                             "All Other Genes"))) %>%
  relocate(Gene_Group, .before = everything()) %>%
  unique()

# How many genes in GTEX data?
gtex %>% 
  group_by(Gene_Group) %>% 
  select(Human_ID) %>% 
  unique() %>%
  count()
# C: 577, W: 112, Other: 55525


## LOAD IN GTEX METADATA ----------
# Participant information
pt <- read_delim(paste0(gtex_dir, "participant.tsv"), delim = "\t") %>%
  select(c(1, 2, 8, 9)) %>%
  rename(Participant = "entity:participant_id", Age = "age", 
                Has_RNAseq = "has_rnaseq", Sex = "sex")

# Sample Attributes
att <- read_delim(
  paste0(gtex_dir, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"), 
  delim = "\t") %>%
  select(c(1, 7)) %>%
  rename(Sample_ID = "SAMPID", Tissue = "SMTSD") %>%
  filter(str_detect(Tissue, "Brain") | Tissue == "Pituitary") %>%
  mutate(Sample_ID_Temp = Sample_ID) %>%
  separate(col = Sample_ID_Temp, into = c("Temp1", "Temp2"), sep = "-") %>%
  unite(col = Participant, c(Temp1, Temp2), sep = "-")

# Merge dataframes
pt_att <- inner_join(pt, att) %>%
  mutate(Group = as.logical(Sex == "Female"),
         Tissue = case_when(
           Tissue == "Brain - Frontal Cortex (BA9)" ~ "Frontal_Cortex",
           Tissue == "Brain - Cortex" ~ "Cortex",
           Tissue == "Brain - Cerebellum" ~ "Cerebellum",
           Tissue == "Brain - Cerebellar Hemisphere" ~ "Cerebellar_Hemisphere",
           Tissue == "Brain - Caudate (basal ganglia)" ~ "Caudate",
           Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~ "Nucleus_Accumbens",
           Tissue == "Brain - Putamen (basal ganglia)" ~ "Putamen",
           Tissue == "Brain - Hypothalamus" ~ "Hypothalamus",
           Tissue == "Brain - Spinal cord (cervical c-1)" ~ "Spinal_Cord_C1",
           Tissue == "Brain - Hippocampus" ~ "Hippocampus",
           Tissue == "Brain - Substantia nigra" ~ "Substantia_Nigra",
           Tissue == "Brain - Frontal Cortex (BA9)" ~ "Frontal_Cortex",
           Tissue == "Brain - Anterior cingulate cortex (BA24)" ~ "Anterior_Cingulate",
           Tissue == "Brain - Amygdala" ~ "Amygdala",
           TRUE ~ Tissue)) %>%
  mutate(Sample_ID = str_to_upper(str_replace_all(Sample_ID, "-", "\\.")))

# Number of male/female samples for each tissue before matching
pt_att %>% 
  filter(Age >= 55, Has_RNAseq == TRUE) %>%
  group_by(Tissue, Sex) %>%
  count()



## LOAD IN BY-TISSUE TPM FILES ----------
# Get list of files, then load in all files
subset_files <- list.files(path = paste0(gtex_dir, "count_subsets"), 
                           pattern = "subset*", full.names = TRUE)

# Tissue names in alphabetical order (except Pituitary last), as a vector
tissue_ids <- c("Amygdala", "Anterior_Cingulate", "Caudate", "Cerebellar_Hemisphere", 
                "Cerebellum", "Cortex", "Frontal_Cortex", "Hippocampus", 
                "Hypothalamus", "Nucleus_Accumbens", "Putamen", 
                "Spinal_Cord_C1", "Substantia_Nigra", "Pituitary")

# Load files separately into environment
list2env(envir = .GlobalEnv,
         lapply(setNames(subset_files, make.names(tissue_ids)), 
                read.delim, sep = "\t"))

# Create list of lists for easier transformation into data frame
# Pivot each sub-list so that Sample IDs are in rows
tissue_list <- list("Amygdala" = Amygdala, 
                    "Anterior_Cingulate" = Anterior_Cingulate, 
                    "Caudate" = Caudate, 
                    "Cerebellar_Hemisphere" = Cerebellar_Hemisphere, 
                    "Cerebellum" = Cerebellum, "Cortex" = Cortex, 
                    "Frontal_Cortex" = Frontal_Cortex, 
                    "Hippocampus" = Hippocampus, "Hypothalamus" = Hypothalamus, 
                    "Nucleus_Accumbens" = Nucleus_Accumbens, 
                    "Putamen" = Putamen, 
                    "Spinal_Cord_C1" = Spinal_Cord_C1, 
                    "Substantia_Nigra" = Substantia_Nigra,
                    "Pituitary" = Pituitary) %>%
  lapply(pivot_longer, -c(1:2), names_to = "Sample_ID", values_to = "TPM") %>%
  map2(., names(.), ~cbind(.x, Tissue = .y)) %>%
  lapply(mutate, 
         Sample_ID = str_to_upper(str_replace_all(Sample_ID, "-", "\\.")),
         ID_Temp = as.character(Sample_ID)) %>%
  lapply(separate, col = ID_Temp, into = c("Tmp1", "Tmp2"), sep = "\\.") %>%
  lapply(unite, col = Participant, c(Tmp1, Tmp2), sep = "-")

# Create dataframe from modified tissue lists
tissue_df <- Reduce(function(x, y) 
  full_join(x, y), tissue_list)



## FIX METADATA ----------
# The sample metadata (pt and att dataframes) indicate that there are samples
# that aren't actually present in the data. 
# E.g. the frontal cortex should have the sample GTEX-1117F-0011-R10a-SM-AHZ7F
# corresponding to participant GTEX-1117F, but this sample is not present in 
# the frontal cortex samples (or any other tissues of interest).

# Join metadata and TPM dataframes by inner_join to get rid of samples that
# don't actually exist
all_df <- inner_join(pt_att, tissue_df)

# Get "new" metadata dataframe based only on samples that exist
new_metadata <- all_df %>%
  select(-c(Name, Description, TPM)) %>%
  unique()



## AGE-MATCH SAMPLES ----------
# Age-match brain samples for males and females 55+
matched <- new_metadata %>% 
  filter(Age >= 55) %>%
  matchit(Group ~ Tissue + Age, data = ., 
          method = "nearest", exact = c("Age", "Tissue"))

# New dataframe after matching
matched_df <- match.data(matched)[1:7] %>%
  mutate(Matched = TRUE)

# How many samples age 55+ per tissue and sex? BEFORE MATCHING
before_match <- new_metadata %>% 
  filter(Age >= 55) %>%
  group_by(Tissue, Sex) %>%
  summarize(Before_Match = n())

# How many samples age 55+ per tissue and sex? AFTER MATCHING
after_match <- matched_df %>%
  group_by(Tissue, Sex) %>%
  summarize(After_Match = n())

# Join w/ TPM values
tpm <- matched_df %>% inner_join(., all_df)

# Verify that the number of samples stayed the same
after_match_verify <- tpm %>%
  select(Participant, Sex, Age, Matched, Tissue) %>%
  unique() %>%
  group_by(Tissue, Sex) %>%
  summarize(After_Match_Verify = n())

# Join summary dataframes and verify
match_summary <- Reduce(function(x, y) 
  full_join(x, y), list(before_match, after_match, after_match_verify))
match_summary

# Final TPM dataframe for stats
final_bysex_tpm <- tpm %>%
  rename(Human_ID = "Name", Human_Symbol = "Description") %>%
  mutate(Human_ID = sub("\\.[0-9]+", "", Human_ID)) %>%
  mutate(Tissue = case_when(
    Tissue == "Frontal_Cortex" ~ "BRN-CTX-a",
    Tissue == "Cerebellar_Hemisphere" ~ "BRN-CB-a",
    Tissue == "Caudate" ~ "BRN-CAU",
    Tissue == "Nucleus_Accumbens" ~ "BRN-NAC",
    Tissue == "Putamen" ~ "BRN-PUT",
    Tissue == "Hypothalamus" ~ "BRN-HYP",
    Tissue == "Spinal_Cord_C1" ~ "BRN-SPC",
    Tissue == "Hippocampus" ~ "BRN-HIPP",
    Tissue == "Substantia_Nigra" ~ "BRN-SN",
    Tissue == "Anterior_Cingulate" ~ "BRN-ACC",
    Tissue == "Amygdala" ~ "BRN-AMY",
    Tissue == "Pituitary" ~ "BRN-PTRY",
    TRUE ~ Tissue))


## SEX/TISSUE ANOVAS FOR INDIVIDUAL GENES ----------
all_hs_ids <- c(sig_carp_hs_ids, sig_walk_hs_ids)

# Overall ANOVAs
bysex_tissue_anovas_temp <- final_bysex_tpm %>%
  filter(Tissue != "Cerebellum", Tissue != "Cortex", 
         Human_ID %in% all_hs_ids) %>%
  unique() %>%
  group_by(Human_Symbol) %>%
  anova_test(TPM ~ Sex + Tissue + Sex:Tissue, type = "III") %>%
  as.data.frame() %>%
  rename(Sig_p05 = "p<.05")

bysex_tissue_anovas <- gtex %>%
  select(Gene_Group, Human_ID, Human_Symbol) %>%
  right_join(., bysex_tissue_anovas_temp) %>%
  unique() %>%
  filter(Gene_Group != "All Other Genes")

# Genes with main effects or interactions
genes_sex <- bysex_tissue_anovas %>% 
  filter(Effect == "Sex", Sig_p05 == "*") %>% 
  pull(Human_Symbol) %>%
  unique()
# 179 genes

genes_tissue <- bysex_tissue_anovas %>% 
  filter(Effect == "Tissue", Sig_p05 == "*")  %>%
  pull(Human_Symbol) %>%
  unique()
# 671 genes

genes_inter <- bysex_tissue_anovas %>% 
  filter(Effect == "Sex:Tissue", Sig_p05 == "*") %>% 
  pull(Human_Symbol) %>%
  unique()
# 55 genes

# SEX EFFECT:
# No need for post-hoc tests for Sex (only 2 groups)
# Can still determine which sex has higher expression

# Genes with higher expression in females
female_higher <- final_bysex_tpm %>%
  filter(Human_Symbol %in% genes_sex) %>%
  unique() %>%
  group_by(Human_Symbol, Sex) %>%
  summarize(Mean_TPM = mean(TPM))%>%
  pivot_wider(names_from = Sex, values_from = Mean_TPM) %>%
  mutate(Higher_Sex = case_when(Female > Male ~ "Female", 
                                Male > Female ~ "Male")) %>%
  filter(Higher_Sex == "Female") %>%
  pull(Human_Symbol) %>%
  unique()
length(female_higher) # 39 of 179

# Use that vector of genes to label a summary dataframe
sex_summary_temp <- final_bysex_tpm %>%
  filter(Human_Symbol %in% genes_sex) %>%
  unique() %>%
  group_by(Human_Symbol, Sex) %>%
  summarize(Mean_TPM = mean(TPM), SD_TPM = sd(TPM)) %>%
  mutate(Higher_Sex = case_when(Human_Symbol %in% female_higher ~ "Female",
                                TRUE ~ "Male")) %>%
  pivot_wider(names_from = Sex, values_from = c(Mean_TPM, SD_TPM))

sex_summary <- gtex %>%
  select(Gene_Group, Human_ID, Human_Symbol) %>%
  right_join(., sex_summary_temp) %>%
  unique() %>%
  filter(Gene_Group != "All Other Genes")

# TISSUE EFFECT:
# Post-hoc Tukey tests for each gene with significant main effect of Tissue
tukey_tissue_temp <- final_bysex_tpm %>%
  filter(Tissue != "Cerebellum", Tissue != "Cortex") %>%
  mutate(Tissue = gsub("-", "_", Tissue)) %>%
  filter(Human_Symbol %in% genes_tissue | Human_Symbol %in% genes_inter) %>%
  unique() %>%
  group_by(Human_Symbol) %>%
  tukey_hsd(TPM ~ Tissue) %>%
  relocate(Human_Symbol, .before = term) %>%
  as.data.frame()

# Get human IDs and symbols together in a dataframe
hs_sym_ids <- final_bysex_tpm %>% select(Human_Symbol, Human_ID)

tukey_tissue <- gtex %>%
  select(Gene_Group, Human_ID, Human_Symbol) %>%
  right_join(., tukey_tissue_temp) %>%
  unique() %>%
  filter(Gene_Group != "All Other Genes")

# SEX:TISSUE INTERACTION
# Post-hoc Tukey tests for each gene with significant Sex:Tissue interaction
tukey_sex_tissue_temp <- final_bysex_tpm %>%
  filter(Tissue != "Cerebellum", Tissue != "Cortex") %>%
  mutate(Tissue = gsub("-", "_", Tissue)) %>%
  filter(Human_Symbol %in% genes_inter) %>%
  unique() %>%
  group_by(Human_Symbol, Tissue) %>%
  tukey_hsd(TPM ~ Sex) %>%
  relocate(Human_Symbol, .before = everything()) %>%
  as.data.frame()

tukey_sex_tissue <- gtex %>%
  select(Gene_Group, Human_ID, Human_Symbol) %>%
  right_join(., tukey_sex_tissue_temp) %>%
  unique() %>%
  filter(Gene_Group != "All Other Genes")

# Load specific packages
library(emmeans)
library(multcomp)

# Use vector of genes to label a summary dataframe
sex_tissue_summary <- final_bysex_tpm %>%
  filter(Human_Symbol %in% genes_inter) %>%
  unique() %>%
  mutate(Tissue = gsub("-", "_", Tissue)) %>%
  group_by(Human_Symbol, Tissue, Sex) %>%
  summarize(Mean_TPM = mean(TPM), SD_TPM = sd(TPM)) %>%
  mutate(Higher_Sex = case_when(Human_Symbol %in% female_higher ~ "Female",
                                TRUE ~ "Male")) %>%
  pivot_wider(names_from = Sex, values_from = c(Mean_TPM, SD_TPM)) %>%
  full_join(tukey_sex_tissue) %>%
  filter(p.adj.signif != "ns") %>%
  droplevels() %>%
  relocate(Gene_Group, Human_ID, .before = everything()) %>%
  select(1:9)

# CLDs FOR TISSUE EFFECT
# Get compact letter display (CLD) for tissue data
# First create function
get_cld <- function(gene) {
  final_bysex_tpm %>%
    filter(Human_Symbol == gene,
           Tissue != "Cerebellum", Tissue != "Cortex") %>%
    unique() %>%
    mutate(Tissue = as.factor(Tissue)) %>%
    lm(TPM ~ Tissue, data = .) %>%
    emmeans(., pairwise ~ "Tissue", adjust = "tukey") %>%
    cld(., reversed = TRUE)
}

# Create empty list
cld_list <- list()

# Get CLD for each gene and populate list
for (i in 1:length(genes_tissue)) {
  cld_list[[i]] <- get_cld(genes_tissue[i]) %>%
    as.data.frame() %>%
    arrange(Tissue) %>%
    mutate(Human_Gene_Symbol = genes_tissue[i], .before = everything())
}

# Convert list to dataframe
cld_df_temp <- Reduce(function(x, y) 
  full_join(x, y), cld_list) %>%
  mutate(Group = str_trim(.group), .keep = "unused") %>%
  mutate(Group = gsub("1", "a", Group),
         Group = gsub("2", "b", Group),
         Group = gsub("3", "c", Group),
         Group = gsub("4", "d", Group),
         Group = gsub("5", "e", Group),
         Group = gsub("6", "f", Group),
         Group = gsub("7", "g", Group),
         Group = gsub("8", "h", Group),
         Group = gsub("9", "i", Group))

cld_df <- gtex %>%
  select(Gene_Group, Human_ID, Human_Symbol) %>%
  right_join(., cld_df_temp, by = c("Human_Symbol" = "Human_Gene_Symbol")) %>%
  unique() %>%
  filter(Gene_Group != "All Other Genes")

# Unload packages that don't play nicely with tidyverse
detach("package:multcomp", unload = TRUE)
detach("package:TH.data", unload = TRUE)
detach("package:MASS", unload = TRUE)



## BRAIN SPECIFICITY ----------
# Calculate "Brain Specificity" (log2 FC of CNS vs. Body Totals)
# Can apply expression filters here also, e.g. TPM > 1 in all tissues
brain_spec_firstpass <- gtex %>%
  filter(!is.na(System)) %>%
  group_by(Human_ID) %>%
  mutate(Grand_Tissue_Mean = mean(Mean_TPM)) %>%
  select(Human_ID, Human_Symbol, System, System_Total,
                Grand_Tissue_Mean, Gene_Group) %>%
  unique() %>%
  pivot_wider(names_from = System, values_from = System_Total) %>%
  mutate(BrainSpec = log2((CNS)/(Body))) %>%
  ungroup()

# Filter out genes based on grand mean over all tissues
brain_spec <- brain_spec_firstpass %>% filter(Grand_Tissue_Mean > 0)
filter_genes <- brain_spec %>% pull(Human_ID) %>% unique()

# Count genes before and after filtering
# All genes
brain_spec_firstpass %>% group_by(Human_ID) %>% count() %>% nrow() # 56,200
brain_spec %>% group_by(Human_ID) %>% count() %>% nrow() # 55,874 genes after
# "All Other Genes"
brain_spec_firstpass %>% filter(Gene_Group == "All Other Genes") %>% 
  group_by(Human_ID) %>% count() %>% nrow() # 55,525 before
brain_spec %>% filter(Gene_Group == "All Other Genes") %>% 
  group_by(Human_ID) %>% count() %>% nrow() # 55,199 after
# Candidate genes
brain_spec_firstpass %>% filter(Gene_Group != "All Other Genes") %>% 
  group_by(Human_ID) %>% count() %>% nrow() # 675 before
brain_spec %>% filter(Gene_Group != "All Other Genes") %>% 
  group_by(Human_ID) %>% count() %>% nrow() # 675 after

# Calculate "Within-Brain Percentage":
# Mean TPM for a gene in a CNS tissue vs. that gene's CNS Total
wnbrain_perc <- gtex %>%
  # Only keep the filtered genes
  filter(Human_ID %in% brain_spec_firstpass$Human_ID) %>%
  filter(System == "CNS") %>%
  group_by(Human_ID) %>%
  summarize(Human_ID = Human_ID, 
            Human_Symbol = Human_Symbol,
            Tissue = Tissue,
            Mean_TPM = Mean_TPM,
            CNS_Total = System_Total,
            WithinBrain_Perc = Mean_TPM/CNS_Total * 100) %>%
  ungroup()

# Join dataframes and add new columns for plotting later
gtex_spec <- Reduce(full_join, list(gtex, brain_spec_firstpass, wnbrain_perc)) %>%
  filter(Human_ID %in% filter_genes)

cns_means <- gtex_spec %>% 
  filter(System == "CNS", Gene_Group != "All Other Genes") %>%
  select(Human_ID, Human_Symbol, Tissue, Mean_TPM, Median_TPM) %>%
  group_by(Human_ID) %>%
  mutate(CNS_Mean = mean(Mean_TPM)) %>%
  select(Human_ID, Human_Symbol, CNS_Mean) %>%
  arrange(Human_Symbol) %>%
  unique()

# Number of genes with higher expression in brain than body
gtex_spec %>%
  filter(Gene_Group != "All Other Genes", BrainSpec > 0) %>%
  group_by(Gene_Group) %>%
  select(Human_ID, BrainSpec) %>%
  unique() %>%
  count()
  


### PLOT ----------
# Genes to label
shared_label_list <- c("CARTPT")
general_label_list <- c("MOBP", "GAD2",
                        "GPR101", "PRKCG",
                        "FAM53B", "C1QL2", "CREB1", "DRD3")

# All overlaps in ggrepel
options(ggrepel.max.overlaps = Inf)

bs_plot <- brain_spec %>%
  mutate(Gene_Group = as.character(Gene_Group),
         Gene_Group = fct_relevel(case_when(Gene_Group == "All Other Genes" ~ "All Other Genes",
                                            TRUE ~ Gene_Group),
                                  levels = c("Carpenter", "Walker", "All Other Genes"))) %>%
  select(Human_ID, Human_Symbol, BrainSpec, Gene_Group) %>%
  unique() %>%
  filter(BrainSpec != -Inf, BrainSpec != Inf) %>%
  ggplot(aes(x = Gene_Group, y = BrainSpec, fill = Gene_Group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(seq(-16, -2, 2), seq(2, 12, 2)), 
             linetype = "solid", color = "grey85") +
  geom_hline(yintercept = c(seq(-16, -4, 4), seq(4, 12, 4)), 
             linetype = "solid", color = "grey65") +
  geom_violin(width = 0.8) +  
  geom_boxplot(fill = "white", color = "black", width = 0.18, 
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.6) +
  labs(x = "", y = expression("Brain Specificity (Log "[2]*" Fold Change)"),
       title = "Brain Specificity (GTEx)\n") +
  scale_y_continuous(limits = c(-16, 12), breaks = seq(-16, 12, 4), 
                     expand = c(0, 0.00001)) +
  scale_x_discrete(expand = c(0, 0.55)) +
  scale_fill_manual(values = c("#44AA99", "#C148AD", "grey80")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        plot.title = element_text(face = "bold"),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(Human_Symbol %in% shared_label_list, 
                                           Human_Symbol, NA_character_)),
                   aes(label = label), color = "black", fill = "white",
                   nudge_y = 0.01, nudge_x = 0.01,
                   box.padding = 0.5, fontface = "bold", size = 3, seed = 58) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(Human_Symbol %in% general_label_list, 
                                           Human_Symbol, NA_character_)),
                   aes(label = label), color = "black",
                   box.padding = 0.5, fontface = "bold", size = 3, nudge_y = 1, seed = 1000)
bs_plot

# Save as svg
# svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/BrainSpec_9_25_2022.svg", width = 4.2, height = 5.5)
# bs_plot
# dev.off()

### STATISTICS ----------
# Wilcoxon test between groups
# Are craving genes more specific to the brain than other genes?
gtex_spec %>%
  select(Human_ID, System, BrainSpec, Gene_Group) %>%
  unique() %>%
  filter(System == "CNS", BrainSpec != -Inf, BrainSpec != Inf) %>%
  wilcox_test(BrainSpec ~ Gene_Group) %>%
  filter(group2 == "All Other Genes")

# Summary data for brain specificity by gene group
summary_brainspec <- gtex_spec %>%
  select(Human_Symbol, System, BrainSpec, Gene_Group) %>%
  unique() %>%
  filter(System == "CNS", BrainSpec != -Inf, BrainSpec != Inf) %>%
  group_by(Gene_Group) %>%
  summarize(Median_BrainSpec = median(BrainSpec))
summary_brainspec


## FINAL DATAFRAME (Need to run conservation script first) ----------
# Get CNS means
cns_means_df <- gtex %>%
  filter(!is.na(System)) %>%
  group_by(Human_ID, System) %>%
  mutate(System_Mean = mean(Mean_TPM)) %>%
  select(Human_ID, Human_Symbol, System, System_Total,
                System_Mean, Gene_Group) %>%
  unique() %>%
  pivot_wider(names_from = System, values_from = System_Mean) %>%
  select(-Body) %>%
  rename(CNS_Mean = CNS) %>%
  filter(!is.na(CNS_Mean)) %>%
  ungroup()

# Get regions of highest expression
temp_df <- cld_df %>%
  select(Human_Symbol, Tissue, Group) %>%
  filter(grepl("a", Group)) %>%
  unique() %>%
  pivot_wider(names_from = Tissue, values_from = Group)

for (i in 1:nrow(temp_df)){
  for (cl in colnames(temp_df)[2:ncol(temp_df)]){
    temp_df[i, cl] <- ifelse(is.na(temp_df[i, cl]), NA_character_, cl)
  }
}

temp_df <- temp_df %>%
  unite(col = "All", 2:ncol(.), sep = ", ", na.rm = TRUE)
  
highest_region_df <- cld_df %>%
  dplyr::select(Human_Symbol) %>%
  unique() %>%
  mutate(Highest_CNS_Regions = "")

for (i in 1:nrow(temp_df)) {
  highest_region_df[i, 2] <- temp_df %>%
    .[i, ] %>%
    .[ , apply(., 2, function(x) !any(is.na(x)))] %>%
    unite("Highest_CNS_Regions", 2:ncol(.), remove = TRUE, sep = ", ") %>%
    select(2)
}

# Load conservation data
full_cons_df <- readRDS(file = paste0(main_dir, "post_processing/results/other/full_cons_df.RDS"))

# Final dataframe with everything
final_df <- full_join(brain_spec, cns_means_df) %>%
  mutate(Gene_Group = as.factor(ifelse(Gene_Group == "All Other Genes", 
                                       "All Other Orthologs", 
                                       as.character(Gene_Group)))) %>%
  filter(Gene_Group != "All Other Orthologs") %>%
  full_join(., highest_region_df) %>%
  full_join(full_cons_df, .) %>%
  filter(Gene_Group != "All Other Orthologs") %>%
  select(1:18, 23, 26, 25)
# 862 * 21 dimensions

saveRDS(final_df, file = paste0(main_dir, "post_processing/results/other/final_df_05-20-23.RDS"))


## MAKE GTEX PLOTS ----------
### TISSUE PLOTS ----------
vert_gene_plot <- function(gene_name, x_lim, br, strip_color = "#DE8971", 
                           blanks = "black") {
  
  gtex_spec %>%
    mutate(System = fct_rev(case_when(
      Tissue == "BRN-CB-b" | Tissue == "BRN-CTX-b" ~ NA_character_,
      Tissue == "CELL-FB" | Tissue == "CELL-LYM" ~ NA_character_,
      TRUE ~ as.character(System)))) %>%
    filter(Human_Symbol == gene_name, !is.na(System)) %>%
    ggplot(aes(x = fct_rev(Tissue), fill = System, color = System)) +
    geom_point(aes(y = Mean_TPM), shape = 21, size = 3.25, color = "black",
               alpha = 0.6) +
    geom_point(aes(y = Median_TPM), shape = 4, stroke = 1.5, size = 3) +
    scale_x_discrete(expand = c(0, 1.1)) +
    scale_y_continuous(limits = c(0, x_lim), breaks = seq(0, x_lim, br)) +
    scale_fill_manual(values = c("#9359CD", "grey80")) +
    scale_color_manual(values = c("#9359CD", "grey60")) +
    labs(x = "Tissue", y = "Mean/Median TPM", fill = "System") +
    theme(legend.title.align = 0.5, legend.position = "none",
          legend.margin = margin(-1.25, 0, 0, 0, "cm"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_line(color = "grey60"),
          axis.text.y = element_text(size = 7.25, color = blanks),
          axis.title.y = element_text(color = blanks),
          axis.ticks.y = element_line(color = blanks),
          axis.title.x = element_text(size = 11, color = "white"),
          axis.text.x = element_text(size = 11),
          plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"),
          strip.background = element_rect(fill = strip_color),
          strip.text = element_text(size = 11)) +
    coord_flip() +
    facet_wrap(~ Human_Symbol, ncol = 4, nrow = 1, scales = "free_y")
}


# Interesting candidate genes
svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/GAD2_Tissue_9_25_2022.svg", width = 2.35, height = 5.4)
vert_gene_plot("GAD2", 75, 25, "#44AA99", "black")
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/KIF5A_Tissue_9_25_2022.svg", width = 2.35, height = 5.4)
vert_gene_plot("KIF5A", 1200, 400, "#44AA99", "white")
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/KIF5A_Tissue_9_25_2022.svg", width = 2.35, height = 5.4)
vert_gene_plot("GPR101", 8, 2, "#C148AD", "white")
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/PRKCG_Tissue_9_25_2022.svg", width = 2.35, height = 5.4)
vert_gene_plot("PRKCG", 102, 25, "#C148AD", "white")
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/CARTPT_Tissue_9_25_2022.svg", width = 2.35, height = 5.4)
vert_gene_plot("CARTPT", 450, 150, "white", "white")
dev.off()


# Comparison genes
svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/CREB1_Tissue_9_25_2022.svg", width = 2.5, height = 5.4)
vert_gene_plot("CREB1", 18, 3, "grey80", "black")
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/C1QL2_Tissue_9_25_2022.svg", width = 2.5, height = 5.4)
vert_gene_plot("C1QL2", 24, 4, "grey80", "white")
dev.off()


