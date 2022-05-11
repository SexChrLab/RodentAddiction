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
main_dir <- "C:/Annika/Github Repositories/RodentAddiction/"
de_dir <- paste0(main_dir, "post_processing/results/gene_info/")
gtex_dir <- paste0(main_dir, "post_processing/downloaded_data/gtex/")

# Reduce chances of scientific notation
options(scipen = 9999)

# Update and set ggplot theme globally
new_theme <- theme_bw(base_size = 14) %+replace%
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))
theme_set(new_theme)



## LOAD IN GENE LISTS ----------
# All in one list
all_genes_list <- readRDS(file = paste0(main_dir, "post_processing/results/other/all_genes_list.RDS"))

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
    dplyr::select("Name", "Description", contains(cols_to_convert))
  # Rename columns to match tissues
  names(df) <- c("Human_ID", "Human_Symbol", names(list_of_lists))
  # Pivot longer
  pivot_longer(df, cols = 3:55, names_to = "Tissue", values_to = new_col)
}

# Prepare overall mean/median dataframes
gtex_all_means <- gtex_prep(all_tissue_list, "mean", "Mean_TPM")
gtex_all_medians <- gtex_prep(all_tissue_list, "median", "Median_TPM")

# Merge dataframes & get rid of version numbers on ENSEMBL Gene IDs
gtex_all <- full_join(gtex_all_means, gtex_all_medians) %>%
  mutate(Human_ID = sub("\\.[0-9]+", "", Human_ID))

# Mark experimental and GWAS genes based on Human_ID
gtex <- gtex_all %>%
  mutate(Gene_Group = fct_relevel(
    case_when(Human_ID %in% sd_crave_hs ~ "Craving",
              TRUE ~ "All Other Genes"),
    c("Craving", "All Other Genes"))) %>%
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

# How many craving genes in GTEX data? (32)
gtex %>% 
  group_by(Gene_Group) %>% 
  dplyr::select(Human_ID) %>% 
  unique() %>%
  count()



## LOAD IN GTEX METADATA ----------
# Participant information
pt <- read_delim(paste0(gtex_dir, "participant.tsv"), delim = "\t") %>%
  dplyr::select(c(1, 2, 8, 9)) %>%
  dplyr::rename(Participant = "entity:participant_id", Age = "age", 
                Has_RNAseq = "has_rnaseq", Sex = "sex")

# Sample Attributes
att <- read_delim(
  paste0(gtex_dir, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"), 
  delim = "\t") %>%
  dplyr::select(c(1, 7)) %>%
  dplyr::rename(Sample_ID = "SAMPID", Tissue = "SMTSD") %>%
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
                           pattern = "subset*", full.names = TRUE) %>%
  .[-1]

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
# the frontal cortex samples (or any other tissues).

# Join metadata and TPM dataframes by inner_join to get rid of samples that
# don't actually exist
all_df <- inner_join(pt_att, tissue_df)

# Get "new" metadata dataframe based only on samples that exist
new_metadata <- all_df %>%
  dplyr::select(-c(Name, Description, TPM)) %>%
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
  dplyr::summarize(Before_Match = n())

# How many samples age 55+ per tissue and sex? AFTER MATCHING
after_match <- matched_df %>%
  group_by(Tissue, Sex) %>%
  dplyr::summarize(After_Match = n())

# Join w/ TPM values
tpm <- matched_df %>% inner_join(., all_df)

# Verify that the number of samples stayed the same
after_match_verify <- tpm %>%
  dplyr::select(Participant, Sex, Age, Matched, Tissue) %>%
  unique() %>%
  group_by(Tissue, Sex) %>%
  dplyr::summarize(After_Match_Verify = n())

# Join summary dataframes and verify
match_summary <- Reduce(function(x, y) 
  full_join(x, y), list(before_match, after_match, after_match_verify))
match_summary

# Final TPM dataframe for stats
final_bysex_tpm <- tpm %>%
  dplyr::rename(Human_ID = "Name", Human_Symbol = "Description") %>%
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
# Overall ANOVAs
bysex_tissue_anovas <- final_bysex_tpm %>%
  filter(Tissue != "Cerebellum", Tissue != "Cortex", 
         Human_ID %in% sd_crave_hs) %>%
  group_by(Human_Symbol) %>%
  anova_test(TPM ~ Sex + Tissue + Sex:Tissue, type = "III") %>%
  as.data.frame() %>%
  dplyr::rename(Sig_p05 = "p<.05")

# Genes with main effects or interactions
genes_sex <- bysex_tissue_anovas %>% 
  filter(Effect == "Sex", Sig_p05 == "*") %>% 
  pull(Human_Symbol)
# [1] "AMZ1"    "B2M"     "NTS"     "PITPNM3" "USP46"  
genes_tissue <- bysex_tissue_anovas %>% 
  filter(Effect == "Tissue", Sig_p05 == "*") %>% 
  pull(Human_Symbol)
# [1] "AGK"     "AMZ1"    "B2M"     "BCAS1"   "BTG1"    "CACYBP"  "CARTPT"  "CCDC88C"
# [9] "CELF6"   "EGR2"    "FABP7"   "FKBP4"   "FTH1"    "GPD1"    "GUCY1A3" "HAPLN2" 
# [17] "HSPA8"   "IRS2"    "KIF5A"   "LYPD1"   "MBP"     "MOBP"    "NTS"     "PHLDA1" 
# [25] "PITPNM3" "RGS5"    "RPS6KA2" "SOX17"   "TIPARP"  "TTLL1"   "USP46"   "VIM"    
genes_inter <- bysex_tissue_anovas %>% 
  filter(Effect == "Sex:Tissue", Sig_p05 == "*") %>% 
  pull(Human_Symbol)
# [1] "CARTPT"  "CELF6"   "HAPLN2"  "NTS"     "RPS6KA2"

# SEX EFFECT:
# No need for post-hoc tests for Sex (only 2 groups)
# Can still determine which sex has higher expression

# Genes with higher expression in females (just B2M)
female_higher <- final_bysex_tpm %>%
  filter(Human_Symbol %in% genes_sex) %>%
  group_by(Human_Symbol, Sex) %>%
  summarize(Mean_TPM = mean(TPM))%>%
  pivot_wider(names_from = Sex, values_from = Mean_TPM) %>%
  mutate(Higher_Sex = case_when(Female > Male ~ "Female", 
                                Male > Female ~ "Male")) %>%
  filter(Higher_Sex == "Female") %>%
  pull(Human_Symbol)


# Use that vector of genes to label a summary dataframe
sex_summary <- final_bysex_tpm %>%
  filter(Human_Symbol %in% genes_sex) %>%
  group_by(Human_Symbol, Sex) %>%
  summarize(Mean_TPM = mean(TPM), SD_TPM = sd(TPM)) %>%
  mutate(Higher_Sex = case_when(Human_Symbol %in% female_higher ~ "Female",
                                TRUE ~ "Male"))

# TISSUE EFFECT:
# Post-hoc Tukey tests for each gene with significant main effect of Tissue
tukey_tissue <- final_bysex_tpm %>%
  filter(Tissue != "Cerebellum", Tissue != "Cortex") %>%
  mutate(Tissue = gsub("-", "_", Tissue)) %>%
  filter(Human_Symbol %in% genes_tissue | Human_Symbol %in% genes_inter) %>%
  group_by(Human_Symbol) %>%
  tukey_hsd(TPM ~ Tissue) %>%
  as.data.frame()

# SEX:TISSUE INTERACTION
# Post-hoc Tukey tests for each gene with significant Sex:Tissue interaction
tukey_sex_tissue <- final_bysex_tpm %>%
  filter(Tissue != "Cerebellum", Tissue != "Cortex") %>%
  mutate(Tissue = gsub("-", "_", Tissue)) %>%
  filter(Human_Symbol %in% genes_inter) %>%
  group_by(Human_Symbol, Tissue) %>%
  tukey_hsd(TPM ~ Sex) %>%
  as.data.frame()

# Load specific packages
library(emmeans)
library(multcomp)

# Get compact letter display (CLD) for tissue data
# First create function
get_cld <- function(gene) {
  final_bysex_tpm %>%
    filter(Human_Symbol == gene,
           Tissue != "Cerebellum", Tissue != "Cortex") %>%
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
cld_df <- Reduce(function(x, y) 
  full_join(x, y), cld_list) %>%
  mutate(Group = str_trim(.group), .keep = "unused") %>%
  mutate(Group = gsub("1", "a", Group),
         Group = gsub("2", "b", Group),
         Group = gsub("3", "c", Group),
         Group = gsub("4", "d", Group),
         Group = gsub("5", "e", Group),
         Group = gsub("6", "f", Group),
         Group = gsub("7", "g", Group),
         Group = gsub("8", "h", Group))

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
  dplyr::select(Human_ID, Human_Symbol, System, System_Total,
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
brain_spec %>% nrow() # 55,873 genes after
# "All Other Genes"
brain_spec_firstpass %>% filter(Gene_Group == "All Other Genes") %>% 
  group_by(Human_ID) %>% count() %>% nrow() # 56,184 before
brain_spec %>% filter(Gene_Group == "All Other Genes") %>% nrow() # 55,857 genes
# CE, CW, and GWAS genes after
brain_spec %>% filter(Gene_Group == "Craving") %>% nrow() # 16/16 Craving genes after

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
  filter(System == "CNS", Gene_Group == "Craving") %>%
  dplyr::select(Human_ID, Human_Symbol, Tissue, Mean_TPM, Median_TPM) %>%
  group_by(Human_ID) %>%
  mutate(CNS_Mean = mean(Mean_TPM)) %>%
  dplyr::select(Human_ID, Human_Symbol, CNS_Mean) %>%
  arrange(Human_Symbol) %>%
  unique()

# Genes to label
label_list <- c("MOBP",
                "FAM53B", "C1QL2", "DRD3", "CREB1", "NR4A1")

# All overlaps in ggrepel
options(ggrepel.max.overlaps = Inf)

# Brain specificity plot
bs_plot <- brain_spec %>%
  dplyr::select(Human_ID, Human_Symbol, BrainSpec, Gene_Group) %>%
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
               outlier.fill = "black", outlier.alpha = 0.9) +
  labs(x = "", y = expression("Brain Specificity (Log "[2]*" Fold Change)"),
       title = "Specificity of Craving Genes to the Brain (GTEx)\n") +
  scale_y_continuous(limits = c(-16, 12), breaks = seq(-16, 12, 4), 
                     expand = c(0, 0.00001)) +
  scale_x_discrete(expand = c(0, 0.55)) +
  scale_fill_manual(values = c("#DE8971", "grey80")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        plot.title = element_blank()) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(Human_Symbol %in% label_list, 
                                           Human_Symbol, NA_character_)),
                   aes(label = label), color = "black",
                   box.padding = 0.8, fontface = "bold", size = 3.5) +
  stat_compare_means(comparisons = list(c("Craving", "All Other Genes")), method = "wilcox.test", 
                     hide.ns = FALSE, label.y = c(9.3), tip.length = c(0, 0.05))

# Save as pdf
# svglite("G:/Users/Annika/Documents/Figures for Gene Conservation Project/BrainSpec_1_24_2022.svg", width = 4, height = 5.5)
bs_plot
# dev.off()

# One-sample Wilcoxon tests
# Function
one_sample_wilcox <- function(group) {
  gtex_spec %>%
    dplyr::select(Human_ID, System, BrainSpec, Gene_Group) %>%
    unique() %>%
    filter(System == "CNS", Gene_Group == group,
           BrainSpec != -Inf, BrainSpec != Inf) %>%
    wilcox_test(BrainSpec ~ 1, mu = 0)
}

# Are craving genes highly expressed in the brain?
one_sample_wilcox("Craving")
one_sample_wilcox("All Other Genes")

# Wilcoxon test between groups
# Are craving genes more specific to the brain than other genes?
gtex_spec %>%
  dplyr::select(Human_ID, System, BrainSpec, Gene_Group) %>%
  unique() %>%
  filter(System == "CNS", BrainSpec != -Inf, BrainSpec != Inf) %>%
  wilcox_test(BrainSpec ~ Gene_Group)
# .y.       group1  group2             n1    n2 statistic        p
# 1 BrainSpec Craving All Other Genes    32 55305   1233008 0.000116

# Summary data for brain specificity by gene group
summary_brainspec <- gtex_spec %>%
  dplyr::select(Human_Symbol, System, BrainSpec, Gene_Group) %>%
  unique() %>%
  filter(System == "CNS", BrainSpec != -Inf, BrainSpec != Inf) %>%
  group_by(Gene_Group) %>%
  summarize(Mean_BrainSpec = mean(BrainSpec), 
            Median_BrainSpec = median(BrainSpec))
summary_brainspec

# This can be sorted to look for interesting targets
crav_gtex_info <- gtex_spec %>%
  filter(System == "CNS", Gene_Group == "Craving") %>%
  dplyr::select(Human_ID, Human_Symbol, Tissue, BrainSpec, WithinBrain_Perc) %>%
  unique()

# 11 genes have higher expression in brain than the rest of the body
crav_gtex_info %>% 
  dplyr::select(-c(Tissue, WithinBrain_Perc)) %>%
  filter(BrainSpec > 0) %>%
  unique() %>%
  arrange(desc(BrainSpec)) %>%
  pull(Human_Symbol)
# [1] "MOBP"   "KIF5A"  "MBP"    "HAPLN2" "BCAS1"  "FABP7" 
# [7] "LYPD1"  "CARTPT" "AMZ1"   "NTS"

# Rank brain tissues by highest expression
rank_brain <- gtex_spec %>% 
  group_by(Human_ID) %>%
  filter(System == "CNS") %>%
  mutate(Rank = dense_rank(dplyr::desc(WithinBrain_Perc)), .before = Tissue) %>%
  dplyr::select(Human_ID, Human_Symbol, BrainSpec, Rank, Tissue, WithinBrain_Perc,
                Mean_TPM, Gene_Group)

# Short tables w/ top 5 regions for publication
# Function
wnbrain_table <- function(group) {
  rank_brain %>%
    filter(Gene_Group == group, Rank <= 5) %>%
    dplyr::select(-Gene_Group) %>%
    arrange(Human_Symbol, Rank)
}

# Create table
rank_crave_tbl <- wnbrain_table("Craving")

# Get list only for craving genes
rank_crave <- rank_brain %>% filter(Gene_Group == "Craving")

# Get genes with high expression (>= 10% of brain expression) in reward Systems
reward_rank <- rank_crave %>%
  dplyr::select(-Gene_Group) %>%
  filter(WithinBrain_Perc >= 0.10) %>%
  filter(!Tissue %in% c("BRN-CB-a", "BRN-HYP", "BRN-SPC", "BRN-PTRY")) %>%
  arrange(Human_Symbol, Rank)



## COPY DATA INTO SPREADSHEETS ----------
write_clip(cld_df)

tukey_sex_tissue %>%
  relocate(Human_Symbol, .before = everything()) %>%
  write_clip()

tukey_tissue %>% write_clip()

rank_crave %>%
  relocate(Human_Symbol, .before = everything()) %>%
  dplyr::select(1:3) %>%
  unique() %>%
  arrange(Human_Symbol) %>%
  write_clip()

bysex_tissue_anovas %>% write_clip()


## MAKE GTEX PLOTS ----------
# Function for creating a 2-panel plot with by-system and by-sex data for a
# single gene


# Plot means by tissue for specific genes (all samples)


horiz_gene_plot <- function(gene, high_limit, break_unit, rect = "#DE8971",
                      adjust = 0) {
  # Create labeller for plots
  gene_label <- c("All Samples")
  names(gene_label) <- c(gene)
  
  gtex_spec %>%
    mutate(System = fct_rev(case_when(
      Tissue == "BRN-CB-b" | Tissue == "BRN-CTX-b" ~ "CNS",
      Tissue == "CELL-FB" | Tissue == "CELL-LYM" ~ "Body",
      TRUE ~ as.character(System)))) %>%
    filter(Human_Symbol == gene) %>%
    ggplot(aes(x = Tissue, fill = System, color = System)) +
    geom_point(aes(y = Mean_TPM), shape = 21, size = 5, color = "black",
               alpha = 0.6) +
    geom_point(aes(y = Median_TPM), shape = 4, stroke = 1.5, size = 5) +
    scale_y_continuous(limits = c(0, high_limit),
                       breaks = seq(0, high_limit, break_unit)) +
    scale_fill_manual(values = c("#44AA99", "#332288")) +
    scale_color_manual(values = c("#44AA99", "#332288")) +
    labs(x = "", y = "Mean/Median TPM", fill = "System") +
    theme(legend.title.align = 0.5, legend.position = "bottom",
          legend.margin = margin(-1.25, 0, 0, 0, "cm"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(color = "grey70"),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          plot.margin = margin(0.25, 0.25, 0.75, 0.25, "cm"),
          strip.background = element_rect(fill = rect)) +
    facet_wrap(~ Human_Symbol, labeller = labeller(Human_Symbol = gene_label))
}

sex_plot <- function(gene, high_limit, break_unit, rect = "#DE8971",
                            adjust = 0, pos = high_limit-0.5*break_unit) {
  
  # Compact letter display
  cld_info <- cld_df %>% 
    filter(Human_Gene_Symbol == gene) %>% 
    dplyr::select(Tissue, Group)
  
  # Plot means by tissue for specific genes (by sex)
  final_bysex_tpm %>%
    filter(Human_Symbol == gene, Tissue != "Cerebellum", Tissue != "Cortex") %>%
    ggplot(aes(x = Tissue)) +
    stat_boxplot(aes(y = TPM, fill = Sex), geom = "errorbar", size = 1, color = "grey10",
                 position = position_dodge(width = 0.8), width = 0.6) +
    geom_boxplot(aes(y = TPM, fill = Sex), color = "grey10", outlier.color = "black", outlier.alpha = 0.8,
                 outlier.shape = 21, outlier.size = 3.75, width = 0.8,
                 position = position_dodge(width = 0.8)) +
    stat_summary(aes(y = TPM, fill = Sex), fun.data = "mean_se", geom = "point", shape = 4, stroke = 1.5, 
                 position = position_dodge(width = 0.8)) +
    labs(x = "") +
    scale_fill_manual(values = c("#CC6677", "#88CCEE")) +
    scale_y_continuous(limits = c(-0.1, high_limit), 
                       breaks = seq(0, high_limit, break_unit)) +
    theme(legend.title.align = 0.5, legend.position = "bottom",
          legend.margin = margin(-0.75, 0, 0, 0, "cm"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(color = "grey70"),
          axis.text.x = element_text(size = 11),
          axis.ticks.x = element_blank(),
          plot.margin = margin(0.25, 0.25, 0.75, 0.25, "cm"),
          strip.background = element_rect(fill = rect),
          legend.text = element_text(size = 14),
          legend.key.size = unit(1.5, "cm"),
          strip.text.x = element_text(size = 14)) +
    facet_wrap(~ Human_Symbol) +
    geom_text(data = cld_info, aes(x = Tissue, y = pos, label = Group), size = 5)
}
# 1200, 425 is 4 per page


svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/AGK_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("AGK", 60, 10, pos = c(38, 23, 26, 55, 38, 28, 35, 28, 33, 25, 26, 28))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/AMZ1_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("AMZ1", 8, 2, pos = c(5, 7, 5.7, 3.2, 7.5, 5, 5.5, 7.2, 5, 5.5, 5, 3.6))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/B2M_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("B2M", 8200, 1000, pos = c(6700, 8200, 1500, 1600, 5300, 5300, 7300, 5400, 2700, 1300, 1600, 2400))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/BCAS1_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("BCAS1", 2100, 400, pos = c(900, 950, 500, 320, 1500, 900, 950, 1050, 200, 1300, 1100, 2100))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/BTG1_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("BTG1", 565, 50, pos = c(185, 230, 70, 335, 85, 145, 170, 120, 230, 70, 70, 255))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/CACYBP_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("CACYBP", 565, 50, pos = c(220, 185, 240, 435, 270, 180, 220, 420, 380, 150, 140, 250))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/CARTPT_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("CARTPT", 5150, 1000, pos = c(400, 600, 400, 400, 700, 500, 5150, 800, 600, 400, 400, 500))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/CCDC88C_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("CCDC88C", 50, 10, pos = c(15, 19, 42, 6, 8, 15, 17, 49, 25, 35, 6, 12))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/CELF6_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("CELF6", 120, 20, pos = c(45, 30, 25, 65, 43, 35, 101, 105, 117, 25, 38, 35))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/EGR2_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("EGR2", 205, 20, pos = c(30, 35, 55, 65, 28, 50, 37, 42, 205, 46, 26, 50))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/FABP7_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("FABP7", 265, 50, pos = c(130, 90, 80, 220, 90, 85, 265, 85, 50, 45, 55, 94))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/FKBP4_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("FKBP4", 600, 100, pos = c(240, 190, 320, 560, 380, 220, 240, 430, 530, 280, 190, 350))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/FTH1_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("FTH1", 5200, 1000, pos = c(3800, 3300, 2400, 1500, 4600, 3500, 3000, 3300, 5200, 4000, 2900, 5200))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/GPD1_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("GPD1", 375, 40, pos = c(50, 55, 30, 45, 30, 85, 85, 45, 35, 65, 85, 375))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/GUCY1A3_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("GUCY1A3", 35, 5, pos = c(15, 14, 24, 7, 16, 9, 15, 34.5, 30, 24, 8, 24))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/HAPLN2_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("HAPLN2", 2400, 400, pos = c(820, 360, 360, 320, 320, 740, 670, 420, 200, 950, 930, 2300))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/HSPA8_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("HSPA8", 4650, 750, pos = c(1200, 1000, 1800, 4650, 2000, 1100, 1400, 2450, 3950, 1000, 1200, 1400))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/IRS2_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("IRS2", 105, 20, pos = c(56, 46, 35, 75, 47, 58, 88, 84, 105, 27, 27, 45))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/KIF5A_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("KIF5A", 2575, 500, pos = c(1900, 1400, 1100, 1600, 2575, 1600, 800, 1100, 300, 950, 400, 650))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/LYPD1_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("LYPD1", 72, 10, pos = c(45, 58, 72, 20, 35, 58, 65, 68, 72, 49, 45, 42))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/MBP_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("MBP", 25000, 5000, pos = c(9000, 11000, 4100, 3400, 14000, 9000, 10000, 8000, 2000, 11000, 11000, 24800))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/MOBP_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("MOBP", 3500, 500, pos = c(1250, 1350, 450, 400, 1650, 1300, 1650, 850, 330, 1350, 1530, 3300))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/NTS_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("NTS", 565, 50, pos = c(40, 120, 125, 80, 40, 270, 220, 70, 565, 40, 320, 220))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/PHLDA1_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("PHLDA1", 120, 20, pos = c(32, 35, 32, 31, 40, 43, 52, 45, 115, 32, 32, 32))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/PITPNM3_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("PITPNM3", 150, 25, pos = c(90, 102, 80, 148, 90, 70, 60, 95, 20, 65, 25, 35))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/RGS5_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("RGS5", 365, 50, pos = c(180, 145, 210, 170, 200, 190, 255, 248, 365, 225, 220, 220))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/RPS6KA2_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("RPS6KA2", 210, 50, pos = c(95, 75, 92, 55, 80, 60, 95, 90, 90, 75, 110, 210))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/SOX17_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("SOX17", 30.5, 5, pos = c(11, 20, 9.5, 7, 9, 8.5, 11, 7, 22, 9, 11, 30.5))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/TIPARP_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("TIPARP", 120, 20, pos = c(43, 46, 38, 44, 50, 55, 75, 73, 80, 30, 27, 120))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/TTLL1_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("TTLL1", 50, 10, pos = c(29, 38, 25, 48, 32, 33, 42, 32, 38, 25, 24, 32))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/USP46_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("USP46", 62, 10, pos = c(42, 42, 37, 62, 46, 45, 36, 42, 28, 32, 33, 44))
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Supp_Fig_5_Genes/VIM_bySex_4_1_2022.svg", width = 10.25, height = 5.25)
sex_plot("VIM", 5000, 1000, pos = c(1000, 1200, 1200, 600, 850, 1300, 1100, 3400, 4900, 1200, 700, 1200))
dev.off()

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
#    geom_point(aes(y = Median_TPM), shape = 4, stroke = 1.5, size = 3) +
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

# svglite("G:/Users/Annika/Documents/Figures for Gene Conservation Project/MOBP_Tissue_1_24_2022.svg", width = 2.6, height = 5.4)
# vert_gene_plot("MOBP", 1500, 500, "#DE8971", "black")
# dev.off()
# svglite("G:/Users/Annika/Documents/Figures for Gene Conservation Project/KIF5A_Tissue_1_24_2022.svg", width = 2.6, height = 5.4)
# vert_gene_plot("KIF5A", 1200, 400, "#DE8971", "white")
# dev.off()
# svglite("G:/Users/Annika/Documents/Figures for Gene Conservation Project/MBP_Tissue_1_24_2022.svg", width = 2.6, height = 5.4)
# vert_gene_plot("MBP", 10000, 2000, "#DE8971", "white")
# dev.off()
# svglite("G:/Users/Annika/Documents/Figures for Gene Conservation Project/HAPLN2_Tissue_1_24_2022.svg", width = 2.6, height = 5.4)
# vert_gene_plot("HAPLN2", 800, 200, "#DE8971", "white")
# dev.off()


# Comparison genes
# svglite("G:/Users/Annika/Documents/Figures for Gene Conservation Project/CREB1_Tissue_1_24_2022.svg", width = 2.5, height = 5.4)
# vert_gene_plot("CREB1", 18, 3, "grey80", "black")
# dev.off()
# svglite("G:/Users/Annika/Documents/Figures for Gene Conservation Project/C1QL2_Tissue_1_24_2022.svg", width = 2.5, height = 5.4)
# vert_gene_plot("C1QL2", 24, 4, "grey80", "white")
# dev.off()
# svglite("G:/Users/Annika/Documents/Figures for Gene Conservation Project/DRD3_Tissue_1_24_2022.svg", width = 2.5, height = 5.4)
# vert_gene_plot("DRD3", 4, 1, "grey80", "black")
# dev.off()

# Craving genes
# svglite("G:/Users/Annika/Documents/Figures for Gene Conservation Project/SYT1_Tissue.svg", width = 2.65, height = 5.4)
# vert_gene_plot("SYT1", 350, 75, "#DE8971", "black")
# dev.off()

