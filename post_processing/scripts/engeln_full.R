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

# ENSEMBL
library(biomaRt)
library(org.Mm.eg.db)
library(org.Rn.eg.db)

# Matching samples
library(MatchIt)

# Other
library(clipr)


## SET PARAMETERS ----
# Change the tidyverse functions
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
arrange <- dplyr::arrange

options(dplyr.summarise.inform = FALSE)

# Reduce chances of scientific notation
options(scipen = 9999)

# Update and set ggplot theme globally
new_theme <- theme_bw(base_size = 14) %+replace%
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))
theme_set(new_theme)

# File locations
cm_file <- paste0("D:/GitHub Repositories/RodentAddiction/external_downloads/1-s2.0-S2211124720312973-mmc4.xlsx")


## GET MARTS ----
# Human (GRCh38.p13)
grch38 <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", 
                     host = "https://apr2022.archive.ensembl.org")
# Mouse (GRCm39)
grcm39 <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl",
                     host = "https://apr2022.archive.ensembl.org")
# Rat (mRatBN7.2)
mRatBN7.2 <- useEnsembl("ensembl", dataset = "rnorvegicus_gene_ensembl", 
                        host = "https://apr2022.archive.ensembl.org")


### OLDER MART FOR dN/dS (JANUARY 2020) ----------
# Human (Ensembl release 99)
hs_ens99 <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                    host = "https://jan2020.archive.ensembl.org")


## IMPORT DATA ----
engeln <- read_excel("C:/Users/Annika/Downloads/41380_2022_1668_MOESM2_ESM.xlsx", 
                        sheet = 1) %>% filter(FDR < 0.05)
mouse_gene_ids <- str_to_upper(engeln$Feature.ID)



### GET HOMOLOGS ----
# Get dataframe with all relevant mouse or rat info
rodent_bm <- getBM(attributes = c("external_gene_name", # Gene symbol
                                  "ensembl_gene_id", # Mouse ENSEMBL
                                  "rnorvegicus_homolog_ensembl_gene", # Rat
                                  "hsapiens_homolog_ensembl_gene"), # Human
                   filter = "ensembl_gene_id", 
                   values = mouse_gene_ids,
                   mart = grcm39) %>%
  rename(BM_Gene_Name = external_gene_name,
         Mouse_ID = ensembl_gene_id,
         Rat_ID = rnorvegicus_homolog_ensembl_gene,
         Human_ID = hsapiens_homolog_ensembl_gene) %>%
  mutate(across(everything(), function(x) na_if(x, "")),
         Has_Rat_Hom = case_when(!is.na(Rat_ID) ~ "Yes",
                                 TRUE ~ "No"),
         Has_Human_Hom = case_when(!is.na(Human_ID) ~ "Yes",
                                   TRUE ~ "No"))
rodent_ids <- unique(rodent_bm$Mouse_ID)

# Extract human and rat orthologs
human_ids <- rodent_bm %>% filter(!is.na(Human_ID)) %>% pull(Human_ID) %>% unique()
rat_ids <- rodent_bm %>% filter(!is.na(Rat_ID)) %>% pull(Rat_ID) %>% unique()

# Rename mouse
mouse_ids <- rodent_ids

# How many genes do or do not have human or rat orthologs?
# Human
has_human_homolog <- length(rodent_bm %>% filter(Has_Human_Hom == "Yes") %>% pull(Mouse_ID))
no_human_homolog <- length(rodent_bm %>% filter(Has_Human_Hom == "No") %>% pull(Mouse_ID))
message(paste("Mouse genes with a human ortholog:", has_human_homolog))
message(paste("Mouse genes without a human ortholog:", no_human_homolog))
# Rat
has_rat_homolog <- length(rodent_bm %>% filter(Has_Rat_Hom == "Yes") %>% pull(Mouse_ID))
no_rat_homolog <- length(rodent_bm %>% filter(Has_Rat_Hom == "No") %>% pull(Mouse_ID))
message(paste("Mouse genes with a rat ortholog:", has_rat_homolog))
message(paste("Mouse genes without a rat ortholog:", no_rat_homolog))


## SEQUENCE SIMILARITY ----
# Get sequence similarity data from newer version of Ensembl using human IDs
new_ortho <- getBM(attributes = c("external_gene_name", # Gene symbol
                                  "ensembl_gene_id", # Human ENSEMBL
                                  "mmusculus_homolog_ensembl_gene", # Mouse ID
                                  "rnorvegicus_homolog_ensembl_gene", # Rat ID
                                  # Mouse homology
                                  "mmusculus_homolog_subtype",
                                  "mmusculus_homolog_orthology_type",
                                  "mmusculus_homolog_perc_id",
                                  "mmusculus_homolog_perc_id_r1",
                                  "mmusculus_homolog_orthology_confidence",
                                  # Rat homology
                                  "rnorvegicus_homolog_subtype",
                                  "rnorvegicus_homolog_orthology_type",
                                  "rnorvegicus_homolog_perc_id",
                                  "rnorvegicus_homolog_perc_id_r1",
                                  "rnorvegicus_homolog_orthology_confidence"),
                   mart = grch38) %>%
  rename(BM_Human_Gene_Name = external_gene_name,
         Human_ID = ensembl_gene_id,
         Mouse_ID = mmusculus_homolog_ensembl_gene,
         Rat_ID = rnorvegicus_homolog_ensembl_gene,
         Mouse_Ortho_Subtype = mmusculus_homolog_subtype,
         Mouse_Ortho_Type = mmusculus_homolog_orthology_type,
         MouseToHuman_PercMatch = mmusculus_homolog_perc_id,
         HumanToMouse_PercMatch = mmusculus_homolog_perc_id_r1,
         Mouse_Ortho_Conf = mmusculus_homolog_orthology_confidence,
         Rat_Ortho_Subtype = rnorvegicus_homolog_subtype,
         Rat_Ortho_Type = rnorvegicus_homolog_orthology_type,
         RatToHuman_PercMatch = rnorvegicus_homolog_perc_id,
         HumanToRat_PercMatch = rnorvegicus_homolog_perc_id_r1,
         Rat_Ortho_Conf = rnorvegicus_homolog_orthology_confidence) %>%
  mutate(across(everything(), function(x) na_if(x, "NaN"))) %>%
  mutate(across(everything(), function(x) na_if(x, "")))

# Separate into mouse information
new_mouse_ortho <- new_ortho %>%
  select(-contains("Rat")) %>%
  unique()

# Separate by mouse and rat
new_mouse_ortho <- new_ortho %>%
  select(-contains("Rat")) %>%
  unique() %>%
  mutate(Gene_Group = fct_relevel(ifelse(
    Human_ID %in% human_ids & Mouse_ID %in% mouse_ids, 
    "My Genes", "All Other Orthologs"),
    levels = c("My Genes", "All Other Orthologs")))

new_rat_ortho <- new_ortho %>%
  select(-contains("Mouse")) %>%
  unique() %>%
  mutate(Gene_Group = fct_relevel(ifelse(
    Human_ID %in% human_ids & Rat_ID %in% rat_ids, 
    "My Genes", "All Other Orthologs"),
    levels = c("My Genes", "All Other Orthologs")))

# Join dataframes for reference later
full_new_ortho <- Reduce(function(x, y) 
  full_join(x, y), list(new_mouse_ortho, new_rat_ortho)) %>%
  relocate(Gene_Group, .before = everything()) %>%
  unique() %>%
  select(-Mouse_Ortho_Subtype, -Rat_Ortho_Subtype)


#### STATISTICS ----
# Stats comparing each dataset to the ones in no datasets
stat_mouse_perc <- new_mouse_ortho %>%
  filter(!is.na(MouseToHuman_PercMatch), Mouse_Ortho_Type == "ortholog_one2one") %>%
  wilcox_test(MouseToHuman_PercMatch ~ Gene_Group) %>%
  filter(group2 == "All Other Orthologs") %>%
  rename(Comparison = ".y.")
stat_rat_perc <- new_rat_ortho %>%
  filter(!is.na(RatToHuman_PercMatch), Rat_Ortho_Type == "ortholog_one2one") %>%
  wilcox_test(RatToHuman_PercMatch ~ Gene_Group) %>%
  filter(group2 == "All Other Orthologs") %>%
  rename(Comparison = ".y.")

# Medians
mouse_perc_medians <- new_mouse_ortho %>%
  filter(!is.na(MouseToHuman_PercMatch), Mouse_Ortho_Type == "ortholog_one2one") %>%
  group_by(Gene_Group) %>%
  summarize(MouseToHuman_PercMatch_Median = median(MouseToHuman_PercMatch))
rat_perc_medians <- new_rat_ortho %>%
  filter(!is.na(RatToHuman_PercMatch), Rat_Ortho_Type == "ortholog_one2one") %>%
  group_by(Gene_Group) %>%
  summarize(RatToHuman_PercMatch_Median = median(RatToHuman_PercMatch))

# View stats and medians
stat_perc <- full_join(stat_mouse_perc, stat_rat_perc)
medians_perc <- full_join(mouse_perc_medians, rat_perc_medians)
stat_perc
medians_perc



#### CHECK OUTLIERS ----
lowest_sim_mouse <- full_new_ortho %>% 
  filter(!is.na(MouseToHuman_PercMatch), Mouse_Ortho_Type == "ortholog_one2one") %>%
  filter(Gene_Group == "My Genes") %>%
  arrange(MouseToHuman_PercMatch) %>%
  select(BM_Human_Gene_Name, Mouse_ID, MouseToHuman_PercMatch) %>%
  head()
lowest_sim_mouse

lowest_sim_rat <- full_new_ortho %>% 
  filter(!is.na(RatToHuman_PercMatch), Rat_Ortho_Type == "ortholog_one2one") %>%
  filter(Gene_Group == "My Genes") %>%
  arrange(RatToHuman_PercMatch) %>%
  select(BM_Human_Gene_Name, Rat_ID, RatToHuman_PercMatch) %>%
  head()
lowest_sim_rat


#### PLOT ----
mouse_perc_label_list <- c("ERVW-1", "CDR1", "CCDC187", "RESP18", "ZDBF2")
rat_perc_label_list <- c("GPR179", "ERVW-1", "CDR1", "CCDC187", "RESP18")

new_mouse_ortho_for_plot <- new_mouse_ortho %>% 
  filter(!is.na(MouseToHuman_PercMatch), Mouse_Ortho_Type == "ortholog_one2one")
mouse_perc_plot <- new_mouse_ortho_for_plot %>%
  ggplot(aes(x = Gene_Group, y = MouseToHuman_PercMatch)) +
  geom_violin(aes(fill = Gene_Group), width = 0.8) +  
  geom_boxplot(data = new_mouse_ortho_for_plot %>% filter(Gene_Group == "My Genes"),
               fill = "white", color = "black", width = 0.15,
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  geom_boxplot(data = new_mouse_ortho_for_plot %>% filter(Gene_Group != "My Genes"),
               fill = "white", color = "black", width = 0.08,
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  labs(x = "", y = "Sequence Match to Human (%)", title = "Mouse-Human\n") +
  scale_x_discrete(labels= c("Engeln", "All Other Orthologs")) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_fill_manual(values = c("#296D98", "grey80")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey85"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14)) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(BM_Human_Gene_Name %in% mouse_perc_label_list &
                                             Gene_Group != "All Other Orthologs",
                                           BM_Human_Gene_Name, NA_character_)),
                   aes(label = label), color = "black", fill = "#358dc5",
                   nudge_y = 0.01, nudge_x = 0.01,
                   box.padding = 0.5, fontface = "bold", size = 4, seed = 58)

# Rat-Human
new_rat_ortho_for_plot <- new_rat_ortho %>% 
  filter(!is.na(RatToHuman_PercMatch), Rat_Ortho_Type == "ortholog_one2one")
rat_perc_plot <- new_rat_ortho_for_plot %>%
  filter(!is.na(RatToHuman_PercMatch), Rat_Ortho_Type == "ortholog_one2one") %>%
  ggplot(aes(x = Gene_Group, y = RatToHuman_PercMatch)) +
  geom_violin(aes(fill = Gene_Group), width = 0.8) +  
  geom_boxplot(data = new_rat_ortho_for_plot %>% filter(Gene_Group == "My Genes"),
               fill = "white", color = "black", width = 0.15,
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  geom_boxplot(data = new_rat_ortho_for_plot %>% filter(Gene_Group != "My Genes"),
               fill = "white", color = "black", width = 0.08,
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  labs(x = "", y = "Sequence Match to Human (%)", title = "Rat-Human\n") +
  scale_x_discrete(labels= c("Engeln", "All Other Orthologs")) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_fill_manual(values = c("#296D98", "grey80")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey85"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        axis.title.y = element_text(size = 14, color = "white")) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(BM_Human_Gene_Name %in% rat_perc_label_list &
                                             Gene_Group != "All Other Orthologs",
                                           BM_Human_Gene_Name, NA_character_)),
                   aes(label = label), color = "black", fill = "#358dc5",
                   nudge_y = 0.01, nudge_x = 0.01,
                   box.padding = 0.5, fontface = "bold", size = 4, seed = 58)

perc_plot <- ggarrange(mouse_perc_plot, rat_perc_plot)
perc_plot


## DN/DS ----
old_ortho <- getBM(attributes = c("external_gene_name", # Gene symbol
                                  "ensembl_gene_id", # Human ENSEMBL
                                  "mmusculus_homolog_ensembl_gene", # Mouse ID
                                  "rnorvegicus_homolog_ensembl_gene", # Rat ID
                                  # Mouse homology
                                  "mmusculus_homolog_dn",
                                  "mmusculus_homolog_ds",
                                  "mmusculus_homolog_subtype",
                                  "mmusculus_homolog_orthology_type",
                                  "mmusculus_homolog_perc_id",
                                  "mmusculus_homolog_perc_id_r1",
                                  "mmusculus_homolog_orthology_confidence",
                                  # Rat homology
                                  "rnorvegicus_homolog_dn",
                                  "rnorvegicus_homolog_ds",
                                  "rnorvegicus_homolog_subtype",
                                  "rnorvegicus_homolog_orthology_type",
                                  "rnorvegicus_homolog_perc_id",
                                  "rnorvegicus_homolog_perc_id_r1",
                                  "rnorvegicus_homolog_orthology_confidence"),
                   mart = hs_ens99) %>%
  rename(BM_Human_Gene_Name = external_gene_name,
         Human_ID = ensembl_gene_id,
         Mouse_ID = mmusculus_homolog_ensembl_gene,
         Rat_ID = rnorvegicus_homolog_ensembl_gene,
         Mouse_DN = mmusculus_homolog_dn,
         Mouse_DS = mmusculus_homolog_ds,
         Mouse_Ortho_Subtype = mmusculus_homolog_subtype,
         Mouse_Ortho_Type = mmusculus_homolog_orthology_type,
         MouseToHuman_PercMatch = mmusculus_homolog_perc_id,
         HumanToMouse_PercMatch = mmusculus_homolog_perc_id_r1,
         Mouse_Ortho_Conf = mmusculus_homolog_orthology_confidence,
         Rat_DN = rnorvegicus_homolog_dn,
         Rat_DS = rnorvegicus_homolog_ds,
         Rat_Ortho_Subtype = rnorvegicus_homolog_subtype,
         Rat_Ortho_Type = rnorvegicus_homolog_orthology_type,
         RatToHuman_PercMatch = rnorvegicus_homolog_perc_id,
         HumanToRat_PercMatch = rnorvegicus_homolog_perc_id_r1,
         Rat_Ortho_Conf= rnorvegicus_homolog_orthology_confidence) %>%
  mutate(Mouse_DNDS = Mouse_DN/Mouse_DS,
         Rat_DNDS = Rat_DN/Rat_DS) %>%
  mutate(across(everything(), function(x) na_if(x, "NaN")))

# Separate by mouse and rat
old_mouse_ortho <- old_ortho %>%
  select(-contains("Rat")) %>%
  unique() %>%
  mutate(Gene_Group = fct_relevel(ifelse(
    Human_ID %in% human_ids & Mouse_ID %in% mouse_ids, 
    "My Genes", "All Other Orthologs"),
    levels = c("My Genes", "All Other Orthologs")))

old_rat_ortho <- old_ortho %>%
  select(-contains("Mouse")) %>%
  unique() %>%
  mutate(Gene_Group = fct_relevel(ifelse(
    Human_ID %in% human_ids & Rat_ID %in% rat_ids, 
    "My Genes", "All Other Orthologs"),
    levels = c("My Genes", "All Other Orthologs")))

# Join dataframes for reference later
full_old_ortho <- Reduce(function(x, y) 
  full_join(x, y), list(old_mouse_ortho, old_rat_ortho)) %>%
  relocate(Gene_Group, .before = everything()) %>%
  unique() %>%
  select(-Mouse_Ortho_Subtype, -Rat_Ortho_Subtype)


#### STATISTICS ----
# Stats comparing dataset to other genes
stat_mouse_dnds <- old_mouse_ortho %>%
  filter(!is.na(Mouse_DNDS), Mouse_DNDS != Inf, 
         Mouse_Ortho_Type == "ortholog_one2one") %>%
  wilcox_test(Mouse_DNDS ~ Gene_Group)
stat_rat_dnds <- old_rat_ortho %>%
  filter(!is.na(Rat_DNDS), Rat_DNDS != Inf, 
         Rat_Ortho_Type == "ortholog_one2one") %>%
  wilcox_test(Rat_DNDS ~ Gene_Group)

# Also calculate mouse statistics removing CDR1
stat_mouse_dnds2 <- old_mouse_ortho %>%
  filter(!is.na(Mouse_DNDS), Mouse_DNDS != Inf, 
         Mouse_Ortho_Type == "ortholog_one2one",
         BM_Human_Gene_Name != "CDR1") %>%
  wilcox_test(Mouse_DNDS ~ Gene_Group)
mouse_dnds_medians2 <- old_mouse_ortho %>%
  filter(!is.na(Mouse_DNDS), Mouse_DNDS != Inf, 
         Mouse_Ortho_Type == "ortholog_one2one",
         BM_Human_Gene_Name != "CDR1") %>%
  group_by(Gene_Group) %>%
  summarize(DNDS_withMouse_Median = median(Mouse_DNDS))

# Also calculate mouse statistics changing CDR1 category
stat_mouse_dnds3 <- old_mouse_ortho %>%
  filter(!is.na(Mouse_DNDS), Mouse_DNDS != Inf, 
         Mouse_Ortho_Type == "ortholog_one2one") %>%
  mutate(Gene_Group = fct_relevel(case_when(BM_Human_Gene_Name == "CDR1" ~ "My Genes",
                                            TRUE ~ as.character(Gene_Group)),
                                  levels = c("My Genes", "All Other Orthologs"))) %>%
  wilcox_test(Mouse_DNDS ~ Gene_Group)
mouse_dnds_medians3 <- old_mouse_ortho %>%
  filter(!is.na(Mouse_DNDS), Mouse_DNDS != Inf, 
         Mouse_Ortho_Type == "ortholog_one2one") %>%
  mutate(Gene_Group = fct_relevel(case_when(BM_Human_Gene_Name == "CDR1" ~ "My Genes",
                                            TRUE ~ as.character(Gene_Group)),
         levels = c("My Genes", "All Other Orthologs"))) %>%
  group_by(Gene_Group) %>%
  summarize(DNDS_withMouse_Median = median(Mouse_DNDS))

other <- stat_mouse_dnds
eng <- stat_mouse_dnds3
neither <- stat_mouse_dnds2

other <- mouse_dnds_medians
eng <- mouse_dnds_medians3
neither <- mouse_dnds_medians2

# Medians
mouse_dnds_medians <- old_mouse_ortho %>%
  filter(!is.na(Mouse_DNDS), Mouse_DNDS != Inf, Mouse_Ortho_Type == "ortholog_one2one") %>%
  group_by(Gene_Group) %>%
  summarize(DNDS_withMouse_Median = median(Mouse_DNDS))
rat_dnds_medians <- old_rat_ortho %>%
  filter(!is.na(Rat_DNDS), Rat_DNDS != Inf, Rat_Ortho_Type == "ortholog_one2one") %>%
  group_by(Gene_Group) %>%
  summarize(DNDS_withRat_Median = median(Rat_DNDS))

# View stats and medians
stat_dnds <- full_join(stat_mouse_dnds, stat_rat_dnds)
medians_dnds <- full_join(mouse_dnds_medians, rat_dnds_medians)
stat_dnds
medians_dnds


#### CHECK OUTLIERS ----
highest_dnds_mouse <- full_old_ortho %>% 
  filter(!is.na(Mouse_DNDS), Mouse_Ortho_Type == "ortholog_one2one") %>%
  filter(Gene_Group == "My Genes") %>%
  arrange(desc(Mouse_DNDS)) %>%
  select(BM_Human_Gene_Name, Mouse_ID, Mouse_DNDS) %>%
  mutate(outlier = is_outlier(Mouse_DNDS))
highest_dnds_mouse

highest_dnds_rat <- full_old_ortho %>% 
  filter(!is.na(Rat_DNDS), Rat_Ortho_Type == "ortholog_one2one") %>%
  filter(Gene_Group == "My Genes") %>%
  arrange(desc(Rat_DNDS)) %>%
  select(BM_Human_Gene_Name, Rat_ID, Rat_DNDS) %>%
  mutate(outlier = is_outlier(Rat_DNDS))
highest_dnds_rat


#### PLOT ----
mouse_dnds_label_list <- c("SPINK8", "ZDBF2", "RESP18")
rat_dnds_label_list <- c("SPINK8", "ZDBF2", "RESP18")

old_mouse_ortho_for_plot <- old_mouse_ortho %>% 
  filter(!is.na(Mouse_DNDS), Mouse_Ortho_Type == "ortholog_one2one")
mouse_dnds_plot <- old_mouse_ortho_for_plot %>%
  ggplot(aes(x = Gene_Group, y = Mouse_DNDS)) +
  geom_violin(aes(fill = Gene_Group), width = 0.8) +  
  geom_boxplot(data = old_mouse_ortho_for_plot %>% filter(Gene_Group == "My Genes"),
               fill = "white", color = "black", width = 0.15,
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  geom_boxplot(data = old_mouse_ortho_for_plot %>% filter(Gene_Group != "My Genes"),
               fill = "white", color = "black", width = 0.08,
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  labs(x = "", y = "dN/dS with Human Ortholog\n",
       title = "Mouse-Human\n") +
  scale_x_discrete(expand = c(0, 0.55), labels= c("Engeln", "All Other Orthologs")) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5), expand = c(0, 0)) +
  scale_fill_manual(values = c("#296D98", "grey80")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey85"),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 13)) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(BM_Human_Gene_Name %in% mouse_dnds_label_list, 
                                           BM_Human_Gene_Name, NA_character_)),
                   aes(label = label), color = "black", fill = "#358dc5",
                   nudge_y = 0.01, nudge_x = 0.01,
                   box.padding = 0.5, fontface = "bold", size = 4, seed = 19010)

old_rat_ortho_for_plot <- old_rat_ortho %>% 
  filter(!is.na(Rat_DNDS), Rat_Ortho_Type == "ortholog_one2one")
rat_dnds_plot <- old_rat_ortho_for_plot %>%
  ggplot(aes(x = Gene_Group, y = Rat_DNDS)) +
  geom_violin(aes(fill = Gene_Group), width = 0.75) +  
  geom_boxplot(data = old_rat_ortho_for_plot %>% filter(Gene_Group == "My Genes"),
               fill = "white", color = "black", width = 0.15,
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  geom_boxplot(data = old_rat_ortho_for_plot %>% filter(Gene_Group != "My Genes"),
               fill = "white", color = "black", width = 0.08,
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  labs(x = "", y = "dN/dS with Human Ortholog\n",
       title = "Rat-Human\n") +
  scale_x_discrete(expand = c(0, 0.55), labels= c("Engeln", "All Other Orthologs")) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5), expand = c(0, 0)) +
  scale_fill_manual(values = c("#296D98", "grey80")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey85"),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 13, color = "white")) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(BM_Human_Gene_Name %in% rat_dnds_label_list, 
                                           BM_Human_Gene_Name, NA_character_)),
                   aes(label = label), color = "black", fill = "#358dc5",
                   nudge_y = 0.01, nudge_x = 0.01,
                   box.padding = 0.5, fontface = "bold", size = 4, seed = 10)

# Set up plots
dnds_plot <- ggarrange(mouse_dnds_plot, rat_dnds_plot)
dnds_plot



## DEVELOPMENTAL CONSERVATION ----
# "Developmental gene expression differences between humans & mammalian models"
# Cardoso-Moriera et al. (2020)
# Table S3A. Comparison of brain temporal trajectories

# Includes human, mouse, rat, rabbit, rhesus
cons_orig <- read_excel(cm_file, sheet = 1, skip = 2)

# Tidy and create new columns
# Note that with the new column "Conservation" we DO NOT have enough information
# to determine whether a gene's temporal developmental expression is shared
# between mouse and rat.
cons2020 <- cons_orig %>%
  select(-c(starts_with("Prob"), Rabbit, Rhesus)) %>%
  # Categorize by conservation (see above note)
  mutate(Conservation = fct_relevel(case_when(
    (Mouse == "Same" & Rat == "Same") ~ "HMR",
    (Mouse == "Same" & Rat == "Different") ~ "HM",
    (Mouse == "Different" & Rat == "Same") ~ "HR",
    (Mouse == "Different" & Rat == "Different") ~ "H",
    TRUE ~ NA_character_),
    c("HMR", "HM", "HR", "H"))) %>%
  mutate(across(everything(), function(x) na_if(x, "NA")))

cons2020_with_incompletes <- cons_orig %>%
  select(-c(starts_with("Prob"), Rabbit, Rhesus)) %>%
  # Categorize by conservation (see above note)
  mutate(Conservation = fct_relevel(case_when(
    (Mouse == "Same" & Rat == "Same") ~ "HMR",
    (Mouse == "Same" & Rat == "Different") ~ "HM",
    (Mouse == "Different" & Rat == "Same") ~ "HR",
    (Mouse == "Different" & Rat == "Different") ~ "H",
    # NA genes were not tested because they didn't have complete data.
    (Mouse == "Same" & Rat == "NA") ~ "HM*",
    (Mouse == "Different" & Rat == "NA") ~ "H*",
    (Mouse == "NA" & Rat == "Same") ~ "HR*",  
    (Mouse == "NA" & Rat == "Different") ~ "H*",    
    TRUE ~ NA_character_),
    c("HMR", "HM", "HR", "H", "HM*", "HR*", "H*"))) %>%
  mutate(across(everything(), function(x) na_if(x, "NA")))

cons2020 <- cons2020 %>%
  mutate(Gene_Group = fct_relevel(ifelse(
    Human_ID %in% human_ids, "My Genes", "All Other Orthologs"),
    levels = c("My Genes", "All Other Orthologs")))

cons2020_with_incompletes <- cons2020_with_incompletes %>%
  mutate(Gene_Group = fct_relevel(ifelse(
    Human_ID %in% human_ids, "My Genes", "All Other Orthologs"),
    levels = c("My Genes", "All Other Orthologs")))

# How many genes are in each Gene_Group in the cons2020 dataset?
cons2020 %>%
  group_by(Gene_Group) %>%
  count()

# Genes not in dataset
human_ids[human_ids %in% (cons2020 %>% 
                            filter(Gene_Group == "My Genes") %>% 
                            pull(Human_ID))]

# Which genes are HM, HR, or H?
cons2020 %>% filter(Gene_Group == "My Genes", Conservation != "HMR")


#### STATISTICS ----
# Fisher's Exact Test
fisher_table <- table(cons2020$Gene_Group, cons2020$Conservation)
fisher_test <- fisher.test(fisher_table, simulate.p.value = FALSE)
fisher_table
fisher_test


#### PLOT ----
cm_plot <- cons2020 %>%
  filter(!is.na(Conservation)) %>%
  ggplot(aes(x = Conservation, fill = Gene_Group)) +
  geom_bar(aes(group = 1), show.legend = FALSE,
           position = position_dodge2(preserve = "single"),
           color = "black") +
  scale_fill_manual(values = c("#296D98", "grey80")) +
  labs(y = "Count", 
       title = "Developmental Trajectories in the Forebrain\n") +
  theme(axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(face = "bold")) +
  facet_wrap(~ Gene_Group, scales = "free_y", nrow = 1, 
             labeller = as_labeller(c(`My Genes` = "Engeln",
                                      `All Other Orthologs` = "All Other Orthologs")))

cm_gtable <- ggplot_gtable(ggplot_build(cm_plot))
striprt <- which(grepl('strip-r', cm_gtable$layout$name) | grepl('strip-t', cm_gtable$layout$name))
fills <- rep(c("#358dc5", "grey80"), 2)
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', cm_gtable$grobs[[i]]$grobs[[1]]$childrenOrder))
  cm_gtable$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid::grid.draw(cm_gtable)


#### GET CONSERVATION DATAFRAME ----
# Get dataframes with just candidate genes
full_old_ortho_candidates <- full_old_ortho %>% filter(Gene_Group == "My Genes")
full_new_ortho_candidates <- full_new_ortho %>% filter(Gene_Group == "My Genes")

# For some genes, the RAT orthologs changed between the new and old Ensembl 
# versions. The mouse orthologs are all the same (Any_Mouse_Same == TRUE for all
# genes). Have separate columns for new and old rat orthologs.
full_all_cons_df <- full_old_ortho_candidates %>%
  select(-contains(c("Perc", "Ortho_Conf")), -BM_Human_Gene_Name) %>%
  rename(Ensembl_99_Mouse_Ortho_Type = Mouse_Ortho_Type,
         Ensembl_99_Rat_Ortho_Type = Rat_Ortho_Type,
         Ensembl_99_Rat_ID = Rat_ID,
         Ensembl_99_Mouse_ID = Mouse_ID) %>%
  mutate(across(everything(), function(x) na_if(x, ""))) %>%
  full_join(full_new_ortho_candidates) %>%
  rename(Ensembl_106_Mouse_Ortho_Type = Mouse_Ortho_Type,
         Ensembl_106_Rat_Ortho_Type = Rat_Ortho_Type,
         Ensembl_106_Rat_ID = Rat_ID,
         Ensembl_106_Mouse_ID = Mouse_ID) %>%
  select(-contains("Ortho_Conf")) %>%
  select(-c(Mouse_DN, Mouse_DS, Rat_DN, Rat_DS)) %>%
  relocate(Ensembl_99_Rat_Ortho_Type, .after = Mouse_DNDS) %>%
  relocate(BM_Human_Gene_Name, .after = Gene_Group) %>%
  mutate(Mouse_Same_ID = ifelse(Ensembl_99_Mouse_ID == Ensembl_106_Mouse_ID, 
                                TRUE, FALSE),
         Rat_Same_ID = ifelse(Ensembl_99_Rat_ID == Ensembl_106_Rat_ID, 
                              TRUE, FALSE),
         Mouse_Same_Ortho = case_when(Ensembl_99_Mouse_Ortho_Type == 
                                        Ensembl_106_Mouse_Ortho_Type ~ TRUE, 
                                      Ensembl_99_Mouse_Ortho_Type != 
                                        Ensembl_106_Mouse_Ortho_Type |
                                        is.na(Ensembl_99_Mouse_Ortho_Type) |
                                        is.na(Ensembl_106_Mouse_Ortho_Type) 
                                      ~ FALSE),
         Rat_Same_Ortho = case_when(Ensembl_99_Rat_Ortho_Type == 
                                      Ensembl_106_Rat_Ortho_Type ~ TRUE, 
                                    Ensembl_99_Rat_Ortho_Type != 
                                      Ensembl_106_Rat_Ortho_Type |
                                      is.na(Ensembl_99_Rat_Ortho_Type) |
                                      is.na(Ensembl_106_Rat_Ortho_Type) 
                                    ~ FALSE)) %>%
  group_by(Human_ID) %>%
  mutate(Any_Mouse_Same = ifelse(sum(Mouse_Same_ID) >= 1, TRUE, FALSE),
         Any_Rat_Same = ifelse(sum(Rat_Same_ID) >= 1, TRUE, FALSE),
         Any_Mouse_Ortho_Same = ifelse(sum(Mouse_Same_Ortho) >= 1, TRUE, FALSE),
         Any_Rat_Ortho_Same = ifelse(sum(Rat_Same_Ortho) >= 1, TRUE, FALSE),
         Ensembl_99_Mouse_ID = "Same ID as Ensembl release 106",
         Ensembl_99_Rat_ID = case_when(Rat_Same_ID == 1 ~ "Same ID as Ensembl release 106",
                                       Any_Rat_Same == TRUE ~ Ensembl_106_Rat_ID,
                                       TRUE ~ NA_character_)) %>%
  ungroup() %>%
  filter(Rat_Same_ID == TRUE | (Rat_Same_ID == FALSE & Any_Rat_Same == FALSE) |
           is.na(Any_Rat_Same)) %>%
  mutate(across(everything(), function(x) 
    gsub("ortholog_one2one", "One-to-one ortholog", x))) %>%
  mutate(across(everything(), function(x) 
    gsub("ortholog_many2many", "Many-to-many orthologs", x))) %>%
  mutate(across(contains("Mouse_Ortho_Type"), function(x) 
    gsub("ortholog_one2many", "Many mouse to one human ortholog", x))) %>%
  mutate(across(contains("Rat_Ortho_Type"), function(x) 
    gsub("ortholog_one2many", "Many rat to one human ortholog", x))) %>%
  mutate(across(everything(), function(x) na_if(x, ""))) %>%
  left_join(cons2020_with_incompletes) %>%
  droplevels() %>%
  select(c(Gene_Group, BM_Human_Gene_Name, Human_ID, Conservation, 
           Ensembl_106_Mouse_ID, Ensembl_106_Mouse_Ortho_Type, 
           MouseToHuman_PercMatch, HumanToMouse_PercMatch, Ensembl_99_Mouse_ID, 
           Ensembl_99_Mouse_Ortho_Type, Mouse_DNDS, Any_Mouse_Same, 
           Ensembl_106_Rat_ID, Ensembl_106_Rat_Ortho_Type, RatToHuman_PercMatch, 
           HumanToRat_PercMatch, Ensembl_99_Rat_ID, Ensembl_99_Rat_Ortho_Type, 
           Rat_DNDS, Any_Rat_Same))



## LOAD IN GTEX DATA (OVERALL MEANS/MEDIANS) -----------
# Get list of files
main_dir <- "G:/Annika/Github Repositories/RodentAddiction/"
gtex_dir <- paste0(main_dir, "post_processing/downloaded_data/gtex/")
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

# Get dataframe for dataset
my_hs_ids <- full_new_ortho %>%
  filter(Gene_Group == "My Genes") %>%
  pull(Human_ID)
my_gtex <- gtex_temp %>%
  filter(Human_ID %in% my_hs_ids) %>%
  mutate(Gene_Group = "My Genes")

# Join for supplementary material
gtex <- Reduce(function(x, y) 
  full_join(x, y), list(my_gtex, gtex_temp)) %>%
  mutate(Gene_Group = fct_relevel(ifelse(is.na(Gene_Group), 
                                         "All Other Genes", Gene_Group),
                                  levels = c("My Genes",  "All Other Genes"))) %>%
  relocate(Gene_Group, .before = everything()) %>%
  unique()

# How many genes in GTEX data?
gtex %>% 
  group_by(Gene_Group) %>% 
  select(Human_ID) %>% 
  unique() %>%
  count()


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
subset_files <- list.files(path = paste0(gtex_dir, "count_subsets_new"), 
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
all_hs_ids <- my_hs_ids

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
genes_sex

genes_tissue <- bysex_tissue_anovas %>% 
  filter(Effect == "Tissue", Sig_p05 == "*")  %>%
  pull(Human_Symbol) %>%
  unique()
genes_tissue

genes_inter <- bysex_tissue_anovas %>% 
  filter(Effect == "Sex:Tissue", Sig_p05 == "*") %>% 
  pull(Human_Symbol) %>%
  unique()
genes_inter


# SEX EFFECT:
# Check male vs. female, which is higher
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
print(paste("Of the", length(genes_sex), "with a significant main effect of Sex,", 
            length(female_higher), "have higher expression in females"))

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
tukey_tissue

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
tukey_sex_tissue

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
  select(Human_Symbol, BrainSpec) %>%
  arrange(desc(BrainSpec)) %>%
  unique() %>%
  pull(Human_Symbol)



### PLOT ----------
# Genes to label
general_label_list <- c("GFAP", "AVP", "C1QL2")
shared_label_list <- c("FAM53B", "CREB1", "DRD3")

# All overlaps in ggrepel
options(ggrepel.max.overlaps = Inf)

bs_plot <- brain_spec %>%
  mutate(Gene_Group = as.character(Gene_Group),
         Gene_Group = fct_relevel(case_when(Gene_Group == "All Other Genes" ~ "All Other Genes",
                                            TRUE ~ Gene_Group),
                                  levels = c("My Genes", "All Other Genes"))) %>%
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
  labs(x = "", y = expression("Brain Specificity (Log "[2]*" FC)"),
       title = "Brain Specificity (GTEx)\n") +
  scale_y_continuous(limits = c(-16, 12), breaks = seq(-16, 12, 4), 
                     expand = c(0, 0.00001)) +
  scale_x_discrete(expand = c(0, 0.55), labels= c("Engeln", "All Other Genes")) +
  scale_fill_manual(values = c("#296D98", "grey80")) +
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
                   box.padding = 0.5, fontface = "bold", size = 4, seed = 58) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(Human_Symbol %in% general_label_list, 
                                           Human_Symbol, NA_character_)),
                   aes(label = label), color = "black", fill = "#358dc5",
                   box.padding = 0.5, fontface = "bold", size = 4, nudge_y = 1, seed = 1000)
bs_plot


#### STATISTICS ----------
# Wilcoxon test between groups
# Are my genes more specific to the brain than other genes?
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


## FINAL DATAFRAME ----
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


## WRITE FILES ----
full_csv_file <- "C:/Users/Annika/Documents/Manuscript - RNA-seq Comparison/full_csv_file_052023.csv"
anova_file <- "C:/Users/Annika/Documents/Manuscript - RNA-seq Comparison/anova_052023.csv"
tissue_effect_file <- "C:/Users/Annika/Documents/Manuscript - RNA-seq Comparison/tissue_effect_052023.csv"
tissue_cld_file <- "C:/Users/Annika/Documents/Manuscript - RNA-seq Comparison/tissue_cld_052023.csv"
highest_cld_file <- "C:/Users/Annika/Documents/Manuscript - RNA-seq Comparison/highest_cld_052023.csv"
interaction_file <- "C:/Users/Annika/Documents/Manuscript - RNA-seq Comparison/interaction_052023.csv"

full_final_df <- full_all_cons_df %>%
  rename(Human_Symbol = "BM_Human_Gene_Name") %>%
  full_join(brain_spec %>% filter(Gene_Group == "My Genes")) %>%
  select(Gene_Group, Human_Symbol, Human_ID, Ensembl_99_Mouse_ID, Ensembl_99_Mouse_Ortho_Type,
         Mouse_DNDS, Ensembl_99_Rat_ID, Ensembl_99_Rat_Ortho_Type, Rat_DNDS,
         Ensembl_106_Mouse_ID, Ensembl_106_Mouse_Ortho_Type, MouseToHuman_PercMatch,
         HumanToMouse_PercMatch, Ensembl_106_Rat_ID, Ensembl_106_Rat_Ortho_Type,
         RatToHuman_PercMatch, HumanToRat_PercMatch, Conservation, BrainSpec, 
         Grand_Tissue_Mean) %>% 
  full_join(highest_region_df) %>%
  unique()

write.csv(full_final_df, full_csv_file)
write.csv(bysex_tissue_anovas, anova_file)
write.csv(tukey_tissue, tissue_effect_file)
write.csv(cld_df, tissue_cld_file)
write.csv(highest_region_df, highest_cld_file)
write.csv(tukey_sex_tissue, interaction_file)



## SAVE PLOTS ----
perc_plot_loc <- "C:/Users/Annika/Documents/new_perc_052023.svg"
dnds_plot_loc <- "C:/Users/Annika/Documents/new_dnds_052023.svg"
cm_plot_loc <- "C:/Users/Annika/Documents/new_cm_052023.svg"
bs_plot_loc <- "C:/Users/Annika/Documents/new_bs_052023.svg"

svglite(perc_plot_loc, width = 9, height = 3.5)
perc_plot
dev.off()

svglite(dnds_plot_loc, width = 9, height = 3.5)
dnds_plot
dev.off()

svglite(cm_plot_loc, width = 4.9, height = 3.5)
grid::grid.draw(cm_gtable)
dev.off()

svglite(bs_plot_loc, width = 3.5, height = 3.5)
bs_plot
dev.off()

# Save gene names
my_genes_list <- list(human_ids, mouse_ids)
names(my_genes_list) <- c("my_hs_ids", "my_mm_ids")
saveRDS(my_genes_list, file = paste0(main_dir, "post_processing/results/other/my_genes_list.RDS"))



