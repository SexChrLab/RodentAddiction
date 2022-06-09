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

# Gene annotations
library(biomaRt)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(limma)

# Statistics
library(rstatix)



## BEFORE USING THIS SCRIPT ----------
# User should alter these variables as necessary
main_dir <- "G:/Noctsol/GitHub Repositories/RodentAddiction/"
de_dir <- paste0(main_dir, "results/gene_info/")
cm_dir <- paste0(main_dir, "post_processing/downloaded_data/cm_2020/")

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



## GET HUMAN MARTS ---------
### NEWEST MART (APRIL 2022) ----------
# Human (GRCh38.p13)
grch38 <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", 
                     host = "https://apr2022.archive.ensembl.org")

### OLDER MART FOR dN/dS (JANUARY 2020) ----------
# Last available mart with DN and DS values: Ensembl Release 99
hs_ens99 <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                    host = "https://jan2020.archive.ensembl.org")


## SEQUENCE SIMILARITY ----------
new_hom <- getBM(attributes = c("external_gene_name", # Gene symbol
                                 "ensembl_gene_id", # Human ENSEMBL
                                 "mmusculus_homolog_ensembl_gene", # Mouse ID
                                 "rnorvegicus_homolog_ensembl_gene", # Rat ID
                                 # Mouse Homology
                                 "mmusculus_homolog_subtype",
                                 "mmusculus_homolog_orthology_type",
                                 "mmusculus_homolog_perc_id",
                                 "mmusculus_homolog_perc_id_r1",
                                 "mmusculus_homolog_orthology_confidence",
                                 # Rat Homology
                                 "rnorvegicus_homolog_subtype",
                                 "rnorvegicus_homolog_orthology_type",
                                 "rnorvegicus_homolog_perc_id",
                                 "rnorvegicus_homolog_perc_id_r1",
                                 "rnorvegicus_homolog_orthology_confidence"),
                  mart = grch38) %>%
  dplyr::rename(BM_Human_Gene_Name = external_gene_name,
                Human_ID = ensembl_gene_id,
                Mouse_ID = mmusculus_homolog_ensembl_gene,
                Rat_ID = rnorvegicus_homolog_ensembl_gene,
                Mouse_Hom_Subtype = mmusculus_homolog_subtype,
                Mouse_Ortho_Type = mmusculus_homolog_orthology_type,
                MouseToHuman_PercMatch = mmusculus_homolog_perc_id,
                HumanToMouse_PercMatch = mmusculus_homolog_perc_id_r1,
                Mouse_Ortho_Conf = mmusculus_homolog_orthology_confidence,
                Rat_Hom_Subtype = rnorvegicus_homolog_subtype,
                Rat_Ortho_Type = rnorvegicus_homolog_orthology_type,
                RatToHuman_PercMatch = rnorvegicus_homolog_perc_id,
                HumanToRat_PercMatch = rnorvegicus_homolog_perc_id_r1,
                Rat_Ortho_Conf = rnorvegicus_homolog_orthology_confidence) %>%
  mutate(across(everything(), function(x) na_if(x, "NaN"))) %>%
  mutate(across(everything(), function(x) na_if(x, "")))

# Separate by mouse and rat
new_mouse_hom <- new_hom %>%
  dplyr::dplyr::select(-contains("Rat")) %>%
  unique()
  
new_rat_hom <- new_hom %>%
  dplyr::dplyr::select(-contains("Mouse")) %>%
  unique()

# Get general datasets for supplementary material
carp_new_hom <- new_hom %>%
  filter(Human_ID %in% sig_carp_hs_ids, Mouse_ID %in% sig_carp_mm_ids) %>%
  mutate(Gene_Group = "Carpenter")
walk_new_hom <- new_hom %>%
  filter(Human_ID %in% sig_walk_hs_ids, Mouse_ID %in% sig_walk_mm_ids) %>%
  mutate(Gene_Group = "Walker")
pow_new_hom <- new_hom %>%
  filter(Human_ID %in% sig_pow_hs_ids, Rat_ID %in% sig_pow_rn_ids) %>%
  mutate(Gene_Group = "Powell")

# Join for supplementary material
full_new_hom <- Reduce(function(x, y) 
  full_join(x, y), list(carp_new_hom, walk_new_hom, pow_new_hom, new_hom)) %>%
  mutate(Gene_Group = fct_relevel(ifelse(is.na(Gene_Group), 
                                         "All Other Homologs", Gene_Group),
                                  levels = c("Carpenter", "Walker", "Powell", "All Other Homologs"))) %>%
  relocate(Gene_Group, .before = everything()) %>%
  unique() %>%
  dplyr::select(-Mouse_Hom_Subtype, -Rat_Hom_Subtype)
dim(full_new_hom) # 194771     13

# Split into mouse and rat
new_mouse_hom <- full_new_hom %>%
  dplyr::select(-contains("Rat")) %>%
  unique()

new_rat_hom <- full_new_hom %>%
  dplyr::select(-contains("Mouse")) %>%
  unique()
  
# Get dataframes for analysis
# Mouse
carp_new_mouse_hom <- new_hom %>%
  filter(Human_ID %in% sig_carp_hs_ids, Mouse_ID %in% sig_carp_mm_ids,
         Mouse_Ortho_Type == "ortholog_one2one") %>%
  dplyr::select(-contains("Rat")) %>%
  unique() %>%
  mutate(Gene_Group = "Carpenter")
walk_new_mouse_hom <- new_hom %>%
  filter(Human_ID %in% sig_walk_hs_ids, Mouse_ID %in% sig_walk_mm_ids,
         Mouse_Ortho_Type == "ortholog_one2one") %>%
  dplyr::select(-contains("Rat")) %>%
  unique() %>%
  mutate(Gene_Group = "Walker")
pow_new_mouse_hom<- new_hom %>%
  filter(Human_ID %in% sig_pow_hs_ids, Rat_ID %in% sig_pow_rn_ids,
         Mouse_Ortho_Type == "ortholog_one2one") %>%
  dplyr::select(-contains("Rat")) %>%
  unique() %>%
  mutate(Gene_Group = "Powell")

# Rat
carp_new_rat_hom <- new_hom %>%
  filter(Human_ID %in% sig_carp_hs_ids, Mouse_ID %in% sig_carp_mm_ids,
         Rat_Ortho_Type == "ortholog_one2one") %>%
  dplyr::select(-contains("Mouse")) %>%
  unique() %>%
  mutate(Gene_Group = "Carpenter")
walk_new_rat_hom <- new_hom %>%
  filter(Human_ID %in% sig_walk_hs_ids, Mouse_ID %in% sig_walk_mm_ids,
         Rat_Ortho_Type == "ortholog_one2one") %>%
  dplyr::select(-contains("Mouse")) %>%
  unique() %>%
  mutate(Gene_Group = "Walker")
pow_new_rat_hom <- new_hom %>%
  filter(Human_ID %in% sig_pow_hs_ids, Rat_ID %in% sig_pow_rn_ids,
         Rat_Ortho_Type == "ortholog_one2one") %>%
  dplyr::select(-contains("Mouse")) %>%
  unique() %>%
  mutate(Gene_Group = "Powell")

# Join together; all genes not DEGs in any dataset are "other" homologs  
full_new_mouse_hom <- Reduce(function(x, y) 
  full_join(x, y), list(carp_new_mouse_hom, walk_new_mouse_hom, pow_new_mouse_hom, new_mouse_hom)) %>%
  mutate(Gene_Group = fct_relevel(ifelse(is.na(Gene_Group), 
                                         "All Other\nHomologs", Gene_Group),
                                  levels = c("Carpenter", "Walker", "Powell", "All Other\nHomologs"))) %>%
  filter(Mouse_Ortho_Type == "ortholog_one2one") %>%
  unique()

full_new_rat_hom <- Reduce(function(x, y) 
  full_join(x, y), list(carp_new_rat_hom, walk_new_rat_hom, pow_new_rat_hom, new_rat_hom)) %>%
  mutate(Gene_Group = fct_relevel(ifelse(is.na(Gene_Group), 
                                         "All Other\nHomologs", Gene_Group),
                                  levels = c("Carpenter", "Walker", "Powell", "All Other Homologs"))) %>%
  filter(Rat_Ortho_Type == "ortholog_one2one") %>%
  unique()
  

### STATISTICS ----------
# Stats comparing each dataset to the ones in no datasets
stat_mouse_perc <- full_new_mouse_hom %>%
  filter(!is.na(MouseToHuman_PercMatch)) %>%
  wilcox_test(MouseToHuman_PercMatch ~ Gene_Group, 
              p.adjust.method = "bonferroni") %>%
  filter(group2 == "All Other Homologs")

stat_rat_perc <- full_new_rat_hom %>%
  filter(!is.na(RatToHuman_PercMatch)) %>%
  wilcox_test(RatToHuman_PercMatch ~ Gene_Group,
              p.adjust.method = "bonferroni") %>%
  filter(group2 == "All Other Homologs")

# Medians
full_new_mouse_hom %>%
  filter(!is.na(MouseToHuman_PercMatch)) %>%
  group_by(Gene_Group) %>%
  summarize(Median = median(MouseToHuman_PercMatch))

full_new_rat_hom %>%
  filter(!is.na(RatToHuman_PercMatch)) %>%
  group_by(Gene_Group) %>%
  summarize(Median = median(RatToHuman_PercMatch))


### PLOTS ----------
# Label dplyr::select outliers
shared_mouse_perc_labels <- c("BCAS1")
general_mouse_perc_labels <- c("MUC3A", "XAF1", "ZDBF2",
                               "BST2", "PPP1R15A",
                               "CABYR", "SPINT3")

mouse_perc_plot <- full_new_mouse_hom %>%
  mutate(Gene_Group = as.character(Gene_Group),
         Gene_Group = fct_relevel(case_when(Gene_Group == "All Other Homologs" ~ "Other",
                                            TRUE ~ Gene_Group),
                                  levels = c("Carpenter", "Walker", "Powell", "Other"))) %>%
  ggplot(aes(x = Gene_Group, y = MouseToHuman_PercMatch)) +
  geom_violin(aes(fill = Gene_Group), width = 0.8) +  
  geom_boxplot(fill = "white", color = "black", width = 0.12, 
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  labs(x = "", y = "Sequence Match to Human (%)", title = "Mouse-Human") +
  scale_y_continuous(limits = c(0, 102), expand = c(0, 0)) +
  scale_fill_manual(values = c("#44AA99", "#C148AD", "#DDCC77", "grey60")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey85"),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 13.5),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(BM_Human_Gene_Name %in% shared_mouse_perc_labels, 
                                           BM_Human_Gene_Name, NA_character_)),
                   aes(label = label), color = "black", fill = "white",
                   nudge_x = c(0, 0.2), nudge_y = c(0, 2),
                   box.padding = 0.5, fontface = "bold", size = 3, seed = 31) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(BM_Human_Gene_Name %in% general_mouse_perc_labels, 
                                           BM_Human_Gene_Name, NA_character_)),
                   aes(label = label, fill = Gene_Group), color = "black",
                   box.padding = 0.5, fontface = "bold", size = 3, seed = 45)

# Label dplyr::select outliers
general_rat_perc_labels <- c("BST2", "PPP1R15A",
                             "MUC3A",
                             "RYR3", "CABYR")

rat_perc_plot <- full_new_rat_hom %>%
  mutate(Gene_Group = as.character(Gene_Group),
         Gene_Group = fct_relevel(case_when(Gene_Group == "All Other Homologs" ~ "Other",
                                            TRUE ~ Gene_Group),
                                  levels = c("Carpenter", "Walker", "Powell", "Other"))) %>%
  ggplot(aes(x = Gene_Group, y = RatToHuman_PercMatch)) +
  geom_violin(aes(fill = Gene_Group), width = 0.8) +  
  geom_boxplot(fill = "white", color = "black", width = 0.12, 
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  labs(x = "", y = "Sequence Match to Human (%)", title = "Rat-Human") +
  scale_y_continuous(limits = c(0, 102), expand = c(0, 0)) +
  scale_fill_manual(values = c("#44AA99", "#C148AD", "#DDCC77", "grey60")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey85"),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 13.5, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(BM_Human_Gene_Name %in% general_rat_perc_labels, 
                                           BM_Human_Gene_Name, NA_character_)),
                   aes(label = label, fill = Gene_Group), color = "black",
                   box.padding = 0.35, fontface = "bold", size = 3, seed = 495)

# Save plots
# svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Cons_Perc_5_31_2022.svg", width = 9, height = 4)
ggarrange(mouse_perc_plot, rat_perc_plot)
# dev.off()


## DN/DS VALUES ----------
old_hom <- getBM(attributes = c("external_gene_name", # Gene symbol
                                "ensembl_gene_id", # Human ENSEMBL
                                "mmusculus_homolog_ensembl_gene", # Mouse ID
                                "rnorvegicus_homolog_ensembl_gene", # Rat ID
                                # Mouse Homology
                                "mmusculus_homolog_dn",
                                "mmusculus_homolog_ds",
                                "mmusculus_homolog_subtype",
                                "mmusculus_homolog_orthology_type",
                                "mmusculus_homolog_perc_id",
                                "mmusculus_homolog_perc_id_r1",
                                "mmusculus_homolog_orthology_confidence",
                                # Rat Homology
                                "rnorvegicus_homolog_dn",
                                "rnorvegicus_homolog_ds",
                                "rnorvegicus_homolog_subtype",
                                "rnorvegicus_homolog_orthology_type",
                                "rnorvegicus_homolog_perc_id",
                                "rnorvegicus_homolog_perc_id_r1",
                                "rnorvegicus_homolog_orthology_confidence"),
                 mart = hs_ens99) %>%
  dplyr::rename(BM_Human_Gene_Name = external_gene_name,
                Human_ID = ensembl_gene_id,
                Mouse_ID = mmusculus_homolog_ensembl_gene,
                Rat_ID = rnorvegicus_homolog_ensembl_gene,
                Mouse_DN = mmusculus_homolog_dn,
                Mouse_DS = mmusculus_homolog_ds,
                Mouse_Hom_Subtype = mmusculus_homolog_subtype,
                Mouse_Ortho_Type = mmusculus_homolog_orthology_type,
                MouseToHuman_PercMatch = mmusculus_homolog_perc_id,
                HumanToMouse_PercMatch = mmusculus_homolog_perc_id_r1,
                Mouse_Ortho_Conf = mmusculus_homolog_orthology_confidence,
                Rat_DN = rnorvegicus_homolog_dn,
                Rat_DS = rnorvegicus_homolog_ds,
                Rat_Hom_Subtype = rnorvegicus_homolog_subtype,
                Rat_Ortho_Type = rnorvegicus_homolog_orthology_type,
                RatToHuman_PercMatch = rnorvegicus_homolog_perc_id,
                HumanToRat_PercMatch = rnorvegicus_homolog_perc_id_r1,
                Rat_Ortho_Conf= rnorvegicus_homolog_orthology_confidence) %>%
  mutate(Mouse_DNDS = Mouse_DN/Mouse_DS,
         Rat_DNDS = Rat_DN/Rat_DS) %>%
  mutate(across(everything(), function(x) na_if(x, "NaN")))

# Separate by mouse and rat
old_mouse_hom <- old_hom %>%
  dplyr::select(-contains("Rat")) %>%
  unique()

old_rat_hom <- old_hom %>%
  dplyr::select(-contains("Mouse")) %>%
  unique()

# Get general datasets for supplementary material
carp_old_hom <- old_hom %>%
  filter(Human_ID %in% sig_carp_hs_ids, Mouse_ID %in% sig_carp_mm_ids) %>%
  mutate(Gene_Group = "Carpenter")
walk_old_hom <- old_hom %>%
  filter(Human_ID %in% sig_walk_hs_ids, Mouse_ID %in% sig_walk_mm_ids) %>%
  mutate(Gene_Group = "Walker")
pow_old_hom <- old_hom %>%
  filter(Human_ID %in% sig_pow_hs_ids, Rat_ID %in% sig_pow_rn_ids) %>%
  mutate(Gene_Group = "Powell")

# Join for supplementary material
full_old_hom <- Reduce(function(x, y) 
  full_join(x, y), list(carp_old_hom, walk_old_hom, pow_old_hom, old_hom)) %>%
  mutate(Gene_Group = fct_relevel(ifelse(is.na(Gene_Group), 
                                         "All Other Homologs", Gene_Group),
                                  levels = c("Carpenter", "Walker", "Powell", "All Other Homologs"))) %>%
  relocate(Gene_Group, .before = everything()) %>%
  unique() %>%
  dplyr::select(-Mouse_Hom_Subtype, -Rat_Hom_Subtype)

# Get dataframes for analysis
# Mouse
carp_old_mouse_hom <- old_hom %>%
  filter(Human_ID %in% sig_carp_hs_ids, Mouse_ID %in% sig_carp_mm_ids,
         Mouse_Ortho_Type == "ortholog_one2one") %>%
  dplyr::select(-contains("Rat")) %>%
  unique() %>%
  mutate(Gene_Group = "Carpenter")
walk_old_mouse_hom <- old_hom %>%
  filter(Human_ID %in% sig_walk_hs_ids, Mouse_ID %in% sig_walk_mm_ids,
         Mouse_Ortho_Type == "ortholog_one2one") %>%
  dplyr::select(-contains("Rat")) %>%
  unique() %>%
  mutate(Gene_Group = "Walker")
pow_old_mouse_hom <- old_hom %>%
  filter(Human_ID %in% sig_pow_hs_ids, Rat_ID %in% sig_pow_rn_ids,
         Mouse_Ortho_Type == "ortholog_one2one") %>%
  dplyr::select(-contains("Rat")) %>%
  unique() %>%
  mutate(Gene_Group = "Powell")

# Rat
carp_old_rat_hom <- old_hom %>%
  filter(Human_ID %in% sig_carp_hs_ids, Mouse_ID %in% sig_carp_mm_ids,
         Rat_Ortho_Type == "ortholog_one2one") %>%
  dplyr::select(-contains("Mouse")) %>%
  unique() %>%
  mutate(Gene_Group = "Carpenter")
walk_old_rat_hom <- old_hom %>%
  filter(Human_ID %in% sig_walk_hs_ids, Mouse_ID %in% sig_walk_mm_ids,
         Rat_Ortho_Type == "ortholog_one2one") %>%
  dplyr::select(-contains("Mouse")) %>%
  unique() %>%
  mutate(Gene_Group = "Walker")
pow_old_rat_hom <- old_hom %>%
  filter(Human_ID %in% sig_pow_hs_ids, Rat_ID %in% sig_pow_rn_ids,
         Rat_Ortho_Type == "ortholog_one2one") %>%
  dplyr::select(-contains("Mouse")) %>%
  unique() %>%
  mutate(Gene_Group = "Powell")

# Join together; all genes not DEGs in any dataset are "other" homologs  
full_old_mouse_hom <- Reduce(function(x, y) 
  full_join(x, y), list(carp_old_mouse_hom, walk_old_mouse_hom, pow_old_mouse_hom, old_mouse_hom)) %>%
  mutate(Gene_Group = fct_relevel(ifelse(is.na(Gene_Group), 
                                         "All Other\nHomologs", Gene_Group),
                                  levels = c("Carpenter", "Walker", "Powell", "All Other\nHomologs"))) %>%
  filter(Mouse_Ortho_Type == "ortholog_one2one") %>%
  unique()

full_old_rat_hom <- Reduce(function(x, y) 
  full_join(x, y), list(carp_old_rat_hom, walk_old_rat_hom, pow_old_rat_hom, old_rat_hom)) %>%
  mutate(Gene_Group = fct_relevel(ifelse(is.na(Gene_Group), 
                                         "All Other\nHomologs", Gene_Group),
                                  levels = c("Carpenter", "Walker", "Powell", "All Other\nHomologs"))) %>%
  filter(Rat_Ortho_Type == "ortholog_one2one") %>%
  unique()

### STATISTICS ----------
# Stats comparing each dataset to the ones in no datasets
stat_mouse_dnds <- full_old_mouse_hom %>%
  filter(!is.na(Mouse_DNDS), Mouse_DNDS != Inf) %>%
  wilcox_test(Mouse_DNDS ~ Gene_Group, p.adjust.method = "bonferroni") %>%
  filter(group2 == "All Other\nHomologs")

# Medians
full_old_mouse_hom %>%
  filter(!is.na(Mouse_DNDS), Mouse_DNDS != Inf) %>%
  group_by(Gene_Group) %>%
  summarize(Median = median(Mouse_DNDS))

# NOTE: Rat analysis is at the end of the script, so that only those
# rat homologs that were consistent between Ensembl versions were included in
# the statistical analysis.


### PLOTS ----------
# Label select outliers
# Genes to label
shared_mouse_dnds_label_list <- c("BCAS1")
general_mouse_dnds_label_list <- c("XAF1", "ARMCX4", "ZDBF2",
                                   "PPP1R15A", "CFAP74",
                                   "CABYR", "GPRIN3", "SPINT3")

mouse_dnds_plot <- full_old_mouse_hom %>%
  mutate(Gene_Group = as.character(Gene_Group),
         Gene_Group = fct_relevel(case_when(Gene_Group == "All Other\nHomologs" ~ "Other",
                                            TRUE ~ Gene_Group),
                                  levels = c("Carpenter", "Walker", "Powell", "Other"))) %>%
  filter(!is.na(Mouse_DNDS), Mouse_DNDS != Inf) %>%
  ggplot(aes(x = Gene_Group, y = Mouse_DNDS)) +
  geom_violin(aes(fill = Gene_Group), width = 0.8) +  
  geom_boxplot(fill = "white", color = "black", width = 0.18, 
               outlier.size = 2, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.6) +
  labs(x = "", y = "dN/dS with Human Homolog",
       title = "Mouse-Human") +
  scale_x_discrete(expand = c(0, 0.55)) +
  scale_y_continuous(limits = c(0, 2.45), breaks = seq(0, 2.4, 0.4), expand = c(0, 0)) +
  scale_fill_manual(values = c("#44AA99", "#C148AD", "#DDCC77", "grey60")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey85"),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 13.5),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(BM_Human_Gene_Name %in% shared_mouse_dnds_label_list, 
                                           BM_Human_Gene_Name, NA_character_)),
                   aes(label = label), color = "black", fill = "white",
                   nudge_x = c(0.45, 0.25), nudge_y = c(-0.2, 0.2), 
                   fontface = "bold", size = 3, seed = 41) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(BM_Human_Gene_Name %in% general_mouse_dnds_label_list, 
                                           BM_Human_Gene_Name, NA_character_)),
                   aes(label = label, fill = Gene_Group), color = "black",
                   box.padding = 0.5, fontface = "bold", size = 3, seed = 43)
# NOTE THAT 1 "ALL OTHER" OUTLIER HAS BEEN REMOVED FROM PLOT FOR CLARITY
# Outlier: CDR1 (98.98)

shared_rat_dnds_label_list <- c("BCAS1")
general_rat_dnds_label_list <- c("XAF1", "ARMCX4", "ZDBF2",
                                 "PPP1R15A",
                                 "CABYR", "GPRIN3")

rat_dnds_plot <- full_old_rat_hom %>%
  mutate(Gene_Group = as.character(Gene_Group),
         Gene_Group = fct_relevel(case_when(Gene_Group == "All Other\nHomologs" ~ "Other",
                                            TRUE ~ Gene_Group),
                                  levels = c("Carpenter", "Walker", "Powell", "Other"))) %>%
  ggplot(aes(x = Gene_Group, y = Rat_DNDS)) +
  geom_violin(aes(fill = Gene_Group), width = 0.75) +  
  geom_boxplot(fill = "white", color = "black", width = 0.18, 
               outlier.size = 2, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  labs(x = "", y = "dN/dS with Human Homolog",
       title = "Rat-Human") +
  scale_y_continuous(limits = c(0, 1.63), breaks = seq(0, 1.6, 0.4), expand = c(0, 0)) +
  scale_fill_manual(values = c("#44AA99", "#C148AD", "#DDCC77", "grey60")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey85"),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 13.5, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_label_repel(data = . %>%
                   mutate(label = ifelse(BM_Human_Gene_Name %in% shared_rat_dnds_label_list, 
                                           BM_Human_Gene_Name, NA_character_)),
                 aes(label = label), color = "black", fill = "white",
                 box.padding = 0.5, fontface = "bold", size = 3, seed = 41) +
  geom_label_repel(data = . %>%
                     mutate(label = ifelse(BM_Human_Gene_Name %in% general_rat_dnds_label_list, 
                                           BM_Human_Gene_Name, NA_character_)),
                   aes(label = label, fill = Gene_Group), color = "black",
                   box.padding = 0.4, fontface = "bold", size = 3, seed = 80)
# NOTE THAT 2 "ALL OTHER" OUTLIERS HAVE BEEN REMOVED FROM PLOT FOR CLARITY
# Outliers: WDR83OS (6.01); LPA (3.49)

# Save plots
svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Cons_DNDS_5_31_2022.svg", width = 9, height = 4)
ggarrange(mouse_dnds_plot, rat_dnds_plot)
dev.off()


# COMPARE TO DEVELOPMENTAL EXPRESSION ACROSS SPECIES -----------
# "Developmental gene expression differences between humans & mammalian models"
# Cardoso-Moriera et al. (2020)
# Table S3A. Comparison of brain temporal trajectories

# Includes human, mouse, rat, rabbit, rhesus
cons_orig <- read_excel(paste0(cm_dir, "1-s2.0-S2211124720312973-mmc4.xlsx"), 
                        sheet = 1, skip = 2)

# Tidy and create new columns
# Note that with the new column "Conservation" we DO NOT have enough information
# to determine whether a gene's temporal developmental expression is shared
# between mouse and rat.
cons2020 <- cons_orig %>%
  dplyr::select(-c(starts_with("Prob"), Rabbit, Rhesus)) %>%
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
  dplyr::select(-c(starts_with("Prob"), Rabbit, Rhesus)) %>%
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

# Get general datasets for supplementary material
carp_cons_incomp <- cons2020_with_incompletes %>%
  filter(Human_ID %in% sig_carp_hs_ids) %>%
  mutate(Gene_Group = "Carpenter")
walk_cons_incomp  <- cons2020_with_incompletes %>%
  filter(Human_ID %in% sig_walk_hs_ids) %>%
  mutate(Gene_Group = "Walker")
pow_cons_incomp  <- cons2020_with_incompletes %>%
  filter(Human_ID %in% sig_pow_hs_ids) %>%
  mutate(Gene_Group = "Powell")

# Join datasets for supplementary material later
full_cons2020_with_incompletes <- Reduce(function(x, y) 
  full_join(x, y), list(carp_cons_incomp, walk_cons_incomp, pow_cons_incomp, 
                        cons2020_with_incompletes)) %>%
  mutate(Gene_Group = fct_relevel(ifelse(is.na(Gene_Group), 
                                         "All Other Homologs", Gene_Group),
                                  levels = c("Carpenter", "Walker", "Powell", 
                                             "All Other Homologs"))) %>%
  unique()

# Get general datasets again for analysis
carp_cons <- cons2020 %>%
  filter(Human_ID %in% sig_carp_hs_ids) %>%
  mutate(Gene_Group = "Carpenter")
walk_cons  <- cons2020 %>%
  filter(Human_ID %in% sig_walk_hs_ids) %>%
  mutate(Gene_Group = "Walker")
pow_cons  <- cons2020 %>%
  filter(Human_ID %in% sig_pow_hs_ids) %>%
  mutate(Gene_Group = "Powell")

# Join datasets for analysis with only genes with complete data
full_cons2020 <- Reduce(function(x, y) 
  full_join(x, y), list(carp_cons, walk_cons, pow_cons, cons2020)) %>%
  mutate(Gene_Group = fct_relevel(ifelse(is.na(Gene_Group), 
                                         "All Other Homologs", Gene_Group),
                                  levels = c("Carpenter", "Walker", "Powell",
                                             "All Other Homologs"))) %>%
  unique()

# How many genes are in each Gene_Group in the cons2020 dataset?
full_cons2020 %>%
  group_by(Gene_Group) %>%
  count()
# C: 431 (of 585 human genes), W: 65 (of 116 human genes)
# P: 391 (of 537 human genes); All Others: 9191

# How many have Conservation scores for both H-M and H-R?
full_cons2020 %>%
  group_by(Gene_Group) %>%
  filter(!is.na(Conservation)) %>%
  count()
# C: 330 (of 585 human genes), W: 47 (of 116 human genes)
# P: 323 (of 537 human genes); All Others: 9191

### PLOT ----------
# Plot counts on flexible y scales so that they are more comparable
cm_plot <- full_cons2020 %>%
  filter(!is.na(Conservation)) %>%
  mutate(ymax = case_when(Gene_Group == "Carpenter" ~ 306,
                          Gene_Group == "Walker" ~ 40.8,
                          Gene_Group == "Powell" ~ 306,
                          Gene_Group == "All Other Homologs" ~ 5100,
                          TRUE ~ NA_real_)) %>%
  mutate(Gene_Group = as.character(Gene_Group),
         Gene_Group = fct_relevel(case_when(Gene_Group == "All Other Homologs" ~ "Other",
                                            TRUE ~ Gene_Group),
                                  levels = c("Carpenter", "Walker", "Powell", "Other"))) %>%
  ggplot(aes(x = Conservation, fill = Gene_Group)) +
  geom_bar(aes(group = 1), show.legend = FALSE,
           position = position_dodge2(preserve = "single"),
           color = "black") +
  scale_fill_manual(values = c("#44AA99", "#C148AD", "#DDCC77", "grey60")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Count", 
       title = "Developmental Trajectories in the Forebrain\n") +
  theme(axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 13.5),
        panel.border = element_blank(),
        strip.text.x = element_text(size = 14),
        plot.title = element_text(face = "bold"),
        axis.line = element_line(color = "black")) +
  facet_wrap(~ Gene_Group, scales = "free_y", nrow = 1) +
  geom_blank(aes(y = ymax))
cm_plot

cm_gtable <- ggplot_gtable(ggplot_build(cm_plot))
striprt <- which(grepl('strip-r', cm_gtable$layout$name) | grepl('strip-t', cm_gtable$layout$name))
fills <- rep(c("#44AA99", "#C148AD", "#DDCC77", "grey60"), 4)
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', cm_gtable$grobs[[i]]$grobs[[1]]$childrenOrder))
  cm_gtable$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid::grid.draw(cm_gtable)

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/CM_Dev_Conservation_5_31_2022.svg", width = 8.7, height = 3.4)
grid::grid.draw(cm_gtable)
dev.off()


### STATISTICS ----------
# Fisher's Exact Test (and Pearson's Chi-Squared Test)
# Separate by each dataset
carp_table_setup <- full_cons2020 %>%
  filter(Gene_Group == "Carpenter" | Gene_Group == "All Other Homologs") %>%
  droplevels()
carp_table <- table(carp_table_setup$Gene_Group, carp_table_setup$Conservation)
fisher.test(carp_table, simulate.p.value = FALSE) # p-value = 0.1327

walk_table_setup <- full_cons2020 %>%
  filter(Gene_Group == "Walker" | Gene_Group == "All Other Homologs") %>%
  droplevels()
walk_table <- table(walk_table_setup$Gene_Group, walk_table_setup$Conservation)
fisher.test(walk_table, simulate.p.value = FALSE) # p-value = 0.9191

pow_table_setup <- full_cons2020 %>%
  filter(Gene_Group == "Powell" | Gene_Group == "All Other Homologs") %>%
  droplevels()
pow_table <- table(pow_table_setup$Gene_Group, pow_table_setup$Conservation)
fisher.test(pow_table, simulate.p.value = FALSE) # p-value = 0.03519 *
# Post-hoc test
# Separate into fully conserved vs. not
pow_table_setup_posthoc <- full_cons2020 %>%
  filter(Gene_Group == "Powell" | Gene_Group == "All Other Homologs") %>%
  mutate(Conservation = ifelse(Conservation == "HMR", "HMR", "Other")) %>%
  droplevels()
pow_table_posthoc <- table(pow_table_setup_posthoc$Gene_Group, 
                           pow_table_setup_posthoc$Conservation)
fisher.test(pow_table_posthoc, simulate.p.value = FALSE) # p-value = 0.004173 *

# Names of genes that were not fully conserved despite full data
full_cons2020 %>%
  filter(Gene_Group == "Carpenter",
         !is.na(Conservation), Conservation != "HMR") %>%
  pull(Gene_symbol)

full_cons2020 %>%
  filter(Gene_Group == "Walker",
         !is.na(Conservation), Conservation != "HMR") %>%
  pull(Gene_symbol)

full_cons2020 %>%
  filter(Gene_Group == "Powell",
         !is.na(Conservation), Conservation != "HMR") %>%
  pull(Gene_symbol)


## FINAL CONSERVATION DATAFRAME ----------
full_all_cons_df <- full_old_hom %>%
  dplyr::select(-contains(c("Perc", "Ortho_Conf")), -BM_Human_Gene_Name) %>%
  rename(Mouse_Ortho_Type = "Ensembl_99_Mouse_Ortho_Type",
         Rat_Ortho_Type = "Ensembl_99_Rat_Ortho_Type",
         Rat_ID = "Ensembl_99_Rat_ID") %>%
  mutate(across(everything(), function(x) na_if(x, ""))) %>%
  full_join(full_new_hom) %>%
  rename(Mouse_Ortho_Type = "Ensembl_106_Mouse_Ortho_Type",
         Rat_Ortho_Type = "Ensembl_106_Rat_Ortho_Type",
         Rat_ID = "Ensembl_106_Rat_ID") %>%
  dplyr::select(-contains("Ortho_Conf")) %>%
  dplyr::select(-c(Mouse_DN, Mouse_DS, Rat_DN, Rat_DS)) %>%
  relocate(Ensembl_99_Rat_Ortho_Type, .after = Mouse_DNDS) %>%
  relocate(BM_Human_Gene_Name, .after = Gene_Group) %>%
  mutate(Rat_Same_ID = ifelse(Ensembl_99_Rat_ID == Ensembl_106_Rat_ID, 
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
  mutate(Any_Rat_Same = ifelse(sum(Rat_Same_ID) >= 1, TRUE, FALSE),
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
  left_join(full_cons2020_with_incompletes) %>%
  droplevels() %>%
  dplyr::select(1:4, 6:7, 5, 8:9, 23, 11:13, 10, 14:16, 27)

# Subset to only candidates
full_cons_df <- full_all_cons_df %>% filter(Gene_Group != "All Other Homologs")
dim(full_cons_df) # 1424, 18

# Save as RDS object
saveRDS(full_cons_df, file = paste0(main_dir, "post_processing/results/other/full_cons_df.RDS"))


## STATISTICAL ANALYSIS FOR RAT DN/DS ----------
# NOTE: Only those rat homologs that were consistent between Ensembl versions 
# were included in the statistical analysis.

# All homologs
stat_rat_dnds <- full_old_rat_hom %>%
  filter(!is.na(Rat_DNDS), Rat_DNDS != Inf) %>%
  wilcox_test(Rat_DNDS ~ Gene_Group, p.adjust.method = "bonferroni") %>%
  filter(group2 == "All Other\nHomologs")

# Medians
full_old_rat_hom %>%
  filter(!is.na(Rat_DNDS), Rat_DNDS != Inf) %>%
  group_by(Gene_Group) %>%
  summarize(Median = median(Rat_DNDS))

# Only homologs consistent between versions
consistent_hs_ids <- full_all_cons_df %>% 
  filter(Ensembl_99_Rat_ID == "Same ID as Ensembl release 106") %>%
  pull(Human_ID) %>%
  unique()
stat_rat_dnds_consistent <- full_old_rat_hom %>%
  filter(!is.na(Rat_DNDS), Rat_DNDS != Inf, Human_ID %in% consistent_hs_ids) %>%
  wilcox_test(Rat_DNDS ~ Gene_Group, p.adjust.method = "bonferroni") %>%
  filter(group2 == "All Other\nHomologs")

# Medians
full_old_rat_hom %>%
  filter(!is.na(Rat_DNDS), Rat_DNDS != Inf, Human_ID %in% consistent_hs_ids) %>%
  group_by(Gene_Group) %>%
  summarize(Median = median(Rat_DNDS))

