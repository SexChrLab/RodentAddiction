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
main_dir <- "C:/Annika/GitHub Repositories/RodentAddiction/"
de_dir <- paste0(main_dir, "results/gene_info/")
cm_dir <- paste0(main_dir, "downloaded_data/cm_2020/")

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

# IDs to keep
keep_mm_ids <- sd_crave_mm
keep_rn_ids <- sd_crave_rn
keep_hs_ids <- sd_crave_hs



## GET NEWEST MARTS ----------
# Mouse (GRCm39)
grcm39 <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl",
                     host = "https://dec2021.archive.ensembl.org")
# Rat (mRatBN7.2)
mRatBN7.2 <- useEnsembl("ensembl", dataset = "rnorvegicus_gene_ensembl", 
                        host = "https://dec2021.archive.ensembl.org")
# Human (GRCh38.p13)
grch38 <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", 
                     host = "https://dec2021.archive.ensembl.org")



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
                Mouse_Ortho_Conf= mmusculus_homolog_orthology_confidence,
                Rat_Hom_Subtype = rnorvegicus_homolog_subtype,
                Rat_Ortho_Type = rnorvegicus_homolog_orthology_type,
                RatToHuman_PercMatch = rnorvegicus_homolog_perc_id,
                HumanToRat_PercMatch = rnorvegicus_homolog_perc_id_r1,
                Rat_Ortho_Conf = rnorvegicus_homolog_orthology_confidence) %>%
  mutate(across(everything(), function(x) na_if(x, "NaN"))) %>%
  mutate(across(everything(), function(x) na_if(x, "")))




test <- getBM(attributes = c("external_gene_name", # Gene symbol
                                "ensembl_gene_id", # Human ENSEMBL
                                "rnorvegicus_homolog_ensembl_gene", # Rat ID
                                # Rat Homology
                                "rnorvegicus_homolog_subtype",
                                "rnorvegicus_homolog_orthology_type",
                                "rnorvegicus_homolog_perc_id",
                                "rnorvegicus_homolog_perc_id_r1",
                                "rnorvegicus_homolog_orthology_confidence"),
                 filter = "ensembl_gene_id",
                 values = c("ENSG00000167996"),
                 mart = grch38)

test <- getBM(attributes = c("external_gene_name", # Gene symbol
                             "ensembl_gene_id", # Human ENSEMBL
                             "rnorvegicus_homolog_ensembl_gene", # Rat ID
                             # Rat Homology
                             "rnorvegicus_homolog_subtype",
                             "rnorvegicus_homolog_orthology_type",
                             "rnorvegicus_homolog_perc_id",
                             "rnorvegicus_homolog_perc_id_r1",
                             "rnorvegicus_homolog_orthology_confidence"),
              filter = "ensembl_gene_id",
              values = c("ENSMUSG00000024661"),
              mart = grcm39)

ENSRNOG00000022619

# Separate by mouse and rat
new_mouse_hom <- new_hom %>%
  mutate(Gene_Group = fct_relevel(
    case_when(Human_ID %in% keep_hs_ids & Mouse_ID %in% keep_mm_ids ~ "Craving",
              TRUE ~ "All Other Homologs"),
    c("Craving", "All Other Homologs"))) %>%
  filter(!is.na(MouseToHuman_PercMatch)) %>%
  dplyr::select(-contains("Rat")) %>%
  unique()
  
new_rat_hom <- new_hom %>%
  mutate(Gene_Group = fct_relevel(
    case_when(Human_ID %in% keep_hs_ids & Rat_ID %in% keep_rn_ids ~ "Craving",
              TRUE ~ "All Other Homologs"),
    c("Craving", "All Other Homologs"))) %>%
  filter(!is.na(RatToHuman_PercMatch)) %>%
  dplyr::select(-contains("Mouse")) %>%
  unique()

# Are any homologous pairs NOT one-to-one?
new_mouse_hom %>%
  filter(Gene_Group == "Craving", Mouse_Ortho_Type != "ortholog_one2one")
# NO for mouse
new_rat_hom %>%
  filter(Gene_Group == "Craving", Rat_Ortho_Type != "ortholog_one2one")
# YES for rat; HSPA8, CELF6 - this will be excluded

# Get stats for plots
stat_mouse_perc <- new_mouse_hom %>% 
  filter(Mouse_Ortho_Type == "ortholog_one2one") %>%
  wilcox_test(MouseToHuman_PercMatch ~ Gene_Group) %>%
  add_xy_position(x = "Gene_Group") %>%
  mutate(xmax = 1.8)

stat_rat_perc <- new_rat_hom %>% 
  filter(Rat_Ortho_Type == "ortholog_one2one") %>%
  wilcox_test(RatToHuman_PercMatch ~ Gene_Group) %>%
  add_xy_position(x = "Gene_Group") %>%
  mutate(xmax = 1.8)

# Create plots
perc_labels_mouse_crave <- c("BCAS1")
mouse_perc_plot <- new_mouse_hom %>%
  ggplot(aes(x = Gene_Group, y = MouseToHuman_PercMatch)) +
  geom_violin(aes(fill = Gene_Group), width = 0.8) +  
  geom_boxplot(fill = "white", color = "black", width = 0.12, 
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  labs(x = "", y = "Sequence Match to Human(%)", title = "Human-Mouse\n") +
  scale_x_discrete(expand = c(0, 0.55)) +
  scale_fill_manual(values = c("#DE8971", "grey80")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey85"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14)) +
  geom_label_repel(
    data = . %>%
      mutate(label = ifelse(BM_Human_Gene_Name %in% perc_labels_mouse_crave &
                              Gene_Group == "Craving",
                            BM_Human_Gene_Name, NA_character_)),
                   aes(label = label), color = "black", fill = "#DE8971",
                   box.padding = 0.5, fontface = "bold", size = 3) +
  stat_pvalue_manual(stat_mouse_perc, label = "p", y.position = 43, 
                     tip.length = c(-0.13, 0), vjust = 3)

perc_labels_rat_crave <- c("B2M", "BCAS1")
rat_perc_plot <- new_rat_hom %>%
  ggplot(aes(x = Gene_Group, y = RatToHuman_PercMatch)) +
  geom_violin(aes(fill = Gene_Group), width = 0.75) +  
  geom_boxplot(fill = "white", color = "black", width = 0.12, 
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  labs(x = "", y = "Sequence Match to Human(%)", title = "Human-Rat\n") +
  scale_x_discrete(expand = c(0, 0.55)) +
  scale_fill_manual(values = c("#DE8971", "grey80")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey85"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        axis.title.y = element_text(color = "white")) +
  geom_label_repel(
    data = . %>%
      mutate(label = ifelse(BM_Human_Gene_Name %in% perc_labels_rat_crave &
                              Gene_Group == "Craving",
                            BM_Human_Gene_Name, NA_character_)),
    aes(label = label, fill = Gene_Group), color = "black",
    box.padding = 0.5, fontface = "bold", size = 3) +
  stat_pvalue_manual(stat_rat_perc, label = "p", y.position = 40, 
                     tip.length = c(-0.08, 0), vjust = 1.4)

# Save plots
svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Cons_percentage_1_24_2022.svg", width = 9, height = 3.75)
ggarrange(mouse_perc_plot, rat_perc_plot)
dev.off()

# Summaries of means and medians for mouse and rat homologs
new_mouse_hom %>% group_by(Gene_Group) %>%
  filter(Mouse_Ortho_Type == "ortholog_one2one") %>%
  summarize(Mean_SeqSim = mean(MouseToHuman_PercMatch, na.rm = TRUE), 
            Median_SeqSim = median(MouseToHuman_PercMatch, na.rm = TRUE))
new_rat_hom %>% group_by(Gene_Group) %>%
  filter(Rat_Ortho_Type == "ortholog_one2one") %>%
  summarize(Mean_SeqSim = mean(RatToHuman_PercMatch, na.rm = TRUE), 
            Median_SeqSim = median(RatToHuman_PercMatch, na.rm = TRUE))


## DN/DS VALUES ----------
# Last available mart with DN and DS values: Ensembl Release 99
hs_ens99 <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                    host = "https://jan2020.archive.ensembl.org")

# Get info
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
                                "mmusculus_homolog_wga_coverage",
                                "mmusculus_homolog_orthology_confidence",
                                # Rat Homology
                                "rnorvegicus_homolog_dn",
                                "rnorvegicus_homolog_ds",
                                "rnorvegicus_homolog_subtype",
                                "rnorvegicus_homolog_orthology_type",
                                "rnorvegicus_homolog_perc_id",
                                "rnorvegicus_homolog_perc_id_r1",
                                "rnorvegicus_homolog_wga_coverage",
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
                Mouse_WGA_Cov = mmusculus_homolog_wga_coverage,
                Mouse_Ortho_Conf= mmusculus_homolog_orthology_confidence,
                Rat_DN = rnorvegicus_homolog_dn,
                Rat_DS = rnorvegicus_homolog_ds,
                Rat_Hom_Subtype = rnorvegicus_homolog_subtype,
                Rat_Ortho_Type = rnorvegicus_homolog_orthology_type,
                RatToHuman_PercMatch = rnorvegicus_homolog_perc_id,
                HumanToRat_PercMatch = rnorvegicus_homolog_perc_id_r1,
                Rat_WGA_Cov = rnorvegicus_homolog_wga_coverage,
                Rat_Ortho_Conf= rnorvegicus_homolog_orthology_confidence) %>%
  mutate(Mouse_DNDS = Mouse_DN/Mouse_DS,
         Rat_DNDS = Rat_DN/Rat_DS) %>%
  mutate(across(everything(), function(x) na_if(x, "NaN")))

# Separate by mouse and rat
old_mouse_hom <- old_hom %>%
  mutate(Gene_Group = fct_relevel(
    case_when(Human_ID %in% keep_hs_ids & Mouse_ID %in% keep_mm_ids ~ "Craving",
              TRUE ~ "All Other Homologs"),
    c("Craving", "All Other Homologs"))) %>%
  dplyr::select(-contains("Rat")) %>%
  unique()

old_rat_hom <- old_hom %>%
  mutate(Gene_Group = fct_relevel(
    case_when(Human_ID %in% keep_hs_ids & Rat_ID %in% keep_rn_ids ~ "Craving",
              TRUE ~ "All Other Homologs"),
    c("Craving", "All Other Homologs"))) %>%
  dplyr::select(-contains("Mouse")) %>%
  filter(BM_Human_Gene_Name != "FTH1") %>%
  unique()

# Are any homologous pairs NOT one-to-one?
old_mouse_hom %>%
  filter(Gene_Group == "Craving", Mouse_Ortho_Type != "ortholog_one2one")
# YES for mouse; BTG1, CELF6 - these will be excluded 
old_rat_hom %>%
  filter(Gene_Group == "Craving", Rat_Ortho_Type != "ortholog_one2one")
# YES for rat; BTG1, CELF6, HSPA8 - these will be excluded

# Get stats for plots
stat_mouse_dnds <- old_mouse_hom %>% 
  filter(Mouse_Ortho_Type == "ortholog_one2one", !is.na(Mouse_DNDS)) %>%
  filter(Mouse_DNDS != Inf) %>%
  wilcox_test(Mouse_DNDS ~ Gene_Group) %>%
  add_xy_position(x = "Gene_Group")

stat_rat_dnds <- old_rat_hom %>% 
  filter(Rat_Ortho_Type == "ortholog_one2one", !is.na(Rat_DNDS)) %>%
  wilcox_test(Rat_DNDS ~ Gene_Group) %>%
  add_xy_position(x = "Gene_Group")

# Create plots
dnds_labels_mouse_crave <- c("BCAS1", "MBP", "FABP7")
mouse_dnds_plot <- old_mouse_hom %>%
  ggplot(aes(x = Gene_Group, y = Mouse_DNDS)) +
  geom_violin(aes(fill = Gene_Group), width = 0.8) +  
  geom_boxplot(fill = "white", color = "black", width = 0.18, 
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  labs(x = "", y = "dN/dS Compared to Human Homolog\n",
       title = "Human-Mouse\n") +
  scale_x_discrete(expand = c(0, 0.55)) +
  scale_y_continuous(limits = c(0, 3)) +
  scale_fill_manual(values = c("#DE8971", "grey80")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey85"),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 13)) +
  geom_label_repel(
    data = . %>%
      mutate(label = ifelse(BM_Human_Gene_Name %in% dnds_labels_mouse_crave &
                              Gene_Group == "Craving",
                            BM_Human_Gene_Name, NA_character_)),
  aes(label = label, fill = Gene_Group), color = "black",
  box.padding = 0.5, fontface = "bold", size = 3) +
  stat_pvalue_manual(stat_mouse_dnds, label = "p", y.position = 2.8,
                     tip.length = c(0.0175, 0.0005))

dnds_labels_rat_crave <- c("BCAS1", "MBP", "FABP7")
rat_dnds_plot <- old_rat_hom %>%
  ggplot(aes(x = Gene_Group, y = Rat_DNDS)) +
  geom_violin(aes(fill = Gene_Group), width = 0.75) +  
  geom_boxplot(fill = "white", color = "black", width = 0.18, 
               outlier.size = 3, outlier.shape = 21, outlier.color = "grey80", 
               outlier.fill = "black", outlier.alpha = 0.9) +
  labs(x = "", y = "dN/dS Compared to Human Homolog\n",
       title = "Human-Rat\n") +
  scale_x_discrete(expand = c(0, 0.55)) +
  scale_y_continuous(limits = c(0, 2.6), breaks = seq(0, 2.5, 0.5)) +
  scale_fill_manual(values = c("#DE8971", "grey80")) +
  theme(legend.position = "none", axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey85"),
        axis.text = element_text(size = 13),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 13, color = "white")) +
  geom_label_repel(
    data = . %>%
      mutate(label = ifelse(BM_Human_Gene_Name %in% dnds_labels_rat_crave &
                              Gene_Group == "Craving",
                            BM_Human_Gene_Name, NA_character_)),
    aes(label = label, fill = Gene_Group), color = "black",
    box.padding = 0.5, fontface = "bold", size = 3) +
  stat_pvalue_manual(stat_rat_dnds, label = "p", y.position = 2.55,
                     tip.length = c(0.02, 0.001))

# Save plots
svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Cons_DNDS_1_24_2022.svg", width = 9, height = 3.75)
ggarrange(mouse_dnds_plot, rat_dnds_plot)
dev.off()

# Summaries of means and medians for mouse and rat homologs
old_mouse_hom %>% group_by(Gene_Group) %>%
  filter(Mouse_Ortho_Type == "ortholog_one2one") %>%
  summarize(Mean_DNDS = mean(Mouse_DNDS, na.rm = TRUE), 
            Median_DNDS = median(Mouse_DNDS, na.rm = TRUE))
old_rat_hom %>% group_by(Gene_Group) %>%
  filter(Rat_Ortho_Type == "ortholog_one2one") %>%
  summarize(Mean_DNDS = mean(Rat_DNDS, na.rm = TRUE), 
            Median_DNDS = median(Rat_DNDS, na.rm = TRUE))




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
  mutate(Human_ID = factor(Human_ID),
         Human_Symbol = factor(Gene_symbol),
         Gene_Group = fct_relevel(
           case_when(Human_ID %in% sd_crave_hs ~ "Craving",
                     TRUE ~ "All Other Genes"),
           c("Craving", "All Other Genes")),
         .before = everything(),
         .keep = "unused") %>%
  # Categorize by conservation (see above note)
  mutate(Conservation = fct_relevel(case_when(
    (Mouse == "Same" & Rat == "Same") ~ "HMR",
    (Mouse == "Same" & Rat == "Different") ~ "HM",
    (Mouse == "Different" & Rat == "Same") ~ "HR",
    (Mouse == "Different" & Rat == "Different") ~ "H",
    # Assume NA = different
    # NA genes were not tested because they didn't have different temporal 
    # expression in a particular species.
    TRUE ~ NA_character_),
    c("HMR", "HM", "HR", "H"))) %>%
  mutate(across(everything(), function(x) na_if(x, "NA")))

# How many genes are in each Gene_Group in the cons2020 dataset?
cons2020 %>%
  group_by(Gene_Group) %>%
  count()
# Craving: 26/33; All Other Genes: 9824


# Which genes are conserved?
# Function
cons_type_genes <- function(cons_type) {
  cons2020 %>%
    filter(Gene_Group == "Craving", Conservation == cons_type) %>%
    pull(Human_Symbol) %>%
    droplevels() %>%
    sort()
}

# 18 Conserved in all 3 species
cons_type_genes("HMR")
# [1] AGK     AMZ1    B2M     BCAS1   BTG1    CARTPT  CCDC88C
# [8] EGR2    FABP7   GUCY1A3 HAPLN2  IRS2    KIF5A   MBP    
# [15] MOBP    PITPNM3 TTLL1   VIM   

# 2 Conserved in Human-Mouse
cons_type_genes("HM")
# N/A

# 1 Conserved in Human-Rat
cons_type_genes("HR")
# [1] NTS  RGS5

# 4 Not conserved in Humans-Rodents
cons_type_genes("H")
# [1] CACYBP

# Genes that aren't testable (not in cons2020 dataset)
sd_crave_hs[!sd_crave_hs %in% (cons2020 %>% pull(Human_ID))]
# [1] "ENSG00000167996" "ENSG00000150551" "ENSG00000004478"
# [4] "ENSG00000109189" "ENSG00000163659"
# Corresponds to: FTH1, LYPD1, FKBP4, USP46, TIPARP

# Plot counts on flexible y scales so that they are more comparable
cm_plot <- cons2020 %>%
  filter(!is.na(Conservation)) %>%
  mutate(ymax = case_when(Gene_Group == "Craving" ~ 10,
                          Gene_Group == "All Other Genes" ~ 8000,
                          TRUE ~ NA_real_)) %>%
  ggplot(aes(x = Conservation, fill = Gene_Group)) +
  geom_bar(aes(group = 1), show.legend = FALSE,
           position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#DE8971", "grey80")) +
  labs(y = "Count", 
       title = "Developmental Trajectories in the Forebrain\n") +
  theme(axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14)) +
  facet_grid(~ Gene_Group, scales = "free_y") +
  facet_wrap(~ Gene_Group, scales = "free_y") +
  geom_blank(aes(y = ymax))
# 655, 280

cm_plot

cm_gtable <- ggplot_gtable(ggplot_build(cm_plot))
striprt <- which(grepl('strip-r', cm_gtable$layout$name) | grepl('strip-t', cm_gtable$layout$name))
fills <- rep(c("#DE8971", "grey80"), 2)
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', cm_gtable$grobs[[i]]$grobs[[1]]$childrenOrder))
  cm_gtable$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid::grid.draw(cm_gtable)

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/CM_Dev_Conservation_1_24_2022.svg", width = 8.1, height = 3.8)
grid::grid.draw(cm_gtable)
dev.off()

# Fisher's Exact Test (and Pearson's Chi-Squared Test)
cons2020_table <- table(cons2020$Gene_Group, cons2020$Conservation) # Only test Craving vs. All Other Genes
fisher.test(cons2020_table, simulate.p.value = FALSE) # p-value = 0.553
# Craving genes DO NOT have different distributions for conservation of developmental trajectories

# Look into the temporal expression patterns for the genes that were not
# conserved in all 3 species.
cons_check_genes <- append(append(cons_type_genes("HM"), cons_type_genes("HR")), 
                           cons_type_genes("H"))

new_hom %>% filter(Human_ID %in% keep_hs_ids) %>% write_clip()

# Final dataframe
last_hom_df <- old_hom %>% 
  dplyr::select(!c(contains("Perc"), contains("WGA"))) %>%
  full_join(new_hom) %>%
  filter(Human_ID %in% keep_hs_ids) %>%
  dplyr::select(c(1:4, 15:17, 19)) %>%
  left_join(cons2020) %>%
  dplyr::select(-c(Human_Symbol, Mouse, Rat, Gene_Group)) %>%
  relocate(Conservation, .after = Human_ID) %>%
  relocate(Rat_ID, .before = Rat_DNDS) %>%
  relocate(MouseToHuman_PercMatch, .after = Mouse_ID) %>%
  relocate(Rat_DNDS, .after = everything()) %>%
  arrange(BM_Human_Gene_Name)

last_hom_df %>% write_clip

old_hom %>%
  filter(Human_ID %in% keep_hs_ids) %>%
  View()











