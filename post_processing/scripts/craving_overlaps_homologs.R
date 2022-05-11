## LOAD PACKAGES ----------
# Data manipulation
library(tidyverse)

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



## BEFORE USING THIS SCRIPT ----------
# User should alter these variables as necessary
main_dir <- "C:/Annika/GitHub Repositories/RodentAddiction/"
de_dir <- paste0(main_dir, "post_processing/results/gene_info/")

# Reduce chances of scientific notation
options(scipen = 9999)

# Update and set ggplot theme globally
new_theme <- theme_bw(base_size = 14) %+replace%
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))
theme_set(new_theme)



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



## LOAD IN DIFFERENTIAL EXPRESSION DATA ----------
# Load in 3 relevant dataframes together as a list of lists

# Get list of files
DE_files <- list.files(path = de_dir, pattern = "DE_FINAL", full.names = TRUE)
DE_files

# Verify they are in the correct order with next line of code before renaming
DE_names <- list.files(path = de_dir, pattern = "DE_FINAL", full.names = FALSE)
DE_names
DE_listnames <- c("Carpenter", "Powell", "Walker")

# Load in list of lists and rename
DE_lists <- lapply(DE_files, read.delim, sep = " ", header = T)
names(DE_lists) <- DE_listnames

# Check column names of each list; should all be the same
lapply(DE_lists, colnames)
# [1] "Gene_ID"    "Gene.Name"  "Chromosome"
# [4] "logFC"      "AveExpr"    "t"         
# [7] "P.Value"    "adj.P.Val"  "B"     

# Convert list of lists to a dataframe
full_df_nonames <- DE_lists %>%
  map2(., DE_listnames, ~ cbind(.x, Study = .y)) %>%
  Reduce(function(x, y) 
    full_join(x, y), .) %>%
  mutate(Study = fct_relevel(Study, 
                             levels = c("Carpenter", "Powell", "Walker"))) %>%
  dplyr::rename(Original_Gene_ID = Gene_ID)



## ADD GENE NAMES TO DATAFRAME ----------
# Get mouse and rat IDs
all_mm_ids_carp <- full_df_nonames %>% 
  filter(Study == "Carpenter") %>% 
  pull(Original_Gene_ID) %>% 
  unique()

all_mm_ids_walk <- full_df_nonames %>% 
  filter(Study == "Walker") %>% 
  pull(Original_Gene_ID) %>% 
  unique()

all_rn_ids <- full_df_nonames %>% 
  filter(Study == "Powell") %>% 
  pull(Original_Gene_ID) %>% 
  unique()

# Use IDs to get gene names
mart_mm_names_carp <- getBM(attributes = c("ensembl_gene_id", 
                                           "external_gene_name"),
                            filter = "ensembl_gene_id",
                            values = all_mm_ids_carp, mart = grcm39) %>%
  setNames(c("Original_Gene_ID", "Gene.Name"))

mart_mm_names_walk <- getBM(attributes = c("ensembl_gene_id", 
                                           "external_gene_name"),
                            filter = "ensembl_gene_id", 
                            values = all_mm_ids_walk, mart = grcm39) %>%
  setNames(c("Original_Gene_ID", "Gene.Name"))

mart_rn_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                       filter = "ensembl_gene_id", values = all_rn_ids,
                       mart = mRatBN7.2) %>%
  setNames(c("Original_Gene_ID", "Gene.Name"))

# Check that dimensions match number of genes
nrow(mart_mm_names_carp) # 10026
nrow(mart_mm_names_walk) # 9948
nrow(mart_rn_names) # 7104

# New dataframe with gene names, NA if missing
full_df <- Reduce(function(x, y) 
  full_join(x, y), list(full_df_nonames, mart_mm_names_carp, 
                          mart_mm_names_walk, mart_rn_names)) %>%
  relocate(Gene.Name, .after = Original_Gene_ID) %>%
  mutate(across(everything(), function(x) na_if(x, ""))) %>%
  filter(!is.na(Study))

# Get Powell genes
full_df %>%
  filter(Study == "Powell") %>%
  dplyr::select(-Gene.Name) %>%
  inner_join(mart_rn_names) %>%
  relocate(Gene.Name, .after = Original_Gene_ID) %>%
  write_clip()


## GET SIGNIFICANT GENES ----------
# Number of total genes per study after filtering
full_df %>%
  group_by(Study) %>%
  count()
# Study         n
# 1 Carpenter 10026
# 2 Powell    7104
# 3 Walker    9948


# Get summary stats about fold changes for 1) all genes, and 2) DEGs only
full_df %>%
  group_by(Study) %>%
  mutate(logFC = abs(logFC), FC = 2^(abs(logFC))) %>%
  summarize(min_logFC = min(logFC), mean_logFC = mean(logFC), 
            median_logFC = median(logFC), max_logFC = max(logFC), 
            min_FC = min(FC), mean_FC = mean(FC), median_FC = median(FC), 
            max_FC = max(FC))

# Significant with p-value only
sig_df_ponly <- full_df %>% filter(P.Value < 0.05)
sig_df_ponly %>%
  group_by(Study) %>%
  mutate(logFC = abs(logFC), FC = 2^(abs(logFC))) %>%
  summarize(min_logFC = min(logFC), mean_logFC = mean(logFC), 
            median_logFC = median(logFC), max_logFC = max(logFC), 
            min_FC = min(FC), mean_FC = mean(FC), median_FC = median(FC), 
            max_FC = max(FC))

# Sig with p-value and FC
sig_df <- full_df %>% filter(P.Value < 0.05, abs(logFC) > log2(1.2))
sig_df %>%
  group_by(Study) %>%
  mutate(logFC = abs(logFC), FC = 2^(abs(logFC))) %>%
  summarize(min_logFC = min(logFC), mean_logFC = mean(logFC), 
            median_logFC = median(logFC), max_logFC = max(logFC), 
            min_FC = min(FC), mean_FC = mean(FC), median_FC = median(FC), 
            max_FC = max(FC))

# Significant genes only
sig_df_fc <- full_df %>% filter(P.Value < 0.05, abs(logFC) > log2(1.2))
sig_df <- full_df %>% filter(P.Value < 0.05)


## VOLCANO PLOTS ----------
# How many DE genes for each study/comparison?
sig_df_ponly %>% group_by(Study) %>% count()
# Study         n
# 1 Carpenter   633
# 2 Powell      580
# 3 Walker      139

sig_df_fc %>% group_by(Study) %>% count()
# Study         n
# 1 Carpenter   136
# 2 Powell      580
# 3 Walker      109

# Number of DE genes for each study, as variables
carp_count <- sig_df %>% filter(Study == "Carpenter") %>% count() %>% pull(n)
pow_count <- sig_df %>% filter(Study == "Powell") %>% count() %>% pull(n)
walk_count <- sig_df %>% filter(Study == "Walker") %>% count() %>% pull(n)

# How many DE genes for each study/comparison that are up or down-regulated?
sig_df %>% mutate(Direction = ifelse(logFC > 0, "Up", "Down")) %>% group_by(Study, Direction) %>% count()
# Study     Direction     n
# 1 Carpenter Down        329
# 2 Carpenter Up          304
# 3 Powell    Down        291
# 4 Powell    Up          289
# 5 Walker    Down        100
# 6 Walker    Up           39

# Function for creating volcano plots
volcano_plot <- function(df, study, rect_fill, x_lim = NULL, x_breaks = waiver(),
                         y_lim = 5, x_color = "black", y_color = "black",  
                         x_axis_lab = waiver()) {
  
  # Add labels to dataframe for significant genes
  new_df <- df %>%
    mutate(DE = fct_relevel(case_when((logFC > log2(1.2) & P.Value < 0.05 ~ "Up"),
                                      (logFC < -log2(1.2) & P.Value < 0.05 ~ "Down"),
                                      TRUE ~ "None"),
                            levels = c("None", "Down", "Up")),
           DE2 = fct_relevel(case_when((logFC < -log2(1.2) & adj.P.Val < 0.05 ~ "Down"),
                                      TRUE ~ "None"),
                            levels = c("None", "Down")),
           Abs_logFC = abs(logFC),
           Study = case_when(str_detect(Study, "Powell") ~ "Powell",
                             str_detect(Study, "Walker") ~ "Walker",
                             str_detect(Study, "Carpenter") ~ "Carpenter")) %>%
    group_by(Study, DE) %>%
    arrange(-Abs_logFC, .by_group = TRUE) %>%
    mutate(DE_label = case_when(DE2 != "None" ~ Original_Gene_ID,
                                TRUE ~ NA_character_))
  
  # Volcano plot w/ labels for top 10 up- and down-regulated genes
  new_df %>%
    filter(Study == study) %>%
    ggplot(aes(x = logFC, y = -log10(P.Value), fill = DE, label = DE_label)) +
    geom_point(alpha = 0.35, color = "black", shape = 21, 
               size = 4, stroke = 1, show.legend = FALSE) +
    # geom_vline(xintercept = c(-log2(1.2), log2(1.2)), col = "grey10", 
    #            lty = "dashed", lwd = 1.25) + # log2FC threshold
    geom_hline(yintercept = -log10(0.05), col = "grey10", 
               lty = "dashed", lwd = 1.25) + # P-value threshold
    geom_label_repel(color = "black", size = 3.5, show.legend = FALSE,
                     min.segment.length = unit(0, "pt"), segment.size = 1,
                     box.padding = unit(5, "pt")) +
    scale_fill_manual(values = c("grey70", "#2EAECF", "#CA1A18")) +
    scale_x_continuous(limits = x_lim, breaks = x_breaks, labels = x_axis_lab) +
    scale_y_continuous(limits = c(0, y_lim), breaks = seq(0, 5, 1)) +
    labs(x = expression("Log"[2]* " Fold Change"), 
         y = expression("-Log"[10]* " P-Value")) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = c(rect_fill, "white")),
          strip.text.x = element_text(size = 14),
          axis.title.x = element_text(color = x_color),
          axis.title.y = element_text(color = y_color)) +
    facet_wrap(~ Study, nrow = 1)
}

vo1 <- volcano_plot(full_df, "Carpenter", "#44AA99", c(-1.5, 1.5), seq(-1.5, 1.5, 0.5), 4.2, "white", "black")
vo2 <- volcano_plot(full_df, "Walker", "#C148AD", c(-2, 2), seq(-2, 2, 1), 5, "black", "white")
vo3 <- volcano_plot(full_df, "Powell", "#DDCC77", c(-10, 4), seq(-10, 4, 2), 5.2)
# 390, 420

# svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Carpenter_volcano_1_24_2022.svg", width = 4.25, height = 4.25)
# vo1
# dev.off()
# svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Walker_volcano__1_24_2022.svg", width = 4.25, height = 4.25)
# vo2
# dev.off()
# svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Powell_volcano__1_24_2022.svg", width = 4.25, height = 4.25)
# vo3
# dev.off()


library(VennDiagram)
draw.pairwise.venn(area1 = 633, area2 = 139, cross.area = 15)
draw.pairwise.venn(area1 = 139, area2 = 548, cross.area = 5)
draw.pairwise.venn(area1 = 633, area2 = 548, cross.area = 45)




## GET SIGNIFICANT MOUSE GENES ----------
# Get ENSEMBL IDs for significant mouse genes
sig_mm_ids <- sig_df %>% 
  filter(Study == "Carpenter" | Study == "Walker") %>% 
  pull(Original_Gene_ID)
length(sig_mm_ids) # 772 (245 w/ LFC)

# Get shared gene IDs
shared_mm_ids <- sig_mm_ids[duplicated(sig_mm_ids)]
length(shared_mm_ids) # 15 (8 w/ LFC)

# Get dataframe with all relevant mouse info
sig_mm_bm <- getBM(attributes = c("external_gene_name", # Gene symbol
                                  "ensembl_gene_id", # Mouse ENSEMBL
                                  "rnorvegicus_homolog_ensembl_gene", # Rat
                                  "hsapiens_homolog_ensembl_gene"), # Human
                   filter = "ensembl_gene_id", 
                   values = sig_mm_ids,
                   mart = grcm39) %>%
  dplyr::rename(BM_Gene_Name = external_gene_name,
                Mouse_ID = ensembl_gene_id,
                Rat_ID = rnorvegicus_homolog_ensembl_gene,
                Human_ID = hsapiens_homolog_ensembl_gene) %>%
  mutate(across(everything(), function(x) na_if(x, "")),
         Has_Rat_Hom = case_when(!is.na(Rat_ID) ~ "Yes",
                                 TRUE ~ "No"),
         Has_Human_Hom = case_when(!is.na(Human_ID) ~ "Yes",
                                   TRUE ~ "No"))

# Join with other data for more complete information
all_mouse_sig_info <- sig_df %>%
  #  filter(Original_Gene_ID %in% shared_mm_id) %>%
  full_join(sig_mm_bm, by = c("Original_Gene_ID" = "Mouse_ID")) %>%
  dplyr::rename(ORIG_Gene_ID = Original_Gene_ID,
                ORIG_Gene_Name = Gene.Name) %>%
  dplyr::select(c(ORIG_Gene_ID, ORIG_Gene_Name, BM_Gene_Name, Study, 
                  Has_Rat_Hom, Has_Human_Hom, logFC, P.Value)) %>%
  unique()

# Create dataframe of directional info
mouse_reg <- all_mouse_sig_info %>%
  mutate(Direction = case_when(logFC > 0 ~ "Up", logFC < 0 ~ "Down")) %>%
  mutate(BM_Gene_Name = case_when(is.na(BM_Gene_Name) ~ ORIG_Gene_ID,
                                  TRUE ~ BM_Gene_Name)) %>%
  dplyr::select(-c(logFC, P.Value)) %>%
  pivot_wider(names_from = Study, values_from = Direction) %>%
  dplyr::rename(Carpenter_Direction = "Carpenter",
                Walker_Direction = "Walker") %>%
  mutate(Is_Same_Direction = 
           case_when(Walker_Direction == "Up" & Carpenter_Direction == "Up" ~ "Both_Up",
                     Walker_Direction == "Down" & Carpenter_Direction == "Down" ~ "Both_Down",
                     Walker_Direction == "Up" & Carpenter_Direction == "Down" ~ "WUp_CDown",
                     Walker_Direction == "Down" & Carpenter_Direction == "Up" ~ "WDown_CUp"))

# Get genes that are in the same direction for both studies
mouse_same_reg <- mouse_reg %>% filter(str_detect(Is_Same_Direction, "Both"))

# Pull Mouse IDs
mouse_id_same_reg <- mouse_same_reg %>% pull(ORIG_Gene_ID)

# How many genes are regulated in the same direction?
length(mouse_id_same_reg) # 10/15 (5/8 w/ LFC)

# Do the genes have human and rat homologs? 
# This is repetitive with earlier code to double-check results.
final_mouse_bm <- getBM(attributes = c("external_gene_name", # Gene symbol
                                       "ensembl_gene_id", # Mouse ENSEMBL
                                       "rnorvegicus_homolog_ensembl_gene", # Rat
                                       "hsapiens_homolog_ensembl_gene"), # Human
                        filter = "ensembl_gene_id", 
                        values = shared_mm_ids,
                        mart = grcm39) %>%
  dplyr::rename(BM_Gene_Name = external_gene_name,
                Mouse_ID = ensembl_gene_id,
                Rat_ID = rnorvegicus_homolog_ensembl_gene,
                Human_ID = hsapiens_homolog_ensembl_gene) %>%
  mutate(across(everything(), function(x) na_if(x, "")))

# Pull unique rat and human IDs that have mouse homologs
mouse_hs_ids <- final_mouse_bm %>% filter(!is.na(Human_ID), Mouse_ID %in% mouse_id_same_reg) %>% pull(Human_ID) %>% unique()
mouse_rn_ids <- final_mouse_bm %>% filter(!is.na(Rat_ID), Mouse_ID %in% mouse_id_same_reg) %>% pull(Rat_ID) %>% unique()

# How many unique rat and human genes?
length(mouse_hs_ids) # 9 (4 w/ LFC)
length(mouse_rn_ids) # 11 (4 w/ LFC)

# Get vectors of shared Carpenter/Walker genes, w/ and w/o direction
# All shared
carp_walk_mm <- shared_mm_ids
carp_walk_rn <- final_mouse_bm %>% filter(Mouse_ID %in% shared_mm_ids, !is.na(Rat_ID)) %>% pull(Rat_ID) %>% unique()
carp_walk_hs <- final_mouse_bm %>% filter(Mouse_ID %in% shared_mm_ids, !is.na(Human_ID)) %>% pull(Human_ID) %>% unique()
# Shared direction
carp_walk_mm_same <- mouse_id_same_reg
carp_walk_rn_same <- mouse_rn_ids
carp_walk_hs_same <- mouse_hs_ids



## OVERLAPS: RAT ----------
# Significant Powell genes, rat IDs
sig_powell_rat <- sig_df %>%
  filter(Study == "Powell") %>%
  pull(Original_Gene_ID)

# Get mouse homologs for Powell rat genes
powell_bm <- getBM(attributes = c("external_gene_name", # Gene symbol
                                  "ensembl_gene_id", # Mouse ENSEMBL
                                  "mmusculus_homolog_ensembl_gene", # Mouse
                                  "hsapiens_homolog_ensembl_gene"), # Human
                   filter = "ensembl_gene_id", 
                   values = sig_powell_rat,
                   mart = mRatBN7.2) %>%
  dplyr::rename(BM_Gene_Name = external_gene_name,
                Rat_ID = ensembl_gene_id,
                Mouse_ID = mmusculus_homolog_ensembl_gene,
                Human_ID = hsapiens_homolog_ensembl_gene) %>%
  mutate(across(everything(), function(x) na_if(x, "")))

# Signficant Powell genes, mouse homologs
# As a dataframe
sig_powell_mouse_hom_df <- sig_df %>%
  filter(Study == "Powell") %>%
  full_join(powell_bm, by = c("Original_Gene_ID" = "Rat_ID")) %>%
  mutate(across(everything(), function(x) na_if(x, ""))) %>%
  filter(!is.na(Mouse_ID))

# As a vector
sig_powell_mouse_hom <- sig_powell_mouse_hom_df %>%
  pull(Mouse_ID) %>%
  unique()
length(sig_powell_mouse_hom) # 548 genes

# How many mouse homologs up and downregulated? (268 + 280 = 548)
sig_powell_mouse_hom_df %>%
  mutate(Direction = case_when(logFC > 0 ~ "Up", logFC < 0 ~ "Down")) %>%
  dplyr::select(Mouse_ID, Direction) %>%
  unique() %>%
  group_by(Direction) %>%
  count()

# How many rat genes? (534 genes)
sig_powell_mouse_hom_df %>%
  dplyr::select(Original_Gene_ID) %>%
  unique() %>%
  count()

# How many rat genes up and downregulated? (264 + 270 = 534 genes)
sig_powell_mouse_hom_df %>%
  mutate(Direction = ifelse(logFC > 0, "Up", "Down")) %>%
  dplyr::select(Original_Gene_ID, Direction) %>%
  unique() %>%
  group_by(Direction) %>%
  count()

# Cross-reference with mouse dataframe
powell_mouse_id <- sig_mm_bm %>% filter(Mouse_ID %in% sig_powell_mouse_hom) %>% pull(Mouse_ID) %>% unique()
powell_rat_id <- sig_mm_bm %>% filter(Rat_ID %in% sig_powell_rat) %>% pull(Rat_ID) %>% unique()

# Check length - 46 genes (16 w/ LFC)
powell_mouse_id %>% length()
powell_rat_id %>% length()

# Shared genes and directionality for Carpenter-Walker is shown above

# Shared genes and directionality for Carpenter-Powell
# Shared genes
carp_pow_shared_mouse <- sig_df %>%
  filter(Study == "Carpenter", Original_Gene_ID %in% powell_mouse_id) %>%
  pull(Original_Gene_ID)

carp_ids <- sig_df %>%
  filter(Study == "Carpenter") %>%
  pull(Original_Gene_ID)
carp_pow_shared_rat <- sig_powell_mouse_hom_df %>%
  filter(Mouse_ID %in% carp_ids) %>%
  pull(Original_Gene_ID)

length(carp_pow_shared_mouse) # 45 genes (less b/c it doesn't include other dataset)
length(carp_pow_shared_rat) # 47 genes

# Direction
carp_pow_sig_info  <- sig_df %>%
  filter(Study != "Walker") %>%
  dplyr::rename(ORIG_Gene_ID = Original_Gene_ID,
                ORIG_Gene_Name = Gene.Name) %>%
  full_join(sig_mm_bm, by = c("ORIG_Gene_ID" = "Mouse_ID")) %>%
  filter(ORIG_Gene_ID %in% carp_pow_shared_mouse & Rat_ID %in% carp_pow_shared_rat | 
           ORIG_Gene_ID %in% carp_pow_shared_rat)

carp_pow_reg_info <- carp_pow_sig_info %>%
  mutate(Rat_ID = case_when(Study == "Powell" ~ ORIG_Gene_ID,
                            TRUE ~ Rat_ID)) %>%
  dplyr::select(logFC, Study, Rat_ID) %>%
  mutate(Direction = case_when(logFC > 0 ~ "Up", logFC < 0 ~ "Down")) %>%
  dplyr::select(-logFC) %>%
  pivot_wider(names_from = Study, values_from = Direction) %>%
  dplyr::rename(Powell_Direction = "Powell",
                Carpenter_Direction = "Carpenter") %>%
  filter(Carpenter_Direction != "NULL" & Powell_Direction != "NULL") %>%
  mutate(Is_Same_Direction = 
           case_when(Carpenter_Direction == "Up" & Powell_Direction == "Up" ~ "Both_Up",
                     Carpenter_Direction == "Down" & Powell_Direction == "Down" ~ "Both_Down",
                     Carpenter_Direction == "Up" & Powell_Direction == "Down" ~ "CUp_PDown",
                     Carpenter_Direction == "Down" & Powell_Direction == "Up" ~ "CDown_PUp",
                     TRUE ~ "Check_Manually"))

# How many are in the same direction? (20/45 genes)
rat_carp_same_reg <- carp_pow_reg_info %>% 
  filter(str_detect(Is_Same_Direction, "Both")) %>%
  pull(Rat_ID)
length(rat_carp_same_reg)

# Any duplicates? NO
sig_mm_bm %>%
  filter(Rat_ID %in% rat_carp_same_reg) %>%
  pull(Rat_ID) %>%
  duplicated() %>%
  sum()

# Any that need to be checked manually? YES
carp_pow_reg_info %>% filter(Is_Same_Direction == "Check_Manually") %>% View()

# Add ENSRNOG00000052224 to the shared genes by modifying dataframe
carp_pow_reg_info2 <- carp_pow_reg_info %>%
  mutate(Is_Same_Direction = case_when(Rat_ID == "ENSRNOG00000052224" ~ "Both_Down",
                                       Rat_ID == "ENSRNOG00000031855" ~ "Both_Up",
                                       TRUE ~ Is_Same_Direction))

# How many are in the same direction? (22/45 genes)
rat_carp_same_reg <- carp_pow_reg_info2 %>% 
  filter(str_detect(Is_Same_Direction, "Both")) %>%
  pull(Rat_ID)
length(rat_carp_same_reg)


# Shared genes and directionality for Powell-Walker
# Shared genes
walk_pow_shared_mouse <- sig_df %>%
  filter(Study == "Walker", Original_Gene_ID %in% powell_mouse_id) %>%
  pull(Original_Gene_ID)

walk_ids <- sig_df %>%
  filter(Study == "Walker") %>%
  pull(Original_Gene_ID)
walk_pow_shared_rat <- sig_powell_mouse_hom_df %>%
  filter(Mouse_ID %in% walk_ids) %>%
  pull(Original_Gene_ID)

length(walk_pow_shared_mouse) # 5 genes (less b/c it doesn't include other dataset)
length(walk_pow_shared_rat) # 5 genes

# Direction
walk_pow_sig_info  <- sig_df %>%
  filter(Study != "Carpenter") %>%
  dplyr::rename(ORIG_Gene_ID = Original_Gene_ID,
                ORIG_Gene_Name = Gene.Name) %>%
  full_join(sig_mm_bm, by = c("ORIG_Gene_ID" = "Mouse_ID")) %>%
  filter(ORIG_Gene_ID %in% walk_pow_shared_mouse & Rat_ID %in% walk_pow_shared_rat | 
           ORIG_Gene_ID %in% walk_pow_shared_rat)

walk_pow_reg_info <- walk_pow_sig_info %>%
  mutate(Rat_ID = case_when(Study == "Powell" ~ ORIG_Gene_ID,
                            TRUE ~ Rat_ID)) %>%
  dplyr::select(logFC, Study, Rat_ID) %>%
  mutate(Direction = case_when(logFC > 0 ~ "Up", logFC < 0 ~ "Down")) %>%
  dplyr::select(-logFC) %>%
  pivot_wider(names_from = Study, values_from = Direction) %>%
  dplyr::rename(Powell_Direction = "Powell",
                Walker_Direction = "Walker") %>%
  filter(Walker_Direction != "NULL" & Powell_Direction != "NULL") %>%
  mutate(Is_Same_Direction = 
           case_when(Walker_Direction == "Up" & Powell_Direction == "Up" ~ "Both_Up",
                     Walker_Direction == "Down" & Powell_Direction == "Down" ~ "Both_Down",
                     Walker_Direction == "Up" & Powell_Direction == "Down" ~ "CUp_PDown",
                     Walker_Direction == "Down" & Powell_Direction == "Up" ~ "CDown_PUp",
                     TRUE ~ "Check_Manually"))

# How many are in the same direction? (3/5 genes)
rat_walk_same_reg <- walk_pow_reg_info %>% 
  filter(str_detect(Is_Same_Direction, "Both")) %>%
  pull(Rat_ID)
length(rat_walk_same_reg)

# Any duplicates? NO
sig_mm_bm %>%
  filter(Rat_ID %in% rat_walk_same_reg) %>%
  pull(Rat_ID) %>%
  duplicated() %>%
  sum()

carp_only <- sig_df %>%
  full_join(sig_mm_bm, by = c("Original_Gene_ID" = "Mouse_ID", "Gene.Name" = "BM_Gene_Name")) %>%
  filter(Study == "Carpenter") %>%
  mutate(Direction = case_when(logFC > 0 ~ "Up", logFC < 0 ~ "Down")) %>%
  dplyr::select(Original_Gene_ID, Study, Direction)

carp_up <- carp_only %>% filter(Direction == "Up") %>% pull(Original_Gene_ID)
carp_down <- carp_only %>% filter(Direction == "Down") %>% pull(Original_Gene_ID)

walk_only <- sig_df %>%
  full_join(sig_mm_bm, by = c("Original_Gene_ID" = "Mouse_ID", "Gene.Name" = "BM_Gene_Name")) %>%
  filter(Study == "Walker") %>%
  mutate(Direction = case_when(logFC > 0 ~ "Up", logFC < 0 ~ "Down")) %>%
  dplyr::select(Original_Gene_ID, Study, Direction)

walk_up <- walk_only %>% filter(Direction == "Up") %>% pull(Original_Gene_ID)
walk_down <- walk_only %>% filter(Direction == "Down") %>% pull(Original_Gene_ID)

# Get vectors of shared Carpenter/Powell genes, w/ and w/o direction
# All shared (45 genes)
carp_pow_mm <- carp_pow_shared_mouse
carp_pow_rn <- sig_mm_bm %>% filter(Mouse_ID %in% carp_pow_shared_mouse, Rat_ID %in% carp_pow_shared_rat) %>% pull(Rat_ID) %>% unique()
carp_pow_hs <- sig_mm_bm %>% filter(Mouse_ID %in% carp_pow_shared_mouse) %>% filter(Human_ID != "NA") %>% pull(Human_ID) %>% unique()
# Shared direction (22 genes, 21 in human)
carp_pow_mm_same <- sig_mm_bm %>% filter(Rat_ID %in% rat_carp_same_reg) %>% pull(Mouse_ID) %>% unique()
carp_pow_rn_same <- rat_carp_same_reg
carp_pow_hs_same <- sig_mm_bm %>% filter(Rat_ID %in% rat_carp_same_reg) %>% filter(Human_ID != "NA") %>% pull(Human_ID) %>% unique()

# Get vectors of shared Walker/Powell genes, w/ and w/o direction
# All shared (5 genes)
walk_pow_mm <- walk_pow_shared_mouse
walk_pow_rn <- sig_mm_bm %>% filter(Mouse_ID %in% walk_pow_shared_mouse, Rat_ID %in% walk_pow_shared_rat) %>% pull(Rat_ID) %>% unique()
walk_pow_hs <- sig_mm_bm %>% filter(Mouse_ID %in% walk_pow_shared_mouse) %>% pull(Human_ID) %>% unique()
# Shared direction (3 genes)
walk_pow_mm_same <- sig_mm_bm %>% filter(Rat_ID %in% rat_walk_same_reg) %>% pull(Mouse_ID) %>% unique()
walk_pow_rn_same <- rat_walk_same_reg
walk_pow_hs_same <- sig_mm_bm %>% filter(Rat_ID %in% rat_walk_same_reg) %>% pull(Human_ID) %>% unique()

# Consolidate vectors
# All Craving
crave_mm <- c(carp_walk_mm, carp_pow_mm, walk_pow_mm) %>% unique() # 63
crave_rn <- c(carp_walk_rn, carp_pow_rn, walk_pow_rn) %>% unique() # 75
crave_hs <- c(carp_walk_hs, carp_pow_hs, walk_pow_hs) %>% unique() # 60
# Same direction
sd_crave_mm <- c(carp_walk_mm_same, carp_pow_mm_same, walk_pow_mm_same) %>% unique() # 34
sd_crave_rn <- c(carp_walk_rn_same, carp_pow_rn_same, walk_pow_rn_same) %>% unique() # 34
sd_crave_hs <- c(carp_walk_hs_same, carp_pow_hs_same, walk_pow_hs_same) %>% unique() # 33

sd_crave_hs %>% sort()


# All in one list
all_genes_list <- list(carp_pow_mm, carp_pow_rn, carp_pow_hs, carp_pow_mm_same, 
                       carp_pow_rn_same, carp_pow_hs_same, walk_pow_mm, 
                       walk_pow_rn, walk_pow_hs, walk_pow_mm_same, 
                       walk_pow_rn_same, walk_pow_hs_same, carp_walk_mm, 
                       carp_walk_rn, carp_walk_hs, carp_walk_mm_same, 
                       carp_walk_rn_same, carp_walk_hs_same, crave_mm, crave_rn, 
                       crave_hs, sd_crave_mm, sd_crave_rn, sd_crave_hs) %>%
  setNames(c("carp_pow_mm", "carp_pow_rn", "carp_pow_hs", "carp_pow_mm_same", 
             "carp_pow_rn_same", "carp_pow_hs_same", "walk_pow_mm",
             "walk_pow_rn", "walk_pow_hs", "walk_pow_mm_same", 
             "walk_pow_rn_same",  "walk_pow_hs_same", "carp_walk_mm", 
             "carp_walk_rn", "carp_walk_hs", "carp_walk_mm_same", 
             "carp_walk_rn_same", "carp_walk_hs_same", "crave_mm", "crave_rn", 
             "crave_hs", "sd_crave_mm", "sd_crave_rn", "sd_crave_hs"))

# Save big list as an R object
saveRDS(all_genes_list, file = paste0(main_dir, "post_processing/results/other/all_genes_list.RDS"))












