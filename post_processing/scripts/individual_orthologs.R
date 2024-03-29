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

# Other
library(clipr)



## BEFORE USING THIS SCRIPT ----------
# User should alter these variables as necessary
main_dir <- "G:/Annika/GitHub Repositories/RodentAddiction/"
de_dir <- paste0(main_dir, "post_processing/results/gene_info/")

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



## GET NEWEST MARTS ----------
# Mouse (GRCm39)
grcm39 <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl",
                     host = "https://apr2022.archive.ensembl.org")
# Rat (mRatBN7.2)
mRatBN7.2 <- useEnsembl("ensembl", dataset = "rnorvegicus_gene_ensembl", 
                        host = "https://apr2022.archive.ensembl.org")
# Human (GRCh38.p13)
grch38 <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", 
                     host = "https://apr2022.archive.ensembl.org")



## LOAD IN DIFFERENTIAL EXPRESSION DATA ----------
# Load in 3 relevant dataframes together as a list of lists

# Get list of files
DE_files <- list.files(path = de_dir, pattern = "DE_FINAL", full.names = TRUE)
DE_files

# Verify they are in the correct order with next line of code before renaming
DE_names <- list.files(path = de_dir, pattern = "DE_FINAL", full.names = FALSE)
DE_names
DE_listnames <- c("Carpenter", "Walker")

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
  Reduce(function(x, y) full_join(x, y), .) %>%
  mutate(Study = fct_relevel(Study, levels = c("Carpenter", "Walker"))) %>%
  rename(Original_Gene_ID = Gene_ID)



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

# Check that dimensions match number of genes (all should be TRUE)
nrow(mart_mm_names_carp) == 10026 # 10026
nrow(mart_mm_names_walk) == 9948 # 9948

# New dataframe with gene names, NA if missing
full_df <- Reduce(function(x, y) full_join(x, y),
                  list(full_df_nonames, mart_mm_names_carp, 
                       mart_mm_names_walk)) %>%
  relocate(Gene.Name, .after = Original_Gene_ID) %>%
  mutate(across(everything(), function(x) na_if(x, ""))) %>%
  filter(!is.na(Study))


## GET SIGNIFICANT GENES ----------
# Number of total genes per study after filtering
full_df %>% group_by(Study) %>% count()
# Study         n
# 1 Carpenter 10026
# 2 Walker    9948

# Get summary stats about fold changes for 1) all genes, and 2) DEGs only
full_df %>%
  group_by(Study) %>%
  mutate(logFC = abs(logFC), FC = 2^(abs(logFC))) %>%
  summarize(min_logFC = min(logFC), mean_logFC = mean(logFC), 
            median_logFC = median(logFC), max_logFC = max(logFC), 
            min_FC = min(FC), mean_FC = mean(FC), median_FC = median(FC), 
            max_FC = max(FC))

# Significant with uncorrected p-value only
sig_df <- full_df %>% filter(P.Value < 0.05)
sig_df %>%
  group_by(Study) %>%
  mutate(logFC = abs(logFC), FC = 2^(abs(logFC))) %>%
  summarize(min_logFC = min(logFC), mean_logFC = mean(logFC), 
            median_logFC = median(logFC), max_logFC = max(logFC), 
            min_FC = min(FC), mean_FC = mean(FC), median_FC = median(FC), 
            max_FC = max(FC))

# Sig with uncorrected p-value and fold change threshold of 2 (FC = 1.25)
sig_p_fc_df <- full_df %>% 
  mutate(logFC = abs(logFC), FC = 2^(abs(logFC))) %>% 
  filter(P.Value < 0.05, FC >= 1.25)
sig_p_fc_df %>% group_by(Study) %>% count()
# Study         n
# 1 Carpenter   94
# 2 Walker      91

# Sig with FDR < 0.1 
sig_df_fdr <- full_df %>% filter(adj.P.Val < 0.1)
nrow(sig_df_fdr) # No significant genes



## VOLCANO PLOTS ----------
# How many DE genes for each study/comparison?
sig_df %>% group_by(Study) %>% count()
# Study         n
# 1 Carpenter   633
# 2 Walker      139

# Number of DE genes for each study, as variables
carp_count <- sig_df %>% filter(Study == "Carpenter") %>% count() %>% pull(n)
walk_count <- sig_df %>% filter(Study == "Walker") %>% count() %>% pull(n)

# How many DE genes for each study/comparison that are up or down-regulated?
sig_df %>% 
  mutate(Direction = ifelse(logFC > 0, "Up", "Down")) %>% 
  group_by(Study, Direction) %>% count()
# Study     Direction     n
# 1 Carpenter Down        329
# 2 Carpenter Up          304
# 3 Walker    Down        100
# 4 Walker    Up           39

# Function for creating volcano plots
volcano_plot <- function(df, study, rect_fill, x_lim = NULL, x_breaks = waiver(),
                         y_lim = 5, x_color = "black", y_color = "black",  
                         x_axis_lab = waiver()) {
  
  # Add labels to dataframe for significant genes
  new_df <- df %>%
    mutate(DE = fct_relevel(case_when((logFC > log2(1) & P.Value < 0.05 ~ "Up"),
                                      (logFC < -log2(1) & P.Value < 0.05 ~ "Down"),
                                      TRUE ~ "None"),
                            levels = c("None", "Down", "Up")),
           Abs_logFC = abs(logFC),
           Study = case_when(str_detect(Study, "Walker") ~ "Walker",
                             str_detect(Study, "Carpenter") ~ "Carpenter")) %>%
    group_by(Study, DE) %>%
    arrange(-Abs_logFC, .by_group = TRUE)
  
  # Volcano plot w/ labels for top 10 up- and down-regulated genes
  new_df %>%
    filter(Study == study) %>%
    ggplot(aes(x = logFC, y = -log10(P.Value), fill = DE)) +
    geom_point(alpha = 0.35, color = "black", shape = 21, 
               size = 4, stroke = 1, show.legend = FALSE) +
    geom_hline(yintercept = -log10(0.05), col = "grey10", 
               lty = "dashed", lwd = 1.25) + # P-value threshold
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

vo1 <- volcano_plot(full_df, "Carpenter", "#44AA99", c(-1.5, 1.5), seq(-1.5, 1.5, 0.5), 4.2, "black", "black")
vo2 <- volcano_plot(full_df, "Walker", "#C148AD", c(-2, 2), seq(-2, 2, 1), 5, "black", "black")

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Carp_Walk_volcano_09-18-2022.svg", width = 8, height = 4.25)
ggarrange(vo1, vo2)
dev.off()

pdf("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Carp_Walk_volcano_09-18-2022.pdf", width = 8, height = 4.25)
ggarrange(vo1, vo2)
dev.off()



## HUMAN AND RODENT HOMOLOGS ----------
### CARPENTER -----
sig_carp_ids <- sig_df %>% 
  filter(Study == "Carpenter") %>% 
  pull(Original_Gene_ID)
length(sig_carp_ids) # 633

# Get dataframe with all relevant orthologs
sig_carp_bm <- getBM(attributes = c("external_gene_name", # Gene symbol
                                    "ensembl_gene_id", # Mouse ENSEMBL
                                    "rnorvegicus_homolog_ensembl_gene", # Rat
                                    "hsapiens_homolog_ensembl_gene"), # Human
                     filter = "ensembl_gene_id", 
                     values = sig_carp_ids,
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

# Pull human, mouse, and rat IDs
sig_carp_mm_ids <- sig_carp_bm %>% filter(!is.na(Mouse_ID)) %>% pull(Mouse_ID) %>% unique() # 633
sig_carp_rn_ids <- sig_carp_bm %>% filter(!is.na(Rat_ID)) %>% pull(Rat_ID) %>% unique() # 590
sig_carp_hs_ids <- sig_carp_bm %>% filter(!is.na(Human_ID)) %>% pull(Human_ID) %>% unique() # 585


### WALKER -----
sig_walk_ids <- sig_df %>% 
  filter(Study == "Walker") %>% 
  pull(Original_Gene_ID)
length(sig_walk_ids) # 139

# Get dataframe with all relevant orthologs
sig_walk_bm <- getBM(attributes = c("external_gene_name", # Gene symbol
                                    "ensembl_gene_id", # Mouse ENSEMBL
                                    "rnorvegicus_homolog_ensembl_gene", # Rat
                                    "hsapiens_homolog_ensembl_gene"), # Human
                     filter = "ensembl_gene_id", 
                     values = sig_walk_ids,
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

# Pull human, mouse, and rat IDs
sig_walk_mm_ids <- sig_walk_bm %>% filter(!is.na(Mouse_ID)) %>% pull(Mouse_ID) %>% unique() # 139
sig_walk_rn_ids <- sig_walk_bm %>% filter(!is.na(Rat_ID)) %>% pull(Rat_ID) %>% unique() # 123
sig_walk_hs_ids <- sig_walk_bm %>% filter(!is.na(Human_ID)) %>% pull(Human_ID) %>% unique() # 116



## ANY GENES SHARED? ----------
# Shared in both datasets (mouse gene/ortholog)
shared_genes <- data.frame(Gene = c(sig_carp_mm_ids, sig_walk_mm_ids)) %>%
  group_by(Gene) %>%
  mutate(Count = length(Gene)) %>%
  arrange(desc(Count)) %>%
  unique() %>%
  filter(Count == 2) %>%
  pull(Gene)
shared_genes

# Create dataframe of directional info
sig_df_direction <- sig_df %>%
  mutate(Direction = case_when(logFC > 0 ~ "Up", logFC < 0 ~ "Down")) %>%
  select(Original_Gene_ID, Gene.Name, Study, Direction) %>%
  pivot_wider(names_from = Study, values_from = Direction) %>%
  rename(Carpenter_Direction = "Carpenter", Walker_Direction = "Walker") %>%
  mutate(Is_Same_Direction = 
           case_when(Walker_Direction == "Up" & Carpenter_Direction == "Up" ~ "BothUp",
                     Walker_Direction == "Down" & Carpenter_Direction == "Down" ~ "BothDown",
                     Walker_Direction == "Up" & Carpenter_Direction == "Down" ~ "WUp_CDown",
                     Walker_Direction == "Down" & Carpenter_Direction == "Up" ~ "WDown_CUp")) %>%
  select(Original_Gene_ID, Gene.Name, ends_with("Direction"))

# Get genes that are in the same or different direction for both studies
same_reg <- sig_df_direction %>% filter(str_detect(Is_Same_Direction, "Both"))
diff_reg <- sig_df_direction %>% filter(str_detect(Is_Same_Direction, "_"))

# Pull Mouse IDs
same_reg_ids <- same_reg %>% select(Original_Gene_ID, Gene.Name)
diff_reg_ids <- diff_reg %>% select(Original_Gene_ID, Gene.Name)

same_reg_ids
diff_reg_ids

# How many genes are regulated in the same direction?
nrow(same_reg_ids) # 10/15
nrow(diff_reg_ids) # 5/15



## CREATE RDS OBJECT ----------
# All in one list
all_genes_list <- list(sig_carp_mm_ids, sig_carp_rn_ids, sig_carp_hs_ids,
                       sig_walk_mm_ids, sig_walk_rn_ids, sig_walk_hs_ids) %>%
  setNames(c("sig_carp_mm_ids", "sig_carp_rn_ids", "sig_carp_hs_ids",
             "sig_walk_mm_ids", "sig_walk_rn_ids", "sig_walk_hs_ids"))

# Save big list as an R object
saveRDS(all_genes_list, 
        file = file.path(main_dir, "post_processing/results/other/all_genes_list.RDS"))





