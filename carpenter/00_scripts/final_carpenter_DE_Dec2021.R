
# PACKAGES ----------
# Data Manipulation
library(tidyverse)
library(reshape2)

# Modeling & Differential Expression
library(edgeR)
library(limma)
library(variancePartition)

# Plots
library(ggpubr)
library(pals)
library(RColorBrewer)
library(wesanderson)
library(svglite)



# BEFORE USING THIS SCRIPT ----------
# User should alter these variables as necessary
main_dir <- "G:/Noctsol/GitHub Repositories/temp/"
pheno_dir <- paste0(main_dir, "carpenter/")
gene_dir <- paste0(main_dir, "carpenter/genes/")
setwd(paste0(main_dir, "results/"))



# LOAD IN DATA & RENAME COLUMNS ----------
# NOTE: Sample order must be maintained for accuracy.
# Sample order is always numerical from low (SRS5770730) to high (SRS5770741).

# Load in phenotype file
pheno <- read_csv(paste0(pheno_dir, "carpenter_pheno.csv")) %>% 
  filter(Abstinence_Length == "28")

# Get sample IDs
sample_ids <- pheno$Sample_Treatment

# Load in gene and transcript counts files
g_count <- read_csv(paste0(gene_dir, "carpenter_gene_count_matrix_28abs.csv"))

# Separate initial gene_id column into Gene.ID & Gene.Name
# E.g. ENSMUSG00000102348|Gm10568 becomes ENSMUSG00000102348 & Gm10568
fix_gene_ids <- function(df) {separate(df, col = "gene_id", 
                                       into = c("Gene.ID", "Gene.Name"), 
                                       sep = "\\|")
}

g_count <- fix_gene_ids(g_count)

# The gene matrix with raw counts doesn't have gene locations and lengths, but 
# the original StringTie files do. Get these and merge into expression matrix.

# Get list of files, then load in all files
g_files <- list.files(path = gene_dir, 
                      pattern = "*_genes.txt", full.names = TRUE)

# Only select columns; not reading in Strand, Coverage, or TPM
# Keep FPKM to use downstream
list2env(envir = .GlobalEnv,
         lapply(setNames(g_files, make.names(sample_ids)), 
                read.delim, sep = "\t", 
                colClasses = c("character", "character", "factor", 
                               "NULL", "numeric", "numeric", 
                               "NULL", "numeric", "NULL")))

# Create list of lists for easier transformation into data frame
sample_list <- list("SRS5770730_C28abs" = SRS5770730_C28abs,
                    "SRS5770731_C28abs" = SRS5770731_C28abs,
                    "SRS5770732_C28abs" = SRS5770732_C28abs,
                    "SRS5770733_C28abs" = SRS5770733_C28abs,
                    "SRS5770734_C28abs" = SRS5770734_C28abs,
                    "SRS5770735_C28abs" = SRS5770735_C28abs,
                    "SRS5770736_S28abs" = SRS5770736_S28abs,
                    "SRS5770737_S28abs" = SRS5770737_S28abs,
                    "SRS5770738_S28abs" = SRS5770738_S28abs,
                    "SRS5770739_S28abs" = SRS5770739_S28abs,
                    "SRS5770740_S28abs" = SRS5770740_S28abs,
                    "SRS5770741_S28abs" = SRS5770741_S28abs)

# Sanity check: Sample IDs in sample_list & pheno file are in the same order
identical(names(sample_list), sample_ids) # Should be TRUE

# Convert list of lists to single data frame
g_original <- Reduce(function(x, y) 
  full_join(x, y, by = c("Gene.ID", "Gene.Name", "Reference",
                         "Start", "End")), sample_list)

# Prepare data: 1) Rename columns with sample IDs, 2) Calculate gene length, 
# 3) Add length as a suffix so duplicates can be differentiated
gene_fpkm <- g_original %>% 
  rename_with(function(x) sample_ids, .cols = c(6:ncol(.))) %>%
  mutate(Length = End - Start + 1, .after = End) %>%
  unite(Gene.ID.Length, Gene.ID, Length, sep = "-len", remove = FALSE)



# FIND DUPLICATE GENES ----------
# NOTE: Duplicates (same Gene.ID & Gene.Name) can be tracked w/ Gene.ID.Length
# Before fixing duplicates, there are 53,718 genes
dim(gene_fpkm) # 53718, 19

# Find duplicate genes
# Sometimes genes have duplicated names but different IDs, so check by Gene.ID
dup_genes_info <- gene_fpkm %>% 
  group_by(across(Gene.ID)) %>% 
  filter(n() > 1) %>% 
  ungroup()

# Number of duplicates
nrow(dup_genes_info) # 26

# Get *unique* duplicated genes for each data set
# Same ID and chromosome (NOT paralogs)
distinct_dups <- dup_genes_info %>%
  distinct(Gene.ID, Gene.Name, Reference) %>%
  arrange(Gene.Name)

# Number of unique duplicates
nrow(distinct_dups) # 8

# Sanity check: Duplicated genes should be the only genes with NAs in the data
# Identify genes with NAs in 1 or more samples (FPKM)
na_genes_info <- gene_fpkm[!complete.cases(gene_fpkm), ]

# Number of NAs
nrow(na_genes_info) # 23

# Get *unique* NAs
distinct_nas <- na_genes_info %>%
  distinct(Gene.ID, Gene.Name, Reference) %>%
  arrange(Gene.Name)

# Number of unique NAs
nrow(distinct_nas) # 7

# The duplicate genes and NA genes are not the same
all_equal(distinct_dups, distinct_nas) # "Different number of rows"



## COMBINE "DUPLICATES" (SUM EXPRESSION) ----------
# Sum duplicates
dup_fixed <- dup_genes_info %>%
  group_by(Gene.ID, Gene.Name, Reference) %>%
  summarize_at(vars("SRS5770730_C28abs":"SRS5770741_S28abs"),
               sum, na.rm = TRUE) %>%
  mutate(Start = NA, End = NA, Length = NA, .after = Reference)

# Remove duplicates from original dataset, then add back 
gene_fpkm_nodup <- gene_fpkm[
  !gene_fpkm$Gene.ID %in% distinct_dups$Gene.ID, ] %>%
  full_join(dup_fixed)

# Check that this was done correctly (# of genes by FPKM matches raw counts)
nrow(gene_fpkm_nodup) == nrow(g_count) # Should be TRUE



## FILTER LOWLY-EXPRESSED GENES (FPKM) ----------
# Number of genes before filtering
dim(gene_fpkm_nodup) # 53,700 genes, 19 columns (12 samples)

# Get mean FPKM per treatment per gene
C28abs_means <- gene_fpkm_nodup %>%
  select(ends_with("C28abs")) %>%
  apply(1, mean)

S28abs_means <- gene_fpkm_nodup %>%
  select(ends_with("S28abs")) %>%
  apply(1, mean)

# Create data frames of FPKM means
fpkm_means <- data.frame(
  cbind("Gene.ID" = gene_fpkm_nodup$Gene.ID,
        "Reference" = gene_fpkm_nodup$Reference,
        "C28abs" = C28abs_means, "S28abs" = S28abs_means))

# Only keep genes that have a mean of > 0.5 FPKM in 1 treatment group
keep_genes <- fpkm_means %>% filter(C28abs > 0.5 | S28abs > 0.5)

# Remove genes & verify that this has been done correctly
keep_gene_fpkm <- gene_fpkm_nodup %>%
  semi_join(keep_genes, by = "Gene.ID") %>%
  select(-Gene.ID.Length)

nrow(keep_genes) == nrow(keep_gene_fpkm) # Should be TRUE

# Filter by raw read count, at least 6 in at least 6 samples (size of smallest treatment group)
keep_count <- rowSums(g_count >= 6) >= 6
keep_count <- g_count[keep_count, ]

# Get full list of "keep" genes from raw counts data frame
final_genes <- semi_join(keep_count, keep_gene_fpkm, by = "Gene.ID")
nrow(final_genes) # Total genes after filtering out low expression: 15,837



## WRITE GENE INFO & COUNT FILES ----------
# Separate into data frames:
# 1) raw counts after filtering by FPKM, 2) gene info (ID, name, chromosome)

# Raw counts after filtering
final_gene_matrix <- final_genes[, -(1:2)]

# Gene info
gene_info <- keep_gene_fpkm %>% 
  filter(Gene.ID %in% (final_genes %>% pull(Gene.ID))) %>%
  select("Gene.ID", "Gene.Name", "Reference") %>%
  dplyr::rename("Chromosome" = "Reference")

# Write files (Note that 28abs dataframe & count matrix were created earlier)
write.table(final_gene_matrix, file = "carpenter_FINAL_gene_FPKM.txt")
write.table(gene_info, file = "carpenter_FINAL_gene_info.txt")



## CREATE & MODIFY DGE LIST ----------
# Create initial DGE List object using filtered counts (not FPKM)
gene_dge <- DGEList(counts = final_gene_matrix, genes = gene_info)

# Add pheno info to samples
gene_dge$samples$treatment <- factor(pheno$Treatment_Group, 
                                     levels = c("C28abs", "S28abs"))

# Convert to log2 CPM (not normalized)
gene_lcpm <- cpm(gene_dge, log = TRUE)



## MULTI-DIMENSIONAL SCALING PLOTS ----------
# Use MDS plots to check for outliers

# Get treatment variables
treatment <- gene_dge$samples$treatment %>%
  gsub("abs", "", .) %>%
  fct_rev()

# Get and set colors for variables
color_treatment <- treatment
levels(color_treatment) <- c("#377EB8", "#E41A1C")
color_treatment <- as.character(color_treatment)

# Simplify sample names
simple_sample <- rownames(gene_dge$samples) %>%
  gsub("_.{6}", "", .) %>%
  gsub("SRS.{4}", "SRS-", .)

# MDS plots
par(mfrow = c(2, 2))

plotMDS(gene_lcpm, main = "Treatment Group",
        label = treatment, col = color_treatment, 
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm, main = "Samples", label = simple_sample,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm, main = "Treatment Group",
        label = treatment, col = color_treatment, 
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)
plotMDS(gene_lcpm, main = "Samples", label = simple_sample,
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)


## VOOM TRANSFORMATION ----------
# Design: Treatment (S28abs, C28abs)
# Note that this is the same as using Drug (Saline, Cocaine)
design <- model.matrix( ~ 0 + gene_dge$samples$treatment)
colnames(design) <- gsub("gene_dge\\$samples\\$treatment", "", 
                         colnames(design))

# Get TMM-normalization factors before using voom
par(mfrow = c(1, 1))
gene_dge <- calcNormFactors(gene_dge, normalize = "TMM")
v <- voom(gene_dge, design, plot = TRUE)

# MDS Plots based on voom models
par(mfrow = c(1, 2))

# Treatment
plotMDS(v, main = "Treatment Group", pch = 16, cex = 2,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(v, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)

## VARIANCE PARTITION -----------
# No need to do variance partition, since all parameters (Sequencing Lane and 
# Instrument, Batch, Abstinence Length, etc.) are the same other than Drug.



## DIFFERENTIAL EXPRESSION ANALYSIS ---
# Create contrast matrices
cm <- makeContrasts(Drug_28dAbs = S28abs - C28abs,
                    levels = levels(gene_dge$samples$treatment))

# Fit the model using limma
fm <- lmFit(v, design)

# Fit the contrasts
fc <- contrasts.fit(fm, contrasts = cm)

# Run contrasts
results <- decideTests(fc)

# Calculate the t-statistics for the contrasts
# Robustified against outlier sample variances
ebayes <- eBayes(fc)

# Final figure
svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Carpenter_RNAseq.svg", width = 6.5, height = 10)
par(mfrow = c(3, 2))
plotMDS(gene_lcpm,
        pch = 18, cex = 3, col = color_treatment, cex.axis = 1.5, cex.lab = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm,
        pch = 18, cex = 3, col = color_treatment, cex.axis = 1.5, cex.lab = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)
# Original voom data
voom_info <- voom(gene_dge, design, plot = FALSE, save.plot = TRUE)
plot(voom_info$voom.xy, pch = 19, cex = 0.5, cex.main = 1.5, cex.axis = 1.5, 
     cex.lab = 1.5, main = "voom: Mean-variance trend", 
     xlab = "log2(count size + 0.5)", ylab = "Sqrt(standard deviation)")
lines(voom_info$voom.line, col = "red")
# Final voom model
plotSA(ebayes, main = "Final model: Mean-variance trend", 
       cex = 0.5, pch = 19, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
plotMDS(v, pch = 18, cex = 3, cex.axis = 1.5, cex.lab = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(v, pch = 18, cex = 3, cex.axis = 1.5, cex.lab = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)
dev.off()

# Final figure for presentation
svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Carpenter_RNAseq_simple_presentation.svg", width = 4.5, height = 4.5)
par(mfrow = c(1, 1))
plotMDS(gene_lcpm,
        pch = 16, cex = 3, col = color_treatment, cex.axis = 1.5, cex.lab = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
dev.off()


# Get stats for all genes, in same order as they were originally
stats <- topTable(ebayes, number = nrow(ebayes), sort.by = "none")

stats %>% filter(P.Value < 0.05) %>% nrow() # 1016 genes
stats %>% filter(P.Value < 0.05, abs(logFC) > log2(1.2)) %>% nrow() # 374 genes
stats %>% filter(adj.P.Val < 0.1) %>% nrow() # None

# Write files
write.table(stats, file = "carpenter_DE_FINAL.txt")



## CREATE OBJECT FOR POWER ANALYSIS ----------
carpenter_dge <- gene_dge # Already calculated norm. factors for TMM method
carpenter_stats <- stats
carpenter_summary <- carpenter_stats %>%
  filter(P.Value < 0.05) %>%
  summarize(Up = sum(logFC > 0), Down = sum(logFC < 0), Total = Up + Down)
carpenter_gene_mat <- final_gene_matrix
carpenter_gene_info <- gene_info
carpenter_fc <- carpenter_stats %>%
  filter(P.Value < 0.05) %>%
  mutate(Abs_FC = 2^abs(logFC)) %>%
  dplyr::select(Abs_FC) %>%
  summary()
carpenter_pheno <- pheno
carpenter_treat <- fct_relevel(pheno$Treatment_Group, levels = c("S28abs", "C28abs"))
carpenter_control <- levels(carpenter_treat)[1]

carpenter_obj <- list(DGE = carpenter_dge, Stats = carpenter_stats, 
                      Stats_Summary = carpenter_summary, 
                      Gene_Matrix = carpenter_gene_mat, 
                      Gene_Info = carpenter_gene_info, 
                      FC = carpenter_fc, Pheno = carpenter_pheno, 
                      Treat = carpenter_treat, 
                      Control = carpenter_control)

saveRDS(carpenter_obj, file = "carpenter_obj_4power.RDS")

