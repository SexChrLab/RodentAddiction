
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
main_dir <- "C:/Annika/GitHub Repositories/RodentAddiction/"
pheno_dir <- paste0(main_dir, "carpenter/00_misc_and_scripts/miscellaneous/")
gene_dir <- paste0(main_dir, "carpenter/00_misc_and_scripts/gene_counts/")
setwd(paste0(main_dir, "results/"))

# Scientific notation is a pain
options(scipen = 99999)

# Update and set ggplot theme globally
new_theme <- theme_bw(base_size = 14) %+replace%
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))
theme_set(new_theme)



# LOAD IN DATA & RENAME COLUMNS ----------
# NOTE: Sample order must be maintained for accuracy.
# Sample order is always numerical from low (SRS5770694) to high (SRS5770741).

# Load in phenotype file
pheno <- read_csv(paste0(pheno_dir, "carpenter_pheno.csv")) %>%
  filter(Abstinence_Length == "28")
sample_ids <- pheno %>% 
  filter(Abstinence_Length == "28") %>% 
  pull(Sample_Treatment) # Vector of sample IDs


# Get list of files, then load in all files
g_files <- list.files(path = gene_dir, 
                      pattern = "*_genes_s2.txt", full.names = TRUE)

# Only select columns; not reading in Strand, Coverage, or TPM
# Keep FPKM to use downstream
list2env(envir = .GlobalEnv,
         lapply(setNames(g_files, make.names(sample_ids)), 
                read.delim, sep = "\t",
                colClasses = c("character", "character", "character", 
                               "character", "character", "numeric", 
                               "numeric"),
                col.names = c("Gene_ID", "Chr", "Start", "End", "Strand", 
                              "Length", "Count"),
                skip = 1))

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
                    "SRS5770741_S28abs" = SRS5770741_S28abs) %>%
  map2(., names(.), ~cbind(.x, Sample = .y))

# Sanity check: Sample IDs in sample_list & pheno file are in the same order
identical(names(sample_list), sample_ids) # Should be TRUE

# Convert list of lists to single data frame
g_original <- Reduce(function(x, y) 
  full_join(x, y), sample_list) %>%
  pivot_wider(names_from = Sample, values_from = Count)

# Are there any discrepancies in the dataset?
# Search for NAs - No NAs
g_original %>% mutate(across(everything(), as.factor)) %>% summary()
# Search for duplicate genes - No duplicates (n = 1, not higher)
g_original %>% group_by(Gene_ID) %>% count() %>% arrange(desc(n)) %>% head()



## FILTER LOWLY-EXPRESSED GENES (FPKM) ----------
# To create DGE object, split matrix from gene info and raw counts
gene_matrix <- g_original %>% select(starts_with("SRS"))
gene_info <- g_original %>% select(!starts_with("SRS"))
gene_dge <- DGEList(counts = gene_matrix, genes = gene_info)

# Add phenotype info - treatment group
gene_dge$samples$treatment <- factor(pheno$Treatment_Group, 
                                     levels = c("S28abs", "C28abs"))

# Check number of genes
dim(gene_dge) # 55,414 genes, 12 samples/columns

# How many genes aren't expressed in any samples?
table(rowSums(gene_dge$counts==0)==12)
# TRUE (unexpressed) = 25,731

# Convert raw counts to FPKM using edgeR and create new dataframe
norm_dge <- calcNormFactors(gene_dge, method = "TMM")
fpkm <- rpkm(norm_dge, gene.length = norm_dge$genes$Length)
fpkm_df <- cbind(gene_info, fpkm)

# Verify that the dimensions of the new dataframe are the same as the old one
all.equal(dim(fpkm_df), dim(g_original)) # TRUE

# Get mean FPKM per treatment per gene
S28abs_means <- as.data.frame(fpkm) %>%
  dplyr::select(ends_with("S28abs")) %>%
  apply(1, mean, na.rm = TRUE) %>%
  mutate(as.data.frame(fpkm_df), S28abs_mean = .)

C28abs_means <- as.data.frame(fpkm) %>%
  dplyr::select(ends_with("C28abs")) %>%
  apply(1, mean, na.rm = TRUE) %>%
  mutate(as.data.frame(fpkm_df), C28abs_mean = .)

# Create single data frame of FPKM means
fpkm_means <- full_join(S28abs_means, C28abs_means)

# Check that row number is still the same
dim(fpkm_means)

# Only keep genes that are > 0.5 FPKM in 1 treatment group
keep_genes_fpkm <- fpkm_means %>% 
  filter(S28abs_mean > 0.5 | C28abs_mean > 0.5) %>%
  pull(Gene_ID)

# Must also have a read count of at least 6 for 2 (min. sample size) samples
keep_genes_count <- g_original %>%
  mutate(across(starts_with("SRS"), ~ ifelse(.x >=6, TRUE, FALSE))) %>%
  group_by(Gene_ID) %>%
  mutate(Qualifying_Samples = sum(across(starts_with("SRS")))) %>%
  ungroup() %>%
  filter(Qualifying_Samples >= 6) %>%
  pull(Gene_ID)

# Remove genes from count dataframe & verify that this has been done correctly
g_filtered <- g_original %>% filter(Gene_ID %in% keep_genes_fpkm &
                                      Gene_ID %in% keep_genes_count)

# Looking at dataframes before and after, it looks like we did get rid of 0's
g_original
g_filtered

# Number of genes after filtering by FPKM
nrow(g_filtered) # 10,026 genes




## MAKE GENE INFO & COUNT MATRICIES ----------
# To create DGE object, split matrix from gene info and raw counts
filtered_gene_matrix <- g_filtered %>% select(starts_with("SRS"))
filtered_gene_info <- g_filtered %>% select(!starts_with("SRS"))
filtered_gene_dge <- DGEList(counts = filtered_gene_matrix, 
                             genes = filtered_gene_info)



## CREATE & MODIFY DGE LIST ----------
# Create initial DGE List object using filtered counts (not FPKM)
filtered_gene_dge <- DGEList(counts = filtered_gene_matrix,
                             genes = filtered_gene_info)

# Add treatment group to samples
filtered_gene_dge$samples$treatment <- factor(pheno$Treatment_Group, 
                                              levels = c("S28abs", "C28abs"))


## CHECK FILTERING ----------
# Convert to log2 CPM (not normalized)
gene_lcpm <- cpm(gene_dge, log = TRUE)
filtered_gene_lcpm <- cpm(filtered_gene_dge, log = TRUE)

# Create dataframes from the original and filtered log2 CPM
og_df <- as.data.frame(gene_lcpm) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "log2_CPM")
fgl_df <- as.data.frame(filtered_gene_lcpm) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "log2_CPM")

# Look at plots
og_df %>%
  ggplot(aes(x = log2_CPM, color = Sample)) +
  geom_density()
fgl_df %>%
  ggplot(aes(x = log2_CPM, color = Sample)) +
  geom_density() +
  scale_x_continuous(limits = c(-5, 20))
# It looks like filtering got rid of the lowly-expressed samples and made things
# more normally distributed



## MULTI-DIMENSIONAL SCALING PLOTS ----------
# Use MDS plots to check for outliers

# Get treatment variables
treatment <- filtered_gene_dge$samples$treatment %>%
  gsub("abs", "", .) %>%
  fct_rev()

# Get and set colors for variables
color_treatment <- treatment
levels(color_treatment) <- c("#377EB8", "#E41A1C")
color_treatment <- as.character(color_treatment)

# Simplify sample names
simple_sample <- rownames(filtered_gene_dge$samples) %>%
  gsub("_.{6}", "", .) %>%
  gsub("SRS.{4}", "SRS-", .)

# MDS plots
par(mfrow = c(2, 2))

# No outliers
plotMDS(filtered_gene_lcpm, main = "Treatment Group",
        label = treatment, col = color_treatment, 
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(filtered_gene_lcpm, main = "Samples", label = simple_sample,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)

plotMDS(filtered_gene_lcpm, main = "Treatment Group",
        label = treatment, col = color_treatment, 
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)
plotMDS(filtered_gene_lcpm, main = "Samples", label = simple_sample,
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)

# Create final dge and lcpm
final_gene_dge <- filtered_gene_dge
final_gene_lcpm <- filtered_gene_lcpm



## MAKE GENE INFO & COUNT MATRICIES ----------
# To create DGE object, split matrix from gene info and raw counts
final_gene_matrix <- filtered_gene_matrix
final_gene_info <- filtered_gene_info

# Create final dge and lcpm
final_gene_dge <- filtered_gene_dge
final_gene_lcpm <- filtered_gene_lcpm


## VOOM TRANSFORMATION ----------
# Design: Treatment (S28abs, C28abs)
# Note that this is the same as using Drug (Saline, Cocaine)
design <- model.matrix( ~ 0 + final_gene_dge$samples$treatment)
colnames(design) <- gsub("final_gene_dge\\$samples\\$treatment", "", 
                         colnames(design))

# Get TMM-normalization factors before using voom
par(mfrow = c(1, 1))
final_gene_dge <- calcNormFactors(final_gene_dge, normalize = "TMM")
v <- voom(final_gene_dge, design, plot = TRUE)

# MDS Plots based on voom models
par(mfrow = c(1, 2))

# Treatment
plotMDS(v, main = "Treatment Group", pch = 16, cex = 2,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(v, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)


## VARIANCE PARTITION -----------
# No need to do variance partition, since all parameters (Sequencing Lane and 
# Instrument, Batch, Abstinence Length, etc.) are the same other than Drug.


## DIFFERENTIAL EXPRESSION ANALYSIS ---
# Create contrast matrices
cm <- makeContrasts(Drug_28dAbs = S28abs - C28abs,
                    levels = levels(final_gene_dge$samples$treatment))

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
# svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Carpenter_RNAseq_1_24_2022.svg", width = 6.5, height = 10)
par(mfrow = c(3, 2))
plotMDS(final_gene_lcpm,
        pch = 18, cex = 3, col = color_treatment, cex.axis = 1.5, cex.lab = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(final_gene_lcpm,
        pch = 18, cex = 3, col = color_treatment, cex.axis = 1.5, cex.lab = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)
# Original voom data
voom_info <- voom(final_gene_dge, design, plot = FALSE, save.plot = TRUE)
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
# dev.off()

# Final figure for presentation
# svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Carpenter_RNAseq_simple_presentation_1_24_2022.svg", width = 4.6, height = 4.5)
par(mfrow = c(1, 1))
plotMDS(final_gene_lcpm,
        pch = 16, cex = 3, col = color_treatment, cex.axis = 1.5, cex.lab = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
# dev.off()


# Get stats for all genes, in same order as they were originally
stats <- topTable(ebayes, number = nrow(ebayes), sort.by = "none")

stats %>% filter(P.Value < 0.05) %>% nrow() # 633 genes
stats %>% filter(P.Value < 0.05, abs(logFC) > log2(1.2)) %>% nrow() # 136 genes
stats %>% filter(adj.P.Val < 0.1) %>% nrow() # None

# Write files
write.table(stats, file = paste0(main_dir, "post_processing/results/gene_info/carpenter_DE_FINAL.txt"))


## CREATE OBJECT FOR POWER ANALYSIS ----------
carpenter_dge <- final_gene_dge # Already calculated norm. factors for TMM method
carpenter_stats <- stats
carpenter_summary <- carpenter_stats %>%
  filter(P.Value < 0.05) %>%
  summarize(Up = sum(logFC > 0), Down = sum(logFC < 0), Total = Up + Down)
carpenter_gene_mat <- final_gene_matrix
carpenter_gene_info <- final_gene_info
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

saveRDS(carpenter_obj, file = paste0(main_dir, "post_processing/results/power_objects/carpenter_obj_4power_1_24_2022.RDS"))

