# PACKAGES ----------
# Data Manipulation
library(tidyverse)
library(reshape2)

# Modeling & Differential Expression
library(limma)
library(edgeR)
library(variancePartition)

# Plots
library(ggpubr)
library(pals)
library(RColorBrewer)
library(wesanderson)
library(svglite)

# Other
library(doParallel)



# BEFORE USING THIS SCRIPT ----------
# User should alter these variables as necessary
main_dir <- "C:/Annika/GitHub Repositories/RodentAddiction/"
pheno_dir <- paste0(main_dir, "powell/00_misc_and_scripts/miscellaneous/")
gene_dir <- paste0(main_dir, "powell/00_misc_and_scripts/gene_counts/")
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
# Sample order is always numerical from low (SRS6085261) to high (SRS6085272).

# Load in phenotype file
pheno <- read_csv(paste0(pheno_dir, "powell_pheno.csv", sep = "")) %>%
  filter(Housing_Condition == "Isolation")
sample_ids <- pheno %>% 
  filter(Housing_Condition == "Isolation") %>% 
  pull(Sample_Treatment) # Vector of sample IDs

# Get list of files, then load in all files
g_files <- list.files(path = gene_dir, 
                      pattern = "*_genes.txt", full.names = TRUE)

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
sample_list <- list("SRS6085261_CI1cue" = SRS6085261_CI1cue,
                    "SRS6085262_CI21cue" = SRS6085262_CI21cue,
                    "SRS6085267_CI1cue" = SRS6085267_CI1cue,
                    "SRS6085268_CI1cue" = SRS6085268_CI1cue,
                    "SRS6085269_CI21cue" = SRS6085269_CI21cue,
                    "SRS6085272_CI21cue" = SRS6085272_CI21cue) %>%
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
# Search for duplicate genes - No duplicates
g_original %>% group_by(Gene_ID) %>% count() %>% arrange(desc(n)) %>% head()



## FILTER LOWLY-EXPRESSED GENES (FPKM) ----------
# To create DGE object, split matrix from gene info and raw counts
gene_matrix <- g_original %>% select(starts_with("SRS"))
gene_info <- g_original %>% select(!starts_with("SRS"))
gene_dge <- DGEList(counts = gene_matrix, genes = gene_info)

# Add phenotype info - treatment group
gene_dge$samples$treatment <- factor(pheno$Treatment_Group, 
                                     levels = c("CI1cue", "CI21cue"))
gene_dge$samples$lane <- factor(pheno$Lane, levels = c("25", "26"))

# Check number of genes
dim(gene_dge) # 30,560 genes, 6 samples/columns

# How many genes aren't expressed in any samples?
table(rowSums(gene_dge$counts==0)==6) # TRUE (unexpressed) = 9,426

# Convert raw counts to FPKM using edgeR and create new dataframe
norm_dge <- calcNormFactors(gene_dge, method = "TMM")
fpkm <- rpkm(norm_dge, gene.length = norm_dge$genes$Length)
fpkm_df <- cbind(gene_info, fpkm)

# Verify that the dimensions of the new dataframe are the same as the old one
all.equal(dim(fpkm_df), dim(g_original)) # TRUE

# Get mean FPKM per treatment per gene
CI1cue_means <- as.data.frame(fpkm) %>%
  dplyr::select(ends_with("CI1cue")) %>%
  apply(1, mean, na.rm = TRUE) %>%
  mutate(as.data.frame(fpkm_df), CI1cue_mean = .)

CI21cue_means <- as.data.frame(fpkm) %>%
  dplyr::select(ends_with("CI21cue")) %>%
  apply(1, mean, na.rm = TRUE) %>%
  mutate(as.data.frame(fpkm_df), CI21cue_mean = .)

# Create single data frame of FPKM means
fpkm_means <- full_join(CI1cue_means, CI21cue_means)

# Check that row number is still the same
dim(fpkm_means)

# Only keep genes that are > 0.5 FPKM in 1 treatment group
keep_genes_fpkm <- fpkm_means %>% 
  filter(CI1cue_mean > 0.5 | CI21cue_mean > 0.5) %>%
  pull(Gene_ID)

# Must also have a read count of at least 6 for 3 (min. sample size) samples
keep_genes_count <- g_original %>%
  mutate(across(starts_with("SRS"), ~ ifelse(.x >=6, TRUE, FALSE))) %>%
  group_by(Gene_ID) %>%
  mutate(Qualifying_Samples = sum(across(starts_with("SRS")))) %>%
  ungroup() %>%
  filter(Qualifying_Samples >= 3) %>%
  pull(Gene_ID)

# Remove genes from count dataframe & verify that this has been done correctly
g_filtered <- g_original %>% filter(Gene_ID %in% keep_genes_fpkm &
                                      Gene_ID %in% keep_genes_count)

# Number of genes after filtering by FPKM
nrow(g_filtered) # 7,224 genes



## WRITE GENE INFO & COUNT FILES ----------
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
                                              levels = c("CI1cue", "CI21cue"))
filtered_gene_dge$samples$lane <- factor(pheno$Lane, levels = c("25", "26"))

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



## MULTI-DIMENSIONAL SCALING PLOTS (PART I, NO OUTLIERS REMOVED) ----------
# Use MDS plots to check for outliers

# Get treatment variable
treatment <- filtered_gene_dge$samples$treatment %>%
  gsub("cue", "", .) %>%
  as.factor()
lane <- filtered_gene_dge$samples$lane %>% as.factor()

# Get and set colors for variable
get_set_color <- function(var, brewer.pal = "Set1") {
  levels(var) <- brewer.pal(nlevels(var), brewer.pal)
  colors <- as.character(var)
  colors
}

color_treatment1 <- get_set_color(treatment)
color_lane1 <- get_set_color(lane)

# Simplify sample names
simple_sample <- rownames(filtered_gene_dge$samples) %>%
  gsub("_.{4,7}", "", .) %>%
  gsub("SRS.{5}", "SRS-", .)

# MDS plots
par(mfrow = c(1, 3))

# PCs 1 & 2
# Outlier: SRS6085267_CI1cue and SRS6085269_CI21cue
plotMDS(filtered_gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment1, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(filtered_gene_lcpm, main = "Lane", label = lane, cex = 1.5, 
        col = color_lane1, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(filtered_gene_lcpm, main = "Samples", label = simple_sample, cex = 1.5, 
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)

# PCs 2 & 3
plotMDS(filtered_gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment1, top = 100, gene.selection = "common", 
        dim.plot = c(2, 3), plot = TRUE)
plotMDS(filtered_gene_lcpm, main = "Lane", label = lane, cex = 1.5, 
        col = color_lane1, top = 100, gene.selection = "common", 
        dim.plot = c(2, 3), plot = TRUE)
plotMDS(filtered_gene_lcpm, main = "Samples", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(2, 3), plot = TRUE)



## REMOVE OUTLIER & RE-FILTER GENES ----------
# Remove outlier samples
g_original2 <- g_original %>% select(-SRS6085267_CI1cue, -SRS6085269_CI21cue)
fpkm2 <- as.data.frame(fpkm) %>% select(-SRS6085267_CI1cue, -SRS6085269_CI21cue)
fpkm_df2 <- fpkm_df %>% select(-SRS6085267_CI1cue, -SRS6085269_CI21cue)
pheno2 <- pheno %>% filter(Sample_Treatment != "SRS6085267_CI1cue" & 
                             Sample_Treatment != "SRS6085269_CI21cue")

# To create DGE object, split matrix from gene info and raw counts
gene_matrix2 <- g_original2 %>% select(starts_with("SRS"))
gene_info2 <- g_original2 %>% select(!starts_with("SRS"))

# How many genes aren't expressed in any samples?
gene_dge2 <- DGEList(counts = gene_matrix2, genes = gene_info)
table(rowSums(gene_dge2$counts==0)==4)
# TRUE (unexpressed) = 9,949

# Check that all are filtered
(ncol(g_original2) - 6) == nrow(pheno2) # Must be TRUE
(ncol(fpkm_df2) - 6) == nrow(pheno2) # Must be TRUE
ncol(fpkm2) == nrow(pheno2) # Must be TRUE

# Number of genes before filtering
dim(fpkm_df2) # 30,560 genes, 11 columns (5 samples)

# Get mean FPKM per treatment per gene
CI1cue_means2 <- as.data.frame(fpkm2) %>%
  dplyr::select(ends_with("CI1cue")) %>%
  apply(1, mean, na.rm = TRUE) %>%
  mutate(as.data.frame(fpkm_df2), CI1cue_mean = .)

CI21cue_means2 <- as.data.frame(fpkm2) %>%
  dplyr::select(ends_with("CI21cue")) %>%
  apply(1, mean, na.rm = TRUE) %>%
  mutate(as.data.frame(fpkm_df2), CI21cue_mean = .)

# Create single data frame of FPKM means
fpkm_means2 <- full_join(CI1cue_means2, CI21cue_means2)

# Check that row number is still the same
dim(fpkm_means2)

# Only keep genes that are > 0.5 FPKM in 1 treatment group
keep_genes_fpkm2 <- fpkm_means2 %>% 
  filter(CI1cue_mean > 0.5 | CI21cue_mean > 0.5) %>%
  pull(Gene_ID)

# Must also have a read count of at least 6 for 2 (min. sample size) samples
keep_genes_count2 <- g_original2 %>%
  mutate(across(starts_with("SRS"), ~ ifelse(.x >=6, TRUE, FALSE))) %>%
  group_by(Gene_ID) %>%
  mutate(Qualifying_Samples = sum(across(starts_with("SRS")))) %>%
  ungroup() %>%
  filter(Qualifying_Samples >= 2) %>%
  pull(Gene_ID)

# Remove genes from count dataframe & verify that this has been done correctly
g_filtered2 <- g_original2 %>% filter(Gene_ID %in% keep_genes_fpkm2 &
                                      Gene_ID %in% keep_genes_count2)

# Number of genes after filtering by FPKM
nrow(g_filtered2) # 7,104 genes


## WRITE GENE INFO & COUNT FILES ----------
# To create DGE object, split matrix from gene info and raw counts
final_gene_matrix <- g_filtered2 %>% select(starts_with("SRS"))
final_gene_info <- g_filtered2 %>% select(!starts_with("SRS"))


## REMAKE DGE OBJECT & MDS PLOTS WITHOUT OUTLIER ----------
# Create DGE List object using filtered counts (not FPKM)
final_gene_dge <- DGEList(counts = final_gene_matrix, 
                          genes = final_gene_info)

# Add pheno info to samples
# Treatment_Group = Housing + Abstinence_Length
final_gene_dge$samples$treatment <- factor(pheno2$Treatment_Group,
                                           levels = c("CI1cue", "CI21cue"))
final_gene_dge$samples$lane <- factor(pheno2$Lane, levels = c("25", "26"))

# Convert to log2 CPM (not normalized)
final_gene_lcpm <- cpm(final_gene_dge, log = TRUE)

# Get variables
treatment <- final_gene_dge$samples$treatment %>%
  gsub("cue", "", .) %>%
  as.factor()
lane <- final_gene_dge$samples$lane

# Set colors for variables
color_treatment <- get_set_color(fct_rev(treatment))
color_lane <- get_set_color(lane)

# Simplify sample names
simple_sample <- rownames(final_gene_dge$samples) %>%
  gsub("_.{4,7}", "", .) %>%
  gsub("SRS.{5}", "SRS-", .)

# MDS plots by each variable
par(mfrow = c(1, 3))

# PCs 1 & 2
plotMDS(final_gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(final_gene_lcpm, main = "Sequencing Lane", label = lane, cex = 1.5, 
        col = color_lane, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(final_gene_lcpm, main = "Sample", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)

# PCs 2 & 3
plotMDS(final_gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(2, 3), plot = TRUE)
plotMDS(final_gene_lcpm, main = "Sequencing Lane", label = lane, cex = 1.5, 
        col = color_lane, top = 100, gene.selection = "common", 
        dim.plot = c(2, 3), plot = TRUE)
plotMDS(final_gene_lcpm, main = "Sample", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)



## VOOM TRANSFORMATIONS & VARIANCE PARTITION ----------
# Get TMM-normalization factors before using voom
final_gene_dge <- calcNormFactors(final_gene_dge, normalize = "TMM")

# Add pheno info to samples
final_gene_dge$samples$treatment <- factor(pheno2$Treatment_Group,
                                           levels = c("CI1cue", "CI21cue"))
final_gene_dge$samples$lane <- factor(pheno2$Lane, levels = c("25", "26"))

# Re-define group variables
treatment <- final_gene_dge$samples$treatment
lane <- final_gene_dge$samples$lane

# Test several designs/models to see how variances partition based on design
# Treatment (CI1_cue, CI21_cue). Note that Treatment = Abstinence Length.
# To control for batch effects, we may need to use the Lane variable.

# Treatment
des_treat <- model.matrix(~ 0 + treatment)
model.matrix(~ treatment + 0)
colnames(des_treat) <- gsub("treatment", "", colnames(des_treat))

# Treatment + Lane
des_treat_lane <- model.matrix(~ 0 + treatment + lane)
colnames(des_treat_lane) <- gsub("treatment", "", colnames(des_treat_lane))
colnames(des_treat_lane) <- gsub("lane", "L", colnames(des_treat_lane))

# Voom transformations
v_treat <- voom(final_gene_dge, des_treat, plot = TRUE)
v_treat_lane <- voom(final_gene_dge, des_treat_lane, plot = TRUE)

# MDS Plots based on voom models - warnings with very low sample size
par(mfrow = c(2, 3))

# Treatment
plotMDS(v_treat, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common",
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(v_treat, main = "Sequencing Lane", label = lane, cex = 1.5,
        col = color_lane, top = 100, gene.selection = "common",
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(v_treat, main = "Sample", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(v_treat, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common",
        dim.plot = c(3, 4), plot = TRUE)
plotMDS(v_treat, main = "Sequencing Lane", label = lane, cex = 1.5,
        col = color_lane, top = 100, gene.selection = "common",
        dim.plot = c(3, 4), plot = TRUE)
plotMDS(v_treat, main = "Sample", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)

# Treatment + Lane
plotMDS(v_treat_lane, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common",
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(v_treat_lane, main = "Sequencing Lane", label = lane, cex = 1.5,
        col = color_lane, top = 100, gene.selection = "common",
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(v_treat_lane, main = "Sample", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(v_treat_lane, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common",
        dim.plot = c(3, 4), plot = TRUE)
plotMDS(v_treat_lane, main = "Sequencing Lane", label = lane, cex = 1.5,
        col = color_lane, top = 100, gene.selection = "common",
        dim.plot = c(3, 4), plot = TRUE)
plotMDS(v_treat_lane, main = "Sample", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)

# Treatment + Lane, batch corrected for Lane
vCorrectLane <- removeBatchEffect(v_treat_lane, batch = v_treat_lane$targets$lane)
plotMDS(vCorrectLane, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common",
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(vCorrectLane, main = "Sequencing Lane", label = lane, cex = 1.5,
        col = color_lane, top = 100, gene.selection = "common",
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(vCorrectLane, main = "Sample", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(vCorrectLane, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common",
        dim.plot = c(3, 4), plot = TRUE)
plotMDS(vCorrectLane, main = "Sequencing Lane", label = lane, cex = 1.5,
        col = color_lane, top = 100, gene.selection = "common",
        dim.plot = c(3, 4), plot = TRUE)
plotMDS(vCorrectLane, main = "Sample", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(2, 3), plot = TRUE)


# Final plots
# Form for variance partition; include all variables
form <- ~ (1|treatment) + (1|lane)
# Only run fitExtractVarPartModel() if you're prepared for a time-consuming analysis
cl <- makeCluster(2)
registerDoParallel(cl)

# Palette
palette <- c(wes_palette("Chevalier1", 4, type = "discrete"), "grey85")[-c(2,3)]

# Treatment model
#vp_treat <- fitExtractVarPartModel(v_treat, form, final_gene_dge$samples)
#saveRDS(vp_treat, file = paste0(main_dir, "post_processing/results/variance_partition/powell/powell_FINAL_treat.RDS")
vp_treat <- readRDS(file = paste0(main_dir, "post_processing/results/variance_partition/powell/powell_FINAL_treat.RDS"))
vp_treat <- sortCols(vp_treat)
bars_treat <- plotPercentBars(vp_treat[1:10, ]) +
  scale_fill_manual(values = palette,
                    labels = c("Lane", "Treatment", "Residuals")) +
  theme(legend.position = "bottom")
varpart_treat <- plotVarPart(vp_treat) +
  scale_x_discrete(labels = c("Lane", "Treatment", "Residuals")) +
  scale_fill_manual(values = palette) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
ggarrange(varpart_treat, bars_treat)

# Treatment + Lane model
#vp_treat_lane <- fitExtractVarPartModel(v_treat_lane, form, final_gene_dge$samples)
#saveRDS(vp_treat_lane, file = paste0(main_dir, "post_processing/results/variance_partition/powell/powell_FINAL_treat_lane.RDS"))
vp_treat_lane <- readRDS(file = paste0(main_dir, "post_processing/results/variance_partition/powell/powell_FINAL_treat_lane.RDS"))
vp_treat_lane <- sortCols(vp_treat_lane)
bars_lane <- plotPercentBars(vp_treat_lane[1:10, ]) +
  scale_fill_manual(values = palette,
                    labels = c("Lane", "Treatment", "Residuals")) +
  theme(legend.position = "bottom")
varpart_lane <- plotVarPart(vp_treat_lane) +
  scale_x_discrete(labels = c("Lane", "Treatment", "Residuals")) +
  scale_fill_manual(values = palette) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
ggarrange(varpart_lane, bars_lane)

# Lane is accounting for a significant amount of variation in the model.
# Interpretation: the variance explained by each variables after correcting for 
# all other variables
stopCluster(cl)


## DIFFERENTIAL EXPRESSION ANALYSIS ---
# Create contrasts matrix
treat_cm <- makeContrasts(Incubation = CI1cue - CI21cue, 
                          levels = colnames(des_treat))
treat_lane_cm <- makeContrasts(Incubation = CI1cue - CI21cue,
                               levels = colnames(des_treat_lane))

# Fit the models and contrasts, then run contrasts
treat_fit <- lmFit(v_treat, des_treat)
treat_fc <- contrasts.fit(treat_fit, contrasts = treat_cm)

treat_lane_fit <- lmFit(v_treat_lane, des_treat_lane)
treat_lane_fc <- contrasts.fit(treat_lane_fit, contrasts = treat_lane_cm)

# Calculate the t-statistics for the contrasts
treat_ebayes <- eBayes(treat_fc)
treat_lane_ebayes <- eBayes(treat_lane_fc)

# SA Plot
par(mfrow = c(1, 2))

v_treat <- voom(final_gene_dge, des_treat, plot = TRUE)
plotSA(treat_ebayes, main = "Final model: mean-variance trend")

v_treat_lane <- voom(final_gene_dge, des_treat_lane, plot = TRUE)
plotSA(treat_lane_ebayes, main = "Final model: mean-variance trend")


# Write files
# Get stats for all genes, in same order as they were originally
treat_stats <- topTable(treat_ebayes, number = nrow(treat_ebayes))
treat_lane_stats <- topTable(treat_lane_ebayes, number = nrow(treat_lane_ebayes))

sig_treat <- treat_stats %>% filter(P.Value < 0.05) %>% pull(Gene_ID) # 244 genes
sig_treat_lane <- treat_lane_stats %>% filter(P.Value < 0.05) %>% pull(Gene_ID) # 286 genes

sig_treat[sig_treat %in% sig_treat_lane] %>% length() # Only 215 shared

treat_lane_stats %>%
  filter(P.Value < 0.05) %>%
  summarize(Up = sum(logFC > 0), Down = sum(logFC < 0))

# Write files - can be used for getting overlaps (UpSetR plots)
# write.table(treat_stats, file = paste0(main_dir, "post_processing/results/gene_info/powell_DE_FINAL_treat.txt"))
# write.table(treat_lane_stats, file = paste0(main_dir, "post_processing/results/gene_info/powell_DE_FINAL_treat_lane.txt"))



## CREATE OBJECT FOR POWER ANALYSIS ----------
powell_dge <- final_gene_dge # Already calculated norm. factors for TMM method
powell_stats <- treat_lane_stats
powell_summary <- powell_stats %>%
  filter(P.Value < 0.05) %>%
  summarize(Up = sum(logFC > 0), Down = sum(logFC < 0), Total = Up + Down)
powell_gene_mat <- final_gene_matrix
powell_gene_info <- final_gene_info
powell_fc <- powell_stats %>%
  filter(P.Value < 0.05) %>%
  mutate(Abs_FC = 2^abs(logFC)) %>%
  dplyr::select(Abs_FC) %>%
  summary()
powell_pheno <- pheno2
powell_treat <- fct_relevel(pheno2$Treatment_Group, levels = c("CI1cue", "CI21cue"))
powell_control <- levels(powell_treat)[1]

powell_obj <- list(DGE = powell_dge, Stats = powell_stats, 
                   Stats_Summary = powell_summary, 
                   Gene_Matrix = powell_gene_mat, Gene_Info = powell_gene_info, 
                   FC = powell_fc, Pheno = powell_pheno, Treat = powell_treat, 
                   Control = powell_control)

saveRDS(powell_obj, file = paste0(main_dir, "post_processing/results/power_objects/powell_obj_4power.RDS"))


## FINAL FIGURES ----------
points <- color_lane %>% 
  str_replace_all("#E41A1C", "16") %>%
  str_replace_all("#377EB8", "18") %>%
  as.numeric()

# Final figure
# svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Powell_RNAseq_1_24_2022.svg", width = 5.4, height = 8)
par(mfrow = c(3, 2))
plotMDS(filtered_gene_lcpm,
        pch = points, cex = 3, col = color_treatment1, cex.axis = 1.25, cex.lab = 1.15,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(filtered_gene_lcpm,
        pch = points, cex = 3, col = color_treatment1, cex.axis = 1.25, cex.lab = 1.15,
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)
# Original voom data
voom_info <- voom(final_gene_dge, des_treat_lane, plot = FALSE, save.plot = TRUE)
plot(voom_info$voom.xy, pch = 19, cex = 0.5, cex.main = 1.25, cex.axis = 1.25, 
     cex.lab = 1.25, main = "voom: Mean-variance trend", 
     xlab = "log2(count size + 0.5)", ylab = "Sqrt(standard deviation)")
lines(voom_info$voom.line, col = "red")
# Final voom model
plotSA(treat_lane_ebayes, main = "Final model: Mean-variance trend", 
       cex = 0.5, pch = 19, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.25)
plotMDS(v_treat_lane, pch = points,  cex = 3, cex.axis = 1.25, cex.lab = 1.15,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(v_treat_lane, pch = points,  cex = 3, cex.axis = 1.25, cex.lab = 1.15,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)
# dev.off()


# Final figure for presentation
# svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Powell_RNAseq_simple_presentation_1_24_2022.svg", width = 4.5, height = 4.5)
par(mfrow = c(1, 1))
plotMDS(filtered_gene_lcpm,
        pch = 16, cex = 3, col = color_treatment, cex.axis = 1.5, cex.lab = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
# dev.off()

# svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Powell_varPart_1_24_2022.svg", width = 6.9, height = 2.6)
bars_lane <- plotPercentBars(vp_treat[1:10, ]) +
  scale_fill_manual(values = palette,
                    labels = c("Lane", "Treatment", "Residuals")) +
  theme(legend.position = "none")
varpart_lane <- plotVarPart(vp_treat_lane) +
  scale_x_discrete(labels = c("Lane", "Treatment", "Residuals")) +
  scale_fill_manual(values = palette, 
                    labels = c("Lane", "Treatment", "Residuals")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())
ggarrange(varpart_lane, bars_lane, common.legend = TRUE, legend = "right")
# dev.off()






