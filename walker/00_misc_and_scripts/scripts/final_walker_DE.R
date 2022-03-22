
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
pheno_dir <- paste0(main_dir, "walker/00_misc_and_scripts/miscellaneous/")
gene_dir <- paste0(main_dir, "walker/00_misc_and_scripts/gene_counts/")
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
# Sample order is always numerical from low (SRS2926555) to high (SRS2926592).

# Load in phenotype file
pheno <- read_csv(paste0(pheno_dir, "walker_pheno.csv")) %>%
  filter(Abstinence_Length == "30")
sample_ids <- pheno %>% 
  filter(Abstinence_Length == "30") %>% 
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
sample_list <- list("SRS2926555_C30sal" = SRS2926555_C30sal, 
                    "SRS2926557_S30sal" = SRS2926557_S30sal, 
                    "SRS2926558_S30sal" = SRS2926558_S30sal, 
                    "SRS2926560_S30sal" = SRS2926560_S30sal,
                    "SRS2926561_C30sal" = SRS2926561_C30sal, 
                    "SRS2926565_C30sal" = SRS2926565_C30sal,
                    "SRS2926566_S30sal" = SRS2926566_S30sal, 
                    "SRS2926568_C30sal" = SRS2926568_C30sal,
                    "SRS2926571_S30sal" = SRS2926571_S30sal, 
                    "SRS2926574_C30sal" = SRS2926574_C30sal,
                    "SRS2926576_S30sal" = SRS2926576_S30sal) %>%
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
                                     levels = c("S30sal", "C30sal"))

# Check number of genes
dim(gene_dge) # 55,414 genes, 11 samples/columns

# How many genes aren't expressed in any samples?
table(rowSums(gene_dge$counts==0)==11)
# TRUE (unexpressed) = 26,879

# Convert raw counts to FPKM using edgeR and create new dataframe
norm_dge <- calcNormFactors(gene_dge, method = "TMM")
fpkm <- rpkm(norm_dge, gene.length = norm_dge$genes$Length)
fpkm_df <- cbind(gene_info, fpkm)

# Verify that the dimensions of the new dataframe are the same as the old one
all.equal(dim(fpkm_df), dim(g_original)) # TRUE

# Get mean FPKM per treatment per gene
S30sal_means <- as.data.frame(fpkm) %>%
  dplyr::select(ends_with("S30sal")) %>%
  apply(1, mean, na.rm = TRUE) %>%
  mutate(as.data.frame(fpkm_df), S30sal_mean = .)

C30sal_means <- as.data.frame(fpkm) %>%
  dplyr::select(ends_with("C30sal")) %>%
  apply(1, mean, na.rm = TRUE) %>%
  mutate(as.data.frame(fpkm_df), C30sal_mean = .)

# Create single data frame of FPKM means
fpkm_means <- full_join(S30sal_means, C30sal_means)

# Check that row number is still the same
dim(fpkm_means)

# Only keep genes that are > 0.5 FPKM in 1 treatment group
keep_genes_fpkm <- fpkm_means %>% 
  filter(S30sal_mean > 0.5 | C30sal_mean > 0.5) %>%
  pull(Gene_ID)

# Must also have a read count of at least 6 for 2 (min. sample size) samples
keep_genes_count <- g_original %>%
  mutate(across(starts_with("SRS"), ~ ifelse(.x >=6, TRUE, FALSE))) %>%
  group_by(Gene_ID) %>%
  mutate(Qualifying_Samples = sum(across(starts_with("SRS")))) %>%
  ungroup() %>%
  filter(Qualifying_Samples >= 5) %>%
  pull(Gene_ID)

# Remove genes from count dataframe & verify that this has been done correctly
g_filtered <- g_original %>% filter(Gene_ID %in% keep_genes_fpkm &
                                      Gene_ID %in% keep_genes_count)

# Number of genes after filtering by FPKM
nrow(g_filtered) # 9,948 genes



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
                                              levels = c("S30sal", "C30sal"))

# # Convert to log2 CPM (not normalized)
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


## MULTI-DIMENSIONAL SCALING PLOTS ----------
# Use MDS plots to check for outliers

# Get treatment variable
treatment <- filtered_gene_dge$samples$treatment

# Get and set colors for variable
get_set_color <- function(var, brewer.pal = "Set1") {
  levels(var) <- brewer.pal(nlevels(var), brewer.pal)
  colors <- as.character(var)
  colors
}

color_treatment <- get_set_color(treatment)

# Simplify sample names
simple_sample <- rownames(filtered_gene_dge$samples) %>%
  gsub("_.{4,6}", "", .) %>%
  gsub("SRS.{5}", "SRS-", .)

# MDS plots
par(mfrow = c(1, 2))

# PCs 1 & 2
# Outlier: SRS2926555_C30sal?
plotMDS(filtered_gene_lcpm, main = "Treatment Group", label = treatment, 
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(filtered_gene_lcpm, main = "Samples", label = simple_sample, top = 100, 
        gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)

# Check on 3rd and 4th PCs - Doesn't look like an outlier
plotMDS(filtered_gene_lcpm, main = "Treatment Group", label = treatment, 
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)
plotMDS(filtered_gene_lcpm, main = "Samples", label = simple_sample, top = 100, 
        gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)

# Check again using top 50 genes instead of 100 - Looks less like an outlier
plotMDS(filtered_gene_lcpm, main = "Treatment Group", label = treatment, 
        col = color_treatment, top = 50, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(filtered_gene_lcpm, main = "Samples", label = simple_sample, top = 50, 
        gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)


## WRITE GENE INFO & COUNT FILES ----------
# To create DGE object, split matrix from gene info and raw counts
final_gene_matrix <- g_filtered %>% select(starts_with("SRS"))
final_gene_info <- g_filtered %>% select(!starts_with("SRS"))
final_gene_dge <- DGEList(counts = final_gene_matrix, 
                          genes = final_gene_info)


## REMAKE DGE OBJECT & MDS PLOTS WITHOUT OUTLIER ----------
# Create DGE List object using filtered counts (not FPKM)
final_gene_dge <- DGEList(counts = final_gene_matrix, 
                          genes = final_gene_info)

# Add pheno info to samples
final_gene_dge$samples$treatment <- factor(pheno$Treatment_Group,
                                           levels = c("S30sal", "C30sal"))
final_gene_dge$samples$batch <- factor(pheno$Lab_Run_ID,
                                       levels = c("Run_714", "Run_715", 
                                                  "Run_716", "Run_717"))
final_gene_dge$samples$instrument <- factor(pheno$Instrument_Name,
                                            levels = c("HISEQ", "HWI-D00743",
                                                       "HWI-D00381"))

# Convert to log2 CPM (not normalized)
final_gene_lcpm <- cpm(final_gene_dge, log = TRUE)

# Get variables
treatment <- final_gene_dge$samples$treatment 
batch <- final_gene_dge$samples$batch
instrument <- final_gene_dge$samples$instrument

# Set colors for variables
color_treatment <- get_set_color(fct_rev(treatment))
color_batch <- get_set_color(batch)
color_instrument <- get_set_color(instrument)

# Simplify sample names
simple_sample <- rownames(final_gene_dge$samples) %>%
  gsub("_.{4,6}", "", .) %>%
  gsub("SRS.{5}", "SRS-", .)

# Simplify Batch and Instrument names
simple_batch <- gsub("Run_", "", final_gene_dge$samples$batch)
simple_instrument <- final_gene_dge$samples$instrument %>%
  gsub("D00", "", .) %>%
  strtrim(7)

# MDS plots
par(mfrow = c(2, 2))

# PCs 1 & 2
plotMDS(final_gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(final_gene_lcpm, main = "Samples", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(final_gene_lcpm, main = "Batch", label = simple_batch, cex = 1.5,
        col = color_batch, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(final_gene_lcpm, main = "Sequencing Instrument", label = simple_instrument, 
        col = color_instrument, top = 100, gene.selection = "common", cex = 1.5,
        dim.plot = c(1, 2), plot = TRUE)

# PCs 2 & 3
plotMDS(final_gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)
plotMDS(final_gene_lcpm, main = "Samples", label = simple_sample, cex = 1.5, 
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)
plotMDS(final_gene_lcpm, main = "Batch", label = simple_batch, cex = 1.5,
        col = color_batch, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)
plotMDS(final_gene_lcpm, main = "Sequencing Instrument", label = simple_instrument, 
        cex = 1.5, col = color_instrument, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)



## VOOM TRANSFORMATIONS & VARIANCE PARTITION ----------
# Get TMM-normalization factors before using voom
final_gene_dge <- calcNormFactors(final_gene_dge, normalize = "TMM")

# Add pheno info to samples
final_gene_dge$samples$treatment <- factor(pheno$Treatment_Group,
                                           levels = c("S30sal", "C30sal"))
final_gene_dge$samples$batch <- factor(pheno$Lab_Run_ID,
                                       levels = c("Run_714", "Run_715", 
                                                  "Run_716", "Run_717"))
final_gene_dge$samples$instrument <- factor(pheno$Instrument_Name,
                                            levels = c("HISEQ", "HWI-D00743",
                                                       "HWI-D00381"))

# Re-define group variables
treatment <- final_gene_dge$samples$treatment
batch <- final_gene_dge$samples$batch
instrument <- final_gene_dge$samples$instrument

# Test several designs/models to see how variances partition based on design
# Treatment (S30sal, C30sal). Note that Treatment = Drug in this case, since
# both abstinence lengths are the same.
# To control for batch effects, we need to include either Batch or Instrument.

# Treatment
des_treat <- model.matrix(~ 0 + treatment)
colnames(des_treat) <- gsub("treatment", "", colnames(des_treat))

# Treatment + Batch
des_treat_batch <- model.matrix(~ 0 + treatment + batch)
colnames(des_treat_batch) <- gsub("treatment", "", colnames(des_treat_batch))

# Treatment + Instrument
des_treat_ins <- model.matrix(~ 0 + treatment + instrument)
colnames(des_treat_ins) <- gsub("treatment", "", colnames(des_treat_ins))

# Treatment + Batch + Instrument
des_treat_batch_ins <- model.matrix(~ 0 + treatment + batch + instrument)
colnames(des_treat_batch_ins) <- gsub("treatment", "", colnames(des_treat_batch_ins))

# Voom transformations
par(mfrow = c(1, 1))
v_treat <- voom(final_gene_dge, des_treat, plot = TRUE)
v_treat_batch <- voom(final_gene_dge, des_treat_batch, plot = TRUE)
v_treat_ins <- voom(final_gene_dge, des_treat_ins, plot = TRUE)
v_treat_batch_ins <- voom(final_gene_dge, des_treat_batch_ins, plot = TRUE)

# Form for variance partition; include all variables
form <- ~ (1|treatment) + (1|batch) + (1|instrument)
# Only run fitExtractVarPartModel() if you're prepared for a time-consuming analysis
cl <- makeCluster(2)
registerDoParallel(cl)

# Palette
palette <- c(wes_palette("Chevalier1", 4, type = "discrete")[-3], "grey85")

# Treatment model
# More variance is explained by Instrument than Batch in this model
#vp_treat <- fitExtractVarPartModel(v_treat, form, final_gene_dge$samples)
#saveRDS(vp_treat, file = paste0(main_dir, "post_processing/results/variance_partition/walker_FINAL_treat.RDS"))
vp_treat <- readRDS(file = paste0(main_dir, "post_processing/results/variance_partition/walker_FINAL_treat.RDS"))
vp_treat <- sortCols(vp_treat)
bars_treat <- plotPercentBars(vp_treat[1:10, ]) +
  scale_fill_manual(values = palette,
                    labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  theme(legend.position = "bottom")
varpart_treat <- plotVarPart(vp_treat) +
  scale_x_discrete(labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  scale_fill_manual(values = palette) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major.x = element_blank())
ggarrange(varpart_treat, bars_treat)

# Treatment + Batch model
# More variance is explained by Batch than Instrument with this model
#vp_treat_batch <- fitExtractVarPartModel(v_treat_batch, form, final_gene_dge$samples)
#saveRDS(vp_treat_batch, file = paste0(main_dir, "post_processing/results/variance_partition/walker_FINAL_treat_batch.RDS"))
vp_treat_batch <- readRDS(file = paste0(main_dir, "post_processing/results/variance_partition/walker_FINAL_treat_batch.RDS"))
vp_treat_batch <- sortCols(vp_treat_batch)
bars_batch <- plotPercentBars(vp_treat_batch[1:10, ]) +
  scale_fill_manual(values = palette,
                    labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  theme(legend.position = "bottom")
varpart_batch <- plotVarPart(vp_treat_batch) +
  scale_x_discrete(labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  scale_fill_manual(values = palette) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major.x = element_blank())
ggarrange(varpart_batch, bars_batch)

# Treatment + Instrument model
# More variance is explained by Instrument than Batch with this model
#vp_treat_ins <- fitExtractVarPartModel(v_treat_ins, form, final_gene_dge$samples)
#saveRDS(vp_treat_ins, file = paste0(main_dir, "post_processing/results/variance_partition/walker_FINAL_treat_ins.RDS"))
vp_treat_ins <- readRDS(file = paste0(main_dir, "post_processing/results/variance_partition/walker_FINAL_treat_ins.RDS"))
vp_treat_ins <- sortCols(vp_treat_ins)
bars_ins <- plotPercentBars(vp_treat_ins[1:10, ]) +
  scale_fill_manual(values = palette,
                    labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  theme(legend.position = "bottom")
varpart_ins <- plotVarPart(vp_treat_ins) +
  scale_x_discrete(labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  scale_fill_manual(values = palette) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major.x = element_blank())
ggarrange(varpart_ins, bars_ins)

# Treatment + Instrument model
# More variance is explained by Instrument than Batch or Treatment with this model
#vp_treat_batch_ins <- fitExtractVarPartModel(v_treat_batch_ins, form, final_gene_dge$samples)
#saveRDS(vp_treat_batch_ins , file = paste0(main_dir, "post_processing/results/variance_partition/walker_FINAL_treat_batch_ins.RDS"))
vp_treat_batch_ins <- readRDS(file = paste0(main_dir, "post_processing/results/variance_partition/walker_FINAL_treat_batch_ins.RDS"))
vp_treat_batch_ins <- sortCols(vp_treat_batch_ins)
bars_batch_ins <- plotPercentBars(vp_treat_batch_ins[1:10, ]) +
  scale_fill_manual(values = palette,
                    labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  theme(legend.position = "bottom")
varpart_batch_ins <- plotVarPart(vp_treat_batch_ins) +
  scale_x_discrete(labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  scale_fill_manual(values = palette) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major.x = element_blank())
ggarrange(varpart_batch_ins, bars_batch_ins)

stopCluster(cl)

# Using both Batch and Instrument in our model could overfit the data (and in
# fact, gives us the same results as just using Batch - see below with voom).
# We will select the Treatment + Batch model because, using this model,
# we are able to account for some of the variation from Batch. Batch is explains
# the most variation on PC2 (other than residuals) - this is the only variable
# on the first 3 PCs.
# In the full model (Treatment + Batch + Instrument), we can see that even when
# the variables are considered independently, there is a larger contribution of
# Batch than the other variables.



## DIFFERENTIAL EXPRESSION ANALYSIS ---
# Create contrast matrices
treat_cm <- makeContrasts(S30sal_C30sal = S30sal - C30sal, 
                          levels = colnames(des_treat))
treat_batch_cm <- makeContrasts(S30sal_C30sal = S30sal - C30sal, 
                                levels = colnames(des_treat_batch))
treat_ins_cm <- makeContrasts(S30sal_C30sal = S30sal - C30sal, 
                              levels = colnames(des_treat_ins))
treat_batch_ins_cm <- makeContrasts(S30sal_C30sal = S30sal - C30sal, 
                                    levels = colnames(des_treat_batch_ins))

# Fit the models and contrasts, then run contrasts
treat_fm <- lmFit(v_treat, des_treat)
treat_fc <- contrasts.fit(treat_fm, contrasts = treat_cm)

treat_batch_fm <- lmFit(v_treat_batch, des_treat_batch)
treat_batch_fc <- contrasts.fit(treat_batch_fm, contrasts = treat_batch_cm)

treat_ins_fm <- lmFit(v_treat_ins, des_treat_ins)
treat_ins_fc <- contrasts.fit(treat_ins_fm, contrasts = treat_ins_cm)

treat_batch_ins_fm <- lmFit(v_treat_batch_ins, des_treat_batch_ins)
treat_batch_ins_fc <- contrasts.fit(treat_batch_ins_fm, contrasts = treat_batch_ins_cm)

# Calculate the t-statistics for the contrasts
# Robustified against outlier sample variances
treat_ebayes <- eBayes(treat_fc)
treat_batch_ebayes <- eBayes(treat_batch_fc)
treat_ins_ebayes <- eBayes(treat_ins_fc)
treat_batch_ins_ebayes <- eBayes(treat_batch_ins_fc)

# SA Plot
par(mfrow = c(1, 2))

v_treat <- voom(final_gene_dge, des_treat, plot = TRUE)
plotSA(treat_ebayes, main = "Final model: mean-variance trend")

v_treat_batch <- voom(final_gene_dge, des_treat_batch, plot = TRUE)
plotSA(treat_batch_ebayes, main = "Final model: mean-variance trend")

v_treat_ins <- voom(final_gene_dge, des_treat_ins, plot = TRUE)
plotSA(treat_ins_ebayes, main = "Final model: mean-variance trend")

v_treat_batch_ins <- voom(final_gene_dge, des_treat_batch_ins, plot = TRUE)
plotSA(treat_batch_ins_ebayes, main = "Final model: mean-variance trend")

# Get stats for all genes, in same order as they were originally
treat_stats <- topTable(treat_ebayes, number = nrow(treat_ebayes))
treat_batch_stats <- topTable(treat_batch_ebayes, number = nrow(treat_batch_ebayes))
treat_ins_stats <- topTable(treat_ins_ebayes, number = nrow(treat_ins_ebayes))
treat_batch_ins_stats <- topTable(treat_batch_ins_ebayes, number = nrow(treat_batch_ins_ebayes))

# Get genes w/ P-value < 0.05
sig_treat <- treat_stats %>% filter(P.Value < 0.05) %>% select(Gene_ID) %>% mutate(Design = "treat")
sig_treat_batch <- treat_batch_stats %>% filter(P.Value < 0.05) %>% select(Gene_ID) %>% mutate(Design = "treat_batch")
sig_treat_ins <- treat_ins_stats %>% filter(P.Value < 0.05) %>% select(Gene_ID) %>% mutate(Design = "treat_ins")
sig_treat_batch_ins <- treat_batch_ins_stats %>% filter(P.Value < 0.05) %>% select(Gene_ID) %>% mutate(Design = "treat_batch_ins")

df <- Reduce(function(x, y) 
  full_join(x, y), list(sig_treat, sig_treat_batch, sig_treat_ins, sig_treat_batch_ins)) %>%
  group_by(Gene_ID) %>%
  mutate(count = length(Gene_ID))

sig_treat <- sig_treat %>% pull(Gene_ID) # 156 genes
sig_treat_batch <- sig_treat_batch %>% pull(Gene_ID) # 139 genes
sig_treat_ins <- sig_treat_ins %>% pull(Gene_ID) # 161 genes
sig_treat_batch_ins <- sig_treat_batch_ins %>% pull(Gene_ID) # 139 genes

sig_treat_batch[sig_treat_batch %in% sig_treat_batch_ins] %>% length() # All the same (139)

sig_treat[sig_treat %in% sig_treat_batch] %>% length() # Mostly the same (110/122 to 139)
sig_treat_ins[sig_treat_ins %in% sig_treat_batch] %>% length() # Mostly the same (126/161 to 139)
sig_treat[sig_treat %in% sig_treat_ins] %>% length() # Mostly the same (126/156 to 161)

# Write files
write.table(treat_stats, file = paste0(main_dir, "post_processing/results/gene_info/walker_DE_FINAL_treat.txt"))
write.table(treat_batch_stats, file = paste0(main_dir, "post_processing/results/gene_info/walker_DE_FINAL_treat_batch.txt"))
write.table(treat_ins_stats, file = paste0(main_dir, "post_processing/results/gene_info/walker_DE_FINAL_treat_ins.txt"))
write.table(treat_batch_ins_stats, file = paste0(main_dir, "post_processing/results/gene_info/walker_DE_FINAL_treat_batch_ins.txt"))



## MDS PLOTS FOR VOOM MODEL ----------
vCorrectBatch <- removeBatchEffect(v_treat_batch, batch = v_treat_batch$targets$batch)
plotMDS(v_treat_batch, main = "Treatment & Batch", label = treatment, cex = 1.5, col = color_batch,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(vCorrectBatch, main = "Treatment & Batch", label = treatment, cex = 1.5, col = color_batch,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)



## CREATE OBJECT FOR POWER ANALYSIS ----------
walker_dge <- final_gene_dge # Already calculated norm. factors for TMM method
walker_stats <- treat_batch_stats
walker_summary <- walker_stats %>%
  filter(P.Value < 0.05) %>%
  summarize(Up = sum(logFC > 0), Down = sum(logFC < 0), Total = Up + Down)
walker_gene_mat <- final_gene_matrix
walker_gene_info <- final_gene_info
walker_fc <- walker_stats %>%
  filter(P.Value < 0.05) %>%
  mutate(Abs_FC = 2^abs(logFC)) %>%
  dplyr::select(Abs_FC) %>%
  summary()
walker_pheno <- pheno
walker_treat <- fct_relevel(pheno$Treatment_Group, levels = c("S30sal", "C30sal"))
walker_control <- levels(walker_treat)[1]

walker_obj <- list(DGE = walker_dge, Stats = walker_stats, 
                   Stats_Summary = walker_summary, 
                   Gene_Matrix = walker_gene_mat, Gene_Info = walker_gene_info, 
                   FC = walker_fc, Pheno = walker_pheno, Treat = walker_treat, 
                   Control = walker_control)

saveRDS(walker_obj, file = paste0(main_dir, "post_processing/results/variance_partition/walker_obj_4power.RDS"))

points <- color_batch %>% 
  str_replace_all("#4DAF4A", "15") %>%
  str_replace_all("#377EB8", "16") %>%
  str_replace_all("#E41A1C", "17") %>%
  str_replace_all("#984EA3", "18") %>%
  as.numeric()

# Final figure
# svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Walker_RNAseq_1_24_2022.svg", width = 5.4, height = 8)
par(mfrow = c(3, 2))
plotMDS(filtered_gene_lcpm,
        pch = points, cex = 3, col = color_treatment, cex.axis = 1.25, cex.lab = 1.15,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(filtered_gene_lcpm,
        pch = points, cex = 3, col = color_treatment, cex.axis = 1.25, cex.lab = 1.15,
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)
# Original voom data
voom_info <- voom(final_gene_dge, des_treat_batch, plot = FALSE, save.plot = TRUE)
plot(voom_info$voom.xy, pch = 19, cex = 0.5, cex.main = 1.25, cex.axis = 1.25, 
     cex.lab = 1.25, main = "voom: Mean-variance trend", 
     xlab = "log2(count size + 0.5)", ylab = "Sqrt(standard deviation)")
lines(voom_info$voom.line, col = "red")
# Final voom model
plotSA(treat_batch_ebayes, main = "Final model: Mean-variance trend", 
       cex = 0.5, pch = 19, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.25)
plotMDS(v_treat_batch, pch = points,  cex = 3, cex.axis = 1.25, cex.lab = 1.15,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(v_treat_batch, pch = points,  cex = 3, cex.axis = 1.25, cex.lab = 1.15,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)
# dev.off()


# Final figure for presentation
# svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Walker_RNAseq_simple_presentation_1_24_2022.svg", width = 4.5, height = 4.5)
par(mfrow = c(1, 1))
plotMDS(final_gene_lcpm,
        pch = 16, cex = 3, col = color_treatment, cex.axis = 1.5, cex.lab = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
# dev.off()


# svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Walker_varPart_1_24_2022.svg", width = 6.9, height = 2.6)
bars_batch <- plotPercentBars(vp_treat[1:10, ]) +
  scale_fill_manual(values = palette,
                    labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  theme(legend.position = "bottom")
varpart_batch <- plotVarPart(vp_treat_batch) +
  scale_x_discrete(labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  scale_fill_manual(values = palette, 
                    labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank())
ggarrange(varpart_batch, bars_batch, common.legend = TRUE, legend = "right")
# dev.off()




