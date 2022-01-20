
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
main_dir <- "G:/Noctsol/GitHub Repositories/temp/"
pheno_dir <- paste0(main_dir, "walker/")
gene_dir <- paste0(main_dir, "walker/genes/")
setwd(paste0(main_dir, "results/"))



# LOAD IN DATA & RENAME COLUMNS ----------
# NOTE: Sample order must be maintained for accuracy.
# Sample order is always numerical from low (SRS2926555) to high (SRS2926592).

# Load in phenotype file
pheno <- read_csv(paste0(pheno_dir, "walker_pheno.csv"))
sample_ids <- pheno$Sample_Treatment # Vector of sample IDs

# Load in gene and transcript counts files
g_count <- read_csv(paste0(gene_dir, "walker_gene_count_matrix.csv"))

# Separate initial gene_id column into Gene.ID & Gene.Name
# E.g. ENSMUSG00000102348|Gm10568 becomes ENSMUSG00000102348 & Gm10568
fix_gene_ids <- function(df) {separate(df, col = "gene_id", 
                                       into = c("Gene.ID", "Gene.Name"), 
                                       sep = "\\|")}

# Fix gene IDs and remove samples that aren't of interest
g_count <- fix_gene_ids(g_count) %>%
  select(-contains(c("_S1", "_C1")))

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
                    "SRS2926576_S30sal" = SRS2926576_S30sal, 
                    "SRS2926577_C1no" = SRS2926577_C1no,
                    "SRS2926578_C1no" = SRS2926578_C1no, 
                    "SRS2926579_S1no" = SRS2926579_S1no,
                    "SRS2926580_C1no" = SRS2926580_C1no, 
                    "SRS2926581_S1no" = SRS2926581_S1no,
                    "SRS2926582_S1no" = SRS2926582_S1no, 
                    "SRS2926583_S1no" = SRS2926583_S1no,
                    "SRS2926584_C1no" = SRS2926584_C1no, 
                    "SRS2926585_C1no" = SRS2926585_C1no,
                    "SRS2926586_S1no" = SRS2926586_S1no, 
                    "SRS2926587_C1no" = SRS2926587_C1no,
                    "SRS2926588_S1no" = SRS2926588_S1no, 
                    "SRS2926589_S1no" = SRS2926589_S1no,
                    "SRS2926590_C1no" = SRS2926590_C1no, 
                    "SRS2926591_C1no" = SRS2926591_C1no,
                    "SRS2926592_S1no" = SRS2926592_S1no)

# Sanity check: Sample IDs in sample_list & pheno file are in the same order
identical(names(sample_list), sample_ids) # Should be TRUE

# Convert list of lists to single data frame
g_original <- Reduce(function(x, y) 
  full_join(x, y, by = c("Gene.ID", "Gene.Name", "Reference",
                         "Start", "End")), sample_list)

# Prepare data: 1) Rename columns with sample IDs, 2) Calculate gene length, 
# 3) Add length as a suffix so duplicates can be differentiated, 4) Remove
# samples that are not of interest
gene_fpkm <- g_original %>% 
  rename_with(function(x) sample_ids, .cols = c(6:ncol(.))) %>%
  mutate(Length = End - Start + 1, .after = End) %>%
  unite(Gene.ID.Length, Gene.ID, Length, sep = "-len", remove = FALSE) %>%
  select(-contains(c("_S1", "_C1")))

# Remove samples from pheno that are not of interest
pheno <- pheno %>% filter(Abstinence_Length == "30")



# FIND DUPLICATE GENES ----------
# NOTE: Duplicates (same Gene.ID & Gene.Name) can be tracked w/ Gene.ID.Length
# There are 53,721 genes before fixing duplicates
dim(gene_fpkm) # 53721, 28

# Find duplicate genes
# Sometimes genes have duplicated names but different IDs, so check by Gene.ID
dup_genes_info <- gene_fpkm %>% 
  group_by(across(Gene.ID)) %>% 
  filter(n() > 1) %>% 
  ungroup()

# Number of duplicates
nrow(dup_genes_info) # 29

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
nrow(distinct_nas) # 8

# The duplicate genes and NA genes are the same
all_equal(distinct_dups, distinct_nas) # TRUE



## COMBINE "DUPLICATES" (SUM EXPRESSION) ----------
# Sum duplicates
dup_fixed <- dup_genes_info %>%
  group_by(Gene.ID, Gene.Name, Reference) %>%
  summarize_at(vars("SRS2926555_C30sal":"SRS2926576_S30sal"), 
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
dim(gene_fpkm_nodup) # 53,700 genes, 18 columns (11 samples)

# Get mean FPKM per treatment per gene
S30sal_means <- gene_fpkm_nodup %>%
  select(ends_with("S30sal")) %>%
  apply(1, mean)

C30sal_means <- gene_fpkm_nodup %>%
  select(ends_with("C30sal")) %>%
  apply(1, mean)

# Create single data frame of FPKM means
fpkm_means <- data.frame(
  cbind("Gene.ID" = gene_fpkm_nodup$Gene.ID,
        "Reference" = gene_fpkm_nodup$Reference,
        "S30sal" = S30sal_means, "C30sal" = C30sal_means))

# Only keep genes that are > 0.5 FPKM at least 1 treatment group and expressed 
keep_genes <- fpkm_means %>% 
  filter(S30sal > 0.5 | C30sal > 0.5)

# Remove genes & verify that this has been done correctly
keep_gene_fpkm <- gene_fpkm_nodup %>%
  semi_join(keep_genes, by = "Gene.ID") %>%
  select(-Gene.ID.Length)

nrow(keep_genes) == nrow(keep_gene_fpkm) # Should be TRUE

# Filter by raw read count, at least 6 in at least 5 samples (size of smallest treatment group)
keep_count <- rowSums(g_count > 6) >= 5
keep_count <- g_count[keep_count, ]

# Get full list of "keep" genes from raw counts data frame
final_genes <- semi_join(keep_count, keep_gene_fpkm, by = "Gene.ID")
nrow(final_genes) # Total genes after filtering out low expression: 16,840



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

# Write files
write.table(final_gene_matrix, file = "walker_original_gene_FPKM.txt")
write.table(gene_info, file = "walker_original_gene_info.txt")



## CREATE & MODIFY DGE LIST ----------
# Create initial DGE List object using filtered counts (not FPKM)
gene_dge <- DGEList(counts = final_gene_matrix, genes = gene_info)

# Add treatment group to samples
gene_dge$samples$treatment <- factor(pheno$Treatment_Group, 
                                     levels = c("S30sal", "C30sal"))

# Convert to log2 CPM (not normalized)
gene_lcpm <- cpm(gene_dge, log = TRUE)



## MULTI-DIMENSIONAL SCALING PLOTS ----------
# Use MDS plots to check for outliers

# Get treatment variable
treatment <- gene_dge$samples$treatment

# Get and set colors for variable
get_set_color <- function(var, brewer.pal = "Set1") {
  levels(var) <- brewer.pal(nlevels(var), brewer.pal)
  colors <- as.character(var)
  colors
}

color_treatment <- get_set_color(treatment)

# Simplify sample names
simple_sample <- rownames(gene_dge$samples) %>%
  gsub("_.{4,6}", "", .) %>%
  gsub("SRS.{5}", "SRS-", .)

# MDS plots
par(mfrow = c(1, 2))

# PCs 1 & 2
# Outlier: SRS2926555_C30sal
plotMDS(gene_lcpm, main = "Treatment Group", label = treatment, 
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm, main = "Samples", label = simple_sample, top = 100, 
        gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)

# PCs 3 & 4
# Does not look like there are any outliers for these PCs
plotMDS(gene_lcpm, main = "Treatment Group", label = treatment, 
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)
plotMDS(gene_lcpm, main = "Samples", label = simple_sample, top = 100, 
        gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)



## REMOVE OUTLIER & RE-FILTER GENES ----------
# Remove outlier samples
g_count <- g_count %>% select(-SRS2926555_C30sal)
gene_fpkm_nodup <- gene_fpkm_nodup %>% select(-SRS2926555_C30sal)
pheno <- pheno %>% filter(Sample_Treatment != "SRS2926555_C30sal")

# Check that all are filtered
(ncol(gene_fpkm_nodup) - 7) == nrow(pheno) # Must be TRUE
(ncol(g_count) - 2) == nrow(pheno) # Must be TRUE

# Number of genes before filtering
dim(gene_fpkm_nodup) # 53,700 genes, 17 columns (10 samples)

# Get mean FPKM per treatment per gene
S30sal_means <- gene_fpkm_nodup %>%
  select(ends_with("S30sal")) %>%
  apply(1, mean)

C30sal_means <- gene_fpkm_nodup %>%
  select(ends_with("C30sal")) %>%
  apply(1, mean)

# Create single data frame of FPKM means
fpkm_means <- data.frame(
  cbind("Gene.ID" = gene_fpkm_nodup$Gene.ID,
        "Reference" = gene_fpkm_nodup$Reference,
        "S30sal" = S30sal_means, "C30sal" = C30sal_means))

# Only keep genes that are > 0.5 FPKM at least 1 treatment group
keep_genes <- fpkm_means %>% 
  filter(S30sal > 0.5 | C30sal > 0.5)

# Remove genes & verify that this has been done correctly
keep_gene_fpkm <- gene_fpkm_nodup %>%
  semi_join(keep_genes, by = "Gene.ID") %>%
  select(-Gene.ID.Length)

nrow(keep_genes) == nrow(keep_gene_fpkm) # Should be TRUE

# Filter by raw read count, at least 6 in at least 5 samples (size of smallest treatment group)
keep_count <- rowSums(g_count > 6) >= 5
keep_count <- g_count[keep_count, ]

# Get full list of "keep" genes from raw counts data frame
final_genes <- semi_join(keep_count, keep_gene_fpkm, by = "Gene.ID")
nrow(final_genes) # Total genes after filtering out low expression: 16,848



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

# Write files
write.table(final_gene_matrix, file = "walker_FINAL_gene_FPKM.txt")
write.table(gene_info, file = "walker_FINAL_gene_info.txt")



## REMAKE DGE OBJECT & MDS PLOTS WITHOUT OUTLIER ----------
# Create DGE List object using filtered counts (not FPKM)
gene_dge <- DGEList(counts = final_gene_matrix, genes = gene_info)

# Add pheno info to samples
gene_dge$samples$treatment <- factor(pheno$Treatment_Group,
                                     levels = c("S30sal", "C30sal"))
gene_dge$samples$batch <- factor(pheno$Lab_Run_ID,
                                 levels = c("Run_714", "Run_715", "Run_716",
                                            "Run_717"))
gene_dge$samples$instrument <- factor(pheno$Instrument_Name,
                                      levels = c("HISEQ", "HWI-D00743",
                                                 "HWI-D00381"))


# Convert to log2 CPM (not normalized)
gene_lcpm <- cpm(gene_dge, log = TRUE)

# Get variables
treatment <- gene_dge$samples$treatment 
batch <- gene_dge$samples$batch
instrument <- gene_dge$samples$instrument

# Set colors for variables
color_treatment <- get_set_color(fct_rev(treatment))
color_batch <- get_set_color(batch)
color_instrument <- get_set_color(instrument)

# Simplify sample names
simple_sample <- rownames(gene_dge$samples) %>%
  gsub("_.{4,6}", "", .) %>%
  gsub("SRS.{5}", "SRS-", .)

# Simplify Batch and Instrument names
simple_batch <- gsub("Run_", "", gene_dge$samples$batch)
simple_instrument <- gene_dge$samples$instrument %>%
  gsub("D00", "", .) %>%
  strtrim(7)

# MDS plots
par(mfrow = c(2, 2))

# PCs 1 & 2: No clear outliers
plotMDS(gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm, main = "Samples", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm, main = "Batch", label = simple_batch, cex = 1.5,
        col = color_batch, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm, main = "Sequencing Instrument", label = simple_instrument, 
        col = color_instrument, top = 100, gene.selection = "common", cex = 1.5,
        dim.plot = c(1, 2), plot = TRUE)

# PCs 2 & 3
plotMDS(gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)
plotMDS(gene_lcpm, main = "Samples", label = simple_sample, cex = 1.5, 
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)
plotMDS(gene_lcpm, main = "Batch", label = simple_batch, cex = 1.5,
        col = color_batch, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)
plotMDS(gene_lcpm, main = "Sequencing Instrument", label = simple_instrument, 
        cex = 1.5, col = color_instrument, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)



## VOOM TRANSFORMATIONS & VARIANCE PARTITION ----------
# Get TMM-normalization factors before using voom
gene_dge <- calcNormFactors(gene_dge, normalize = "TMM")

# Add pheno info to samples
gene_dge$samples$treatment <- factor(pheno$Treatment_Group,
                                     levels = c("S30sal", "C30sal"))
gene_dge$samples$batch <- factor(pheno$Lab_Run_ID,
                                 levels = c("Run_714", "Run_715", 
                                            "Run_716", "Run_717"))
gene_dge$samples$instrument <- factor(pheno$Instrument_Name,
                                      levels = c("HISEQ", "HWI-D00743",
                                                 "HWI-D00381"))

# Re-define group variables
treatment <- gene_dge$samples$treatment
batch <- gene_dge$samples$batch
instrument <- gene_dge$samples$instrument

# Test several designs/models to see how variances partition based on design
# Treatment (C1no, C30sal). Note that Treatment = Abstinence Length.
# To control for batch effects, we need to include either Batch or Instrument.

# Treatment
des_treat <- model.matrix(~ 0 + treatment)
colnames(des_treat) <- gsub("treatment", "", colnames(des_treat))

# Treatment + Batch
des_treat_batch <- model.matrix(~ 0 + treatment + batch)
colnames(des_treat_batch) <- gsub("treatment", "", colnames(des_treat_batch))
colnames(des_treat_batch) <- gsub("batch", "", colnames(des_treat_batch))

# Treatment + Instrument
des_treat_ins <- model.matrix(~ 0 + treatment + instrument)
colnames(des_treat_ins) <- gsub("treatment", "", colnames(des_treat_ins))
colnames(des_treat_ins) <- gsub("instrument", "", colnames(des_treat_ins))
colnames(des_treat_ins) <- gsub("-", "_", colnames(des_treat_ins))

# Treatment + Batch + Instrument
des_treat_batch_ins <- model.matrix(~ 0 + treatment + batch + instrument)
colnames(des_treat_batch_ins) <- gsub("treatment", "", colnames(des_treat_batch_ins))
colnames(des_treat_batch_ins) <- gsub("batch", "", colnames(des_treat_batch_ins))
colnames(des_treat_batch_ins) <- gsub("instrument", "", colnames(des_treat_batch_ins))
colnames(des_treat_batch_ins) <- gsub("-", "_", colnames(des_treat_batch_ins))


# Voom transformations
par(mfrow = c(1, 1))
v_treat <- voom(gene_dge, des_treat, plot = TRUE)
v_treat_batch <- voom(gene_dge, des_treat_batch, plot = TRUE)
v_treat_ins <- voom(gene_dge, des_treat_ins, plot = TRUE)
v_treat_batch_ins <- voom(gene_dge, des_treat_batch_ins, plot = TRUE)

# Form for variance partition; include all variables
form <- ~ (1|treatment) + (1|batch) + (1|instrument)
# Only run fitExtractVarPartModel() if you're prepared for a time-consuming analysis
# cl <- makeCluster(2)
# registerDoParallel(cl)
# stopCluster(cl)

# pdf("Test.pdf", width = 5, height = 5)
# dev.off()

# Palette
palette <- c(wes_palette("Chevalier1", 4, type = "discrete")[-3], "grey85")

# Treatment model
# More variance is explained by Instrument than Batch in this model
# vp_treat <- fitExtractVarPartModel(v_treat, form, gene_dge$samples)
# saveRDS(vp_treat, file = "walker_FINAL_treat.RDS")
vp_treat <- readRDS(file = "walker_FINAL_treat.RDS")
vp_treat <- sortCols(vp_treat)
bars_treat <- plotPercentBars(vp_treat[1:10, ]) +
  scale_fill_manual(values = palette,
                    labels = c("Batch", "Instrument", "Treatment", "Residuals"))
varpart_treat <- plotVarPart(vp_treat) +
  scale_x_discrete(labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  scale_fill_manual(values = palette) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major.x = element_blank())
ggarrange(varpart_treat, bars_treat)

# Treatment + Batch model
# More variance is explained by Batch than Instrument with this model
# vp_treat_batch <- fitExtractVarPartModel(v_treat_batch, form, gene_dge$samples)
# saveRDS(vp_treat_batch, file = "walker_FINAL_treat_batch.RDS")
vp_treat_batch <- readRDS(file = "walker_FINAL_treat_batch.RDS")
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
# vp_treat_ins <- fitExtractVarPartModel(v_treat_ins, form, gene_dge$samples)
# saveRDS(vp_treat_ins, file = "walker_FINAL_treat_ins.RDS")
vp_treat_ins <- readRDS(file = "walker_FINAL_treat_ins.RDS")
vp_treat_ins <- sortCols(vp_treat_ins)
bars_ins <- plotPercentBars(vp_treat_ins[1:10, ]) +
  scale_fill_manual(values = palette,
                    labels = c("Batch", "Instrument", "Treatment", "Residuals"))
varpart_ins <- plotVarPart(vp_treat_ins) +
  scale_x_discrete(labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  scale_fill_manual(values = palette) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major.x = element_blank())
ggarrange(varpart_ins, bars_ins)
# 800, 215 dim

# Using both Batch and Instrument in our model could overfit the data (and in
# fact, gives us the same results - see below with voom).
# We will select the Treatment + Instrument model because, using this model,
# Treatment accounts for more variance than Batch does with the Treatment +
# Batch model. With Treatment + Instrument, we also capture the variation from
# Batch, while we don't capture the variation of Instrument with the Treatment +
# Batch model. There's also less variance explained by residuals for this model.



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

v_treat <- voom(gene_dge, des_treat, plot = TRUE)
plotSA(treat_ebayes, main = "Final model: mean-variance trend")

v_treat_batch <- voom(gene_dge, des_treat_batch, plot = TRUE)
plotSA(treat_batch_ebayes, main = "Final model: mean-variance trend")

v_treat_ins <- voom(gene_dge, des_treat_ins, plot = TRUE)
plotSA(treat_ins_ebayes, main = "Final model: mean-variance trend")

v_treat_batch_ins <- voom(gene_dge, des_treat_batch_ins, plot = TRUE)
plotSA(treat_batch_ins_ebayes, main = "Final model: mean-variance trend")

# Get stats for all genes, in same order as they were originally
treat_stats <- topTable(treat_ebayes, number = nrow(treat_ebayes))
treat_batch_stats <- topTable(treat_batch_ebayes, number = nrow(treat_batch_ebayes))
treat_ins_stats <- topTable(treat_ins_ebayes, number = nrow(treat_ins_ebayes))
treat_batch_ins_stats <- topTable(treat_batch_ins_ebayes, number = nrow(treat_batch_ins_ebayes))

# Get genes w/ P-value < 0.05
sig_treat <- treat_stats %>% filter(P.Value < 0.05) %>% select(Gene.ID) %>% mutate(Design = "treat")
sig_treat_batch <- treat_batch_stats %>% filter(P.Value < 0.05) %>% select(Gene.ID) %>% mutate(Design = "treat_batch")
sig_treat_ins <- treat_ins_stats %>% filter(P.Value < 0.05) %>% select(Gene.ID) %>% mutate(Design = "treat_ins")
sig_treat_batch_ins <- treat_batch_ins_stats %>% filter(P.Value < 0.05) %>% select(Gene.ID) %>% mutate(Design = "treat_batch_ins")

df <- Reduce(function(x, y) 
  full_join(x, y), list(sig_treat, sig_treat_batch, sig_treat_ins, sig_treat_batch_ins)) %>%
  group_by(Gene.ID) %>%
  mutate(count = length(Gene.ID))

sig_treat <- sig_treat %>% pull(Gene.ID) # 236 genes
sig_treat_batch <- sig_treat_batch %>% pull(Gene.ID) # 203 genes
sig_treat_ins <- sig_treat_ins %>% pull(Gene.ID) # 227 genes
sig_treat_batch_ins <- sig_treat_batch_ins %>% pull(Gene.ID) # 203 genes

sig_treat_batch[sig_treat_batch %in% sig_treat_batch_ins] %>% length() # All the same (203)

sig_treat[sig_treat %in% sig_treat_batch] %>% length() # Half the same (132/203-236)
sig_treat_ins[sig_treat_ins %in% sig_treat_batch] %>% length() # Mostly the same (176/203-227)
sig_treat[sig_treat %in% sig_treat_ins] %>% length() # Half the same (147/203-236)

# Write files
write.table(treat_batch_stats, file = "walker_DE_FINAL.txt")



## MDS PLOTS FOR VOOM MODELS ----------
vCorrectBatch <- removeBatchEffect(v_treat_batch, batch = v_treat_batch$targets$batch)
plotMDS(v_treat_batch, main = "Samples", label = treatment, cex = 1.5, col = color_treatment,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(vCorrectBatch, main = "Samples", label = treatment, cex = 1.5, col = color_treatment,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)

vCorrectIns <- removeBatchEffect(v_treat_ins, batch = v_treat_ins$targets$ins)
plotMDS(v_treat_ins, main = "Samples", label = treatment, cex = 1.5, col = color_treatment,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(vCorrectIns, main = "Samples", label = treatment, cex = 1.5, col = color_treatment,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)

vCorrectBatchIns <- removeBatchEffect(v_treat_batch_ins, batch = v_treat_batch_ins$targets$batch)
vCorrectBatchIns <- removeBatchEffect(vCorrectBatchIns, batch = v_treat_batch_ins$targets$ins)
plotMDS(v_treat_batch_ins, main = "Samples", label = treatment, cex = 1.5, col = color_treatment,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(vCorrectBatchIns, main = "Samples", label = treatment, cex = 1.5, col = color_treatment,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)



## CREATE OBJECT FOR POWER ANALYSIS ----------
walker_dge <- gene_dge # Already calculated norm. factors for TMM method
walker_stats <- treat_batch_stats
walker_summary <- walker_stats %>%
  filter(P.Value < 0.05) %>%
  summarize(Up = sum(logFC > 0), Down = sum(logFC < 0), Total = Up + Down)
walker_gene_mat <- final_gene_matrix
walker_gene_info <- gene_info
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

saveRDS(walker_obj, file = "walker_obj_4power.RDS")

points <- color_batch %>% 
  str_replace_all("#4DAF4A", "15") %>%
  str_replace_all("#377EB8", "16") %>%
  str_replace_all("#E41A1C", "17") %>%
  str_replace_all("#984EA3", "18") %>%
  as.numeric()

# Final figure
svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Walker_RNAseq.svg", width = 5.4, height = 8)
par(mfrow = c(3, 2))
plotMDS(gene_lcpm,
        pch = points, cex = 3, col = color_treatment, cex.axis = 1.25, cex.lab = 1.25,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm,
        pch = points, cex = 3, col = color_treatment, cex.axis = 1.25, cex.lab = 1.25,
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)
# Original voom data
voom_info <- voom(gene_dge, des_treat_batch, plot = FALSE, save.plot = TRUE)
plot(voom_info$voom.xy, pch = 19, cex = 0.5, cex.main = 1.25, cex.axis = 1.25, 
     cex.lab = 1.25, main = "voom: Mean-variance trend", 
     xlab = "log2(count size + 0.5)", ylab = "Sqrt(standard deviation)")
lines(voom_info$voom.line, col = "red")
# Final voom model
plotSA(treat_batch_ebayes, main = "Final model: Mean-variance trend", 
       cex = 0.5, pch = 19, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.25)
plotMDS(v_treat_batch, pch = points,  cex = 3, cex.axis = 1.25, cex.lab = 1.25,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(v_treat_batch, pch = points,  cex = 3, cex.axis = 1.25, cex.lab = 1.25,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)
dev.off()


# Final figure for presentation
svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Walker_RNAseq_simple_presentation.svg", width = 4.5, height = 4.5)
par(mfrow = c(1, 1))
plotMDS(gene_lcpm,
        pch = 16, cex = 3, col = color_treatment, cex.axis = 1.5, cex.lab = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
dev.off()


svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Walker_varPart.svg", width = 6.9, height = 2.6)
bars_batch <- plotPercentBars(vp_treat[1:10, ]) +
  scale_fill_manual(values = palette,
                    labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  theme(legend.position = "none")
varpart_batch <- plotVarPart(vp_treat_batch) +
  scale_x_discrete(labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  scale_fill_manual(values = palette, 
                    labels = c("Batch", "Instrument", "Treatment", "Residuals")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())
ggarrange(varpart_batch, bars_batch, common.legend = TRUE, legend = "right")
dev.off()




