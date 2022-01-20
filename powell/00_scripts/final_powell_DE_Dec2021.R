
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
pheno_dir <- paste0(main_dir, "powell/")
gene_dir <- paste0(main_dir, "powell/genes/")
setwd(paste0(main_dir, "results/"))



# LOAD IN DATA & RENAME COLUMNS ----------
# NOTE: Sample order must be maintained for accuracy.
# Sample order is always numerical from low (SRS6085261) to high (SRS6085272).

# Load in phenotype file
pheno <- read_csv(paste0(pheno_dir, "powell_pheno.csv", sep = ""))
sample_ids <- pheno$Sample_Treatment # Vector of sample IDs

# Load in gene counts files
g_count <- read_csv(paste0(gene_dir, "powell_gene_count_matrix.csv"))

# Separate initial gene_id column into Gene.ID & Gene.Name
# E.g. ENSRNOG00000040300|Raet1e becomes EENSRNOG00000040300 & Raet1e
# Note that there will be some NAs for genes with no name.
fix_gene_ids <- function(df) {separate(df, col = "gene_id", 
                                       into = c("Gene.ID", "Gene.Name"), 
                                       sep = "\\|")}

# Fix gene IDs and remove samples that aren't of interest
g_count <- fix_gene_ids(g_count) %>% dplyr::select(-contains("_CE"))

# 396 rows don't have gene names and are converted to NA
# Fix NAs by placing Gene.ID in Gene.Name column
g_count$Gene.Name[is.na(g_count$Gene.Name)] <- g_count$Gene.ID[is.na(g_count$Gene.Name)]

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
sample_list <- list("SRS6085261_CI1cue" = SRS6085261_CI1cue,
                    "SRS6085262_CI21cue" = SRS6085262_CI21cue,
                    "SRS6085263_CE1cue" = SRS6085263_CE1cue,
                    "SRS6085264_CE21cue" = SRS6085264_CE21cue,
                    "SRS6085265_CE21cue" = SRS6085265_CE21cue,
                    "SRS6085266_CE1cue" = SRS6085266_CE1cue,
                    "SRS6085267_CI1cue" = SRS6085267_CI1cue,
                    "SRS6085268_CI1cue" = SRS6085268_CI1cue,
                    "SRS6085269_CI21cue" = SRS6085269_CI21cue,
                    "SRS6085270_CE1cue" = SRS6085270_CE1cue,
                    "SRS6085271_CE21cue" = SRS6085271_CE21cue,
                    "SRS6085272_CI21cue" = SRS6085272_CI21cue)

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
  dplyr::select(-contains("_CE"))


# Fix - and . by placing Gene.ID in Gene.Name column
gene_fpkm$Gene.Name[gene_fpkm$Gene.Name == "-" | 
                      gene_fpkm$Gene.Name == "."] <- gene_fpkm$Gene.ID[gene_fpkm$Gene.Name == "-" |
                                                                         gene_fpkm$Gene.Name == "."]

# Remove samples from pheno that are not of interest
pheno <- pheno %>% filter(Housing_Condition == "Isolation")



# FIND DUPLICATE GENES ----------
# NOTE: Duplicates (same Gene.ID & Gene.Name) can be tracked w/ Gene.ID.Length
# There are 33,030 genes before fixing duplicates
dim(gene_fpkm) # 33030, 13

# Find duplicate genes
# Sometimes genes have duplicated names but different IDs, so check by Gene.ID
dup_genes_info <- gene_fpkm %>% 
  group_by(across(Gene.ID)) %>% 
  filter(n() > 1) %>% 
  ungroup()

# Number of duplicates
nrow(dup_genes_info) # 269

# Get *unique* duplicated genes for each data set
# Same ID and chromosome (NOT paralogs)
distinct_dups <- dup_genes_info %>%
  distinct(Gene.ID, Gene.Name, Reference) %>%
  arrange(Gene.Name)

# Number of unique duplicates
nrow(distinct_dups) # 122

# Sanity check: Are duplicated genes the only ones with NAs in the data?
# Identify genes with NAs in 1 or more samples (FPKM)
na_genes_info <- gene_fpkm[!complete.cases(gene_fpkm), ]

# Number of NAs
nrow(na_genes_info) # 198

# Get *unique* NAs
distinct_nas <- na_genes_info %>%
  distinct(Gene.ID, Gene.Name, Reference) %>%
  arrange(Gene.Name)

# Number of unique NAs
nrow(distinct_nas) # 101

# The duplicate genes and NA genes are NOT the same
all_equal(distinct_dups, distinct_nas) # "Different number of rows"

# There are more duplicates than there are NAs.
# What does expression look like for these duplicates?
exp_dup <- cbind(dup_genes_info[, 1:4], 
                 rowMean = rowMeans(dup_genes_info[, 8:13], na.rm = TRUE))
exp_dup %>% summary()

# Most have low expression (when summed across samples) and are likely to be 
# filtered out in the next step, but some will stay for the long haul.

# List of duplicated genes
dup_gene_names <- exp_dup %>% pull(Gene.Name) %>% unique() %>% sort()



## COMBINE "DUPLICATES" (SUM EXPRESSION) ----------
# Sum duplicates
# Start, End, and Length are no longer reliable and should be converted to NAs
dup_fixed <- dup_genes_info %>%
  group_by(Gene.ID, Gene.Name, Reference) %>%
  summarize_at(vars("SRS6085261_CI1cue":"SRS6085272_CI21cue"), 
               sum, na.rm = TRUE) %>%
  mutate(Start = NA, End = NA, Length = NA, .after = Reference)

# Remove duplicates from original dataset, then add back 
gene_fpkm_nodup <- gene_fpkm[!gene_fpkm$Gene.ID %in% distinct_dups$Gene.ID, ] %>%
  full_join(dup_fixed)

# Check that this was done correctly (# of genes by FPKM matches raw counts)
nrow(gene_fpkm_nodup) == nrow(g_count) # Should be TRUE



## FILTER LOWLY-EXPRESSED GENES (FPKM) ----------
dim(gene_fpkm_nodup) # 32,883 genes, 13 columns (6 samples)

# Get mean FPKM per treatment per gene
CI1cue_means <- gene_fpkm_nodup %>%
  dplyr::select(ends_with("CI1cue")) %>%
  apply(1, mean)

CI21cue_means <- gene_fpkm_nodup %>%
  dplyr::select(ends_with("CI21cue")) %>%
  apply(1, mean)

# Create single data frame of FPKM means
fpkm_means <- data.frame(
  cbind("Gene.ID" = gene_fpkm_nodup$Gene.ID,
        "Reference" = gene_fpkm_nodup$Reference,
        "CI1cue" = CI1cue_means, "CI21cue" = CI21cue_means))

# Only keep genes that are > 0.5 FPKM in 1 treatment group
keep_genes <- fpkm_means %>% 
  filter(CI1cue > 0.5 | CI21cue > 0.5)

# Remove genes & verify that this has been done correctly
keep_gene_fpkm <- gene_fpkm_nodup %>%
  semi_join(keep_genes, by = "Gene.ID") %>%
  dplyr::select(-Gene.ID.Length)

nrow(keep_genes) == nrow(keep_gene_fpkm) # Should be TRUE

# Filter by raw read count, at least 6 in at least 3 samples (size of smallest treatment group)
keep_count <- rowSums(g_count >= 6) >= 3
keep_count <- g_count[keep_count, ]

# Get full list of "keep" genes from raw counts data frame
final_genes <- semi_join(keep_count, keep_gene_fpkm, by = "Gene.ID")
nrow(final_genes) # Total genes after filtering out low expression: 14,872



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
write.table(final_gene_matrix, file = "powell_original_gene_FPKM.txt")
write.table(gene_info, file = "powell_original_gene_info.txt")



## CREATE & MODIFY DGE LIST ----------
# Create initial DGE List object using filtered counts (not FPKM)
gene_dge <- DGEList(counts = final_gene_matrix, genes = gene_info)

# Add treatment group to samples
gene_dge$samples$treatment <- factor(pheno$Treatment_Group, 
                                     levels = c("CI1cue", "CI21cue"))

# Convert to log2 CPM (not normalized)
gene_lcpm <- cpm(gene_dge, log = TRUE)



## MULTI-DIMENSIONAL SCALING PLOTS (PART I, NO OUTLIERS REMOVED) ----------
# Use MDS plots to check for outliers

# Get treatment variable
treatment <- gene_dge$samples$treatment %>%
  gsub("cue", "", .) %>%
  as.factor()

# Get and set colors for variable
get_set_color <- function(var, brewer.pal = "Set1") {
  levels(var) <- brewer.pal(nlevels(var), brewer.pal)
  colors <- as.character(var)
  colors
}

color_treatment <- get_set_color(treatment)

# Simplify sample names
simple_sample <- rownames(gene_dge$samples) %>%
  gsub("_.{4,7}", "", .) %>%
  gsub("SRS.{5}", "SRS-", .)

# MDS plots
par(mfrow = c(1, 2))

# PCs 1 & 2
# Outlier: SRS6085269_CI21cue
plotMDS(gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm, main = "Samples", label = simple_sample, cex = 1.5, 
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)

# PCs 2 & 3
plotMDS(gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(2, 3), plot = TRUE)
plotMDS(gene_lcpm, main = "Samples", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(2, 3), plot = TRUE)



## REMOVE OUTLIER & RE-FILTER GENES ----------
# Remove outlier sample
g_count <- g_count %>% dplyr::select(-SRS6085269_CI21cue)
gene_fpkm_nodup <- gene_fpkm_nodup %>% dplyr::select(-SRS6085269_CI21cue)
pheno <- pheno %>% filter(Sample_Treatment != "SRS6085269_CI21cue")

# Check that all are filtered
(ncol(gene_fpkm_nodup) - 7) == nrow(pheno) # Must be TRUE
(ncol(g_count) - 2) == nrow(pheno) # Must be TRUE

# Number of genes before filtering
dim(gene_fpkm_nodup) # 32,883 genes, 12 columns (5 samples)

# Get mean FPKM per treatment per gene
CI1cue_means <- gene_fpkm_nodup %>%
  dplyr::select(ends_with("CI1cue")) %>%
  apply(1, mean)

CI21cue_means <- gene_fpkm_nodup %>%
  dplyr::select(ends_with("CI21cue")) %>%
  apply(1, mean)

# Create single data frame of FPKM means
fpkm_means <- data.frame(
  cbind("Gene.ID" = gene_fpkm_nodup$Gene.ID,
        "Reference" = gene_fpkm_nodup$Reference,
        "CI1cue" = CI1cue_means, "CI21cue" = CI21cue_means))

# Only keep genes that are > 0.5 FPKM in all 1 treatment group
keep_genes <- fpkm_means %>% 
  filter(CI1cue > 0.5 | CI21cue > 0.5)

# Remove genes & verify that this has been done correctly
keep_gene_fpkm <- gene_fpkm_nodup %>%
  semi_join(keep_genes, by = "Gene.ID") %>%
  dplyr::select(-Gene.ID.Length)

nrow(keep_genes) == nrow(keep_gene_fpkm) # Should be TRUE

# Filter by raw read count, at least 6 in at least 2 samples (size of smallest treatment group)
keep_count <- rowSums(g_count >= 6) >= 2
keep_count <- g_count[keep_count, ]

# Get full list of "keep" genes from raw counts data frame
final_genes <- semi_join(keep_count, keep_gene_fpkm, by = "Gene.ID")
nrow(final_genes) # Total genes after filtering out low expression: 15,812



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
write.table(final_gene_matrix, file = "powell_FINAL_gene_FPKM.txt")
write.table(gene_info, file = "powell_FINAL_gene_info.txt")



## REMAKE DGE OBJECT & MDS PLOTS WITHOUT OUTLIER ----------
# Create DGE List object using filtered counts (not FPKM)
gene_dge <- DGEList(counts = final_gene_matrix, genes = gene_info)

# Add pheno info to samples
# Treatment_Group = Housing + Abstinence_Length
gene_dge$samples$treatment <- factor(pheno$Treatment_Group,
                                     levels = c("CI1cue", "CI21cue"))
gene_dge$samples$lane <- factor(pheno$Lane, levels = c("25", "26"))

# Convert to log2 CPM (not normalized)
gene_lcpm <- cpm(gene_dge, log = TRUE)

# Get variables
treatment <- gene_dge$samples$treatment %>%
  gsub("cue", "", .) %>%
  as.factor()
lane <- gene_dge$samples$lane

# Set colors for variables
color_treatment <- get_set_color(fct_rev(treatment))
color_lane <- get_set_color(lane)

# Simplify sample names
simple_sample <- rownames(gene_dge$samples) %>%
  gsub("_.{4,7}", "", .) %>%
  gsub("SRS.{5}", "SRS-", .)

# MDS plots by each variable
par(mfrow = c(1, 3))

# PCs 1 & 2
plotMDS(gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm, main = "Sequencing Lane", label = lane, cex = 1.5, 
        col = color_lane, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm, main = "Sample", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)

# PCs 2 & 3
plotMDS(gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(2, 3), plot = TRUE)
plotMDS(gene_lcpm, main = "Sequencing Lane", label = lane, cex = 1.5, 
        col = color_lane, top = 100, gene.selection = "common", 
        dim.plot = c(2, 3), plot = TRUE)
plotMDS(gene_lcpm, main = "Sample", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)

# Final plots
par(mfrow = c(2, 2))
plotMDS(gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm, main = "Sequencing Lane", label = lane, cex = 1.5, 
        col = color_lane, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(2, 3), plot = TRUE)
plotMDS(gene_lcpm, main = "Sequencing Lane", label = lane, cex = 1.5, 
        col = color_lane, top = 100, gene.selection = "common", 
        dim.plot = c(2, 3), plot = TRUE)
# 550, 600 dim



## VOOM TRANSFORMATIONS & VARIANCE PARTITION ----------
# Get TMM-normalization factors before using voom
gene_dge <- calcNormFactors(gene_dge, normalize = "TMM")

# Add pheno info to samples
gene_dge$samples$treatment <- factor(pheno$Treatment_Group,
                                     levels = c("CI1cue", "CI21cue"))
gene_dge$samples$lane <- factor(pheno$Lane, levels = c("25", "26"))

# Re-define group variables
treatment <- gene_dge$samples$treatment
lane <- gene_dge$samples$lane

# Test several designs/models to see how variances partition based on design
# Treatment (CI1_cue, CI21_cue). Note that Treatment = Abstinence Length.
# To control for batch effects, we may need to use the Lane variable.

# Treatment
des_treat <- model.matrix(~ 0 + treatment)
colnames(des_treat) <- gsub("treatment", "", colnames(des_treat))

# Treatment + Lane
des_treat_lane <- model.matrix(~ 0 + treatment + lane)
colnames(des_treat_lane) <- gsub("treatment", "", colnames(des_treat_lane))
colnames(des_treat_lane) <- gsub("lane", "L", colnames(des_treat_lane))

# Voom transformations
v_treat <- voom(gene_dge, des_treat, plot = TRUE)
v_treat_lane <- voom(gene_dge, des_treat_lane, plot = TRUE)
# 600, 350 dimensions


# MDS Plots based on voom models
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
        col = color_treatment, 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(v_treat_lane, main = "Sequencing Lane", label = lane, cex = 1.5, 
        col = color_lane, top = 1000, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(v_treat_lane, main = "Sample", label = simple_sample, cex = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(v_treat_lane, main = "Treatment Group", label = treatment, cex = 1.5,
        col = color_treatment, top = 1000, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)
plotMDS(v_treat_lane, main = "Sequencing Lane", label = lane, cex = 1.5, 
        col = color_lane, top = 1000, gene.selection = "common", 
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
# cl <- makeCluster(2)
# registerDoParallel(cl)
# stopCluster(cl)


# pdf("Test.pdf", width = 5, height = 5)
# dev.off()

# Palette
palette <- c(wes_palette("Chevalier1", 4, type = "discrete"), "grey85")[-c(2,3)]

# Treatment model
# vp_treat <- fitExtractVarPartModel(v_treat, form, gene_dge$samples)
# saveRDS(vp_treat, file = "powell_FINAL_treat.RDS")
vp_treat <- readRDS(file = "powell_FINAL_treat.RDS")
vp_treat <- sortCols(vp_treat)
bars_treat <- plotPercentBars(vp_treat[1:10, ]) +
  scale_fill_manual(values = c("#440154FF", "#FDE725FF", "grey85"),
                    labels = c("Lane", "Treatment", "Residuals")) +
  theme(legend.position = "bottom")
varpart_treat <- plotVarPart(vp_treat) +
  scale_x_discrete(labels = c("Lane", "Treatment", "Residuals")) +
  scale_fill_manual(values = c("#440154FF", "#FDE725FF", "grey85")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
#pdf("Test.pdf", width = 10, height = 5)
ggarrange(varpart_treat, bars_treat)
#dev.off()

# Treatment + Lane model
# vp_treat_lane <- fitExtractVarPartModel(v_treat_lane, form, gene_dge$samples)
# saveRDS(vp_treat_lane, file = "powell_FINAL_treat_lane.RDS")
vp_treat_lane <- readRDS(file = "powell_FINAL_treat_lane.RDS")
vp_treat_lane <- sortCols(vp_treat_lane)
bars_lane <- plotPercentBars(vp_treat_lane[1:10, ]) +
  scale_fill_manual(values = palette,
                    labels = c("Lane", "Treatment", "Residuals"))
varpart_lane <- plotVarPart(vp_treat_lane) +
  scale_x_discrete(labels = c("Lane", "Treatment", "Residuals")) +
  scale_fill_manual(values = palette) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
ggarrange(varpart_lane, bars_lane)

# The results are almost the same for both, but use lane in the model.



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

v_treat <- voom(gene_dge, des_treat, plot = TRUE)
plotSA(treat_ebayes, main = "Final model: mean-variance trend")

v_treat_lane <- voom(gene_dge, des_treat_lane, plot = TRUE)
plotSA(treat_lane_ebayes, main = "Final model: mean-variance trend")
# 400, 350 dimensions


# Write files
# Get stats for all genes, in same order as they were originally
treat_stats <- topTable(treat_ebayes, number = nrow(treat_ebayes))

treat_lane_stats <- topTable(treat_lane_ebayes, number = nrow(treat_lane_ebayes))

sig_treat <- treat_stats %>% filter(P.Value < 0.05) %>% pull(Gene.ID) # 902 genes
sig_treat_lane <- treat_lane_stats %>% filter(P.Value < 0.05) %>% pull(Gene.ID) # 823 genes

sig_treat[sig_treat %in% sig_treat_lane] %>% length() # Only 658 shared

treat_lane_stats %>%
  filter(P.Value < 0.05) %>%
  summarize(Up = sum(logFC > 0), Down = sum(logFC < 0))

# Write files - can be used for getting overlaps (UpSetR plots)
write.table(treat_lane_stats, file = "powell_DE_FINAL.txt")



## CREATE OBJECT FOR POWER ANALYSIS ----------
powell_dge <- gene_dge # Already calculated norm. factors for TMM method
powell_stats <- treat_lane_stats
powell_summary <- powell_stats %>%
  filter(P.Value < 0.05) %>%
  summarize(Up = sum(logFC > 0), Down = sum(logFC < 0), Total = Up + Down)
powell_gene_mat <- final_gene_matrix
powell_gene_info <- gene_info
powell_fc <- powell_stats %>%
  filter(P.Value < 0.05) %>%
  mutate(Abs_FC = 2^abs(logFC)) %>%
  dplyr::select(Abs_FC) %>%
  summary()
powell_pheno <- pheno
powell_treat <- fct_relevel(pheno$Treatment_Group, levels = c("CI1cue", "CI21cue"))
powell_control <- levels(powell_treat)[1]

powell_obj <- list(DGE = powell_dge, Stats = powell_stats, 
                   Stats_Summary = powell_summary, 
                   Gene_Matrix = powell_gene_mat, Gene_Info = powell_gene_info, 
                   FC = powell_fc, Pheno = powell_pheno, Treat = powell_treat, 
                   Control = powell_control)

saveRDS(powell_obj, file = "powell_obj_4power.RDS")


## FINAL FIGURES ----------
points <- color_lane %>% 
  str_replace_all("#E41A1C", "16") %>%
  str_replace_all("#377EB8", "18") %>%
  as.numeric()

# Final figure
svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Powell_RNAseq.svg", width = 5.4, height = 8)
par(mfrow = c(3, 2))
plotMDS(gene_lcpm,
        pch = points, cex = 3, col = color_treatment, cex.axis = 1.25, cex.lab = 1.25,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
plotMDS(gene_lcpm,
        pch = points, cex = 3, col = color_treatment, cex.axis = 1.25, cex.lab = 1.25,
        top = 100, gene.selection = "common", dim.plot = c(3, 4), plot = TRUE)
# Original voom data
voom_info <- voom(gene_dge, des_treat_lane, plot = FALSE, save.plot = TRUE)
plot(voom_info$voom.xy, pch = 19, cex = 0.5, cex.main = 1.25, cex.axis = 1.25, 
     cex.lab = 1.25, main = "voom: Mean-variance trend", 
     xlab = "log2(count size + 0.5)", ylab = "Sqrt(standard deviation)")
lines(voom_info$voom.line, col = "red")
# Final voom model
plotSA(treat_lane_ebayes, main = "Final model: Mean-variance trend", 
       cex = 0.5, pch = 19, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.25)
plotMDS(v_treat_lane, pch = points,  cex = 3, cex.axis = 1.25, cex.lab = 1.25,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(1, 2), plot = TRUE)
plotMDS(v_treat_lane, pch = points,  cex = 3, cex.axis = 1.25, cex.lab = 1.25,
        col = color_treatment, top = 100, gene.selection = "common", 
        dim.plot = c(3, 4), plot = TRUE)
dev.off()


# Final figure for presentation
svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Powell_RNAseq_simple_presentation.svg", width = 4.5, height = 4.5)
par(mfrow = c(1, 1))
plotMDS(gene_lcpm,
        pch = 16, cex = 3, col = color_treatment, cex.axis = 1.5, cex.lab = 1.5,
        top = 100, gene.selection = "common", dim.plot = c(1, 2), plot = TRUE)
dev.off()

svglite("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Powell_varPart.svg", width = 6.9, height = 2.6)
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
dev.off()

