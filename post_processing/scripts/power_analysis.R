## LOAD PACKAGES ----------
library(Biobase)
library(edgeR)
library(qvalue)
library(ssizeRNA)
library(viridis)
library(tidyverse)
library(tictoc)


## BEFORE USING THIS SCRIPT ----------
# User should alter these variables as necessary
# Get list of all objects for each dataset
work_dir <- "G:/Annika/GitHub Repositories/RodentAddiction/post_processing/results/power_objects/"
carpenter_obj <- readRDS(file.path(work_dir, "carpenter_obj_4power.RDS"))
walker_obj <- readRDS(file.path(work_dir, "walker_obj_4power.RDS"))


## POWER ANALYSES ---------
# Create list of objects
every_obj <- list(carpenter = carpenter_obj, walker = walker_obj)

# Empty list of variables for each dataset
every_var <- list()

# Put variables for each dataset into the empty all_var list
for (i in 1:length(every_obj)) {
  
  # Extract variables
  treat <- every_obj[[i]][["Treat"]]
  control <- every_obj[[i]][["Control"]]
  nGenes <- nrow(every_obj[[i]][["Gene_Info"]]) # Total number of genes
  # Prop. of up-regulated genes among all DE genes
  up <- every_obj[[i]][["Stats_Summary"]][["Up"]]/every_obj[[i]][["Stats_Summary"]][["Total"]]
  actual_fc <- 2^abs(every_obj[[i]][["Stats"]][["logFC"]])
  
  # Calculate means per gene, for the control group
  mu <- apply(every_obj[[i]][["Gene_Matrix"]][, treat == control], 1, mean)
  
  # Dispersion
  temp_dge <- estimateCommonDisp(every_obj[[i]][["DGE"]])
  temp_dge <- estimateTagwiseDisp(temp_dge)
  disp <- temp_dge$tagwise.dispersion
  
  # Create sub-list and put variables there
  every_var[[i]] <- list()
  every_var[[i]][["disp"]] <- disp
  every_var[[i]][["nGenes"]] <- nGenes
  every_var[[i]][["mu"]] <- mu
  every_var[[i]][["up"]] <- up
  every_var[[i]][["actual_fc"]] <- actual_fc
}

# Fix names
every_var <- setNames(every_var, names(every_obj))

# Check dispersion
as.vector(every_var$carpenter$disp) %>% mean()
as.vector(every_var$carpenter$disp) %>% median()
as.vector(every_var$walker$disp) %>% mean()
as.vector(every_var$walker$disp) %>% median()

# Also get dataframe of variables
get_vars <- every_var %>%
  lapply(., as.data.frame) %>%
  map2(., names(every_var), ~ cbind(.x, Study = .y)) %>%
  Reduce(function(x, y) 
    full_join(x, y), .)

# Variable that will stay constant
m <- 200 # Pseudo sample size for generated data

# The variables pi0, fc, fdr, and power will be varied for the calculations.
# pi0 = Proportion of non-DEGs, e.g. 0.8 means 20% of genes are DE
# fdr = FDR threshold; 0.05 is more stringent, but 0.1 is also appropriate
# power = 1 - beta, a measure of ability to detect true positives; typically 
# 0.8, though 0.7 may also be appropriate

# Other variables are taken directly from the data within the below function
set.seed(042092)
compute_power <- function(dataset) {
  data <- every_var[[dataset]]
  print("Start: Compute Dataset")
  tic("Compute Dataset")
  
  for (i in 1:16) {
    # Create counter
    counter = i
    
    # Set variables depending on where we are in the loop
    if((counter %% 4) == 1) {
      fc_val = 1.15
    } else if((counter %% 4) == 2) {
      fc_val = 1.25
    } else if((counter %% 4) == 3) {
      fc_val = 1.5
    } else if((counter %% 4) == 0) {
      fc_val = 2
    }
    
    if(counter <= 4) {
      pi0_val = 0.8
    } else if(counter > 4 & counter <= 8) {
      pi0_val = 0.9
    } else if(counter > 8 & counter <= 12) {
      pi0_val = 0.95
    } else if(counter > 12 & counter <= 16) {
      pi0_val = 0.99
    }
    
    # Timer
    print("Start: Compute One Sample Size Parameter")
    tic("Compute One Sample Size Parameter")
    
    df_80_05 <- ssizeRNA_vary(nGenes = as.vector(data$nGenes), 
                              mu = as.vector(data$mu), 
                              disp = as.vector(data$disp), 
                              pi0 = pi0_val, fc = fc_val, 
                              m = m,
                              fdr = 0.05, power = 0.8, maxN = 30) %>%
      as.data.frame(.) %>%
      select(c(4, 5, 7)) %>%
      setNames(c("n", "actual_power", "crit_val")) %>%
      mutate(pi0 = pi0_val, fc = fc_val, fdr = 0.05, crit_thresh = 0.8)

    df_80_10 <- ssizeRNA_vary(nGenes = as.vector(data$nGenes), 
                              mu = as.vector(data$mu), 
                              disp = as.vector(data$disp), 
                              pi0 = pi0_val, fc = fc_val, 
                              m = m,
                              fdr = 0.1, power = 0.8, maxN = 30) %>%
      as.data.frame(.) %>%
      select(c(4, 5, 7)) %>%
      setNames(c("n", "actual_power", "crit_val")) %>%
      mutate(pi0 = pi0_val, fc = fc_val, fdr = 0.1, crit_thresh = 0.8)

    # Merge dfs
    complete_df <- list(df_80_05, df_80_10) %>%
      Reduce(function(x, y) 
        full_join(x, y), .)
    
    # Put dfs into temp list (will become sublists)
    power_list[[i]] <- complete_df %>% mutate(Dataset = dataset)
    toc()
  }
  
  # Convert future sub-lists to dataframe
  temp_power_df <- power_list %>% Reduce(function(x, y) full_join(x, y), .)

  # Make that dataframe a sublist in a larger list and return that list
  power_list_of_lists[[dataset]] <- temp_power_df
  power_list_of_lists[[dataset]]
}

# Do the calculations
power_list = list()
power_list_of_lists = list()
power_list[["carpenter"]] <- compute_power("carpenter")
power_list[["walker"]] <- compute_power("walker")

saveRDS(power_list, file = "power_list.RDS")



## START FROM HERE IF POWER_LIST OBJECT HAS BEEN CREATED ----------
power_list <- readRDS(file = "power_list.RDS")

# Merge lists into a dataframe
power_df <- power_list %>%
  Reduce(function(x, y) 
    full_join(x, y), .)

# Plot!
# Update and set ggplot theme globally
new_theme <- theme_bw(base_size = 14) %+replace%
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))
theme_set(new_theme)

p <- power_df %>%
  filter(fdr == 0.1,
         !c(Dataset == "carpenter" & fc == 2 & pi0 > 0.8),
         !c(Dataset == "walker" & fc == 2 & pi0 > 0.8)) %>% 
  mutate(`Proportion DEGs` = as.factor(1 - pi0), `Fold Change` = as.factor(fc), 
         FDR = as.factor(fdr), 
         Dataset = fct_relevel(str_to_sentence(Dataset), 
                               levels = c("Carpenter", "Walker"))) %>%
  ggplot(aes(x = n, y = actual_power, 
             shape = `Proportion DEGs`, color = `Fold Change`)) +
  geom_hline(yintercept = 0.8, lty = "dashed", color = "red") +
  geom_point(stroke = 0.5, size = 2) +
  scale_color_viridis(discrete = TRUE, begin = 0.92, end = 0, option = "G") +
  geom_line() +
  labs(x = "Sample Size per Group", y = "Power") +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 1), 
                     labels = c(0, rep("", 4), "5", rep("", 4), "10",
                                rep("", 4), "15", rep("", 4), "20",
                                rep("", 4), "25", rep("", 4), "30")) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # axis.text = element_text(size = 11),
        # axis.title.y = element_text(size = 12),
        # axis.title.x = element_text(size = 12),
        legend.position = "top", 
        plot.margin = grid::unit(c(0, 2, 2, 2), "mm")
        legend.box = "horizontal",
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 14)) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5),
         shape = guide_legend(title.position = "top", title.hjust = 0.5)) +
  facet_wrap(~ Dataset, scales = "free_x")

g <- ggplot_gtable(ggplot_build(p))
striprt <- which(grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name))
fills <- rep(c("#C148AD", "#44AA99"), 2)
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("C:/Users/Annika/Documents/Figures for Gene Conservation Project/Power_analysis_09-18-22.pdf", width = 8, height = 4)
grid::grid.draw(g)
dev.off()


