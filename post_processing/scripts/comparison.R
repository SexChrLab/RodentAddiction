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

# Statistics
library(rstatix)

# Other
library(clipr)


## SET PARAMETERS ----
# Change the tidyverse functions
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
arrange <- dplyr::arrange

options(dplyr.summarise.inform = FALSE)

# Reduce chances of scientific notation
options(scipen = 9999)

# Update and set ggplot theme globally
new_theme <- theme_bw(base_size = 14) %+replace%
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))
theme_set(new_theme)



## LOAD RESULTS FOR COMPARISON ---
carp_walk <- read_excel("C:/Users/Annika/Desktop/Vannan et al Supplemental Tables 05-20-2023.xlsx", 
                        sheet = 7)
engeln <- read_excel("C:/Users/Annika/Desktop/Vannan et al Supplemental Tables 05-20-2023.xlsx", 
                     sheet = 14)
compare_results <- full_join(carp_walk, engeln)


## LOAD IDS ---
main_dir <- "G:/Annika/Github Repositories/RodentAddiction/"
# Carpenter and Walker
all_genes_list <- readRDS(file = paste0(main_dir, "post_processing/results/other/all_genes_list.RDS"))
eng_genes_list <- readRDS(file = paste0(main_dir, "post_processing/results/other/my_genes_list.RDS"))

# Separate lists into multiple objects
list2env(all_genes_list, envir = .GlobalEnv)
list2env(eng_genes_list, envir = .GlobalEnv)


### ADD PRIORTITY LEVELS TO DATAFRAME ----
compare_results2 <- compare_results %>% 
  # select(1:4, 6, 9, 12, 14, 16, 18:24) %>%
  mutate(`BrainSpec and Region Priority Combination` = fct_relevel(`BrainSpec and Region Priority Combination`,
                                                                   levels = c("Low", "Medium", "High")),
         `CNS Expression Priority` = fct_relevel(case_when(`Mean CNS Expression (TPM)` > 0 & `Mean CNS Expression (TPM)` < 10 ~ "Low",
                                                           `Mean CNS Expression (TPM)` >= 10 & `Mean CNS Expression (TPM)` < 20 ~ "Medium",
                                                           `Mean CNS Expression (TPM)` >= 20 ~ "High"), 
                                                 levels = c("Low", "Medium", "High")),
         `Mouse dNdS Priority` = fct_relevel(`Mouse dNdS Priority`, levels = c("Low", "Medium", "High")),
         `Mouse Sequence Similarity Priority` = fct_relevel(case_when(`Mouse-Human Sequence Similarity` < 80 ~ "Low",
                                                                      `Mouse-Human Sequence Similarity` >= 80 & `Mouse-Human Sequence Similarity` < 90 ~ "Medium",
                                                                      `Mouse-Human Sequence Similarity` >= 90 ~ "High"),
                                                            levels = c("Low", "Medium", "High")),
         `Rat dNdS Priority` = fct_relevel(`Rat dNdS Priority`, levels = c("Low", "Medium", "High")),
         `Rat Sequence Similarity Priority` = fct_relevel(case_when(`Rat-Human Sequence Similarity` < 80 ~ "Low",
                                                                    `Rat-Human Sequence Similarity` >= 80 & `Rat-Human Sequence Similarity` < 90 ~ "Medium",
                                                                    `Rat-Human Sequence Similarity` >= 90 ~ "High"),
                                                          levels = c("Low", "Medium", "High")),
         `Developmental Conservation Priority` = fct_relevel(case_when(`Developmental Conservation` == "H" ~ "Low",
                                                                       `Developmental Conservation` %in% c("HR", "HM") ~ "Medium",
                                                                       `Developmental Conservation` == "HMR" ~ "High"),
                                                             levels = c("Low", "Medium", "High")),
         `Ensembl 99 Mouse Gene ID` = ifelse(`Ensembl 99 Mouse Gene ID` == "Same ID as Ensembl release 106", 
                                             `Ensembl 106 Mouse Gene ID`, `Ensembl 99 Mouse Gene ID`),
         Dataset = fct_relevel(Dataset, levels = c("Carpenter", "Walker", "Engeln"))) %>%
  filter((Dataset == "Carpenter" & `Human Gene ID` %in% sig_carp_hs_ids & `Ensembl 99 Mouse Gene ID` %in% sig_carp_mm_ids) |
           (Dataset == "Walker" & `Human Gene ID` %in% sig_walk_hs_ids & `Ensembl 99 Mouse Gene ID` %in% sig_walk_mm_ids) |
           (Dataset == "Engeln" & `Human Gene ID` %in% my_hs_ids & `Ensembl 99 Mouse Gene ID` %in% my_mm_ids))


#### Mouse-Human Conservation Priority ----
mouse_human_cons_df <- compare_results2 %>%
  select(Dataset, `Human Gene Symbol`, `Ensembl 99 Mouse Gene ID`,
         `Mouse dNdS Priority`, `Mouse Sequence Similarity Priority`) %>%
  unique() %>%
  group_by(Dataset, `Mouse Sequence Similarity Priority`) %>%
  mutate(SS_Cat_Num = length(`Mouse Sequence Similarity Priority`)) %>% 
  group_by(Dataset, `Mouse dNdS Priority`) %>%
  mutate(DNDS_Cat_Num = length(`Mouse dNdS Priority`)) %>% 
  ungroup() %>%
  pivot_longer(c(`Mouse dNdS Priority`, `Mouse Sequence Similarity Priority`), names_to = "Category", values_to = "Priority") %>%
  mutate(Category = fct_relevel(ifelse(Category == "Mouse dNdS Priority",
                                       "dN/dS", "Sequence Similarity"),
                                levels = c("Sequence Similarity", "dN/dS"))) %>%
  filter(Priority != "None", !is.na(Priority))
mouse_comp_plot <- mouse_human_cons_df %>%
  ggplot(aes(x = Dataset, fill = Priority)) +
  geom_bar(position = position_fill(), color = "black") +
  scale_fill_manual(values = c("#dc045d", "#f9e101", "#6db302"),
                    guide = guide_legend(reverse = TRUE)) +
  labs(x = "", y = "\nProportion of Ortholog Pairs", title = "Mouse-Human\nConservation Priority\n") +
  theme(legend.title = element_blank(), legend.position = "bottom", panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 13), plot.title = element_text(size = 14),
        panel.spacing = unit(2, "lines"), axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 13)) +
  facet_wrap(~ Category, ncol = 1, scales = "free")

###### Statistics ----
# Fisher's Exact Test
mouse_fisher_table <- table(mouse_human_cons_df$Dataset, mouse_human_cons_df$`Priority`, mouse_human_cons_df$`Category`)
mouse_seq_fisher_table <- mouse_fisher_table[, , 1]
mouse_dnds_fisher_table <- mouse_fisher_table[, , 2]

mouse_seq_fisher_test <- fisher.test(mouse_seq_fisher_table, simulate.p.value = TRUE)
mouse_seq_fisher_test # Sig

mouse_seq_fisher_test_cw <- fisher.test(mouse_seq_fisher_table[c(1, 2), ], simulate.p.value = FALSE)
mouse_seq_fisher_test_cw # Sig
mouse_seq_fisher_test_ce <- fisher.test(mouse_seq_fisher_table[c(1, 3), ], simulate.p.value = FALSE)
mouse_seq_fisher_test_ce # n.s.
mouse_seq_fisher_test_we <- fisher.test(mouse_seq_fisher_table[c(2, 3), ], simulate.p.value = FALSE)
mouse_seq_fisher_test_we # Sig

mouse_dnds_fisher_test <- fisher.test(mouse_dnds_fisher_table, simulate.p.value = TRUE)
mouse_dnds_fisher_test # n.s.


#### Rat-Human Conservation Priority ----
rat_human_cons_df <- compare_results2 %>%
  select(Dataset, `Human Gene Symbol`, `Ensembl 106 Rat Gene ID`,
         `Rat dNdS Priority`, `Rat Sequence Similarity Priority`) %>%
  unique() %>%
  group_by(Dataset, `Rat Sequence Similarity Priority`) %>%
  mutate(SS_Cat_Num = length(`Rat Sequence Similarity Priority`)) %>% 
  group_by(Dataset, `Rat dNdS Priority`) %>%
  mutate(DNDS_Cat_Num = length(`Rat dNdS Priority`)) %>% 
  ungroup() %>%
  pivot_longer(c(`Rat dNdS Priority`, `Rat Sequence Similarity Priority`), names_to = "Category", values_to = "Priority") %>%
  mutate(Category = fct_relevel(ifelse(Category == "Rat dNdS Priority",
                                       "dN/dS", "Sequence Similarity"),
                                levels = c("Sequence Similarity", "dN/dS"))) %>%
  filter(Priority != "None", !is.na(Priority))
rat_comp_plot <- rat_human_cons_df %>%
  ggplot(aes(x = Dataset, fill = Priority)) +
  geom_bar(position = position_fill(), color = "black") +
  scale_fill_manual(values = c("#dc045d", "#f9e101", "#6db302"),
                    guide = guide_legend(reverse = TRUE)) +
  labs(x = "", y = "\nProportion of Ortholog Pairs", title = "Rat-Human\nConservation Priority\n") +
  theme(legend.title = element_blank(), legend.position = "bottom", panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 13), plot.title = element_text(size = 14),
        panel.spacing = unit(2, "lines"), axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 13)) +
  facet_wrap(~ Category, ncol = 1, scales = "free")

###### Statistics ----
# Fisher's Exact Test
rat_fisher_table <- table(rat_human_cons_df$Dataset, rat_human_cons_df$`Priority`, rat_human_cons_df$`Category`)
rat_seq_fisher_table <- rat_fisher_table[, , 1]
rat_dnds_fisher_table <- rat_fisher_table[, , 2]

rat_seq_fisher_test <- fisher.test(rat_seq_fisher_table, simulate.p.value = TRUE)
rat_seq_fisher_test # Sig

rat_seq_fisher_test_cw <- fisher.test(rat_seq_fisher_table[c(1, 2), ], simulate.p.value = FALSE)
rat_seq_fisher_test_cw # Sig
rat_seq_fisher_test_ce <- fisher.test(rat_seq_fisher_table[c(1, 3), ], simulate.p.value = FALSE)
rat_seq_fisher_test_ce # Sig
rat_seq_fisher_test_we <- fisher.test(rat_seq_fisher_table[c(2, 3), ], simulate.p.value = FALSE)
rat_seq_fisher_test_we # Sig

rat_dnds_fisher_test <- fisher.test(rat_dnds_fisher_table, simulate.p.value = TRUE)
rat_dnds_fisher_test # n.s.



#### Developmental Conservation ----
dev_df <- compare_results2 %>%
  select(Dataset, `Human Gene Symbol`, `Developmental Conservation Priority`, `Developmental Conservation`) %>%
  unique() %>%
  filter(!(`Developmental Conservation` %in% c("H*", "HM*", "HR*")),
         !is.na(`Developmental Conservation`)) %>%
  unique()
comp_dev_plot <-dev_df %>%
  ggplot(aes(x = `Dataset`, fill = `Developmental Conservation Priority`)) +
  geom_bar(position = position_fill(), color = "black") + 
  scale_fill_manual(values = c("#dc045d", "#f9e101", "#6db302"),
                    guide = guide_legend(reverse = TRUE)) +
  labs(x = "", y = "\nProportion of Genes", title = "Developmental Conservation\nPriority\n") +
  theme(legend.title = element_blank(), legend.position = "bottom",
        panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 11), 
        plot.title = element_text(size = 12.5))

###### Statistics ----
# Fisher's Exact Test
dev_fisher_table <- table(dev_df$Dataset, dev_df$`Developmental Conservation Priority`)
dev_fisher_test <- fisher.test(dev_fisher_table, simulate.p.value = FALSE)
dev_fisher_test # n.s.



#### Brain Specificity/Expression Priority ----
brain_df <- compare_results2 %>%
  select(Dataset, `Human Gene Symbol`, `BrainSpec and Region Priority Combination`, `CNS Expression Priority`) %>%
  unique() %>%
  filter(`BrainSpec and Region Priority Combination` != "None") %>%
  pivot_longer(c(`CNS Expression Priority`, `BrainSpec and Region Priority Combination`), 
               names_to = "Category", values_to = "Priority") %>%
  mutate(Category = fct_relevel(ifelse(Category == "BrainSpec and Region Priority Combination", 
                                       "Specificity/Region Enrichment", "Mean CNS Expression"),
                                levels = c("Specificity/Region Enrichment", "Mean CNS Expression")))
comp_brain_plot <- brain_df %>%
  ggplot(aes(x = Dataset, fill = Priority)) +
  geom_bar(position = position_fill(), color = "black") + 
  scale_fill_manual(values = c("#dc045d", "#f9e101", "#6db302"),
                    guide = guide_legend(reverse = TRUE)) +
  labs(x = "", y = "\n\nProportion of Genes", title = "Brain Specificity and Expression Priority\n") +
  theme(legend.title = element_blank(), legend.position = "bottom", panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 11), plot.title = element_text(size = 13)) +
  facet_wrap(~ Category)

###### Statistics ----
# Fisher's Exact Test
brain_fisher_table <- table(brain_df$Dataset, brain_df$`Priority`, brain_df$`Category`)
spec_fisher_table <- brain_fisher_table[, , 1]
mean_fisher_table <- brain_fisher_table[, , 2]

spec_fisher_test <- fisher.test(spec_fisher_table, simulate.p.value = TRUE)
spec_fisher_test # Sig

spec_fisher_test_cw <- fisher.test(spec_fisher_table[c(1, 2), ], simulate.p.value = FALSE)
spec_fisher_test_cw # Sig
spec_fisher_test_ce <- fisher.test(spec_fisher_table[c(1, 3), ], simulate.p.value = FALSE)
spec_fisher_test_ce # Sig
spec_fisher_test_we <- fisher.test(spec_fisher_table[c(2, 3), ], simulate.p.value = FALSE)
spec_fisher_test_we # Sig

mean_fisher_test <- fisher.test(mean_fisher_table, simulate.p.value = TRUE)
mean_fisher_test # Sig

mean_fisher_test_cw <- fisher.test(mean_fisher_table[c(1, 2), ], simulate.p.value = FALSE)
mean_fisher_test_cw # n.s.
mean_fisher_test_ce <- fisher.test(mean_fisher_table[c(1, 3), ], simulate.p.value = FALSE)
mean_fisher_test_ce # Sig
mean_fisher_test_we <- fisher.test(mean_fisher_table[c(2, 3), ], simulate.p.value = FALSE)
mean_fisher_test_we # Sig


## SAVE PLOTS ----
mouse_comp_plot
rat_comp_plot
comp_dev_plot
comp_brain_plot

svglite("C:/Users/Annika/Documents/comp_rat_mouse_cons_052023.svg", width = 9, height = 6.5)
ggarrange(mouse_comp_plot, rat_comp_plot, common.legend = TRUE, legend = "right")
dev.off()

svglite("C:/Users/Annika/Documents/comp_dev_052023.svg", width = 3.2, height = 3.5)
comp_dev_plot + theme(legend.position = "none")
dev.off()

svglite("C:/Users/Annika/Documents/comp_brain_052023.svg", width = 5.5, height = 3.5)
comp_brain_plot + theme(legend.position = "none")
dev.off()



