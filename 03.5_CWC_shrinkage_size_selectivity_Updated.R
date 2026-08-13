# ==============================================================================
# Script Name:  03.5_CWC_shrinkage_size_selectivity_Updated.R
# Description:  Loads the age-integrated master tracked CWC dataset, classifies 
#               colonies into mutually exclusive demographic fates (Total Mortality, 
#               Partial Mortality, Survived/Grew), tests size-selectivity via 
#               Wilcoxon rank-sum tests, and builds four publication-ready 
#               visual alternatives (Boxplot, Density, Logistic Total Mortality, 
#               Logistic Partial Mortality).
# Dependencies: readxl, dplyr, stringr, ggplot2, ragg, patchwork
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(ragg)
library(patchwork)

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & FILE PATHS
# ------------------------------------------------------------------------------
base_dir <- "D:/PhD_Data(Large)/Submission_Dataset/"
input_xlsx_path <- file.path(base_dir, "Master_CWC_Tracked_MDC.xlsx")

output_dir <- file.path(base_dir, "Figures")
out_boxplot         <- file.path(output_dir, "Figure_Size_Selectivity_1_Boxplot.png")
out_density         <- file.path(output_dir, "Figure_Size_Selectivity_2_Density.png")
out_logistic_total   <- file.path(output_dir, "Figure_Size_Selectivity_3_Logistic_TotalMortality.png")
out_logistic_shrink  <- file.path(output_dir, "Figure_Size_Selectivity_4_Logistic_PartialMortality.png")

group_levels <- c("Madrepora oculata", "Desmophyllum pertusum", "Primnoa msp.")

# Master color palette for mutually exclusive fates
fate_colors <- c(
  "Survived (Stable/Grew)"        = "steelblue",
  "Partial Mortality (Shrinkage)" = "darkorange",
  "Total Mortality"               = "darkred"
)

# ------------------------------------------------------------------------------
# 3. LOAD DATASET & CLASSIFY DEMOGRAPHIC FATES (MUTUALLY EXCLUSIVE)
# ------------------------------------------------------------------------------
cat("Loading age-integrated master dataset...\n")
Master_CWC_Tracked_MDC_Age <- read_excel(input_xlsx_path)

standardise_species <- function(x) {
  case_when(
    str_detect(x, "Primnoa") ~ "Primnoa msp.",
    str_detect(x, "pertus")  ~ "Desmophyllum pertusum",
    TRUE                     ~ "Madrepora oculata"
  )
}

# 1. Total mortality: present in 2015, absent in 2022
df_mort <- Master_CWC_Tracked_MDC_Age %>%
  filter(present_2015 == TRUE, present_2022 == FALSE) %>%
  mutate(mortality_type = "Total Mortality")

# 2. Partial mortality (Shrinkage): Survived AND lost area exceeding MDC95
df_shrink <- Master_CWC_Tracked_MDC_Age %>%
  filter(
    present_2015 == TRUE,
    present_2022 == TRUE,
    change_area < 0,
    detectable_change == TRUE
  ) %>%
  mutate(mortality_type = "Partial Mortality (Shrinkage)")

# 3. Survived (Stable/Grew): All surviving colonies excluding those in df_shrink
# (Option A: Uses ID set subtraction to prevent NA evaluation dropouts)
df_survived <- Master_CWC_Tracked_MDC_Age %>%
  filter(
    present_2015 == TRUE,
    present_2022 == TRUE,
    !TagLab.Genet.Id %in% df_shrink$TagLab.Genet.Id
  ) %>%
  mutate(mortality_type = "Survived (Stable/Grew)")

# Combine, drop recruits, standardise species and factor levels
df_combined_mdc <- bind_rows(df_survived, df_shrink, df_mort) %>%
  filter(Species != "Coral Recruit") %>%
  mutate(
    Species = factor(standardise_species(Species), levels = group_levels),
    mortality_type = factor(
      mortality_type,
      levels = c("Total Mortality", "Partial Mortality (Shrinkage)", "Survived (Stable/Grew)")
    )
  )

# Data frame strictly for binary total mortality logistic regression (0 = Survived, 1 = Dead)
df_logistic_total <- Master_CWC_Tracked_MDC_Age %>%
  filter(present_2015 == TRUE, Species != "Coral Recruit") %>%
  mutate(
    Species = factor(standardise_species(Species), levels = group_levels),
    is_dead = if_else(present_2022 == FALSE, 1, 0)
  )

# Data frame strictly for binary partial mortality logistic regression (0 = Stable/Grew, 1 = Shrunk)
df_logistic_shrink <- Master_CWC_Tracked_MDC_Age %>%
  filter(present_2015 == TRUE, present_2022 == TRUE, Species != "Coral Recruit") %>%
  mutate(
    Species   = factor(standardise_species(Species), levels = group_levels),
    is_shrunk = if_else(TagLab.Genet.Id %in% df_shrink$TagLab.Genet.Id, 1, 0)
  )

# ------------------------------------------------------------------------------
# 4. SIZE-SELECTIVITY TESTING (WILCOXON RANK-SUM)
# ------------------------------------------------------------------------------
p_to_stars <- function(p) {
  as.character(cut(p, breaks = c(-Inf, .001, .01, .05, Inf), labels = c("***", "**", "*", "ns")))
}

# --- Test A: Total mortality — Died vs Survived (All survivors combined) ---
cat("\n=== Total Mortality Size-Selectivity (Wilcoxon Rank-Sum) ===\n")
total_mortality_stats <- df_logistic_total %>%
  group_by(Species) %>%
  summarise(
    n_died               = sum(is_dead == 1),
    n_survived           = sum(is_dead == 0),
    median_area_died     = median(Area_2015[is_dead == 1], na.rm = TRUE),
    median_area_survived = median(Area_2015[is_dead == 0], na.rm = TRUE),
    W_stat               = wilcox.test(Area_2015 ~ is_dead)$statistic,
    p_value              = wilcox.test(Area_2015 ~ is_dead)$p.value,
    .groups              = "drop"
  ) %>%
  mutate(sig_stars = p_to_stars(p_value))
print(as.data.frame(total_mortality_stats), row.names = FALSE)

# --- Test B: M. oculata partial mortality / shrinkage — Shrunk vs Stable/Grew ---
df_oculata_persistent <- df_combined_mdc %>%
  filter(Species == "Madrepora oculata", mortality_type != "Total Mortality")

cat("\n=== M. oculata Shrinkage Size-Selectivity ===\n")
wilcox_shrink_oculata <- wilcox.test(Area_2015 ~ mortality_type, data = df_oculata_persistent)
cat("Wilcoxon W =", wilcox_shrink_oculata$statistic, "| p-value =", wilcox_shrink_oculata$p.value, "\n")

oculata_shrink_sig <- p_to_stars(wilcox_shrink_oculata$p.value)
if (oculata_shrink_sig == "ns") oculata_shrink_sig <- ""

# ------------------------------------------------------------------------------
# 5. FIGURE ANNOTATION DATA 
# ------------------------------------------------------------------------------
n_labels <- df_combined_mdc %>%
  group_by(Species, mortality_type) %>%
  summarise(n_count = n(), x_pos = max(Area_2015, na.rm = TRUE) * 1.7, .groups = "drop") %>%
  left_join(total_mortality_stats %>% select(Species, sig_stars), by = "Species") %>%
  mutate(
    sig = case_when(
      mortality_type == "Total Mortality" & sig_stars != "ns"                        ~ sig_stars,
      Species == "Madrepora oculata" & mortality_type == "Partial Mortality (Shrinkage)" ~ oculata_shrink_sig,
      TRUE ~ ""
    ),
    label_text = if_else(sig == "", paste0("n = ", n_count), paste0("n = ", n_count, " ", sig))
  ) %>%
  select(-sig_stars)

# ------------------------------------------------------------------------------
# 6. BUILD VISUAL ALTERNATIVES
# ------------------------------------------------------------------------------
base_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position    = "none",
    strip.text         = element_text(face = "bold.italic", size = 12, color = "black"),
    axis.text          = element_text(size = 10, color = "black"),
    axis.title         = element_text(size = 11, face = "bold", color = "black"),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing      = unit(0.8, "lines"),
    plot.margin        = margin(t = 8, r = 15, b = 8, l = 8, unit = "pt")
  )

x_axis_label <- expression(paste("Initial Planar Area 2015 (cm"^2*"; log scale)"))

# Option 1: Boxplot (Mutually exclusive groups)
cat("Building Alternative 1: Boxplots...\n")
p_boxplot <- ggplot(df_combined_mdc, aes(x = Area_2015, y = mortality_type, fill = mortality_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.5, color = "black", linewidth = 0.4) +
  geom_jitter(aes(color = mortality_type), height = 0.15, width = 0, alpha = 0.4, size = 1.4) +
  geom_text(
    data = n_labels, aes(x = x_pos, y = mortality_type, label = label_text),
    inherit.aes = FALSE, hjust = 0, size = 3.5, color = "black"
  ) +
  facet_wrap(~ Species, ncol = 1, scales = "free_y") +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.30))) +
  scale_fill_manual(values = fate_colors) +
  scale_color_manual(values = fate_colors) +
  base_theme +
  theme(panel.grid.major.y = element_blank(), axis.title.y = element_blank()) +
  labs(x = x_axis_label)

# Option 2: Density Distributions
cat("Building Alternative 2: Density Distributions...\n")
p_density <- ggplot(df_combined_mdc, aes(x = Area_2015, fill = mortality_type, color = mortality_type)) +
  geom_density(alpha = 0.4, linewidth = 0.6) +
  facet_wrap(~ Species, ncol = 1, scales = "free_y") +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.05))) +
  scale_fill_manual(values = fate_colors, name = "Demographic Fate") +
  scale_color_manual(values = fate_colors, name = "Demographic Fate") +
  base_theme +
  theme(legend.position = "bottom") +
  labs(x = x_axis_label, y = "Density Probability")

# Option 3: Logistic Regression Curve (Total Mortality)
cat("Building Alternative 3: Logistic Regression (Total Mortality)...\n")
p_logistic_total <- ggplot(df_logistic_total, aes(x = Area_2015, y = is_dead)) +
  geom_jitter(height = 0.04, width = 0, alpha = 0.25, size = 1.5, color = "darkred") +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), 
              color = "black", fill = "grey50", linewidth = 0.8) +
  facet_wrap(~ Species, ncol = 1) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0% (Survived)", "50%", "100% (Died)")) +
  base_theme +
  labs(x = x_axis_label, y = "Probability of Total Mortality")

# Option 4: Logistic Regression Curve (Partial Mortality / Shrinkage)
cat("Building Alternative 4: Logistic Regression (Partial Mortality)...\n")
p_logistic_shrink <- ggplot(df_logistic_shrink, aes(x = Area_2015, y = is_shrunk)) +
  geom_jitter(height = 0.04, width = 0, alpha = 0.25, size = 1.5, color = "darkorange") +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), 
              color = "black", fill = "grey50", linewidth = 0.8) +
  facet_wrap(~ Species, ncol = 1) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0% (Stable/Grew)", "50%", "100% (Shrunk)")) +
  base_theme +
  labs(x = x_axis_label, y = "Probability of Partial Mortality (Shrinkage)")

# ------------------------------------------------------------------------------
# 7. EXPORT HIGH-RESOLUTION FIGURES
# ------------------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Exporting Alternative 1 (Boxplot)...\n")
ggsave(filename = out_boxplot, plot = p_boxplot, width = 8, height = 9, units = "in", dpi = 600, device = ragg::agg_png, bg = "white")

cat("Exporting Alternative 2 (Density)...\n")
ggsave(filename = out_density, plot = p_density, width = 8, height = 9, units = "in", dpi = 600, device = ragg::agg_png, bg = "white")

cat("Exporting Alternative 3 (Logistic Total Mortality)...\n")
ggsave(filename = out_logistic_total, plot = p_logistic_total, width = 8, height = 9, units = "in", dpi = 600, device = ragg::agg_png, bg = "white")

cat("Exporting Alternative 4 (Logistic Partial Mortality)...\n")
ggsave(filename = out_logistic_shrink, plot = p_logistic_shrink, width = 8, height = 9, units = "in", dpi = 600, device = ragg::agg_png, bg = "white")

cat("\n============================================================\n")
cat("Figures successfully exported to:\n")
cat(" 1.", out_boxplot, "\n")
cat(" 2.", out_density, "\n")
cat(" 3.", out_logistic_total, "\n")
cat(" 4.", out_logistic_shrink, "\n")
cat("============================================================\n")
