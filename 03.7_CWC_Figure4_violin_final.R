# ==============================================================================
# Script Name:  03.7_CWC_Figure4_violin_final.R
# Description:  Figure 4 (size-dependent total and partial mortality), violin
#               version. Direct revision of 03.5_CWC_shrinkage_size_selectivity_
#               Updated.R per reviewer comments AR2/AR3 (curves/violins, not
#               boxplots). NO STATISTICS CHANGED: same fate classification,
#               same Wilcoxon rank-sum tests, same significance-star logic,
#               same n-labels as the original -- this only swaps the display
#               geometry from boxplot to violin (+ inner boxplot for median/IQR,
#               matching the convention already used in your Figure 3D growth
#               panel, so the two figures read consistently).
#
#               The quadratic-GLM "bathtub curve" alternative was tested
#               separately (Script 03.6) and dropped: only M. oculata shrinkage
#               (n=89 events) supported a quadratic term, and even that needs
#               its curvature direction checked before calling it a bathtub.
#               Everything else was flat, monotonic, or too sparse to test. The
#               non-parametric Wilcoxon approach you already have is the more
#               defensible choice given that, so this script keeps it exactly
#               as-is and only changes how the same tested distributions are
#               drawn.
#
#               REVISED (v2): fixed a real bug from the first version -- the
#               fate category order (Survived/Partial/Total) was drifting
#               between species panels because of free y-scales combined with
#               two jitter layers split across filtered subsets. Now uses one
#               shared y-scale and one jitter layer, and zero-count cells
#               (Primnoa Partial Mortality) show explicitly as "n = 0" rather
#               than disappearing from the panel.
# Dependencies: readxl, dplyr, stringr, ggplot2, ragg
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(ragg)

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & FILE PATHS
# ------------------------------------------------------------------------------
base_dir <- "D:/PhD_Data(Large)/Submission_Dataset/"
input_xlsx_path <- file.path(base_dir, "Master_CWC_Tracked_MDC.xlsx")

output_dir <- file.path(base_dir, "Figures")
out_violin <- file.path(output_dir, "Figure4_Size_Dependent_Mortality_Violin.png")
out_audit_csv <- file.path(output_dir, "Figure4_Wilcoxon_Audit_Table.csv")

group_levels <- c("Madrepora oculata", "Desmophyllum pertusum", "Primnoa msp.")

fate_colors <- c(
  "Survived (Stable/Grew)"        = "steelblue",
  "Partial Mortality (Shrinkage)" = "darkorange",
  "Total Mortality"               = "darkred"
)

# EDIT: below this many points, a KDE violin shape carries no real information
# (D. pertusum shrinkage is n=2) -- draws jittered points only, no violin, so
# the figure never implies a distribution shape from 2 colonies. Set to 0 to
# force violins everywhere regardless of n, if you'd rather have visual
# consistency across every cell than this caveat.
MIN_N_FOR_VIOLIN <- 5

# ------------------------------------------------------------------------------
# 3. LOAD DATASET & CLASSIFY DEMOGRAPHIC FATES (unchanged from 03.5)
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

df_mort <- Master_CWC_Tracked_MDC_Age %>%
  filter(present_2015 == TRUE, present_2022 == FALSE) %>%
  mutate(mortality_type = "Total Mortality")

df_shrink <- Master_CWC_Tracked_MDC_Age %>%
  filter(
    present_2015 == TRUE,
    present_2022 == TRUE,
    change_area < 0,
    detectable_change == TRUE
  ) %>%
  mutate(mortality_type = "Partial Mortality (Shrinkage)")

df_survived <- Master_CWC_Tracked_MDC_Age %>%
  filter(
    present_2015 == TRUE,
    present_2022 == TRUE,
    !TagLab.Genet.Id %in% df_shrink$TagLab.Genet.Id
  ) %>%
  mutate(mortality_type = "Survived (Stable/Grew)")

df_combined_mdc <- bind_rows(df_survived, df_shrink, df_mort) %>%
  filter(Species != "Coral Recruit") %>%
  mutate(
    Species = factor(standardise_species(Species), levels = group_levels),
    mortality_type = factor(
      mortality_type,
      levels = c("Total Mortality", "Partial Mortality (Shrinkage)", "Survived (Stable/Grew)")
    )
  )

df_logistic_total <- Master_CWC_Tracked_MDC_Age %>%
  filter(present_2015 == TRUE, Species != "Coral Recruit") %>%
  mutate(
    Species = factor(standardise_species(Species), levels = group_levels),
    is_dead = if_else(present_2022 == FALSE, 1, 0)
  )

# ------------------------------------------------------------------------------
# 4. SIZE-SELECTIVITY TESTING (WILCOXON RANK-SUM) -- unchanged from 03.5
# ------------------------------------------------------------------------------
p_to_stars <- function(p) {
  as.character(cut(p, breaks = c(-Inf, .001, .01, .05, Inf), labels = c("***", "**", "*", "ns")))
}

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

df_oculata_persistent <- df_combined_mdc %>%
  filter(Species == "Madrepora oculata", mortality_type != "Total Mortality")

cat("\n=== M. oculata Shrinkage Size-Selectivity ===\n")
wilcox_shrink_oculata <- wilcox.test(Area_2015 ~ mortality_type, data = df_oculata_persistent)
cat("Wilcoxon W =", wilcox_shrink_oculata$statistic, "| p-value =", wilcox_shrink_oculata$p.value, "\n")

oculata_shrink_sig <- p_to_stars(wilcox_shrink_oculata$p.value)
if (oculata_shrink_sig == "ns") oculata_shrink_sig <- ""

# ------------------------------------------------------------------------------
# 5. FIGURE ANNOTATION DATA -- REVISED from 03.5
# ------------------------------------------------------------------------------
# FIX: .drop = FALSE forces every Species x fate combination to appear,
# including Primnoa x Partial Mortality (Shrinkage), which has zero rows.
# Without this, that combination silently has no label at all and the row
# just vanishes from the figure (same silent-drop problem your Script 03
# comments flag for occlusion counts) -- here it explicitly reads "n = 0".
#
# x_pos is now the SPECIES-level max area (not per-fate-group max), because
# max() on a zero-row group returns -Inf, which would put the "n = 0" label
# at an undefined position and break that panel's layout.
species_max_area <- df_combined_mdc %>%
  group_by(Species) %>%
  summarise(species_max = max(Area_2015, na.rm = TRUE), .groups = "drop")

n_labels <- df_combined_mdc %>%
  group_by(Species, mortality_type, .drop = FALSE) %>%
  summarise(n_count = n(), .groups = "drop") %>%
  left_join(species_max_area, by = "Species") %>%
  mutate(x_pos = species_max * 1.7) %>%
  left_join(total_mortality_stats %>% select(Species, sig_stars), by = "Species") %>%
  mutate(
    sig = case_when(
      mortality_type == "Total Mortality" & !is.na(sig_stars) & sig_stars != "ns"     ~ sig_stars,
      Species == "Madrepora oculata" & mortality_type == "Partial Mortality (Shrinkage)" ~ oculata_shrink_sig,
      TRUE ~ ""
    ),
    label_text = if_else(sig == "", paste0("n = ", n_count), paste0("n = ", n_count, " ", sig))
  ) %>%
  select(-sig_stars, -species_max)

# ------------------------------------------------------------------------------
# 6. BUILD FIGURE 4 -- VIOLIN + INNER BOXPLOT + JITTER
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

# Split by whether each Species x fate cell clears MIN_N_FOR_VIOLIN, so the
# n=2 D. pertusum shrinkage cell never gets a fabricated-looking violin shape.
cell_n <- df_combined_mdc %>% count(Species, mortality_type, name = "n_count")
df_combined_mdc <- df_combined_mdc %>%
  left_join(cell_n, by = c("Species", "mortality_type")) %>%
  mutate(violin_eligible = n_count >= MIN_N_FOR_VIOLIN)

df_violin_part <- df_combined_mdc %>% filter(violin_eligible)
df_points_only <- df_combined_mdc %>% filter(!violin_eligible)

cat("\nCells drawn as points-only (n <", MIN_N_FOR_VIOLIN, ", no violin shape):\n")
print(as.data.frame(cell_n %>% filter(n_count < MIN_N_FOR_VIOLIN, n_count > 0)), row.names = FALSE)

# FIX: scales = "free_y" + two separate jitter layers pulling from different
# filtered subsets is what caused D. pertusum's rows to render in a different
# order (Partial/Survived/Total) than M. oculata's and Primnoa's
# (Survived/Partial/Total) for the SAME factor. Fixed shared scales
# (dropping "free_y" entirely) force one deterministic category order across
# every panel, and a single jitter layer on the full data removes the
# order-dependent split that produced the inconsistency.
p_violin <- ggplot(df_combined_mdc, aes(x = Area_2015, y = mortality_type)) +
  # Violin only for cells with enough points to make a KDE meaningful
  geom_violin(
    data = df_violin_part, aes(fill = mortality_type),
    alpha = 0.55, color = "black", linewidth = 0.4, width = 0.7, trim = TRUE
  ) +
  # Thin inner boxplot for median/IQR -- matches your Figure 3D growth panel convention
  geom_boxplot(
    data = df_violin_part, width = 0.12, fill = "white", outlier.shape = NA,
    alpha = 0.95, color = "black", linewidth = 0.35
  ) +
  # ONE jitter layer, always on the full dataset -- every colony gets a point
  # regardless of whether its cell also gets a violin drawn on top of it.
  geom_jitter(
    aes(color = mortality_type),
    height = 0.15, width = 0, alpha = 0.4, size = 1.4
  ) +
  geom_text(
    data = n_labels, aes(x = x_pos, y = mortality_type, label = label_text),
    inherit.aes = FALSE, hjust = 0, size = 3.5, color = "black"
  ) +
  facet_wrap(~ Species, ncol = 1) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.30))) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_manual(values = fate_colors) +
  scale_color_manual(values = fate_colors) +
  base_theme +
  theme(panel.grid.major.y = element_blank(), axis.title.y = element_blank()) +
  labs(x = x_axis_label)

# ------------------------------------------------------------------------------
# 7. EXPORT
# ------------------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nExporting Figure 4 (violin)...\n")
ggsave(filename = out_violin, plot = p_violin, width = 8, height = 9, units = "in",
       dpi = 600, device = ragg::agg_png, bg = "white")

write.csv(total_mortality_stats, out_audit_csv, row.names = FALSE)

cat("\n======================================================================\n")
cat("Figure 4 exported to:\n ", out_violin, "\n")
cat("Wilcoxon audit table exported to:\n ", out_audit_csv, "\n")
cat("======================================================================\n")
