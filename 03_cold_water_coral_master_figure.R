Here is your complete, updated `03_CWC_master_figure1.R` script, maintained in your exact structure, commenting style, and variable naming format, with the panel clipping and headroom fixes fully incorporated:

```R
# ==============================================================================
# Script Name:  03_CWC_master_figure1.R
# Description:  Loads master tracked cold-water coral dataset (from Script 02),
#               performs Kolmogorov-Smirnov and Kruskal-Wallis/Dunn statistical 
#               testing, computes Wilson score confidence intervals for dynamic 
#               rates, and builds a high-resolution publication-ready 6-panel 
#               composite figure (Panels A-F).
# Dependencies: readxl, dplyr, tidyr, stringr, ggplot2, patchwork, scales, 
#               dunn.test, ggsignif, Hmisc, ragg
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(scales)
library(dunn.test)
library(ggsignif)
library(Hmisc)
library(ragg)

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & FILE PATHS
# ------------------------------------------------------------------------------
# Base directory path
base_dir <- "D:/PhD_Data(Large)/Submission_Dataset/"

input_xlsx_path <- file.path(base_dir, "Master_CWC_Tracked_MDC.xlsx")
output_fig_path <- file.path(base_dir, "Figures/Figure1_Perfect_CWC_Combined.png")

# Global aesthetic & plotting parameters
FONT_FAMILY  <- "sans"
group_levels <- c("Madrepora oculata", "Desmophyllum pertusum", "Primnoa msp.")

# Master color scheme
species_palette <- c(
  "Madrepora oculata"     = "#722082",  # Deep Purple
  "Desmophyllum pertusum" = "#B63679",  # Magenta
  "Primnoa msp."          = "#fb9f3a"   # Coral / Orange
)

c_2015_fill <- "#E5E5E5"
c_2015_line <- "#404040"

# ------------------------------------------------------------------------------
# 3. LOAD TRACKED DATASET & PREPARE ROW 1 (AREA DENSITIES)
# ------------------------------------------------------------------------------
cat("Loading master MDC dataset for statistical modeling and plotting...\n")
df <- read_excel(input_xlsx_path)

# Standardize species groupings across all downstream panels
df_madrepora <- df %>% 
  filter(Species == "Madrepora oculata") %>% 
  mutate(Group = "Madrepora oculata")

df_pertusum  <- df %>% 
  filter(Species %in% c("D. pertusum", "Desmophyllum pertusum")) %>% 
  mutate(Group = "Desmophyllum pertusum")

df_primnoa   <- df %>% 
  filter(Species %in% c("Primnoa msp.5", "Primnoa msp.1")) %>% 
  mutate(Group = "Primnoa msp.")

# Reshape into long-format for density probability distributions
master_long <- bind_rows(df_madrepora, df_pertusum, df_primnoa) %>%
  select(TagLab.Genet.Id, Group, Area_2015, Area_2022) %>%
  pivot_longer(cols = c(Area_2015, Area_2022), names_to = "Year", values_to = "Area") %>%
  mutate(Year = str_remove(Year, "Area_")) %>% 
  filter(Area > 0) %>% 
  mutate(
    Group    = factor(Group, levels = group_levels),
    Log_Area = log(Area)
  )

cat("Calculating Row 1 summary statistics & running Kolmogorov-Smirnov tests...\n")
stats_summary_row1 <- master_long %>%
  group_by(Group) %>%
  summarise(
    n_15    = sum(Year == "2015", na.rm = TRUE),
    n_22    = sum(Year == "2022", na.rm = TRUE),
    med_15  = median(Area[Year == "2015"], na.rm = TRUE),
    med_22  = median(Area[Year == "2022"], na.rm = TRUE),
    
    ks_stat = if(sum(Year == "2015") > 1 & sum(Year == "2022") > 1) {
                ks.test(Area[Year == "2015"], Area[Year == "2022"])$statistic
              } else { NA_real_ },
              
    ks_p    = if(sum(Year == "2015") > 1 & sum(Year == "2022") > 1) {
                ks.test(Area[Year == "2015"], Area[Year == "2022"])$p.value
              } else { NA_real_ },
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
  is.na(ks_p)  ~ "p = NA",
  ks_p < 0.001 ~ "p < 0.001",
  TRUE         ~ paste0("p = ", round(ks_p, 3))
),
    ks_label = if_else(is.na(ks_stat), "D == NA", paste0("D == ", round(ks_stat, 2)))
  )

# ------------------------------------------------------------------------------
# 4. PREPARE ROW 2 DATA & RUN POST-HOC COMPARISONS
# ------------------------------------------------------------------------------
cat("Preparing Row 2 demographic parameters & calculating Wilson CIs...\n")
growth_data <- df %>%
  mutate(Species = case_when(
    Species == "Madrepora oculata" ~ "Madrepora oculata",
    Species %in% c("D. pertusum", "Desmophyllum pertusum") ~ "Desmophyllum pertusum",
    Species %in% c("Primnoa msp.5", "Primnoa msp.1") ~ "Primnoa msp.",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Species)) %>%
  mutate(
    Species          = factor(Species, levels = group_levels),
    state_transition = tolower(state_transition)
  )

df_growth <- growth_data %>% filter(state_transition == "persistence")

# Demographic summaries and Wilson score interval estimation
summary_stats_row2 <- growth_data %>%
  group_by(Species) %>%
  summarise(
    n_2015        = sum(state_transition %in% c("persistence", "mortality")),
    n_2022        = sum(state_transition %in% c("persistence", "recruitment")),
    n_persistence = sum(state_transition == "persistence"),
    n_mortality   = sum(state_transition == "mortality"),
    n_recruitment = sum(state_transition == "recruitment"),
    mean_growth   = mean(Annual_Area_Change[state_transition == "persistence"], na.rm = TRUE),
    mean_mdc      = mean(MDC95_Annual[state_transition == "persistence"], na.rm = TRUE),
    .groups       = "drop"
  )

mort_ci <- Hmisc::binconf(summary_stats_row2$n_mortality,   summary_stats_row2$n_2015, method = "wilson") * 100
rec_ci  <- Hmisc::binconf(summary_stats_row2$n_recruitment, summary_stats_row2$n_2022, method = "wilson") * 100

summary_stats_row2 <- summary_stats_row2 %>%
  mutate(
    mortality_rate   = mort_ci[, "PointEst"],
    mort_ci_lower    = mort_ci[, "Lower"],
    mort_ci_upper    = mort_ci[, "Upper"],
    recruitment_rate = rec_ci[, "PointEst"],
    rec_ci_lower     = rec_ci[, "Lower"],
    rec_ci_upper     = rec_ci[, "Upper"]
  )

mdc_values <- summary_stats_row2 %>% select(Species, mean_mdc)

# Establish sample size dynamic labels
p1_labels <- setNames(paste0("(n = ", summary_stats_row2$n_persistence, ")"), summary_stats_row2$Species)
p2_labels <- setNames(paste0("(n = ", summary_stats_row2$n_2022, ")"), summary_stats_row2$Species) 
p3_labels <- setNames(paste0("(n = ", summary_stats_row2$n_2015, ")"), summary_stats_row2$Species) 

cat("Executing Kruskal-Wallis & Benjamini-Hochberg post-hoc tests...\n")
kw_test  <- kruskal.test(Annual_Area_Change ~ Species, data = df_growth)
dunn_res <- dunn.test::dunn.test(
  x      = df_growth$Annual_Area_Change,
  g      = df_growth$Species,
  method = "bh",
  kw     = TRUE
)

dunn_lookup <- setNames(dunn_res$P.adjusted, dunn_res$comparisons)
comparison_list <- list(
  c("Madrepora oculata", "Desmophyllum pertusum"),
  c("Desmophyllum pertusum", "Primnoa msp."),
  c("Madrepora oculata", "Primnoa msp.")
)

match_dunn_p <- function(pair, lookup) {
  key1 <- paste(pair[1], "-", pair[2])
  key2 <- paste(pair[2], "-", pair[1])
  if (key1 %in% names(lookup)) return(lookup[[key1]])
  if (key2 %in% names(lookup)) return(lookup[[key2]])
  NA_real_
}

dunn_p <- vapply(comparison_list, match_dunn_p, numeric(1), lookup = dunn_lookup)
p_to_stars <- function(p) {
  as.character(cut(p, breaks = c(-Inf, .001, .01, .05, Inf), labels = c("***", "**", "*", "ns")))
}
dunn_labels <- p_to_stars(dunn_p)

# ------------------------------------------------------------------------------
# 5. DEFINE GRAPHICAL THEMES & PLOTTING FUNCTIONS
# ------------------------------------------------------------------------------
strip_labeller <- as_labeller(
  c(
    "Madrepora oculata"     = "bolditalic('Madrepora oculata')",
    "Desmophyllum pertusum" = "bolditalic('Desmophyllum pertusum')",
    "Primnoa msp."          = "bolditalic('Primnoa msp.')"
  ),
  label_parsed
)

theme_horizontal_panel <- function() {
  theme_minimal(base_size = 10, base_family = FONT_FAMILY) +
    theme(
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      axis.text.x      = element_text(color = "black", size = 8.5),
      axis.title.x     = element_text(face = "bold", size = 9.5, margin = margin(t = 8, b = 2)),
      axis.text.y      = element_text(color = "black", size = 8.5),
      axis.title.y     = element_text(face = "bold", size = 9.5, margin = margin(r = 6, b = 2)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
      panel.spacing    = unit(8, "pt"),
      strip.background = element_rect(fill = "white", color = "grey80", linewidth = 0.4),
      strip.text       = element_text(face = "plain", size = 9, margin = margin(t = 4, b = 4)), 
      legend.position  = "none",
      plot.tag         = element_text(face = "bold", size = 12),
      plot.margin      = margin(t = 6, r = 6, b = 8, l = 6)
    )
}

theme_set(theme_classic(base_size = 8, base_family = FONT_FAMILY))

publication_theme <- theme(
  plot.title         = element_blank(), 
  axis.title.x       = element_blank(), 
  axis.title.y       = element_text(margin = margin(r = 6, b = 2), face = "bold", color = "black", size = 8.5),
  axis.text.x        = element_text(color = "black", size = 7.5, hjust = 0.5, vjust = 0.5), 
  axis.text.y        = element_text(color = "black", size = 7.5),
  axis.line          = element_line(linewidth = 0.3, color = "black"),
  axis.ticks         = element_line(linewidth = 0.3, color = "black"),
  panel.grid.major.y = element_line(color = "gray96", linetype = "solid"),
  legend.position    = "none",
  plot.tag           = element_text(face = "bold", size = 11, vjust = 1),
  plot.margin        = margin(t = 8, r = 6, b = 8, l = 6)
)

area_axis_title <- expression(bold(paste("Colony planar area (cm"^2, ", natural-log scale)")))

make_density_plot <- function(grp_name, tag) {
  sub_data   <- master_long %>% filter(Group == grp_name)
  sub_stats  <- stats_summary_row1 %>% filter(Group == grp_name)
  col_22     <- species_palette[[grp_name]]
  
  # Plain text, stacked with \n — no plotmath, no parse=TRUE needed.
  # This is what makes it look "neat": ggplot's native text layout handles
  # line spacing itself instead of us fighting atop() nesting.
  n_text  <- paste0("n: ", sub_stats$n_15, " \u2192 ", sub_stats$n_22)  # \u2192 = "→"
  ks_text <- paste0("KS: D = ", round(sub_stats$ks_stat, 2))
  p_text  <- sub_stats$p_label
  stats_block <- paste(n_text, ks_text, p_text, sep = "\n")

  ggplot(sub_data, aes(x = Log_Area)) +
    geom_density(data = filter(sub_data, Year == "2015"),
                 fill = c_2015_fill, color = c_2015_line, linetype = "dashed", alpha = 0.85, linewidth = 0.5) +
    geom_density(data = filter(sub_data, Year == "2022"),
                 fill = col_22, color = col_22, linetype = "solid", alpha = 0.4, linewidth = 0.5) +
    geom_vline(xintercept = log(sub_stats$med_15), color = c_2015_line, linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = log(sub_stats$med_22), color = col_22,      linetype = "solid",  linewidth = 0.5) +
    # EDIT vjust: 1 = flush against the top of the panel; raise slightly (e.g. 1.1–1.3) for a small margin
    # EDIT lineheight: gap between the 3 lines — 0.9 is close to your reference image
    # EDIT size: font size of the whole block
    annotate("text", x = -Inf, y = Inf, label = stats_block,
             size = 2.3, hjust = -0.1, vjust = 1.15, lineheight = 0.9,
             color = "black", family = FONT_FAMILY) +
    facet_wrap(~Group, labeller = strip_labeller) +
    scale_x_continuous(breaks = log(c(1, 10, 100, 1000)), labels = scales::comma(c(1, 10, 100, 1000)), limits = log(c(0.8, 2500))) +
    coord_cartesian(clip = "off") +
    labs(x = area_axis_title, y = "Density probability", tag = tag) +
    theme_horizontal_panel()
}

# ------------------------------------------------------------------------------
# 6. GENERATE ROW 1 PANELS (PANELS A, B, C)
# ------------------------------------------------------------------------------
cat("Building Row 1 density distribution plots...\n")
# x-axis title is identical across A-C, so it's shown only once, under the middle panel
p1 <- make_density_plot("Madrepora oculata", "A") + 
  theme(axis.title.x = element_blank())
p2 <- make_density_plot("Desmophyllum pertusum", "B") + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())
p3 <- make_density_plot("Primnoa msp.", "C") + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())

# ------------------------------------------------------------------------------
# 7. GENERATE ROW 2 PANELS (PANELS D, E, F)
# ------------------------------------------------------------------------------
cat("Building Row 2 demographic plots...\n")

# Panel D: Growth rate distributions (with headroom & clipping fixes)
p1_growth <- ggplot(df_growth, aes(x = Species, y = Annual_Area_Change)) +
  geom_violin(aes(fill = Species), alpha = 0.08, color = NA) +
  geom_jitter(aes(color = Species, alpha = detectable_change), width = 0.18, size = 0.6) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, alpha = 0.95, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.3) +
  geom_segment(data = mdc_values, aes(x = as.numeric(Species) - 0.22, xend = as.numeric(Species) + 0.22,
                                      y = mean_mdc, yend = mean_mdc), linetype = "dashed", color = "black", linewidth = 0.35) +
  geom_segment(data = mdc_values, aes(x = as.numeric(Species) - 0.22, xend = as.numeric(Species) + 0.22,
                                      y = -mean_mdc, yend = -mean_mdc), linetype = "dashed", color = "black", linewidth = 0.35) +
  geom_signif(
    comparisons = comparison_list, 
    annotations = dunn_labels,
    # EDIT: absolute y-position of each bracket, in cm² yr⁻¹ — order matches comparison_list:
    #   1) Madrepora–Desmophyllum, 2) Desmophyllum–Primnoa, 3) Madrepora–Primnoa
    # Your data currently top out ~63 (Primnoa outlier), so these start just above that
    # and step up by 8 — raise/lower all three together to move the whole stack up/down,
    # or change the gaps between them individually
    y_position  = c(58, 66, 74),
    textsize    = 3.0,
    vjust       = 0.3,
    tip_length  = 0.02, 
    linewidth   = 0.35, 
    color       = "black"
  ) +
  scale_fill_manual(values = species_palette) +
  scale_color_manual(values = species_palette) +
  scale_alpha_manual(values = c("FALSE" = 0.12, "TRUE" = 0.80)) +
  scale_x_discrete(labels = p1_labels) +
  #raise up if clipping observed in the signif levels
  scale_y_continuous(
    breaks = seq(0, 75, by = 25),          # EDIT: keeps gridlines at 0/25/50/75, no stray 100 line
    limits = c(NA, 80),                     # EDIT: hard ceiling — top bracket sits at 84, so 90 leaves ~6 units of clean air above it without ballooning the panel
    expand = expansion(mult = c(0.08, 0))   # top padding now handled by 'limits' above, not expansion
  ) +
  coord_cartesian(clip = "off") +
  labs(
    tag = "D",
    y = expression(bold(paste("Annual growth (", cm^2, " ", yr^-1, ")")))
  ) +
  publication_theme

# Panel E: Recruitment rates
p2_recruitment <- ggplot(summary_stats_row2, aes(x = Species, y = recruitment_rate, fill = Species)) +
  geom_col(width = 0.45, alpha = 0.85, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = rec_ci_lower, ymax = rec_ci_upper), width = 0.10, color = "black", linewidth = 0.4) +
  geom_text(aes(y = rec_ci_upper, label = paste0(round(recruitment_rate, 1), "%")), vjust = -0.6, size = 2.4, fontface = "bold") +
  scale_fill_manual(values = species_palette) +
  scale_x_discrete(labels = p2_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_cartesian(clip = "off") +
  labs(
    tag = "E",
    y = "Recruitment rate (%)"
  ) +
  publication_theme

# Panel F: Mortality rates
p3_mortality <- ggplot(summary_stats_row2, aes(x = Species, y = mortality_rate, fill = Species)) +
  geom_col(width = 0.45, alpha = 0.85, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mort_ci_lower, ymax = mort_ci_upper), width = 0.10, color = "black", linewidth = 0.4) +
  geom_text(aes(y = mort_ci_upper, label = paste0(round(mortality_rate, 1), "%")), vjust = -0.6, size = 2.4, fontface = "bold") +
  scale_fill_manual(values = species_palette) +
  scale_x_discrete(labels = p3_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_cartesian(clip = "off") +
  labs(
    tag = "F",
    y = "Mortality rate (%)"
  ) +
  publication_theme

# ------------------------------------------------------------------------------
# 8. BUILD CUSTOM LEGEND COMPONENTS
# ------------------------------------------------------------------------------
year_key <- ggplot() +
  geom_rect(aes(xmin = 1, xmax = 1.5, ymin = 1, ymax = 1.4), fill = c_2015_fill, color = c_2015_line, linetype = "dashed") +
  annotate("text", x = 1.7, y = 1.2, label = "2015", size = 3.2, fontface = "bold", hjust = 0, family = FONT_FAMILY) +
  geom_rect(aes(xmin = 3.0, xmax = 3.5, ymin = 1, ymax = 1.4), fill = "white", color = "grey30", linetype = "solid") +
  annotate("text", x = 3.7, y = 1.2, label = "2022", size = 3.2, fontface = "bold", hjust = 0, family = FONT_FAMILY) +
  coord_cartesian(xlim = c(0, 6), ylim = c(0.8, 1.6), clip = "off") +
  labs(title = "Survey year") +
  theme_void(base_family = FONT_FAMILY) +
  theme(
    plot.title  = element_text(size = 9.5, face = "bold", hjust = 0.5, margin = margin(b = 2)),
    plot.margin = margin(t = 6, r = 0, b = 6, l = 0)
  )

legend_key <- ggplot() +
  geom_point(aes(x = 1.0, y = 1), color = "#722082", size = 4.0) +
  geom_point(aes(x = 2.7, y = 1), color = "#B63679", size = 4.0) +
  geom_point(aes(x = 5.0, y = 1), color = "#fb9f3a", size = 4.0) +
  annotate("text", x = 1.15, y = 1, label = "Madrepora oculata", 
           fontface = "italic", size = 3.2, hjust = 0, family = FONT_FAMILY) +
  annotate("text", x = 2.85, y = 1, label = "Desmophyllum pertusum", 
           fontface = "italic", size = 3.2, hjust = 0, family = FONT_FAMILY) +
  annotate("text", x = 5.15, y = 1, label = "Primnoa msp.", 
           fontface = "italic", size = 3.2, hjust = 0, family = FONT_FAMILY) +
  coord_cartesian(xlim = c(0.8, 6.2), ylim = c(0.9, 1.1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 6, r = 0, b = 4, l = 0))

# ------------------------------------------------------------------------------
# 9. COMPOSITE LAYOUT ASSEMBLY & EXPORT
# ------------------------------------------------------------------------------
cat("Assembling master composite multi-panel layout...\n")
design_layout <- "
  AAA
  BCD
  EFG
  HHH
"

final_masterpiece <- wrap_plots(
  year_key,                                # A
  p1, p2, p3,                              # B, C, D
  p1_growth, p2_recruitment, p3_mortality, # E, F, G
  legend_key,                              # H
  design = design_layout
) + 
  plot_layout(heights = c(0.12, 1, 1, 0.10))

# Ensure output directory exists before export
dir.create(dirname(output_fig_path), recursive = TRUE, showWarnings = FALSE)

cat("Exporting publication graphic (600 DPI AGG PNG)...\n")
ggsave(
  filename = output_fig_path,
  plot     = final_masterpiece,
  width    = 174, 
  height   = 135, 
  units    = "mm",
  dpi      = 600,
  device   = ragg::agg_png
)

cat("\n======================================================================\n")
cat("Script 03 execution complete!\n")
cat("Master figure exported to:\n ", output_fig_path, "\n")
cat("======================================================================\n")

# ==============================================================================
# 10. REVIEWER AUDIT STATISTICAL SUMMARY TABLES
# ==============================================================================
cat("\n======================================================================\n")
cat("            COMPREHENSIVE REVIEWER AUDIT STATISTICAL SUMMARY            \n")
cat("======================================================================\n\n")

# ------------------------------------------------------------------------------
# A. TABLE 1: SIZE STRUCTURE PERCENTILES & KS TEST STATISTICS (PANELS A-C)
# ------------------------------------------------------------------------------
table1_size_audit <- master_long %>%
  group_by(Group, Year) %>%
  summarise(
    n       = n(),
    p10_cm2 = quantile(Area, 0.10, na.rm = TRUE),
    p50_cm2 = median(Area, na.rm = TRUE),
    p90_cm2 = quantile(Area, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = Year,
    values_from = c(n, p10_cm2, p50_cm2, p90_cm2)
  ) %>%
  left_join(stats_summary_row1 %>% select(Group, ks_stat, ks_p), by = "Group") %>%
  mutate(
    ks_D     = round(ks_stat, 3),
    ks_p_val = if_else(ks_p < 0.001, "< 0.001", as.character(round(ks_p, 4)))
  ) %>%
  select(
    Species = Group, 
    n_2015, n_2022, 
    p10_2015 = p10_cm2_2015, p10_2022 = p10_cm2_2022,
    p50_2015 = p50_cm2_2015, p50_2022 = p50_cm2_2022,
    p90_2015 = p90_cm2_2015, p90_2022 = p90_cm2_2022,
    ks_D, ks_p_val
  )

cat("--- TABLE 1: Size Class Structure Percentiles & KS Tests (Panels A-C) ---\n")
print(as.data.frame(table1_size_audit), row.names = FALSE)
cat("\n")

# ------------------------------------------------------------------------------
# B. TABLE 2: ANNUAL GROWTH DISTRIBUTION & POST-HOC COMPARISONS (PANEL D)
# ------------------------------------------------------------------------------
table2_growth_audit <- df_growth %>%
  group_by(Species) %>%
  summarise(
    n_persistent   = n(),
    mean_growth    = mean(Annual_Area_Change, na.rm = TRUE),
    median_growth  = median(Annual_Area_Change, na.rm = TRUE),
    p25_growth     = quantile(Annual_Area_Change, 0.25, na.rm = TRUE),
    p75_growth     = quantile(Annual_Area_Change, 0.75, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  left_join(mdc_values, by = "Species")

cat("--- TABLE 2A: Growth Rate Distributions (cm² yr⁻¹) & MDC95 (Panel D) ---\n")
print(as.data.frame(table2_growth_audit), row.names = FALSE)

cat(sprintf("\nKruskal-Wallis Global Test: H = %.2f, df = %d, p = %.4e\n\n", 
            unname(kw_test$statistic), unname(kw_test$parameter), unname(kw_test$p.value)))

cat("--- TABLE 2B: Dunn's Post-Hoc Pairwise Comparisons (BH Adjusted) ---\n")
table2_dunn_audit <- data.frame(
  Comparison   = names(dunn_lookup),
  P_Adjusted   = round(unname(dunn_lookup), 5),
  Significance = p_to_stars(unname(dunn_lookup))
)
print(table2_dunn_audit, row.names = FALSE)
cat("\n")

# ------------------------------------------------------------------------------
# C. TABLE 3: DEMOGRAPHIC TURNOVER & WILSON 95% CIs (PANELS E & F)
# ------------------------------------------------------------------------------
table3_demo_audit <- summary_stats_row2 %>%
  mutate(
    net_change_pct = recruitment_rate - mortality_rate
  ) %>%
  select(
    Species, 
    n_2015_total = n_2015, n_mortality, mortality_rate, mort_ci_lower, mort_ci_upper,
    n_2022_total = n_2022, n_recruitment, recruitment_rate, rec_ci_lower, rec_ci_upper,
    net_change_pct
  )

cat("--- TABLE 3: Demographic Rates & 95% Wilson CIs (Panels E-F) ---\n")
print(as.data.frame(table3_demo_audit), row.names = FALSE)
cat("\n======================================================================\n")

# ------------------------------------------------------------------------------
# D. OPTIONAL: EXPORT AUDIT TABLES TO A SINGLE CSV FILE
# ------------------------------------------------------------------------------
audit_csv_path <- file.path(base_dir, "Figures/Statistical_Summary_Table.csv")

suppressWarnings({
  write.table("TABLE 1: SIZE CLASS STRUCTURE & KS TESTS", audit_csv_path, row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(table1_size_audit, audit_csv_path, append = TRUE, row.names = FALSE, sep = ",")
  
  write.table("\nTABLE 2A: GROWTH RATE DISTRIBUTIONS", audit_csv_path, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(table2_growth_audit, audit_csv_path, append = TRUE, row.names = FALSE, sep = ",")
  
  write.table("\nTABLE 2B: DUNN POST-HOC PAIRWISE COMPARISONS", audit_csv_path, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(table2_dunn_audit, audit_csv_path, append = TRUE, row.names = FALSE, sep = ",")
  
  write.table("\nTABLE 3: DEMOGRAPHIC RATES & WILSON 95% CIs", audit_csv_path, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(table3_demo_audit, audit_csv_path, append = TRUE, row.names = FALSE, sep = ",")
})

cat("Statistical audit CSV exported to:\n ", audit_csv_path, "\n")
cat("======================================================================\n")

```
