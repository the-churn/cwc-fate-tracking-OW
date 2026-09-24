# ==============================================================================
# Script Name:  12_recruit_id_size_threshold_SI.R
# Description:  Tests whether the 2015-to-2022 shift in colony size-frequency
#               distributions (ln(Area)) is an artefact of variable taxonomic
#               identification of small recruits, rather than a real
#               demographic signal (addresses manuscript reviewer comment
#               AR2 on the Demographic Characterization section).
#
#               1. Recovers per-year (not coalesced) TagLab.Class.name for
#                  every tracked genet, since master_CWC_tracked.csv only
#                  keeps one coalesced Species column.
#               2. Finds the 90th percentile of planar area among colonies
#                  classed as the generic "Coral Recruit" morphotype, in
#                  EACH year separately, then takes the max of the two as a
#                  single pooling threshold (this script previously used the
#                  strict max, but that made the threshold fully determined
#                  by a single largest colony -- one mislabelled/occluded
#                  "Coral Recruit" annotation could drag the whole pooling
#                  boundary up, as happened with a 57.67 cm^2 outlier found
#                  during review. p90 is robust to that: it takes a genuine
#                  shift in the whole distribution to move the threshold,
#                  not one bad annotation. Manuscript text currently cites
#                  19.01 cm^2 from the old max-based definition -- that
#                  number will need updating to whatever p90 resolves to).
#               3. Pools every colony <= threshold (unidentified recruits +
#                  small D. pertusum / M. oculata) into one scleractinian
#                  cohort, and keeps every colony > threshold in its own
#                  species-specific cohort -- EXCEPT that a colony still
#                  labelled "Coral Recruit" is always pooled regardless of
#                  its size, since by definition it has no reliable species
#                  identity to analyze separately (this matters now that the
#                  threshold is p90, not max: roughly the largest ~10% of
#                  recruit-labelled colonies WILL exceed it, and must stay
#                  pooled rather than being spun off into their own
#                  "Coral Recruit" cohort).
#               4. Runs a two-sample Kolmogorov-Smirnov test (2015 vs 2022,
#                  on ln(Area)) separately within each cohort, exactly as
#                  described in the Demographic Characterization methods.
#               5. Plots the per-cohort size-frequency distributions by year
#                  and exports an SI figure + a results table.
#
#               Inputs are the two files Script 01 already exports --
#               master_data_combined.csv (per-polygon, pre-aggregation, still
#               has per-year class names) and master_CWC_tracked.csv (per-
#               genet, post-cleaning/occlusion-handling). No re-run of
#               Script 01 is required, and Script 01 itself is not modified.
#
# Dependencies: tidyverse (dplyr, readr, stringr, tidyr, purrr, ggplot2)
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# 1. CONFIGURATION & FILE PATHS  (same base_dir convention as Scripts 01-11)
# ------------------------------------------------------------------------------
base_dir <- "D:/PhD_Data(Large)/Submission_Dataset"

master_data_path <- file.path(base_dir, "master_data_combined.csv")  # Script 01 output
tracked_path      <- file.path(base_dir, "master_CWC_tracked.csv")    # Script 01 output

figures_dir <- file.path(base_dir, "figures")
results_dir <- file.path(base_dir, "results")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

output_fig     <- file.path(figures_dir, "FigureS_Recruit_ID_Bias_Check.png")
output_results <- file.path(results_dir, "recruit_id_bias_KS_results.csv")

# Set TRUE to keep Primnoa msp. as its own (always species-specific, never
# pooled) cohort alongside the two scleractinians. Left off by default because
# the manuscript text frames this specifically as a SCLERACTINIAN pooling
# check, and small Primnoa colonies are never labelled "Coral Recruit" in
# Script 01 (they're regex-matched to their own Primnoa class at any size --
# see the `str_detect(TagLab.Class.name, "^Primnoa")` branch), so they were
# never part of the identification-bias concern this test addresses.
include_primnoa <- FALSE

# Must stay in sync with Script 01's exclude_genet_ids -- these genets were
# already dropped before master_data_combined.csv was written in this
# session's run, so re-listing them here is only a safety net if you're
# pointing this script at an older combined file. Copy-pasted verbatim from
# Script 01; update both places together if that list changes.
exclude_genet_ids <- c(
  941, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1008, 1005, 1004,
  1006, 1007, 984, 971, 975, 973, 974, 966, 965, 1086, 989, 993,
  994, 995, 1115
)

# ------------------------------------------------------------------------------
# 2. RECOVER PER-YEAR SPECIES/CLASS LABELS
# ------------------------------------------------------------------------------
# master_CWC_tracked.csv only has one coalesced `Species` column (2015 label
# wins if a genet was present both years), which is exactly wrong for this
# test: a colony first tagged "Coral Recruit" in 2015 and identified to
# species by 2022 must be evaluated against ITS 2022 label when building the
# 2022 population, not its 2015 one. So we rebuild the per-year aggregation
# from master_data_combined.csv, using the identical Join_ID logic as Script
# 01 Section 3 (same filters, same row order -> same IDs), and keep
# Species_2015 / Species_2022 separate instead of coalescing them away.
if (!file.exists(master_data_path)) {
  stop("Missing ", master_data_path, " -- run Script 01 first.", call. = FALSE)
}

master_dataset <- read_csv(master_data_path, show_col_types = FALSE)

regions_cleaned <- master_dataset %>%
  filter(TagLab.Type == "Region") %>%
  filter(!is.na(Image.name) & Image.name != "") %>%
  filter(!TagLab.Genet.Id %in% exclude_genet_ids) %>%
  mutate(
    Join_ID = ifelse(
      TagLab.Genet.Id == 0 | is.na(TagLab.Genet.Id),
      paste0(Year, "_Single_", row_number()),
      as.character(TagLab.Genet.Id)
    )
  )

species_by_year <- regions_cleaned %>%
  group_by(Join_ID, Year) %>%
  summarise(Species = first(TagLab.Class.name), .groups = "drop") %>%
  pivot_wider(names_from = Year, values_from = Species, names_prefix = "Species_")

# ------------------------------------------------------------------------------
# 3. LOAD TRACKED GENETS & ATTACH PER-YEAR SPECIES
# ------------------------------------------------------------------------------
if (!file.exists(tracked_path)) {
  stop("Missing ", tracked_path, " -- run Script 01 first.", call. = FALSE)
}

tracked <- read_csv(tracked_path, show_col_types = FALSE) %>%
  # TagLab.Genet.Id is semantically a string (it can be either a plain numeric
  # genet id or a "<year>_Single_<n>" composite id -- see Script 01 Section
  # 3), but if every genet in THIS file happens to have a plain numeric id
  # (e.g. no untracked singletons survived cleaning), readr's type-guessing
  # will silently infer the column as double on read, which then fails to
  # join against species_by_year$Join_ID (always character, since it's built
  # from paste0()/as.character() below). Force it explicitly so the join
  # works regardless of what readr guessed from this particular file's
  # contents.
  mutate(TagLab.Genet.Id = as.character(TagLab.Genet.Id)) %>%
  left_join(species_by_year, by = c("TagLab.Genet.Id" = "Join_ID"))

n_unmatched <- sum(is.na(tracked$Species_2015) & is.na(tracked$Species_2022))
if (n_unmatched > 0) {
  warning(
    n_unmatched, " genet(s) in master_CWC_tracked.csv had no matching Join_ID in ",
    "master_data_combined.csv -- check that both files came from the same Script 01 run."
  )
}

# Harmonize raw class names the same way Scripts 03/05/06 do (e.g. "D.
# pertusum" -> "Desmophyllum pertusum"), applied separately per year so a
# colony's cohort in 2015 doesn't leak its 2022 identity or vice versa.
harmonize_species <- function(x) {
  case_when(
    x == "Madrepora oculata" ~ "Madrepora oculata",
    x %in% c("D. pertusum", "Desmophyllum pertusum") ~ "Desmophyllum pertusum",
    str_detect(x, "^Primnoa") ~ "Primnoa msp.",
    x == "Coral Recruit" ~ "Coral Recruit",
    TRUE ~ x
  )
}

tracked <- tracked %>%
  mutate(
    Species_2015_clean = harmonize_species(Species_2015),
    Species_2022_clean = harmonize_species(Species_2022)
  )

# ------------------------------------------------------------------------------
# 4. DERIVE THE RECRUIT-ID SIZE THRESHOLD (90th percentile of "Coral Recruit"
#    area in each year, then the max of the two -- i.e. across both
#    timepoints). p90 instead of max() so a single outlier/mislabelled
#    recruit can't single-handedly set the pooling boundary.
# ------------------------------------------------------------------------------
recruit_areas_2015 <- tracked$Area_2015[tracked$Species_2015_clean == "Coral Recruit"]
recruit_areas_2022 <- tracked$Area_2022[tracked$Species_2022_clean == "Coral Recruit"]

n_recruit_2015 <- sum(!is.na(recruit_areas_2015))
n_recruit_2022 <- sum(!is.na(recruit_areas_2022))

max_recruit_2015 <- max(recruit_areas_2015, na.rm = TRUE)   # kept for reference/comparison only
max_recruit_2022 <- max(recruit_areas_2022, na.rm = TRUE)   # kept for reference/comparison only
p90_recruit_2015 <- quantile(recruit_areas_2015, 0.90, na.rm = TRUE, names = FALSE)
p90_recruit_2022 <- quantile(recruit_areas_2022, 0.90, na.rm = TRUE, names = FALSE)

recruit_threshold <- max(p90_recruit_2015, p90_recruit_2022)

cat("======================================================================\n")
cat("RECRUIT IDENTIFICATION-BIAS SIZE THRESHOLD (p90-based)\n")
cat("======================================================================\n")
cat("2015: n =", n_recruit_2015, " | p90 =", round(p90_recruit_2015, 2),
    "cm^2 | max =", round(max_recruit_2015, 2), "cm^2\n")
cat("2022: n =", n_recruit_2022, " | p90 =", round(p90_recruit_2022, 2),
    "cm^2 | max =", round(max_recruit_2022, 2), "cm^2\n")
cat("Pooling threshold (max of the two p90s):", round(recruit_threshold, 2), "cm^2\n")
if (min(n_recruit_2015, n_recruit_2022) < 20) {
  cat("NOTE: fewer than 20 'Coral Recruit' colonies in at least one year --\n")
  cat(" a 90th-percentile estimate is noisy at this sample size; treat it as\n")
  cat(" approximate and sanity-check against the candidate list below.\n")
}
cat("(Manuscript text currently cites 19.01 cm^2 from the old max-based\n")
cat(" definition -- update it to the p90 value above once you're satisfied\n")
cat(" it's not being distorted by an annotation error; see candidates below.)\n\n")

write_excel_csv(
  tibble(
    n_recruit_2015, n_recruit_2022,
    max_recruit_2015, max_recruit_2022,
    p90_recruit_2015, p90_recruit_2022,
    recruit_threshold
  ),
  file.path(results_dir, "recruit_threshold_summary.csv")
)

# ------------------------------------------------------------------------------
# 4b. TOP-N LARGEST "CORAL RECRUIT" CANDIDATES -- with a p90 threshold these
#     are no longer necessarily the colonies SETTING the boundary (p90 by
#     construction ignores the top ~10% tail), but they're still the ones
#     worth a manual look: genuinely large recruits and mislabelled/occluded
#     colonies that should have been ID'd to species both show up here, and
#     only an audit of each one (not the statistic) can tell them apart.
#     Cross-reference TagLab.Genet.Id against
#     data/audit/TagLab_QuickJump_OutlierAudit.xlsx and/or re-open the listed
#     image(s) in TagLab before trusting either 19.01 or the value above.
# ------------------------------------------------------------------------------
top_n_recruit_candidates <- 5

images_by_genet_year <- regions_cleaned %>%
  group_by(Join_ID, Year) %>%
  summarise(Images = paste(unique(Image.name), collapse = "; "), Poly_Count = n(), .groups = "drop")

recruit_candidates <- bind_rows(
  tracked %>%
    filter(Species_2015_clean == "Coral Recruit") %>%
    transmute(TagLab.Genet.Id, Year = 2015L, Area = Area_2015),
  tracked %>%
    filter(Species_2022_clean == "Coral Recruit") %>%
    transmute(TagLab.Genet.Id, Year = 2022L, Area = Area_2022)
) %>%
  arrange(Year, desc(Area)) %>%
  group_by(Year) %>%
  slice_head(n = top_n_recruit_candidates) %>%
  ungroup() %>%
  left_join(images_by_genet_year, by = c("TagLab.Genet.Id" = "Join_ID", "Year"))

cat("Largest 'Coral Recruit'-labelled colonies (top", top_n_recruit_candidates, "per year) --\n")
cat("check these before trusting the threshold above:\n")
print(recruit_candidates, width = Inf, n = Inf)
cat("\n")

write_excel_csv(recruit_candidates, file.path(results_dir, "recruit_threshold_top_candidates.csv"))

# ------------------------------------------------------------------------------
# 5. BUILD PER-YEAR POPULATIONS & ASSIGN COHORTS
# ------------------------------------------------------------------------------
pop_2015 <- tracked %>%
  filter(!is.na(Area_2015), Area_2015 > 0) %>%
  transmute(TagLab.Genet.Id, Year = 2015L, Area = Area_2015, Species = Species_2015_clean)

pop_2022 <- tracked %>%
  filter(!is.na(Area_2022), Area_2022 > 0) %>%
  transmute(TagLab.Genet.Id, Year = 2022L, Area = Area_2022, Species = Species_2022_clean)

demo_pop <- bind_rows(pop_2015, pop_2022) %>%
  mutate(
    scleractinian = Species %in% c("Desmophyllum pertusum", "Madrepora oculata", "Coral Recruit"),
    cohort = case_when(
      # Unconditional: an unidentified "Coral Recruit" is always pooled, at
      # any size -- it has no species identity to analyze separately. This
      # branch must come before the size check below: with a p90 threshold
      # (rather than the old max), roughly the largest ~10% of recruit-
      # labelled colonies in each year WILL have Area > recruit_threshold,
      # and without this branch they'd fall into the next line and get
      # spun off into their own literal "Coral Recruit" cohort/facet.
      Species == "Coral Recruit" ~ "Pooled scleractinian (\u2264 threshold)",
      scleractinian & Area <= recruit_threshold ~ "Pooled scleractinian (\u2264 threshold)",
      scleractinian & Area > recruit_threshold  ~ Species,
      include_primnoa & Species == "Primnoa msp." ~ "Primnoa msp.",
      TRUE ~ NA_character_
    ),
    ln_Area = log(Area)
  ) %>%
  filter(!is.na(cohort))

cat("Colonies per cohort x year going into the K-S tests:\n")
print(demo_pop %>% count(cohort, Year) %>% pivot_wider(names_from = Year, values_from = n, values_fill = 0))
cat("\n")

# ------------------------------------------------------------------------------
# 6. TWO-SAMPLE K-S TESTS, 2015 vs 2022, WITHIN EACH COHORT
# ------------------------------------------------------------------------------
run_ks <- function(df) {
  x2015 <- df$ln_Area[df$Year == 2015]
  x2022 <- df$ln_Area[df$Year == 2022]
  if (length(x2015) < 2 || length(x2022) < 2) {
    return(tibble(n_2015 = length(x2015), n_2022 = length(x2022),
                   D = NA_real_, p_value = NA_real_,
                   note = "Too few colonies for a K-S test (need >=2 per year)"))
  }
  has_ties <- length(unique(x2015)) < length(x2015) || length(unique(x2022)) < length(x2022) ||
    length(intersect(x2015, x2022)) > 0
  res <- suppressWarnings(ks.test(x2015, x2022))
  tibble(
    n_2015 = length(x2015), n_2022 = length(x2022),
    D = unname(res$statistic), p_value = res$p.value,
    note = if (has_ties) "Ties present -- asymptotic p-value (see ?ks.test)" else ""
  )
}

ks_results <- demo_pop %>%
  group_by(cohort) %>%
  group_modify(~ run_ks(.x)) %>%
  ungroup() %>%
  mutate(across(c(D, p_value), ~ round(.x, 4)))

# Consistent 3-decimal display strings for the figure annotations (round() +
# paste0 silently drops trailing zeros, e.g. "0.173" vs "0.2297" side by side)
fmt3 <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.3f", x))

cat("K-S test results (ln(Area), 2015 vs 2022) by cohort:\n")
print(ks_results, width = Inf)

# write_excel_csv() (not write_csv()) adds a UTF-8 BOM so the "<=" character
# in cohort names survives being opened in Excel on Windows -- write_csv()
# alone renders it as mojibake ("â‰¤") there even though the file's bytes are correct UTF-8.
write_excel_csv(ks_results, output_results)
cat("\nResults table saved to:\n ", output_results, "\n\n")

# ------------------------------------------------------------------------------
# 7. SI FIGURE -- SIZE-FREQUENCY DISTRIBUTIONS BY COHORT AND YEAR
# ------------------------------------------------------------------------------
year_colors <- c("2015" = "#3B0F4B", "2022" = "#F69336")  # reuses repo palette hues

# Order facets: pooled cohort first, then species-specific cohorts
cohort_levels <- c("Pooled scleractinian (\u2264 threshold)",
                    setdiff(unique(demo_pop$cohort), "Pooled scleractinian (\u2264 threshold)"))

plot_df <- demo_pop %>%
  mutate(
    cohort = factor(cohort, levels = cohort_levels),
    Year_f = factor(Year)
  )

ann_df <- ks_results %>%
  mutate(
    cohort = factor(cohort, levels = cohort_levels),
    label = ifelse(
      is.na(D), "n too small for K-S",
      paste0("D = ", fmt3(D), "\np ", ifelse(p_value < 0.001, "< 0.001", paste0("= ", fmt3(p_value))))
    )
  )

# Common bin width across panels so histograms are visually comparable
bw <- diff(range(plot_df$ln_Area, na.rm = TRUE)) / 20

si_figure <- ggplot(plot_df, aes(x = ln_Area, fill = Year_f, color = Year_f)) +
  geom_histogram(position = "identity", alpha = 0.55, binwidth = bw, linewidth = 0.2) +
  geom_text(
    data = ann_df, aes(label = label), x = -Inf, y = Inf,
    hjust = -0.08, vjust = 1.15, inherit.aes = FALSE, size = 2.6, lineheight = 0.9
  ) +
  facet_wrap(~cohort, scales = "free_y", ncol = length(cohort_levels)) +
  scale_fill_manual(values = year_colors, name = "Survey year") +
  scale_color_manual(values = year_colors, name = "Survey year") +
  labs(
    x = expression(paste("ln(Colony Planar Area, cm"^2, ")")),
    y = "Number of colonies"
  ) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    strip.background = element_blank(),
    strip.text        = element_text(size = 7, face = "italic"),
    axis.title        = element_text(size = 8, color = "black"),
    axis.text         = element_text(size = 7, color = "black"),
    axis.line         = element_line(linewidth = 0.4, color = "black"),
    axis.ticks        = element_line(linewidth = 0.4, color = "black"),
    legend.position   = "bottom",
    legend.title      = element_text(size = 7),
    legend.text       = element_text(size = 7),
    legend.key.size   = unit(0.35, "cm"),
    plot.margin       = margin(t = 5, r = 8, b = 5, l = 5, unit = "pt")
  )

ggsave(
  filename = output_fig, plot = si_figure,
  width = 178, height = 70, units = "mm", dpi = 600
)

cat("======================================================================\n")
cat("SI figure saved to:\n ", output_fig, "\n")
cat("Script 12 execution complete!\n")
cat("======================================================================\n")
