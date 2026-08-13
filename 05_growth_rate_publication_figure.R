# ==============================================================================
# Script Name:  05_growth_rate_publication_figure.R
# Description:  Loads the fitted heteroscedastic asymptotic (GNLS) growth model
#                from Script 04, extracts biological thresholds dynamically, 
#                and generates an Ecography-styled publication figure mapping 
#                relative growth rate (RGR) across initial sizes.
# Dependencies: readxl, dplyr, nlme, ggplot2, ragg
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(nlme)
library(ggplot2)
library(ragg)

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & FILE PATHS
# ------------------------------------------------------------------------------
base_dir       <- "D:/PhD_Data(Large)/Submission_Dataset"

input_xlsx_path <- file.path(base_dir, "Master_CWC_Tracked_MDC.xlsx")
model_rds_path  <- file.path(base_dir, "models", "m_growth_gnls_SSasymp.rds")
output_fig_path <- file.path(base_dir, "Figures", "Fig_growth_rate_publication.png")

FONT_FAMILY  <- "sans"
group_levels <- c("Madrepora oculata", "Desmophyllum pertusum", "Primnoa msp.")

species_palette <- c(
  "Madrepora oculata"     = "#722082",
  "Desmophyllum pertusum" = "#B63679",
  "Primnoa msp."          = "#fb9f3a"
)

species_labels_expr <- c(
  "Madrepora oculata"     = expression(italic("M. oculata")),
  "Desmophyllum pertusum" = expression(italic("D. pertusum")),
  "Primnoa msp."          = expression(italic("Primnoa") ~ "msp.")
)

# ------------------------------------------------------------------------------
# 3. LOAD DATASET & SAVED GNLS MODEL OBJECT
# ------------------------------------------------------------------------------
cat("Loading cleaned dataset and saved GNLS growth model...\n")
df_raw <- read_excel(input_xlsx_path)

df_clean <- df_raw %>%
  filter(Species != "Coral Recruit") %>%
  mutate(
    Group = case_when(
      Species == "Madrepora oculata" ~ "Madrepora oculata",
      Species %in% c("D. pertusum", "Desmophyllum pertusum") ~ "Desmophyllum pertusum",
      Species %in% c("Primnoa msp.5", "Primnoa msp.1") ~ "Primnoa msp.",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%
  mutate(fSpecies = factor(Group, levels = group_levels))

# Persistent colonies dataset (matching Script 04)
df_growth <- df_clean %>%
  filter(state_transition == "persistence", !is.na(RGR_Planar), !is.na(A1))

# Load model object output from Script 04
if (!file.exists(model_rds_path)) {
  stop("ERROR: Model object not found at ", model_rds_path, "\nPlease run Script 04 first to save the model RDS.")
}
m_growth_gnls <- readRDS(model_rds_path)

cat("Successfully loaded GNLS model object.\n")

# ------------------------------------------------------------------------------
# 4. DYNAMICALLY EXTRACT MODEL PARAMETERS & THRESHOLDS
# ------------------------------------------------------------------------------
cat("Extracting asymptotic parameters and computing 95% size thresholds...\n")

cf <- coef(m_growth_gnls)

# Reference species: Madrepora oculata
asym_mo <- cf["Asym.(Intercept)"]
lrc_mo  <- cf["lrc.(Intercept)"]
rc_mo   <- exp(lrc_mo)

# Desmophyllum pertusum
asym_dp <- cf["Asym.(Intercept)"] + cf["Asym.fSpeciesDesmophyllum pertusum"]
lrc_dp  <- cf["lrc.(Intercept)"]   + cf["lrc.fSpeciesDesmophyllum pertusum"]
rc_dp   <- exp(lrc_dp)

# Primnoa msp.
asym_pr <- cf["Asym.(Intercept)"] + cf["Asym.fSpeciesPrimnoa msp."]
lrc_pr  <- cf["lrc.(Intercept)"]   + cf["lrc.fSpeciesPrimnoa msp."]
rc_pr   <- exp(lrc_pr)

# 95% Asymptote Threshold Calculation
size_95_mo <- -log(0.05) / rc_mo
size_95_dp <- -log(0.05) / rc_dp
size_95_pr <- -log(0.05) / rc_pr

# Summary DataFrames with dynamic text formatting
asym_df <- data.frame(
  fSpecies   = factor(group_levels, levels = group_levels),
  asym_val   = c(asym_mo, asym_dp, asym_pr),
  label_text = sprintf("Asym: %.3f yr⁻¹", c(asym_mo, asym_dp, asym_pr))
)

size_thresholds <- data.frame(
  fSpecies   = factor(group_levels, levels = group_levels),
  asym_size  = c(size_95_mo, size_95_dp, size_95_pr),
  label_text = paste0(round(c(size_95_mo, size_95_dp, size_95_pr), 1), " cm²")
)

# ------------------------------------------------------------------------------
# 6. BUILD ECOGRAPHY PUBLICATION FIGURE
# ------------------------------------------------------------------------------
cat("Assembling publication graphic...\n")

x_breaks <- seq(0, 300, by = 50)
y_breaks <- seq(-0.2, 0.5, by = 0.1)

plot_growth_pub <- ggplot(df_growth, aes(x = A1, y = RGR_Planar, color = fSpecies)) +
  
  # Raw observations
  geom_point(alpha = 0.30, size = 1.2, stroke = 0) +
  
  # Fitted non-linear SSasymp trajectories
  geom_line(data = new_data, aes(y = pred), linewidth = 0.9, show.legend = TRUE) +
  
  # Horizontal baseline asymptotes (Asym plateau)
  geom_hline(
    data = asym_df,
    aes(yintercept = asym_val, color = fSpecies),
    linetype = "dashed", linewidth = 0.4, alpha = 0.8,
    show.legend = FALSE
  ) +
  
  # Vertical 95% threshold size markers
  geom_vline(
    data = size_thresholds,
    aes(xintercept = asym_size, color = fSpecies),
    linetype = "dotted", linewidth = 0.5, alpha = 0.85,
    show.legend = FALSE
  ) +
  
  # Vertical threshold labels (anchored near top)
  geom_text(
    data = size_thresholds,
    aes(x = asym_size, y = 0.48, label = label_text, color = fSpecies),
    angle = 90,
    vjust = -0.4,
    hjust = 1,
    size = 2.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  # Horizontal asymptote labels (anchored near right edge)
  geom_text(
    data = asym_df,
    aes(x = 295, y = asym_val, label = label_text, color = fSpecies),
    hjust = 1,
    vjust = -0.4,
    size = 2.4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  # Scale and Color Mapping
  scale_color_manual(
    values = species_palette,
    labels = species_labels_expr
  ) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks, labels = sprintf("%.1f", y_breaks)) +
  coord_cartesian(xlim = c(0, 300), ylim = c(-0.22, 0.52)) +
  
  # Axis Labels
  labs(
    x = expression(paste("Initial Planar Area ", A[1], " (cm"^2*")")),
    y = expression(paste("Annualized Relative Growth Rate (RGR yr"^-1*")")),
    color = NULL
  ) +
  
  # Ecography Clean Theme
  theme_classic(base_size = 8, base_family = FONT_FAMILY) +
  theme(
    legend.position   = "top",
    legend.key        = element_blank(),
    legend.key.width  = unit(8, "mm"),
    legend.text       = element_text(size = 8, face = "italic"),
    axis.title        = element_text(face = "bold", size = 8.5),
    axis.text         = element_text(color = "black", size = 7.5),
    axis.line         = element_line(linewidth = 0.4),
    axis.ticks        = element_line(linewidth = 0.4)
  )

# ------------------------------------------------------------------------------
# 7. EXPORT HIGH-RESOLUTION GRAPHIC
# ------------------------------------------------------------------------------
dir.create(dirname(output_fig_path), recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = output_fig_path, 
  plot     = plot_growth_pub,
  width    = 110, 
  height   = 90, 
  units    = "mm", 
  dpi      = 600, 
  device   = ragg::agg_png,
  bg       = "white"
)


cat("\n======================================================================\n")
cat("Script 05 Execution Complete!\n")
cat("Publication figure saved to:\n  ", output_fig_path, "\n")
cat("======================================================================\n")
