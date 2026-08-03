# ==============================================================================
# Script Name:  05_growth_rate_publication_figure.R
# Description:  Fits the final heteroscedastic asymptotic (GNLS) growth model,
#               extracts biological thresholds, and generates an Ecography-styled
#               publication figure mapping relative growth rate (RGR) across sizes.
# Dependencies: readxl, dplyr, nlme, ggplot2, ragg
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(nlme)
library(ggplot2)

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & FILE PATHS
# ------------------------------------------------------------------------------
base_dir      <- "D:/PhD_Data(Large)/06_Final_Analysis/Master_Dataset"
input_xlsx    <- file.path(base_dir, "Error_Propagation/Master_CWC_Tracked_MDC_Age.xlsx")
output_fig    <- file.path(base_dir, "Figures/Fig_growth_rate.png")

# ------------------------------------------------------------------------------
# 3. DATA CLEANING & GROWTH SUBSET PREPARATION
# ------------------------------------------------------------------------------
cat("Loading and filtering master dataset...\n")
zuur <- read_excel(input_xlsx)

zuur_clean <- zuur %>%
  filter(state_transition == "persistence" & Species != "Coral Recruit") %>%
  mutate(
    fSpecies = case_when(
      Species == "Madrepora oculata" ~ "Madrepora oculata",
      Species %in% c("D. pertusum", "Desmophyllum pertusum") ~ "D. pertusum",
      Species %in% c("Primnoa msp.5", "Primnoa msp.1") ~ "Primnoa msp.",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(fSpecies)) %>%
  mutate(fSpecies = factor(fSpecies, levels = c("D. pertusum", "Madrepora oculata", "Primnoa msp.")))

# Isolate positive growth subset (RGR >= 0)
zuur_growth <- zuur_clean %>%
  filter(RGR_Planar >= 0)

# ------------------------------------------------------------------------------
# 4. FINAL GNLS MODEL FITTING (VarPower Variance Structure)
# ------------------------------------------------------------------------------
cat("Fitting starting values and final GNLS asymptotic model...\n")

# Self-starting initial estimates via nlsList
m_list <- nlsList(
  RGR_Planar ~ SSasymp(A1, Asym, R0, lrc) | fSpecies, 
  data = zuur_growth,
  na.action = na.omit
)
start_vec <- c(coef(m_list)$Asym, coef(m_list)$R0, coef(m_list)$lrc)

# Final robust GNLS model using the parsimonious VarPower_Sp weighting structure
m_gnls_final <- gnls(
  RGR_Planar ~ SSasymp(A1, Asym, R0, lrc),
  params = list(Asym ~ fSpecies, R0 ~ fSpecies, lrc ~ fSpecies),
  data = zuur_growth,
  weights = varPower(form = ~ A1 | fSpecies),
  start = start_vec,
  na.action = na.omit
)

# ------------------------------------------------------------------------------
# 5. PREDICTION GRID & THRESHOLD SETUP
# ------------------------------------------------------------------------------
asym_df <- data.frame(
  fSpecies = factor(c("D. pertusum", "Madrepora oculata", "Primnoa msp."),
                    levels = c("D. pertusum", "Madrepora oculata", "Primnoa msp.")),
  asym_val = c(0.1112, 0.0328, 0.0699)
)

size_thresholds <- data.frame(
  fSpecies = factor(c("D. pertusum", "Madrepora oculata", "Primnoa msp."),
                    levels = c("D. pertusum", "Madrepora oculata", "Primnoa msp.")),
  asym_size = c(35.9, 59.9, 187.8),
  label_text = c("35.9 cm²", "59.9 cm²", "187.8 cm²")
)

# Generate smooth prediction curve grid bounded up to 300 cm²
new_data <- expand.grid(
  A1 = seq(min(zuur_growth$A1, na.rm = TRUE), 300, length.out = 300),
  fSpecies = levels(zuur_growth$fSpecies)
)
new_data$pred <- predict(m_gnls_final, newdata = new_data)

# ------------------------------------------------------------------------------
# 6. GRAPHICAL STYLING & ECOGRAPHY THEME SETUP
# ------------------------------------------------------------------------------
species_colours <- c(
  "D. pertusum"        = "#B63679",
  "Madrepora oculata"  = "#440154",
  "Primnoa msp."       = "#FB9F3A"
)

species_labels_expr <- c(
  "D. pertusum"        = expression(italic("D. pertusum")),
  "Madrepora oculata"  = expression(italic("M. oculata")),
  "Primnoa msp."       = expression(italic("Primnoa") ~ "msp.")
)

x_breaks <- c(0, 50, 100, 150, 200, 250, 300)
x_labels <- c("0", "50", "100", "150", "200", "250", "300")

y_breaks <- c(0, 0.0328, 0.0699, 0.1112, 0.2, 0.4)
y_labels <- c("0", "0.033", "0.070", "0.111", "0.2", "0.4")

# ------------------------------------------------------------------------------
# 7. BUILD PUBLICATION PLOT
# ------------------------------------------------------------------------------
cat("Assembling final publication figure...\n")

plot_ecography <- ggplot(zuur_growth, aes(x = A1, y = RGR_Planar, color = fSpecies)) +
  
  # Raw observations
  geom_point(alpha = 0.25, size = 1.1, stroke = 0) +
  
  # Fitted non-linear growth trajectories
  geom_line(data = new_data, aes(y = pred), linewidth = 0.9) +
  
  # Horizontal baseline asymptotes
  geom_hline(
    data = asym_df,
    aes(yintercept = asym_val, color = fSpecies),
    linetype = "dashed", linewidth = 0.4, alpha = 0.8
  ) +
  
  # Vertical ontogenetic threshold markers
  geom_vline(
    data = size_thresholds,
    aes(xintercept = asym_size, color = fSpecies),
    linetype = "dotted", linewidth = 0.5, alpha = 0.85
  ) +
  
  # Rotated threshold value labels anchored near top
  geom_text(
    data = size_thresholds,
    aes(x = asym_size, y = 0.48, label = label_text, color = fSpecies),
    angle = 90,
    vjust = -0.4,
    hjust = 1,
    size = 2.4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  # Scales mapping
  scale_color_manual(
    values = species_colours,
    labels = species_labels_expr
  ) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  
  coord_cartesian(xlim = c(0, 300), ylim = c(0, 0.55)) +
  
  labs(
    x = expression("Initial planar area (cm"^2*")"),
    y = expression("Annualized Relative Growth Rate (RGR yr"^-1*")"),
    color = NULL
  ) +
  
  theme_classic(base_size = 8) +
  theme(
    legend.position   = "top",
    legend.key        = element_blank(),
    legend.key.width  = unit(5, "mm"),
    legend.text       = element_text(size = 7),
    axis.title        = element_text(face = "bold", size = 8),
    axis.text         = element_text(color = "black", size = 6.5),
    axis.line         = element_line(linewidth = 0.35),
    axis.ticks        = element_line(linewidth = 0.35)
  )

# ------------------------------------------------------------------------------
# 8. EXPORT HIGH-RESOLUTION PNG
# ------------------------------------------------------------------------------
dir.create(dirname(output_fig), recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = output_fig, 
  plot     = plot_ecography,
  width    = 85, 
  height   = 70, 
  units    = "mm", 
  dpi      = 600, 
  type     = "cairo"
)

cat("Successfully generated and exported:\n ", output_fig, "\n")
