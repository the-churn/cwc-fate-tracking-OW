# ==============================================================================
# Script Name:  04_Model_Selection.R
# Description:  Partitions master cold-water coral dataset into discrete 
#               demographic vital rates (Growth, Shrinkage, Mortality), performs 
#               formal model selection (GNLS, GLS, GLM) following the Zuur protocol,
#               and exports diagnostic tables for open-science archiving.
# Dependencies: readxl, dplyr, tidyr, nlme, bbmle, openxlsx, ggplot2, patchwork, ragg
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(nlme)
library(bbmle)
library(openxlsx)
library(ggplot2)
library(patchwork)
library(ragg)

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & RELATIVE FILE PATHS (Open-Access / RStudio Project Ready)
# ------------------------------------------------------------------------------
# Utilizing relative paths ensures seamless reproducibility for journal archives 
# (e.g., Current Biology supplementary code bundles).
input_xlsx_path <- file.path("data", "Error_Propagation", "Master_CWC_Tracked_MDC_Age.xlsx")
output_tab_path  <- file.path("outputs", "tables", "Script04_Formal_Vital_Rate_Models.xlsx")
output_fig_path  <- file.path("outputs", "figures", "Script04_Vital_Rate_Diagnostics.png")

FONT_FAMILY  <- "sans"
group_levels <- c("Madrepora oculata", "Desmophyllum pertusum", "Primnoa msp.")

species_palette <- c(
  "Madrepora oculata"     = "#722082",
  "Desmophyllum pertusum" = "#B63679",
  "Primnoa msp."          = "#fb9f3a"
)

# ------------------------------------------------------------------------------
# 3. LOAD DATASET & PREPARE VITAL RATE SUBSETS
# ------------------------------------------------------------------------------
cat("Loading master dataset and partitioning vital rate subsets...\n")
df_raw <- read_excel(input_xlsx_path)

# Standardize species factor levels
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

# Partition 1: Positive Growth (Persistence & RGR >= 0)
df_growth <- df_clean %>%
  filter(state_transition == "persistence" & RGR_Planar >= 0)

# Partition 2: Partial Shrinkage (Persistence, RGR < 0, and detectable change)
df_shrink <- df_clean %>%
  filter(state_transition == "persistence" & RGR_Planar < 0 & detectable_change == TRUE) %>%
  mutate(Shrink_Magnitude = abs(RGR_Planar))

# Partition 3: Total Mortality
df_mortality <- df_clean %>%
  filter(fate != "Recruit") %>%
  mutate(Mortality = ifelse(fate == "Total Mortality", 1, 0))

# ------------------------------------------------------------------------------
# 4. GROWTH VITAL RATE & VARIANCE STRUCTURE SELECTION
# ------------------------------------------------------------------------------
cat("Evaluating variance structures for Growth Vital Rate...\n")

# Fit baseline and candidate variance structures for transparent AIC auditing
m_0        <- gls(RGR_Planar ~ A1 * fSpecies, data = df_growth, method = "ML")
m_ident    <- update(m_0, weights = varIdent(form = ~ 1 | fSpecies))
m_power    <- update(m_0, weights = varPower(form = ~ A1))
m_power_sp <- update(m_0, weights = varPower(form = ~ A1 | fSpecies))
m_exp      <- update(m_0, weights = varExp(form = ~ A1))
m_exp_sp   <- update(m_0, weights = varExp(form = ~ A1 | fSpecies))
m_cp       <- update(m_0, weights = varConstPower(form = ~ A1))
m_cp_sp    <- update(m_0, weights = varConstPower(form = ~ A1 | fSpecies))

# Comprehensive AIC Model Comparison Matrix
model_comparison <- AIC(m_0, m_ident, m_power, m_power_sp, 
                        m_exp, m_exp_sp, m_cp, m_cp_sp)

rownames(model_comparison) <- c("Baseline", "VarIdent", "VarPower_Global", 
                                "VarPower_Sp", "VarExp_Global", "VarExp_Sp", 
                                "VarConstPower_Global", "VarConstPower_Sp")

model_comparison <- model_comparison[order(model_comparison$AIC), ]
print(model_comparison)

# ------------------------------------------------------------------------------
# TRANSPARENCY NOTE FOR REVIEWERS:
# Although VarConstPower_Sp yielded a marginally lower AIC (-2171.496) compared 
# to VarPower_Sp (-2169.354), the difference is negligible (delta AIC = 2.1) 
# while consuming 3 additional degrees of freedom (df = 13 vs df = 10). 
# Following Occam's razor and the law of parsimony, VarPower_Sp is selected 
# as the final optimal variance structure to avoid over-parameterization.
# ------------------------------------------------------------------------------

# Fitting the Non-Linear Asymptotic (GNLS) model with VarPower_Sp weighting
m_list <- nlsList(
  RGR_Planar ~ SSasymp(A1, Asym, R0, lrc) | fSpecies, 
  data = df_growth, 
  na.action = na.omit
)
start_vec <- c(coef(m_list)$Asym, coef(m_list)$R0, coef(m_list)$lrc)

m_growth_gnls <- gnls(
  RGR_Planar ~ SSasymp(A1, Asym, R0, lrc),
  params = list(Asym ~ fSpecies, R0 ~ fSpecies, lrc ~ fSpecies),
  data = df_growth,
  weights = varPower(form = ~ A1 | fSpecies),
  start = start_vec,
  na.action = na.omit
)

# ------------------------------------------------------------------------------
# 5. SHRINKAGE & MORTALITY VITAL RATE MODELS
# ------------------------------------------------------------------------------
cat("Fitting Partial Shrinkage and Total Mortality models...\n")

# Shrinkage GLS
m_shrink_ols <- gls(Shrink_Magnitude ~ A1 * fSpecies, data = df_shrink)
m_shrink_exp <- gls(
  Shrink_Magnitude ~ A1 * fSpecies,
  data = df_shrink,
  weights = varExp(form = ~ A1 | fSpecies),
  na.action = na.omit
)
aic_shrink <- bbmle::AICtab(m_shrink_ols, m_shrink_exp, weights = TRUE, logLik = TRUE)
tab_shrink <- as.data.frame(aic_shrink) %>%
  mutate(Model = rownames(.), Process = "Shrinkage (RGR < 0)") %>%
  select(Process, Model, dAIC, df, weight)

# Mortality GLM
m_mort_logit <- glm(Mortality ~ A1 * fSpecies, data = df_mortality, family = binomial(link = "logit"))
m_mort_cloglog <- glm(Mortality ~ A1 * fSpecies, data = df_mortality, family = binomial(link = "cloglog"))
aic_mort <- bbmle::AICtab(m_mort_logit, m_mort_cloglog, weights = TRUE, logLik = TRUE)
tab_mort <- as.data.frame(aic_mort) %>%
  mutate(Model = rownames(.), Process = "Total Mortality") %>%
  select(Process, Model, dAIC, df, weight)

# ------------------------------------------------------------------------------
# 6. EXPORT DIAGNOSTICS & SUMMARY TABLES
# ------------------------------------------------------------------------------
dir.create(dirname(output_tab_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_fig_path), recursive = TRUE, showWarnings = FALSE)

cat("Saving formal model selection summary tables...\n")
wb <- createWorkbook()
addWorksheet(wb, "Growth_Variance_Selection")
addWorksheet(wb, "Shrinkage_Selection")
addWorksheet(wb, "Mortality_Selection")

writeData(wb, "Growth_Variance_Selection", as.data.frame(model_comparison) %>% tibble::rownames_to_column("Model"))
writeData(wb, "Shrinkage_Selection", tab_shrink)
writeData(wb, "Mortality_Selection", tab_mort)

saveWorkbook(wb, output_tab_path, overwrite = TRUE)

# Generate Diagnostic Graphic for Growth GNLS Normalized Residuals
res_df <- data.frame(
  Fitted    = fitted(m_growth_gnls),
  Residuals = residuals(m_growth_gnls, type = "normalized"),
  Species   = df_growth$fSpecies[!is.na(df_growth$RGR_Planar)]
)

p_res <- ggplot(res_df, aes(x = Fitted, y = Residuals, color = Species)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = species_palette) +
  labs(x = "Fitted Values", y = "Normalized Residuals", title = "Growth GNLS Model") +
  theme_classic(base_family = FONT_FAMILY)

p_qq <- ggplot(res_df, aes(sample = Residuals)) +
  stat_qq(color = "#404040", alpha = 0.6) +
  stat_qq_line(color = "red", linetype = "dashed") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = "Q-Q Normal Plot") +
  theme_classic(base_family = FONT_FAMILY)

diag_composite <- p_res + p_qq + plot_layout(guides = "collect")

ggsave(
  filename = output_fig_path,
  plot     = diag_composite,
  width    = 174,
  height   = 85,
  units    = "mm",
  dpi      = 600,
  device   = ragg::agg_png
)

cat("\n======================================================================\n")
cat("Script 04 execution complete successfully!\n")
cat("Summary tables saved to:\n ", output_tab_path, "\n")
cat("Diagnostics saved to:\n ", output_fig_path, "\n")
cat("======================================================================\n")
