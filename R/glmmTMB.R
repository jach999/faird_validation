# ==============================================================================
# glmmTMB.R — Temporal biomass analysis (glmmTMB)
# Replicates Scherber's model: biomass ~ bs(Date,3) * Device_type + (1|Site/Device)
# Produces: Type III ANOVA table + temporal trend plot with CI
# ==============================================================================

library(tidyverse)
library(glmmTMB)
library(car)
library(splines)
library(ggplot2)

source("R/00_load_data.R")
source("R/00_custom_colors.R")

dir.create("outputs/figures",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables",   showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/console",  showWarnings = FALSE, recursive = TRUE)

while (sink.number() > 0) sink()
sink("outputs/console/glmmTMB_output.txt", split = TRUE)

cat("==============================================================================\n")
cat("SCHERBER TEMPORAL BIOMASS ANALYSIS\n")
cat("glmmTMB: biomass ~ bs(Date.num, 3) * Device_type + (1|Site/Device)\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# ==============================================================================
# SECTION 1: DATA PREPARATION
# ==============================================================================

cat("=== DATA PREPARATION ===\n\n")

# FAIRD and ID daily biomass from etraps_daily
etraps_bio <- data_list$etraps_daily %>%
  dplyr::filter(Device_type %in% c("FAIRD", "ID")) %>%
  dplyr::select(Date, Device, Device_type, Site, biomass_mg) %>%
  dplyr::rename(Biomass_mg = biomass_mg)

# AMMOD daily biomass — rename to match structure
ammod_bio <- data_list$ammod_biomass_for_stats %>%
  dplyr::mutate(
    Device     = Device,
    Device_type = "AMMOD",
    Biomass_mg  = Live_mass
  ) %>%
  dplyr::select(Date, Device, Device_type, Site, Biomass_mg)

# Combine into single dataset
bio_combined <- dplyr::bind_rows(etraps_bio, ammod_bio) %>%
  dplyr::arrange(Device_type, Site, Date) %>%
  dplyr::mutate(
    Date.num     = as.numeric(Date),
    Device_type  = factor(Device_type, levels = c("AMMOD", "FAIRD", "ID")),
    Site         = factor(Site),
    Device       = factor(Device)
  )

# Restrict to overlapping date range (AMMOD only has Batch 1: 12 days)
date_range <- range(ammod_bio$Date)
bio_overlap <- bio_combined %>%
  dplyr::filter(Date >= date_range[1] & Date <= date_range[2])

cat(sprintf("Combined dataset: %d observations\n", nrow(bio_combined)))
cat(sprintf("Overlap period (%s to %s): %d observations\n",
            date_range[1], date_range[2], nrow(bio_overlap)))
cat(sprintf("Device types: %s\n", paste(levels(bio_combined$Device_type), collapse = ", ")))

cat("\nObservations per Device_type (overlap period):\n")
print(bio_overlap %>%
        dplyr::group_by(Device_type) %>%
        dplyr::summarise(n = dplyr::n(), 
                         n_devices = dplyr::n_distinct(Device),
                         mean_biomass = round(mean(Biomass_mg, na.rm = TRUE), 1),
                         .groups = "drop"))

cat("\nObservations per Device (overlap period):\n")
print(bio_overlap %>%
        dplyr::group_by(Device_type, Device, Site) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop"))


# ==============================================================================
# SECTION 2: MODEL FITTING
# ==============================================================================

cat("\n=== MODEL FITTING ===\n\n")

# --- Model A: Full model (overlap period only, all 3 device types) ---
cat("--- Model A: Overlap period (AMMOD + FAIRD + ID) ---\n")
cat("Formula: Biomass_mg ~ bs(Date.num, 3) * Device_type + (1 | Site/Device)\n")
cat("Family: Tweedie\n\n")



bio_overlap$Block=factor(ifelse(bio_overlap$Site%in%c("Site1","Site2"),"BLOCK1","BLOCK2"))
bio_overlap


bio_overlap$times=factor(bio_overlap$Date)


#control=glmmTMBControl(optimizer="optim",optArgs=list(method="BFGS")

library(nlme)

var(bio_overlap$Biomass_mg) # total variance is 17.7 Mio mg^2, due to unit of measurement

#standard mixed model using nlme:

BX <- lme(
  Biomass_mg ~ bs(Date.num, 3) * Device_type, random=~Date.num | Block/Site/Device,
  correlation=corAR1(),
  data   = bio_overlap,control=list(opt="optim"))


summary(BX)
library(car)
Anova(BX)

# standard nlme model allows "easier" inspection of random effects,
# including the so-called augmentd predictions:
plot(ranef(BX))
plot(augPred(BX,primary=~Date.num))

# this shows the coefficients, i.e. random effects added on top of the fixed effects
# the so called Best Linear Unbiased Predictions (BLUPs)
coef(BX)

install.packages("ape")
library(ape)

# nice plot of variance components from that model
plot(varcomp(BX))


####################################################

B0 <- glmmTMB(
  Biomass_mg ~ bs(Date.num, 3) * Device_type + (1 | Block/Site/Device)
  ,family="tweedie",
  data   = bio_overlap)

B1 <- glmmTMB(
  Biomass_mg ~ bs(Date.num, 3) * Device_type + (1 | Block/Site/Device)+
    ar1(times+0|Device),family="tweedie",
  data   = bio_overlap)

AIC(B0,B1) #B1 slightly better than B0

summary(B1)

#plot(varcomp(B1)) doesn´t work because model is not of class nlme


#remotes::install_github("SchmidtPaul/BioMathR")

VarCorr(B1)
library(BioMathR)

remotes::install_github("m-clark/mixedup")
library(mixedup)

# the below code gives the variance components (useful figures)
get_varcomp(B0)
get_varcomp(B1)

extract_vc(B1)


install.packages("sjPlot",dependencies=T)
library(sjPlot)


plot_model(B1,type="re",transform=NULL)
plot_model(B0,type="re",transform=NULL)




# Summary
cat("--- Model summary ---\n")
print(summary(B1))

# Type III ANOVA
cat("\n--- Type III ANOVA (car::Anova) ---\n")
anova_B1 <- car::Anova(B1, type = "III")
print(anova_B1)

# Extract key results
cat("\n--- Interpretation ---\n")
p_device   <- anova_B1["Device_type", "Pr(>Chisq)"]
p_temporal <- anova_B1["bs(Date.num, 3)", "Pr(>Chisq)"]
p_interact <- anova_B1["bs(Date.num, 3):Device_type", "Pr(>Chisq)"]

cat(sprintf("  Device_type effect:              p = %.4f %s\n", p_device,
            ifelse(p_device < 0.05, "-> Significant biomass hierarchy", "-> NS")))
cat(sprintf("  Temporal effect (spline):        p = %.4f %s\n", p_temporal,
            ifelse(p_temporal < 0.05, "-> Significant temporal trend", "-> NS")))
cat(sprintf("  Date x Device_type interaction:  p = %.4f %s\n", p_interact,
            ifelse(p_interact < 0.05, "-> Trajectories NOT parallel", "-> Trajectories parallel")))


# --- Model B: FAIRD vs AMMOD only (no ID) ---
cat("\n\n--- Model B: FAIRD vs AMMOD only ---\n")

bio_FA <- bio_overlap %>%
  dplyr::filter(Device_type %in% c("AMMOD", "FAIRD")) %>%
  dplyr::mutate(Device_type = droplevels(Device_type))

B2 <- glmmTMB::glmmTMB(
  Biomass_mg ~ bs(Date.num, 3) * Device_type + (1 | Site/Device),
  family = glmmTMB::tweedie(),
  data   = bio_FA
)

cat("\n--- Type III ANOVA (FAIRD vs AMMOD only) ---\n")
anova_B2 <- car::Anova(B2, type = "III")
print(anova_B2)

p_interact_FA <- anova_B2["bs(Date.num, 3):Device_type", "Pr(>Chisq)"]
cat(sprintf("\n  Date x Device_type (FAIRD vs AMMOD): p = %.4f %s\n", p_interact_FA,
            ifelse(p_interact_FA < 0.05, "-> Trajectories NOT parallel",
                   "-> Trajectories parallel")))


# ==============================================================================
# SECTION 3: TEMPORAL TREND PLOT (replicating Scherber's figure)
# ==============================================================================

cat("\n\n=== GENERATING TEMPORAL PLOT ===\n\n")

# Prediction grid
date_seq <- seq(min(bio_overlap$Date.num), max(bio_overlap$Date.num), length.out = 200)

pred_grid <- expand.grid(
  Date.num    = date_seq,
  Device_type = levels(bio_overlap$Device_type)
) %>%
  dplyr::mutate(
    Site   = NA,
    Device = NA
  )

# Predictions (fixed effects only, no random effects)
pred_grid$fit <- predict(B1, newdata = pred_grid, type = "response",
                         re.form  = NA,  # fixed effects only
                         se.fit   = FALSE)

# SE for confidence intervals
pred_link <- predict(B1, newdata = pred_grid, type = "link",
                     re.form = NA, se.fit = TRUE)

# Back-transform CI from link to response scale
pred_grid$fit_link <- pred_link$fit
pred_grid$se_link  <- pred_link$se.fit
pred_grid$lower    <- B1$modelInfo$family$linkinv(pred_grid$fit_link - 1.96 * pred_grid$se_link)
pred_grid$upper    <- B1$modelInfo$family$linkinv(pred_grid$fit_link + 1.96 * pred_grid$se_link)

# Convert Date.num back to Date
pred_grid$Date <- as.Date(pred_grid$Date.num, origin = "1970-01-01")
bio_overlap$Date_plot <- bio_overlap$Date

# Device type colors
device_colors <- c(
  "AMMOD" = "#E74C3C",   # red
  "FAIRD" = "#D4A017",   # golden/orange
  "ID"    = "#27AE60"    # green
)

# Plot
p_temporal <- ggplot2::ggplot() +
  # Confidence ribbons
  ggplot2::geom_ribbon(
    data = pred_grid,
    ggplot2::aes(x = Date, ymin = lower, ymax = upper, fill = Device_type),
    alpha = 0.15
  ) +
  # Fitted lines
  ggplot2::geom_line(
    data = pred_grid,
    ggplot2::aes(x = Date, y = fit, color = Device_type),
    linewidth = 1
  ) +
  # Raw data points
  ggplot2::geom_point(
    data = bio_overlap,
    ggplot2::aes(x = Date_plot, y = Biomass_mg, color = Device_type),
    alpha = 0.5, size = 1.8
  ) +
  # Rug plot (data density)
  ggplot2::geom_rug(
    data = bio_overlap,
    ggplot2::aes(x = Date_plot),
    sides = "b", alpha = 0.3, length = ggplot2::unit(0.02, "npc")
  ) +
  ggplot2::scale_color_manual(values = device_colors) +
  ggplot2::scale_fill_manual(values = device_colors) +
  ggplot2::scale_y_log10(
    labels = scales::label_comma(),
    breaks = c(10, 50, 100, 250, 500, 1000, 2000, 5000)
  ) +
  ggplot2::labs(
    title    = "Temporal Biomass Trends by Device Type",
    subtitle = sprintf("glmmTMB: Tweedie, bs(Date,3) × Device_type + (1|Site/Device) | Date×Device interaction p = %.4f | log10 y-axis",
                       p_interact),
    x        = "Date",
    y        = "Biomass (mg) — log scale",
    color    = "Device_type",
    fill     = "Device_type"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title        = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle     = ggplot2::element_text(size = 10, color = "grey40"),
    panel.grid.minor  = ggplot2::element_blank(),
    legend.position   = "top"
  )

print(p_temporal)
ggplot2::ggsave("outputs/figures/glmmTMB.png", p_temporal,
                width = 10, height = 7, dpi = 300, bg = "white")
cat("[OK] Saved: outputs/figures/glmmTMB.png\n")


# --- FAIRD vs AMMOD only version ---
pred_FA <- pred_grid %>% dplyr::filter(Device_type %in% c("AMMOD", "FAIRD"))
bio_FA_plot <- bio_overlap %>% dplyr::filter(Device_type %in% c("AMMOD", "FAIRD"))

p_temporal_FA <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data = pred_FA,
    ggplot2::aes(x = Date, ymin = lower, ymax = upper, fill = Device_type),
    alpha = 0.15
  ) +
  ggplot2::geom_line(
    data = pred_FA,
    ggplot2::aes(x = Date, y = fit, color = Device_type),
    linewidth = 1
  ) +
  ggplot2::geom_point(
    data = bio_FA_plot,
    ggplot2::aes(x = Date_plot, y = Biomass_mg, color = Device_type),
    alpha = 0.5, size = 1.8
  ) +
  ggplot2::scale_color_manual(values = device_colors) +
  ggplot2::scale_fill_manual(values = device_colors) +
  ggplot2::scale_y_log10(
    labels = scales::label_comma(),
    breaks = c(10, 50, 100, 250, 500, 1000, 2000, 5000)
  ) +
  ggplot2::labs(
    title    = "Temporal Biomass Trends: FAIRD vs AMMOD",
    subtitle = sprintf("Date×Device interaction p = %.4f | log10 y-axis", p_interact_FA),
    x        = "Date",
    y        = "Biomass (mg) — log scale",
    color    = "Device_type",
    fill     = "Device_type"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title        = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle     = ggplot2::element_text(size = 10, color = "grey40"),
    panel.grid.minor  = ggplot2::element_blank(),
    legend.position   = "top"
  )

print(p_temporal_FA)
ggplot2::ggsave("outputs/figures/glmmTMB_FA_only.png", p_temporal_FA,
                width = 10, height = 7, dpi = 300, bg = "white")
cat("[OK] Saved: outputs/figures/glmmTMB_FA_only.png\n")


# ==============================================================================
# SECTION 4: EXPORT
# ==============================================================================

cat("\n=== EXPORT ===\n\n")

# ANOVA tables
anova_df_B1 <- as.data.frame(anova_B1) %>%
  tibble::rownames_to_column("Term") %>%
  dplyr::mutate(Model = "All_3_devices")

anova_df_B2 <- as.data.frame(anova_B2) %>%
  tibble::rownames_to_column("Term") %>%
  dplyr::mutate(Model = "FAIRD_vs_AMMOD")

anova_export <- dplyr::bind_rows(anova_df_B1, anova_df_B2)
readr::write_csv(anova_export, "outputs/tables/Scherber_glmmTMB_ANOVA_III.csv")
cat("[OK] Saved: outputs/tables/Scherber_glmmTMB_ANOVA_III.csv\n")


p_temporal_nolog <- p_temporal +
  ggplot2::scale_y_continuous(labels = scales::label_comma()) +
  ggplot2::labs(
    y        = "Biomass (mg)",
    subtitle = sprintf("glmmTMB: Tweedie, bs(Date,3) × Device_type + (1|Site/Device) | Date×Device interaction p = %.4f | linear y-axis",
                       p_interact)
  )

print(p_temporal_nolog)
ggplot2::ggsave("outputs/figures/glmmTMB_linear.png", p_temporal_nolog,
                width = 10, height = 7, dpi = 300, bg = "white")

scale_factor <- 10

bio_plot <- bio_overlap %>%
  dplyr::mutate(
    Biomass_plot = dplyr::if_else(Device_type == "AMMOD", 
                                  Biomass_mg / scale_factor, 
                                  Biomass_mg)
  )

pred_plot <- pred_grid %>%
  dplyr::mutate(
    fit_plot   = ifelse(Device_type == "AMMOD", fit / scale_factor, fit),
    lower_plot = ifelse(Device_type == "AMMOD", lower / scale_factor, lower),
    upper_plot = ifelse(Device_type == "AMMOD", upper / scale_factor, upper)
  )

# Find sensible y limit (exclude extreme ID outliers)
y_max <- 1200  # covers FAIRD ~995 and most ID points

p_dual2 <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data = pred_plot,
    ggplot2::aes(x = Date, ymin = lower_plot, ymax = upper_plot, fill = Device_type),
    alpha = 0.15
  ) +
  ggplot2::geom_line(
    data = pred_plot,
    ggplot2::aes(x = Date, y = fit_plot, color = Device_type),
    linewidth = 1
  ) +
  ggplot2::geom_point(
    data = bio_plot,
    ggplot2::aes(x = Date, y = Biomass_plot, color = Device_type),
    alpha = 0.5, size = 1.8
  ) +
  ggplot2::scale_color_manual(values = device_colors) +
  ggplot2::scale_fill_manual(values = device_colors) +
  ggplot2::scale_y_continuous(
    name   = "FAIRD / ID — Biomass (mg)",
    labels = scales::label_comma(),
    limits = c(0, y_max),
    sec.axis = ggplot2::sec_axis(
      ~ . * scale_factor,
      name   = "AMMOD — Biomass (mg)",
      labels = scales::label_comma()
    )
  ) +
  ggplot2::labs(
    title    = "Temporal Biomass Trends by Device Type",
    subtitle = sprintf("glmmTMB: Tweedie | Date×Device interaction p = %.4f | AMMOD on right axis (×10)",
                       p_interact),
    x        = "Date",
    color    = "Device_type",
    fill     = "Device_type"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title            = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle         = ggplot2::element_text(size = 10, color = "grey40"),
    panel.grid.minor      = ggplot2::element_blank(),
    legend.position       = "top",
    axis.title.y.left     = ggplot2::element_text(color = "grey30"),
    axis.title.y.right    = ggplot2::element_text(color = "#E74C3C"),
    axis.text.y.right     = ggplot2::element_text(color = "#E74C3C")
  )

print(p_dual2)
ggplot2::ggsave("outputs/figures/glmmTMB_dualaxis.png", p_dual2,
                width = 10, height = 7, dpi = 300, bg = "white")

cat("\n==============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat(sprintf("Finished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("==============================================================================\n")

while (sink.number() > 0) sink()