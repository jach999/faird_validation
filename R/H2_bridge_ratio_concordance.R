# ==============================================================================
# H2_bridge_ratio_concordance.R
# Exploratory: does the AMMOD/FAIRD biomass ratio move with the AMMOD-FAIRD
# taxonomic concordance (Izsak-Price Delta S, Superfamily level)?
#
# Revision of H2_bridge_exploratory.R. Changes and their reasons:
#   1. Site-day pairing BEFORE aggregation. A unit only uses days on which the
#      site has BOTH AMMOD biomass and FAIRD data, for the ratio AND for the
#      taxonomy. The previous version summed 2 AMMOD days against 3 FAIRD days
#      (Site3, Aug 29-31) and compared AMMOD3 alone with FAIRD3+FAIRD4 at the
#      Meadow and Overall levels.
#   2. Within-site analysis only. The question is temporal (does the ratio
#      change when concordance changes), so between-site differences in level
#      are removed by centring each variable on its site mean. Habitat and
#      Overall levels are dropped: once paired, they leave 4 units, and a
#      Spearman test with n = 4 cannot reach p < 0.05.
#   3. Within-site permutation test. Units are shuffled only within their own
#      site, which is the correct null for a within-site association. Exact
#      enumeration when feasible, Monte Carlo otherwise.
#   4. FAIRD sample size as a covariate. A unit with few FAIRD individuals
#      tends to give a high ratio AND a low Delta S at the same time, so the
#      ratio-concordance association is reported next to the association of
#      each variable with FAIRD sample size.
#   5. Two resolutions: non-overlapping 3-day windows (primary, the H2
#      qualitative unit) and single days (secondary, the resolution of the
#      original question). Both use the same inclusion threshold.
#   6. No hardcoded dates or device exclusions: the window comes from the
#      AMMOD biomass data and the pairing removes AMMOD4 automatically.
#
# Statistics are computed on log(ratio) (multiplicative scale, consistent with
# the log link of the temporal model). Figures use LINEAR axes with the raw
# ratio, following the presentation rule of the project (no log axes).
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(vegan)

source("R/00_load_data.R")

dir.create("outputs/console", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables",  showWarnings = FALSE, recursive = TRUE)

script_name <- "H2_bridge_ratio_concordance"
output_file <- file.path("outputs/console",
                         paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
                                script_name, "_output.txt"))
while (sink.number() > 0) sink()
sink(output_file, split = TRUE)

cat("==============================================================================\n")
cat("H2 BRIDGE: AMMOD/FAIRD BIOMASS RATIO vs TAXONOMIC CONCORDANCE (EXPLORATORY)\n")
cat("Site-day paired units | within-site centring | within-site permutation tests\n")
cat("==============================================================================\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Parameters ---------------------------------------------------------------
nc             <- "#N/C"   # project code for a non-classified taxonomic level
MIN_INDIV      <- 3        # minimum FAIRD individuals with Superfamily per unit
WINDOW_DAYS    <- 3        # window length, same as the H2 qualitative branch
N_MC           <- 9999     # Monte Carlo permutations when exact is not feasible
EXACT_MAX      <- 50000    # enumerate exactly up to this many permutations
SEED           <- 20231001


# ==============================================================================
# 1. SITE-DAY PAIRED BIOMASS
# ==============================================================================
cat("=== 1. SITE-DAY PAIRING ===\n\n")

ammod_day <- data_list$ammod_biomass_for_stats %>%
  dplyr::select(Date, Site, AMMOD = Live_mass)

faird_day <- data_list$etraps_daily %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::select(Date, Site, FAIRD = biomass_mg)

bio_pair <- dplyr::inner_join(ammod_day, faird_day, by = c("Date", "Site")) %>%
  dplyr::filter(!is.na(AMMOD), !is.na(FAIRD)) %>%
  dplyr::arrange(Site, Date)

all_dates <- seq(min(bio_pair$Date), max(bio_pair$Date), by = "day")
window_lookup <- tibble::tibble(
  Date   = all_dates,
  Window = as.integer((as.numeric(all_dates - min(all_dates)) %/% WINDOW_DAYS) + 1)
)
bio_pair <- bio_pair %>% dplyr::left_join(window_lookup, by = "Date")

cat(sprintf("Paired site-days: %d\n", nrow(bio_pair)))
print(bio_pair %>%
        dplyr::group_by(Site) %>%
        dplyr::summarise(n_days = dplyr::n(), from = min(Date), to = max(Date),
                         .groups = "drop"))
cat("\nWindow definition:\n")
print(window_lookup %>%
        dplyr::group_by(Window) %>%
        dplyr::summarise(from = min(Date), to = max(Date), .groups = "drop"))


# ==============================================================================
# 2. TAXONOMY AND IZSAK-PRICE DISTANCE MATRIX (Superfamily)
# ==============================================================================
cat("\n=== 2. TAXONOMY ===\n\n")

clean_sf <- function(x) dplyr::if_else(!is.na(x) & x != nc & x != "", x, NA_character_)

faird_tax <- data_list$etraps %>%
  dplyr::filter(Device_type == "FAIRD") %>%
  dplyr::mutate(Taxon_SF = clean_sf(Superfamily)) %>%
  dplyr::filter(!is.na(Taxon_SF))

ammod_tax <- data_list$ammod_taxonomy %>%
  dplyr::mutate(Taxon_SF = clean_sf(Superfamily)) %>%
  dplyr::filter(!is.na(Taxon_SF))

all_sf_taxa <- dplyr::bind_rows(faird_tax, ammod_tax) %>%
  dplyr::distinct(Taxon_SF, .keep_all = TRUE) %>%
  dplyr::select(Taxon_SF, Class, Order, Suborder, Infraorder, Superfamily) %>%
  as.data.frame()
rownames(all_sf_taxa) <- all_sf_taxa$Taxon_SF
all_sf_taxa$Taxon_SF  <- NULL
dist_matrix_sf <- as.matrix(vegan::taxa2dist(all_sf_taxa, varstep = TRUE, check = TRUE))
cat(sprintf("Superfamilies in distance matrix: %d\n", nrow(dist_matrix_sf)))

# Izsak-Price similarity, same definition as H2_bridge_exploratory.R.
calc_ip <- function(taxa_A, taxa_B, dist_mat) {
  taxa_A <- taxa_A[taxa_A %in% rownames(dist_mat)]
  taxa_B <- taxa_B[taxa_B %in% rownames(dist_mat)]
  if (length(taxa_A) == 0 || length(taxa_B) == 0) return(NA_real_)
  min_A_to_B <- vapply(taxa_A, function(a)
    if (a %in% taxa_B) 0 else min(dist_mat[a, taxa_B], na.rm = TRUE), numeric(1))
  min_B_to_A <- vapply(taxa_B, function(b)
    if (b %in% taxa_A) 0 else min(dist_mat[b, taxa_A], na.rm = TRUE), numeric(1))
  1 - (mean(c(min_A_to_B, min_B_to_A)) / max(dist_mat, na.rm = TRUE))
}


# ==============================================================================
# 3. BUILD ANALYSIS UNITS (window or day)
# ==============================================================================
# For each Site x unit: ratio over the paired days of the unit, taxonomy of
# BOTH systems restricted to exactly those days and that site.
build_units <- function(unit_col) {
  units <- bio_pair %>%
    dplyr::mutate(Unit = .data[[unit_col]]) %>%
    dplyr::group_by(Site, Unit) %>%
    dplyr::summarise(
      first_date = min(Date),
      n_days     = dplyr::n(),
      dates      = list(Date),
      AMMOD_mg   = sum(AMMOD),
      FAIRD_mg   = sum(FAIRD),
      .groups    = "drop"
    )

  rows <- lapply(seq_len(nrow(units)), function(i) {
    u <- units[i, ]
    d <- u$dates[[1]]
    f <- faird_tax %>% dplyr::filter(Site == u$Site, Date %in% d)
    a <- ammod_tax %>% dplyr::filter(Site == u$Site, Date %in% d)
    f_sf <- unique(f$Taxon_SF)
    a_sf <- unique(a$Taxon_SF)
    ok   <- nrow(f) >= MIN_INDIV && length(a_sf) > 0
    tibble::tibble(
      n_FAIRD = nrow(f),
      S_FAIRD = length(f_sf),
      S_AMMOD = length(a_sf),
      DS_pct  = if (ok) 100 * calc_ip(f_sf, a_sf, dist_matrix_sf) else NA_real_
    )
  })

  units %>%
    dplyr::select(-dates) %>%
    dplyr::bind_cols(dplyr::bind_rows(rows)) %>%
    dplyr::mutate(
      ratio    = dplyr::if_else(FAIRD_mg > 0, AMMOD_mg / FAIRD_mg, NA_real_),
      included = !is.na(ratio) & !is.na(DS_pct)
    )
}

# Within-site centring of the analysis variables (included units only).
# Sites with fewer than 2 included units carry no within-site information.
centre_within_site <- function(units) {
  units %>%
    dplyr::filter(included) %>%
    dplyr::group_by(Site) %>%
    dplyr::filter(dplyr::n() >= 2) %>%
    dplyr::mutate(
      c_log_ratio = log(ratio)   - mean(log(ratio)),
      c_DS        = DS_pct       - mean(DS_pct),
      c_log_nF    = log(n_FAIRD) - mean(log(n_FAIRD))
    ) %>%
    dplyr::ungroup()
}


# ==============================================================================
# 4. WITHIN-SITE PERMUTATION TEST
# ==============================================================================
# Statistic: Pearson r between within-site-centred variables. Null: y values
# shuffled within each site (centring is preserved by such shuffles).
perm_matrix <- function(n) {
  if (n == 1) return(matrix(1L, 1, 1))
  sub <- perm_matrix(n - 1)
  do.call(rbind, lapply(seq_len(n), function(i) {
    cbind(i, matrix(setdiff(seq_len(n), i)[sub], nrow = nrow(sub)))
  }))
}

within_site_perm <- function(x, y, site) {
  r_obs <- stats::cor(x, y)
  idx   <- split(seq_along(x), droplevels(factor(site)))
  n_tot <- prod(vapply(idx, function(ix) factorial(length(ix)), numeric(1)))

  if (n_tot <= EXACT_MAX) {
    pm   <- lapply(idx, function(ix) perm_matrix(length(ix)))
    grid <- expand.grid(lapply(pm, function(m) seq_len(nrow(m))))
    r_perm <- apply(grid, 1, function(g) {
      yp <- y
      for (s in seq_along(idx)) yp[idx[[s]]] <- y[idx[[s]]][pm[[s]][g[s], ]]
      stats::cor(x, yp)
    })
    p <- mean(abs(r_perm) >= abs(r_obs) - 1e-12)
    method <- sprintf("exact (%d permutations)", length(r_perm))
  } else {
    set.seed(SEED)
    r_perm <- replicate(N_MC, {
      yp <- y
      for (ix in idx) yp[ix] <- y[ix][sample.int(length(ix))]
      stats::cor(x, yp)
    })
    p <- (sum(abs(r_perm) >= abs(r_obs) - 1e-12) + 1) / (N_MC + 1)
    method <- sprintf("Monte Carlo (%d permutations)", N_MC)
  }
  list(r = r_obs, p = p, method = method)
}

test_pairs <- list(
  c("c_log_ratio", "c_DS",     "log(ratio) vs Delta S"),
  c("c_log_ratio", "c_log_nF", "log(ratio) vs log(n FAIRD)"),
  c("c_DS",        "c_log_nF", "Delta S vs log(n FAIRD)")
)

run_resolution <- function(unit_col, label) {
  cat(sprintf("\n==============================================================================\n"))
  cat(sprintf("RESOLUTION: %s\n", label))
  cat(sprintf("==============================================================================\n\n"))

  units <- build_units(unit_col)
  cat("All units (included = ratio defined and Delta S computable):\n")
  print(units %>%
          dplyr::mutate(dplyr::across(c(AMMOD_mg, FAIRD_mg), ~ round(.x, 1)),
                        ratio = round(ratio, 1), DS_pct = round(DS_pct, 1)),
        n = Inf, width = Inf)

  cen <- centre_within_site(units)
  cat(sprintf("\nIncluded units used in within-site tests: %d (sites: %s)\n",
              nrow(cen), paste(sort(unique(as.character(cen$Site))), collapse = ", ")))
  print(cen %>% dplyr::count(Site, name = "n_units"))

  res <- dplyr::bind_rows(lapply(test_pairs, function(tp) {
    t <- within_site_perm(cen[[tp[1]]], cen[[tp[2]]], cen$Site)
    tibble::tibble(Resolution = label, Association = tp[3], n = nrow(cen),
                   r_within = round(t$r, 3), p_perm = round(t$p, 4),
                   Method = t$method)
  }))
  cat("\nWithin-site associations (Pearson r on site-centred values):\n")
  print(res, width = Inf)

  list(units = units, centred = cen, tests = res)
}

res_win <- run_resolution("Window", sprintf("%d-day windows", WINDOW_DAYS))
res_day <- run_resolution("Date",   "single days")

all_tests <- dplyr::bind_rows(res_win$tests, res_day$tests)
readr::write_csv(all_tests, "outputs/tables/H2_bridge_ratio_concordance_tests.csv")
cat("\n[OK] Saved: outputs/tables/H2_bridge_ratio_concordance_tests.csv\n")

cat("\nHOW TO READ THIS\n")
cat("- Row 1 is the bridge question. Rows 2 and 3 check the shared driver: if\n")
cat("  both ratio and Delta S track FAIRD sample size, an association in row 1\n")
cat("  cannot be read as concordance driving the ratio.\n")
cat("- The expected sign under the bridge hypothesis is negative (lower\n")
cat("  concordance, higher AMMOD/FAIRD ratio). The sampling artefact predicts\n")
cat("  the same sign, which is why rows 2 and 3 are needed.\n")


# ==============================================================================
# 5. FIGURES (linear axes, raw ratio)
# ==============================================================================
plot_bridge <- function(cen, title, label_units) {
  p <- ggplot2::ggplot(cen, ggplot2::aes(x = DS_pct, y = ratio)) +
    ggplot2::geom_path(ggplot2::aes(group = Site), color = "grey60",
                       arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"))) +
    ggplot2::geom_point(ggplot2::aes(color = as.numeric(first_date - min(first_date))),
                        size = 2.8) +
    ggplot2::facet_wrap(~ Site, nrow = 1) +
    ggplot2::scale_color_viridis_c(name = "Days since start") +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::labs(
      title    = title,
      subtitle = "Within each site, points joined in time order. Linear axes.",
      x = "Taxonomic concordance, Izsak-Price Delta S (%)",
      y = "AMMOD / FAIRD biomass ratio",
      caption = "Exploratory. Units use only site-days with both systems present."
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40"),
      plot.caption     = ggplot2::element_text(size = 8, color = "grey50"),
      strip.text       = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top"
    )
  if (label_units) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = format(first_date, "%d %b")),
                                vjust = -1, size = 3, color = "grey30")
  }
  p
}

p_win <- plot_bridge(res_win$centred,
                     sprintf("Biomass ratio vs taxonomic concordance (%d-day windows)", WINDOW_DAYS),
                     label_units = TRUE)
p_day <- plot_bridge(res_day$centred,
                     "Biomass ratio vs taxonomic concordance (single days)",
                     label_units = FALSE)

print(p_win)
ggplot2::ggsave("outputs/figures/H2_bridge_ratio_concordance_windows.png", p_win,
                width = 12, height = 5, dpi = 300, bg = "white")
print(p_day)
ggplot2::ggsave("outputs/figures/H2_bridge_ratio_concordance_days.png", p_day,
                width = 12, height = 5, dpi = 300, bg = "white")
cat("[OK] Saved: outputs/figures/H2_bridge_ratio_concordance_windows.png\n")
cat("[OK] Saved: outputs/figures/H2_bridge_ratio_concordance_days.png\n")

cat("\n==============================================================================\n")
cat("Script completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================================\n")
while (sink.number() > 0) sink()
