# Run the recalibrated model through 1861-2010 and test whether it passes
# through the observed catch and biomass.
#
# WHAT THIS ANSWERS
# 06 refitted the fish to the COBALT resource at a single reference level. That
# only shows the model can sit in the right place when the forcing is frozen.
# The real question is whether it stays there once the forcing moves: does it
# track observed catch and biomass over the years where data exist, without the
# level anchor that 04 applies?
#
# The comparison is deliberately against the SAME observations 05 uses --
# biomass_ts.RDS and new_NAFO_4VW_Catches.RDS -- so the numbers here sit beside
# 05's and can be read against them directly.
#
# Two things are reported separately, because they answer different questions:
#   - all 19 species, which includes the forage catchability multipliers
#     (Herring x50, Sand Lance x3000, Capelin x200). These are about 85% of the
#     observed total biomass, so this comparison is dominated by them.
#   - the 16 species with no multiplier, which is the comparison against data
#     rather than against a calibration assumption.
#
# Run from the project root:  Rscript 07.Dynamics_and_fit_assessment.R

suppressMessages({
  library(mizer)
  library(therMizer)
  library(dplyr)
  library(tidyr)
})

# gaussian_sizeweight() lives here. The Seals gear refers to it by name, so
# upgradeTherParams() fails without it.
source("Ess_Fishmip_data/Useful_functions.R")

region_area_m2 <- 146360555008
first_data_year <- 1961
last_data_year  <- 2010

# The three species whose survey biomass carries a catchability multiplier. They
# are reported apart from the rest for the reason given above.
multiplied_forage_species <- c("Atlantic Herring", "Sand Lance", "Capelin")

calibrated_params_file <- "Ess_Fishmip_data/params_GFDL_steady.RDS"
output_file            <- "Ess_Fishmip_data/sims_GFDL_recalibrated.RDS"

# ---- Inputs ------------------------------------------------------------------

if (!file.exists(calibrated_params_file)) {
  stop("Missing ", calibrated_params_file,
       ". Run 06.Recalibrate_on_GFDL_resource.R first.", call. = FALSE)
}
params <- readRDS(calibrated_params_file)

# The same unanchored resource 06 calibrated against, but now year by year
# rather than held at the reference mean.
resource_log10_by_year <- as.matrix(
  read.table("Ess_Fishmip_data/GFDL_resource_spectra_unanchored.dat",
             header = TRUE, check.names = FALSE))

ocean_temp_table <- read.csv("Ess_Fishmip_data/Temp_by_depth_spin.csv")
rownames(ocean_temp_table) <- as.character(ocean_temp_table$year)
ocean_temp_by_year <- as.matrix(
  ocean_temp_table[, setdiff(names(ocean_temp_table), "year"), drop = FALSE])

vertical_migration_env <- new.env()
load("Ess_Fishmip_data/GFDL_ESS_vertical_migration_array.Rdata",
     envir = vertical_migration_env)
vertical_migration_array <- vertical_migration_env$vertical_migration_array
exposure_array <- as.matrix(
  read.csv("Ess_Fishmip_data/GFDL_ESS_exposure_array.csv", row.names = 1))

effort_table <- read.csv("Ess_Fishmip_data/Effort_spin_1861_2010.csv",
                         check.names = FALSE)
effort_by_year <- as.matrix(effort_table[, setdiff(names(effort_table), "year"),
                                         drop = FALSE])
rownames(effort_by_year) <- as.character(effort_table$year)

# The forcing and the effort must cover the same years in the same order, or the
# run silently pairs the wrong year's temperature with the wrong year's fishing.
stopifnot(
  nrow(resource_log10_by_year) == nrow(ocean_temp_by_year),
  identical(rownames(resource_log10_by_year), rownames(ocean_temp_by_year)),
  identical(rownames(resource_log10_by_year), rownames(effort_by_year)),
  ncol(resource_log10_by_year) == length(params@w_full)
)

# ---- Project through the full 1861-2010 series -------------------------------

resource_time_varying <- resource_log10_by_year
dimnames(resource_time_varying) <- list(time = rownames(resource_log10_by_year),
                                        w = names(params@initial_n_pp))

params_dynamic <- upgradeTherParams(
  params,
  ocean_temp               = ocean_temp_by_year,
  n_pp_array               = resource_time_varying,
  vertical_migration_array = vertical_migration_array,
  exposure_array           = exposure_array,
  aerobic_effect           = TRUE,
  metabolism_effect        = TRUE)

cat("projecting 1861-2010 on the unanchored COBALT resource...\n")
sim <- suppressWarnings(suppressMessages(
  project(params_dynamic, effort = effort_by_year)))
saveRDS(sim, output_file)

# ---- Modelled biomass and catch, per species per year ------------------------

modelled_biomass <- getBiomass(sim)      # years x species, grams
modelled_yield   <- getYield(sim)        # years x species, grams per year

to_long <- function(year_by_species_matrix, value_name) {
  as.data.frame(as.table(year_by_species_matrix)) |>
    setNames(c("Year", "Species", value_name)) |>
    mutate(Year = as.integer(as.character(Year)),
           Species = as.character(Species))
}
modelled_biomass_long <- to_long(modelled_biomass, "modelled")
modelled_yield_long   <- to_long(modelled_yield,   "modelled")

# ---- Observations ------------------------------------------------------------

observed_biomass_long <- readRDS("Ess_Fishmip_data/biomass_ts.RDS") |>
  transmute(Year = as.integer(Year), Species = as.character(Species),
            observed = Biomass)

observed_catch_long <- readRDS("Ess_Fishmip_data/new_NAFO_4VW_Catches.RDS") |>
  pivot_longer(-Year, names_to = "Species", values_to = "observed") |>
  filter(!is.na(observed), observed > 0) |>
  mutate(Year = as.integer(Year), Species = as.character(Species))

# ---- Fit summary -------------------------------------------------------------

# log10(observed / modelled): 0 is perfect, positive means the model sits under
# the observations. Same convention as 05's diagnostics, so the two are directly
# comparable.
summarise_fit <- function(modelled_long, observed_long, label, species_subset) {
  both <- inner_join(modelled_long, observed_long, by = c("Year", "Species")) |>
    filter(Year >= first_data_year, Year <= last_data_year,
           modelled > 0, observed > 0, Species %in% species_subset)
  if (!nrow(both)) { cat(sprintf("%-38s no overlapping data\n", label)); return(invisible(NULL)) }
  annual_totals <- both |> group_by(Year) |>
    summarise(modelled = sum(modelled), observed = sum(observed), .groups = "drop")
  log10_ratio <- median(log10(annual_totals$observed / annual_totals$modelled))
  cat(sprintf("%-38s %+6.2f log10  (x%.2f)  %d years\n",
              label, log10_ratio, 10^log10_ratio, nrow(annual_totals)))
  invisible(both)
}

all_species       <- unique(params@species_params$species)
unmultiplied_only <- setdiff(all_species, multiplied_forage_species)

cat("\n=== fit of the recalibrated model, 1961-2010 ===\n")
cat("log10(observed/modelled); 0 is perfect, positive means the model is under\n\n")
summarise_fit(modelled_biomass_long, observed_biomass_long,
              "biomass, all 19 species", all_species)
summarise_fit(modelled_biomass_long, observed_biomass_long,
              "biomass, 16 without multipliers", unmultiplied_only)
summarise_fit(modelled_yield_long, observed_catch_long,
              "catch, all species with data", all_species)
summarise_fit(modelled_yield_long, observed_catch_long,
              "catch, 16 without multipliers", unmultiplied_only)

# ---- Does it drift? ----------------------------------------------------------

# The failure mode the level anchor was fixing was a model that starts too high
# and climbs through the spin-up. Report the trajectory so that is visible
# rather than hidden inside a median.
biomass_g_per_m2_by_year <- rowSums(modelled_biomass) / region_area_m2
catch_g_per_m2_by_year   <- rowSums(modelled_yield) / region_area_m2
years <- as.integer(rownames(modelled_biomass))
at_year <- function(series, y) series[which(years == y)]

cat("\n=== trajectory (total, g/m2) ===\n")
for (y in c(1861, 1900, 1950, 1961, 1980, 2000, 2010)) {
  cat(sprintf("  %d   biomass %7.3f   catch %7.4f\n",
              y, at_year(biomass_g_per_m2_by_year, y),
              at_year(catch_g_per_m2_by_year, y)))
}
cat(sprintf("\nspin-up drift 1861 -> 1960: x%.2f\n",
            at_year(biomass_g_per_m2_by_year, 1960) /
              at_year(biomass_g_per_m2_by_year, 1861)))

# ---- Per-species biomass fit over the data years ----------------------------

per_species <- inner_join(modelled_biomass_long, observed_biomass_long,
                          by = c("Year", "Species")) |>
  filter(Year >= first_data_year, Year <= last_data_year,
         modelled > 0, observed > 0) |>
  group_by(Species) |>
  summarise(log10_ratio = median(log10(observed / modelled)), .groups = "drop") |>
  arrange(desc(abs(log10_ratio)))

cat("\n=== per-species biomass, worst 8 (log10 obs/mod) ===\n")
for (i in seq_len(min(8, nrow(per_species)))) {
  cat(sprintf("  %-22s %+6.2f  (x%.2f)%s\n",
              per_species$Species[i], per_species$log10_ratio[i],
              10^per_species$log10_ratio[i],
              ifelse(per_species$Species[i] %in% multiplied_forage_species,
                     "   <- multiplied", "")))
}
cat(sprintf("\nwithin 2x of observed: %d/%d species\n",
            sum(abs(per_species$log10_ratio) < log10(2)), nrow(per_species)))
# ---- Skill metrics the FishMIP regional workflow asks for --------------------

# Ortega-Cisneros et al. (2025, Earth's Future), step 6, asks for at least the
# bias and the correlation of the catch time series and, where observations
# exist, the biomass time series, for key species. It also says the catch years
# 2005-2010 must not be used in calibration, because they are held back to test
# predictive skill. The calibration here uses 1978-1981 only, so 2005-2010 is a
# genuine out-of-sample test, and it is reported separately for that reason.
skill_table <- function(modelled_long, observed_long, years) {
  inner_join(modelled_long, observed_long, by = c("Year", "Species")) |>
    filter(Year %in% years, modelled > 0, observed > 0) |>
    group_by(Species) |>
    filter(n() >= 4) |>
    summarise(
      # positive bias means the model is ABOVE the observations
      bias_log10  = median(log10(modelled / observed)),
      correlation = suppressWarnings(cor(modelled, observed)),
      years_used  = n(), .groups = "drop")
}

report_skill <- function(label, modelled_long, observed_long, years) {
  skill <- skill_table(modelled_long, observed_long, years) |>
    filter(!Species %in% multiplied_forage_species)
  if (!nrow(skill)) { cat(sprintf("%-34s no data\n", label)); return(invisible()) }
  cat(sprintf(paste0("%-34s %2d species | median bias %+5.2f log10 (x%.2f) | ",
                     "within 2x: %2d | median correlation %+5.2f | r > 0.5: %2d\n"),
              label, nrow(skill), median(skill$bias_log10),
              10^median(skill$bias_log10),
              sum(abs(skill$bias_log10) < log10(2)),
              median(skill$correlation, na.rm = TRUE),
              sum(skill$correlation > 0.5, na.rm = TRUE)))
  invisible(skill)
}

cat("\n=== bias and correlation by species (species without multipliers) ===\n")
cat("bias is log10(modelled/observed): positive means the model is too high\n\n")
report_skill("catch,   1970-2004 (in sample era)", modelled_yield_long,
             observed_catch_long, 1970:2004)
report_skill("catch,   2005-2010 (held out)", modelled_yield_long,
             observed_catch_long, 2005:2010)
report_skill("biomass, 1970-2004 (in sample era)", modelled_biomass_long,
             observed_biomass_long, 1970:2004)
report_skill("biomass, 2005-2010 (held out)", modelled_biomass_long,
             observed_biomass_long, 2005:2010)

cat(sprintf("\nwrote %s\n", output_file))
