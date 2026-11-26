# Why recalibrating against the COBALT resource with steady() does not work.
#
# THE RESULT THIS DEMONSTRATES
# 06 recalibrates the fish against the raw COBALT resource and reports an
# excellent fit: 16/19 species within 2x of observed biomass, median ratio
# 1.019, total 6.872 g/m2 against an observed ~6.9. 07 then projects that same
# object and it drifts up about 4x.
#
# Those two results are not in conflict, and the drift is not a bug in 07. The
# state steady() reports is simply not a fixed point of the dynamics that
# project() runs. This script shows that directly: it takes 06's output and
# projects it under forcing that does not move at all -- the same reference
# resource and temperature steady() was given, and initial_effort -- so there is
# nothing left that could push it. It still leaves 6.872 g/m2 and settles at
# about 21.7 g/m2 within ten years.
#
# THE CAUSE
# therMizer drives the resource through plankton_forcing(), and mizer has no
# balance_plankton_forcing() to go with it. steady() says so on every call:
#
#   "There is no `balance_plankton_forcing()` function, so the resource capacity
#    could not be rebalanced and the preserved resource abundance need not be a
#    steady state of `plankton_forcing()`."
#
# That warning is the whole finding. steady() holds the resource abundance
# fixed and solves the fish around it, but under plankton_forcing() the resource
# is not free to adjust, so the fish equilibrium it finds is not where the
# forced system actually goes.
#
# WHAT IT MEANS FOR THE PROJECT
# Calibrating this model at the raw COBALT resource via steady() cannot work,
# however the reference window or the matching steps are chosen, because the
# target state is one the model does not stay in. Two routes remain:
#   1. Keep the level anchor in 04. It scales the forcing down to the resource
#      the fish were fitted to, and as of this writing it is what actually fits
#      the observations -- modelled catch within 15%.
#   2. Calibrate on the projected dynamics rather than on steady(): project the
#      full series, compare with the observed years, adjust, repeat. More
#      expensive, and not attempted here.
#
# It is NOT evidence that the COBALT resource is wrong. It says this model, at
# that resource level, has an attractor about 3x the observed biomass.
#
# Run from the project root:  Rscript 08.Why_steady_is_not_the_attractor.R

suppressMessages({
  library(mizer)
  library(therMizer)
})

# The Seals gear refers to gaussian_sizeweight() by name.
source("Ess_Fishmip_data/Useful_functions.R")

region_area_m2  <- 146360555008
reference_years <- 1978:1981      # the window 06 calibrated on
projection_years <- 300           # long enough to show where it settles

calibrated_params_file <- "Ess_Fishmip_data/params_GFDL_steady.RDS"
if (!file.exists(calibrated_params_file)) {
  stop("Missing ", calibrated_params_file,
       ". Run 06.Recalibrate_on_GFDL_resource.R first.", call. = FALSE)
}
params <- readRDS(calibrated_params_file)

biomass_g_per_m2 <- function(x) sum(getBiomass(x)) / region_area_m2
cat(sprintf("state 06 saved: %.3f g/m2\n", biomass_g_per_m2(params)))

# mizer's own residual, as one number. It is near zero, which is exactly the
# trap: by its own measure this IS a steady state. The residual is computed
# against a resource that is held fixed, so it cannot see the problem.
cat(sprintf("max |steady residual| reported by mizer: %.3g\n",
            max(abs(getSteadyResidual(params)))))

# ---- Forcing that cannot move -----------------------------------------------

resource_log10_by_year <- as.matrix(
  read.table("Ess_Fishmip_data/GFDL_resource_spectra_unanchored.dat",
             header = TRUE, check.names = FALSE))
ocean_temp_table <- read.csv("Ess_Fishmip_data/Temp_by_depth_spin.csv")
rownames(ocean_temp_table) <- as.character(ocean_temp_table$year)
ocean_temp_by_year <- as.matrix(
  ocean_temp_table[, setdiff(names(ocean_temp_table), "year"), drop = FALSE])

reference_years_chr <- as.character(reference_years)
resource_at_reference <- colMeans(
  resource_log10_by_year[reference_years_chr, , drop = FALSE])
temperature_at_reference <- colMeans(
  ocean_temp_by_year[reference_years_chr, , drop = FALSE])

time_labels <- as.character(seq_len(projection_years))
resource_held_constant <- matrix(
  rep(resource_at_reference, each = projection_years), nrow = projection_years,
  dimnames = list(time = time_labels, w = names(params@initial_n_pp)))
temperature_held_constant <- matrix(
  rep(temperature_at_reference, each = projection_years), nrow = projection_years,
  dimnames = list(time = time_labels, realm = colnames(ocean_temp_by_year)))

vertical_migration_env <- new.env()
load("Ess_Fishmip_data/GFDL_ESS_vertical_migration_array.Rdata",
     envir = vertical_migration_env)
exposure_array <- as.matrix(
  read.csv("Ess_Fishmip_data/GFDL_ESS_exposure_array.csv", row.names = 1))

params_constant_forcing <- suppressWarnings(suppressMessages(upgradeTherParams(
  params,
  ocean_temp               = temperature_held_constant,
  n_pp_array               = resource_held_constant,
  vertical_migration_array = vertical_migration_env$vertical_migration_array,
  exposure_array           = exposure_array,
  aerobic_effect           = TRUE,
  metabolism_effect        = TRUE)))

# Effort held at the value the calibration itself assumed, so fishing cannot be
# the explanation either.
effort_held_constant <- matrix(
  rep(params@initial_effort, each = projection_years), nrow = projection_years,
  dimnames = list(time = time_labels, gear = names(params@initial_effort)))

# ---- Project ----------------------------------------------------------------

sim <- suppressWarnings(suppressMessages(
  project(params_constant_forcing, effort = effort_held_constant)))
biomass_by_year <- rowSums(getBiomass(sim)) / region_area_m2

cat("\nprojected under constant reference forcing and constant effort:\n")
for (year in c(1, 5, 10, 25, 50, 100, 200, projection_years)) {
  cat(sprintf("  year %3d   %7.3f g/m2\n", year, biomass_by_year[year]))
}
cat(sprintf("\ndrift from the state steady() reported: x%.2f\n",
            biomass_by_year[projection_years] / biomass_by_year[1]))
cat(paste0("\nNothing in this run varies, so the drift is not the time-varying\n",
           "forcing, not the spin-up effort series, and not 07. The state\n",
           "steady() reported was never a fixed point of these dynamics.\n"))
