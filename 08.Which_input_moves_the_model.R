# Which input moves the recalibrated model away from its calibrated state?
#
# 06 fits the fish to the COBALT resource with resource, temperature and effort
# held at their 1978-1981 averages. 07 then runs the real 1861-2010 series. A
# dynamic run is not expected to sit at the calibrated state, but it is worth
# knowing WHICH input moves it, because that says where to look when the fit to
# data is poor.
#
# This script starts every run from 06's calibrated state and releases one
# input at a time, holding the other two at their 1978-1981 mean:
#
#   A  everything constant        -> does the calibrated state hold at all?
#   B  only effort varies
#   C  only temperature varies
#   D  only the resource varies
#
# RESULT WHEN LAST RUN (total biomass, g/m2)
#                                   1861    1870    1961    2010
#   A  all constant                 6.45    6.44    6.42    6.42
#   B  only effort varies           6.45    8.16    9.58    1.74
#   C  only temperature varies      6.45    6.45    6.71    6.44
#   D  only resource varies         6.45    6.43    6.40    6.41
#
# (These are with the estimator-based adult mortality in 06. Before it was
# added, case B ended at 6.04 in 2010 rather than 1.74: with realistic adult
# mortality the community no longer rides out the rise in seal effort.)
#
# So the calibrated state is a genuine steady state of the model (A), the
# year-to-year COBALT resource and temperature barely move it (C, D), and
# nearly all of the movement in the dynamic run comes from effort (B).
#
# WHY THE CALIBRATION WINDOW MATTERS SO MUCH
# The same test on a calibration that used 1970-2010 averages gave 125 g/m2 at
# 1961 in case B instead of 9.95. The Seals "gear" averages 0.278 over
# 1970-2010 but only 0.023 in 1978-1981 and 0.004 over the spin-up. Balancing
# the fish against the 1970-2010 mean builds in seal mortality that did not
# exist in the early years, so the community balloons when it is taken away.
# 1978-1981 is close to the spin-up conditions, so it does not.
#
# Takes roughly ten minutes: four 150-year projections.
# Run from the project root:  Rscript 08.Which_input_moves_the_model.R

suppressMessages({
  library(mizer)
  library(therMizer)
})

# The Seals gear refers to gaussian_sizeweight() by name.
source("Ess_Fishmip_data/Useful_functions.R")

region_area_m2     <- 146360555008
averaging_window   <- as.character(1978:1981)   # the window 06 calibrated on

params <- readRDS("Ess_Fishmip_data/params_GFDL_steady.RDS")

resource_by_year <- as.matrix(
  read.table("Ess_Fishmip_data/GFDL_resource_spectra_unanchored.dat",
             header = TRUE, check.names = FALSE))

temperature_table <- read.csv("Ess_Fishmip_data/Temp_by_depth_spin.csv")
rownames(temperature_table) <- as.character(temperature_table$year)
temperature_by_year <- as.matrix(
  temperature_table[, setdiff(names(temperature_table), "year"), drop = FALSE])

effort_table <- read.csv("Ess_Fishmip_data/Effort_spin_1861_2010.csv",
                         check.names = FALSE)
effort_by_year <- as.matrix(effort_table[, setdiff(names(effort_table), "year")])
rownames(effort_by_year) <- as.character(effort_table$year)
effort_by_year <- effort_by_year[, names(params@initial_effort)]

vertical_migration_env <- new.env()
load("Ess_Fishmip_data/GFDL_ESS_vertical_migration_array.Rdata",
     envir = vertical_migration_env)
exposure_array <- as.matrix(
  read.csv("Ess_Fishmip_data/GFDL_ESS_exposure_array.csv", row.names = 1))

all_years <- rownames(resource_by_year)
stopifnot(identical(all_years, rownames(temperature_by_year)),
          identical(all_years, rownames(effort_by_year)))

# Replace every year of a year-by-something table with its mean over the
# averaging window, keeping the same shape and labels.
held_at_time_average <- function(year_by_column) {
  window_mean <- colMeans(year_by_column[averaging_window, , drop = FALSE])
  constant <- matrix(rep(window_mean, each = nrow(year_by_column)),
                     nrow = nrow(year_by_column))
  dimnames(constant) <- dimnames(year_by_column)
  constant
}

run_one_case <- function(label, resource, temperature, effort) {
  dimnames(resource) <- list(time = all_years, w = names(params@initial_n_pp))
  params_for_run <- suppressWarnings(suppressMessages(upgradeTherParams(
    params,
    ocean_temp               = temperature,
    n_pp_array               = resource,
    vertical_migration_array = vertical_migration_env$vertical_migration_array,
    exposure_array           = exposure_array,
    aerobic_effect           = TRUE,
    metabolism_effect        = TRUE)))
  sim <- suppressWarnings(suppressMessages(
    project(params_for_run, effort = effort)))
  biomass <- rowSums(getBiomass(sim)) / region_area_m2
  cat(sprintf("%-32s %7.2f %7.2f %7.2f %7.2f\n", label,
              biomass["1861"], biomass["1870"], biomass["1961"], biomass["2010"]))
}

cat("total biomass, g/m2\n")
cat(sprintf("%-32s %7s %7s %7s %7s\n", "", "1861", "1870", "1961", "2010"))
run_one_case("A  all constant",
             held_at_time_average(resource_by_year),
             held_at_time_average(temperature_by_year),
             held_at_time_average(effort_by_year))
run_one_case("B  only effort varies",
             held_at_time_average(resource_by_year),
             held_at_time_average(temperature_by_year),
             effort_by_year)
run_one_case("C  only temperature varies",
             held_at_time_average(resource_by_year),
             temperature_by_year,
             held_at_time_average(effort_by_year))
run_one_case("D  only resource varies",
             resource_by_year,
             held_at_time_average(temperature_by_year),
             held_at_time_average(effort_by_year))

cat("\neffort by gear: the calibration window against the spin-up\n")
print(round(rbind(
  calibration_window    = colMeans(effort_by_year[averaging_window, ]),
  spinup_1861_1960      = colMeans(effort_by_year[as.character(1861:1960), ])), 3))
