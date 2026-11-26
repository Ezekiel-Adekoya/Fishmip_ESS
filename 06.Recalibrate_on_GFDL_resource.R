# Recalibrate the fish community against the GFDL/COBALT resource.
#
# WHY THIS EXISTS
# The Stage-1 model in params_stage1_steady.RDS was calibrated against a tuned
# semichemostat resource, then script 05 fed it the COBALT-derived resource,
# which is about 5x larger over the size range the model actually sees. The fish
# were therefore eating a resource they were never fitted to, which inflated
# biomass about 5.7x over the spin-up and modelled catch about 3x.
#
# 04 works around that by scaling the forcing DOWN to the calibrated resource
# (the "level anchor"). This script takes the other route: leave the resource as
# COBALT reports it and fit the fish to THAT. If this works, the anchor in 04 is
# no longer needed.
#
# WHAT IS HELD FIXED, AND WHY
# This starts from params_stage1_steady.RDS rather than rebuilding from
# ESS_species_params.csv, so the species list, gear setup and interaction matrix
# are unchanged and the ONLY thing that differs is the resource the fish are
# calibrated against. A rebuild from scratch would also have to redo the
# interactive tuneParams() step in ESS_Restructured/New_steady.Rmd, which cannot
# be reproduced from a script, and would confound several changes at once.
#
# Calibration needs forcing that does not move, so the resource and temperature
# are held at their reference-period means here. The time-varying run is 07.
#
# Run from the project root:  Rscript 06.Recalibrate_on_GFDL_resource.R

suppressMessages({
  library(mizer)
  library(mizerExperimental)
  library(therMizer)
})

# gaussian_sizeweight() lives here. It is the seal prey size-selectivity kernel
# and mizer calls it by name from the Seals gear, so matchYield() fails with
# "object 'gaussian_sizeweight' of mode 'function' was not found" without it.
source("Ess_Fishmip_data/Useful_functions.R")

region_area_m2 <- 146360555008      # Eastern Scotian Shelf

# steady() and matchBiomasses() need inputs that are constant in time, and those
# constants must be the TIME AVERAGE over the same period the calibration targets
# are averaged over. The observed biomass in species_params is a 1970-2010 mean
# (it matches the 1970-2010 survey average at a ratio of 0.994), so the forcing
# and the effort are averaged over 1970-2010 as well.
#
# This used to be 1978:1981, borrowed from the window 04 uses for its anchor.
# That fitted time-averaged targets against a four-year forcing window, which is
# inconsistent: 1978-81 is 0.1 C warmer than the 1970-2010 mean, and a long way
# from the 1961-1970 conditions the spin-up repeats (0.5 C cooler, 19% more
# effort).
reference_years <- 1970:2010

# steady() needs longer than its default here: a 50-year test run was still
# moving. The convergence tolerance is left at mizer's default (tol = 0.1 * dt),
# because the default distance function keeps measuring motion in the resource
# that the fish biomasses do not actually see -- therMizer has no
# balance_plankton_forcing(), so the resource is never rebalanced and that
# residual never goes to zero. A non-converged return is therefore expected;
# what matters is whether the biomasses have stopped moving, which the
# diagnostics below report.
steady_years <- 300

# preserve = "erepro" keeps each species' reproductive efficiency where Stage 1
# had it and lets R_max absorb the change instead. mizer's default preserves
# reproduction_level, which on this model pushed Spiny Dogfish to erepro ~39 --
# meaning it would need 39 times every egg it produces, which is not physically
# admissible. Preserving erepro keeps the fit inside the admissible range.
steady_preserve <- "erepro"

# How many times to repeat matchBiomasses() -> steady(). Each pass moves R_max
# toward the observed biomass and then re-solves; a few passes is normal.
match_biomass_passes <- 4

output_file <- "Ess_Fishmip_data/params_GFDL_steady.RDS"

# ---- Inputs ------------------------------------------------------------------

params <- validParams(readRDS("Ess_Fishmip_data/params_stage1_steady.RDS"))

# Thermal tolerances. These rows are assumed to line up positionally with
# species_params, so assert it: a silent reorder would give each species someone
# else's temperature limits.
thermal_data <- read.csv("Ess_Fishmip_data/species_params_with_temp.csv")
stopifnot(identical(as.character(thermal_data$Common_Name),
                    as.character(params@species_params$species)))
params@species_params$temp_min <- thermal_data$TempPrefMin
params@species_params$temp_max <- thermal_data$TempPrefMax

# The resource as COBALT reports it, with NO level anchor applied. This is the
# whole point of the script, so read the unanchored file and fail loudly if 04
# has not been rerun to produce it.
unanchored_file <- "Ess_Fishmip_data/GFDL_resource_spectra_unanchored.dat"
if (!file.exists(unanchored_file)) {
  stop("Missing ", unanchored_file, ".\n",
       "Rerun 04.Resource_spectra_prep.Rmd, which now writes it alongside the ",
       "anchored spectra.", call. = FALSE)
}
resource_log10_by_year <- as.matrix(read.table(unanchored_file, header = TRUE,
                                               check.names = FALSE))

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

reference_years_chr <- as.character(reference_years)
stopifnot(
  ncol(resource_log10_by_year) == length(params@w_full),
  all(reference_years_chr %in% rownames(resource_log10_by_year)),
  all(reference_years_chr %in% rownames(ocean_temp_by_year))
)

# ---- Constant reference forcing ---------------------------------------------

# Mean resource and temperature over the reference window. These are the values
# the fish get fitted to.
resource_at_reference <- colMeans(
  resource_log10_by_year[reference_years_chr, , drop = FALSE])
temperature_at_reference <- colMeans(
  ocean_temp_by_year[reference_years_chr, , drop = FALSE])

# therMizer looks the forcing up by year, and steady() projects forward
# internally, so the arrays are made long enough that it cannot run off the end.
# Every row is identical, so the forcing is constant in time either way.
calibration_years <- 2000
calibration_time_labels <- as.character(seq_len(calibration_years))

resource_held_constant <- matrix(
  rep(resource_at_reference, each = calibration_years),
  nrow = calibration_years,
  dimnames = list(time = calibration_time_labels,
                  w = names(params@initial_n_pp)))

temperature_held_constant <- matrix(
  rep(temperature_at_reference, each = calibration_years),
  nrow = calibration_years,
  dimnames = list(time = calibration_time_labels,
                  realm = colnames(ocean_temp_by_year)))

# Effort is held at its time average over the same window, gear by gear, for the
# same reason. Stage 1's initial_effort turns out to be exactly the 1978-1981
# mean, so leaving it alone would quietly reintroduce the four-year window.
effort_table <- read.csv("Ess_Fishmip_data/Effort_spin_1861_2010.csv",
                         check.names = FALSE)
effort_in_window <- effort_table[effort_table$year %in% reference_years,
                                 setdiff(names(effort_table), "year"),
                                 drop = FALSE]
time_averaged_effort <- colMeans(effort_in_window)
stopifnot(setequal(names(time_averaged_effort), names(params@initial_effort)))
cat(sprintf("effort: Stage-1 initial_effort %.4f  ->  %d-%d mean %.4f\n",
            sum(params@initial_effort), min(reference_years),
            max(reference_years), sum(time_averaged_effort)))
initial_effort(params) <- time_averaged_effort[names(params@initial_effort)]

params_thermizer <- upgradeTherParams(
  params,
  ocean_temp               = temperature_held_constant,
  n_pp_array               = resource_held_constant,
  vertical_migration_array = vertical_migration_array,
  exposure_array           = exposure_array,
  aerobic_effect           = TRUE,
  metabolism_effect        = TRUE)

# Put the model's own resource state AT the time-averaged COBALT level and hold
# it there while calibrating.
#
# This step matters more than anything else in the script. Attaching n_pp_array
# only tells therMizer what to feed the model during project(); it does not
# change the resource the params object currently holds, and steady() solves
# the fish around whatever resource it is handed. Without this, steady() fitted
# the fish to Stage 1's old resource (about a fifth of COBALT), and the moment
# project() swapped the COBALT resource in the fish had four times the food and
# biomass tripled within ten years. That looked like the model refusing to stay
# at its calibrated state; it was a calibration done on the wrong food level.
time_averaged_resource_density <- 10^resource_at_reference / params@dw_full
time_averaged_resource_density[!is.finite(time_averaged_resource_density)] <- 0
# therMizer's plankton_forcing() zeroes the resource at and above the cutoff, so
# the calibration resource has to do the same or the two will not agree.
time_averaged_resource_density[
  params@w_full >= params@resource_params$w_pp_cutoff] <- 0
names(time_averaged_resource_density) <- names(params@initial_n_pp)

initialNResource(params_thermizer) <- time_averaged_resource_density
# Constant in time, as steady() and matchBiomasses() require. 07 switches back
# to the year-by-year forcing for the dynamic run.
resource_dynamics(params_thermizer) <- "resource_constant"

# How much bigger the COBALT resource is than the one Stage 1 was fitted to,
# over the bins the model sees. This is the gap the recalibration has to absorb.
bins_the_model_sees <- params@w_full < params@resource_params$w_pp_cutoff
resource_biomass_g <- function(log10_numbers_per_bin) {
  numbers_per_gram <- 10^log10_numbers_per_bin / params@dw_full
  numbers_per_gram[!is.finite(numbers_per_gram)] <- 0
  sum((numbers_per_gram * params@w_full * params@dw_full)[bins_the_model_sees])
}
cat(sprintf(
  "COBALT resource vs the Stage-1 calibrated resource (w < %g g): x%.2f\n",
  params@resource_params$w_pp_cutoff,
  resource_biomass_g(resource_at_reference) /
    sum((params@initial_n_pp * params@w_full * params@dw_full)[bins_the_model_sees])))
cat(sprintf("reference window %d-%d, temperature %.2f / %.2f C\n\n",
            min(reference_years), max(reference_years),
            temperature_at_reference[1], temperature_at_reference[2]))

# ---- Calibrate ---------------------------------------------------------------

# Quietly, because therMizer warns on every call that it cannot rebalance the
# resource (it has no balance_plankton_forcing()), which is expected here and
# would otherwise bury the diagnostics.
solve_steady <- function(p) {
  suppressWarnings(suppressMessages(
    steady(p, t_max = steady_years, progress_bar = FALSE,
           preserve = steady_preserve)))
}

report <- function(label, p) {
  biomass_g_per_m2 <- sum(getBiomass(p)) / region_area_m2
  worst_erepro <- max(p@species_params$erepro)
  n_over_one <- sum(p@species_params$erepro > 1)
  cat(sprintf("%-34s biomass %7.3f g/m2   max erepro %8.3g   erepro>1: %d\n",
              label, biomass_g_per_m2, worst_erepro, n_over_one))
  invisible(p)
}

cat("--- calibration ---\n")
report("start (Stage-1 fish, COBALT food)", params_thermizer)

# calibrateBiomass() is deliberately NOT used. It rescales the whole system,
# resource included, to match total observed biomass -- which would move the
# resource off the COBALT level this script exists to keep. Only the fish are
# allowed to adjust here.
params_calibrated <- report("after steady", solve_steady(params_thermizer))

# matchBiomasses moves each species' R_max toward its own observed biomass;
# steady then re-solves the community around the new values.
for (pass in seq_len(match_biomass_passes)) {
  params_calibrated <- matchBiomasses(params_calibrated)
  params_calibrated <- solve_steady(params_calibrated)
  report(sprintf("after matchBiomasses pass %d", pass), params_calibrated)
}

# Yield depends on both biomass and fishing mortality, so it is matched last.
#
# This step is attempted but NOT trusted. On this model matchYield() diverges:
# it drove total biomass to about 3e12 g/m2 and erepro to 1e15 across seven
# species, which is a broken fit rather than a worse one. So the result is kept
# only if it is actually an improvement and stays physically admissible, and the
# biomass-matched fit is kept otherwise. Yield comes out close without being
# fitted anyway, because matching biomass under the observed effort largely
# determines it.
fit_quality <- function(p) {
  sp <- p@species_params
  modelled <- as.numeric(getBiomass(p))
  has_target <- !is.na(sp$biomass_observed) & sp$biomass_observed > 0
  if (!any(has_target)) return(list(score = Inf, admissible = FALSE))
  log_ratio <- log10(modelled[has_target] / sp$biomass_observed[has_target])
  list(score      = median(abs(log_ratio[is.finite(log_ratio)])),
       admissible = all(is.finite(modelled)) && max(sp$erepro) <= 1)
}

quality_before <- fit_quality(params_calibrated)
params_with_yield <- try(
  suppressWarnings(suppressMessages(matchYield(params_calibrated,
                                               keep = "biomass"))),
  silent = TRUE)
if (inherits(params_with_yield, "try-error")) {
  cat("matchYield errored, keeping the biomass-matched fit:\n  ",
      conditionMessage(attr(params_with_yield, "condition")), "\n")
} else {
  params_with_yield <- solve_steady(matchBiomasses(params_with_yield))
  report("after matchYield + matchBiomasses", params_with_yield)
  quality_after <- fit_quality(params_with_yield)
  if (quality_after$score < quality_before$score &&
      (quality_after$admissible || !quality_before$admissible)) {
    cat("  matchYield improved the fit, keeping it\n")
    params_calibrated <- params_with_yield
  } else {
    cat(sprintf(paste0("  REJECTED matchYield: median |log10| biomass error ",
                       "%.3f -> %.3f, admissible %s -> %s\n",
                       "  keeping the biomass-matched fit\n"),
                quality_before$score, quality_after$score,
                quality_before$admissible, quality_after$admissible))
  }
}

params_final <- params_calibrated

# ---- Checks ------------------------------------------------------------------

# erepro > 1 means a species needs more than every egg produced to survive, so
# the fit is not physically admissible for it. Report rather than stop, because
# the Stage-1 object has the same issue for Spiny Dogfish and the comparison is
# still informative.
species_params_final <- params_final@species_params
over_one <- species_params_final$species[species_params_final$erepro > 1]
cat("\n--- checks ---\n")
if (length(over_one)) {
  cat("erepro > 1 (not physically admissible):",
      paste(over_one, collapse = ", "), "\n")
} else {
  cat("all erepro < 1\n")
}

observed_biomass <- species_params_final$biomass_observed
modelled_biomass <- as.numeric(getBiomass(params_final))
has_target <- !is.na(observed_biomass) & observed_biomass > 0
biomass_ratio <- modelled_biomass[has_target] / observed_biomass[has_target]
cat(sprintf("biomass: %d/%d species within 2x of observed, median ratio %.3f\n",
            sum(biomass_ratio > 0.5 & biomass_ratio < 2), sum(has_target),
            median(biomass_ratio)))

observed_yield <- species_params_final$yield_observed
modelled_yield <- as.numeric(getYield(params_final))
has_yield_target <- !is.na(observed_yield) & observed_yield > 0
yield_ratio <- modelled_yield[has_yield_target] / observed_yield[has_yield_target]
cat(sprintf("yield:   %d/%d species within 2x of observed, median ratio %.3f\n",
            sum(yield_ratio > 0.5 & yield_ratio < 2), sum(has_yield_target),
            median(yield_ratio)))

worst <- order(abs(log10(biomass_ratio)), decreasing = TRUE)[1:5]
cat("furthest from observed biomass:\n")
for (i in worst) {
  cat(sprintf("  %-22s x%.2f\n",
              species_params_final$species[has_target][i], biomass_ratio[i]))
}

# The resource must still be exactly the time-averaged COBALT spectrum. If any
# calibration step rescaled it, the fit is against the wrong food level again.
resource_shift <- max(abs(params_final@initial_n_pp -
                          time_averaged_resource_density)) /
  max(time_averaged_resource_density)
cat(sprintf("resource still at the COBALT level: max relative shift %.2e\n",
            resource_shift))
stopifnot(resource_shift < 1e-8)

saveRDS(params_final, output_file)
cat(sprintf("\nwrote %s\n", output_file))
cat(sprintf("  total biomass %.3f g/m2 at the COBALT resource, no level anchor\n",
            sum(getBiomass(params_final)) / region_area_m2))
