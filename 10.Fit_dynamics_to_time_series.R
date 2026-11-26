# Fit a small number of dynamic parameters to the catch and biomass time series.
#
# WHAT THIS IS, AND WHAT IT IS NOT
# 06 calibrates the model at a steady state (1978-1981 averages). Nothing in 06
# is tuned against how catch and biomass CHANGE through time. This script does
# that next step: it estimates five parameters that only act after the
# calibration period, by minimising the mismatch with the observed series.
#
# The calibrated steady state is not touched. Every parameter below is defined
# so that it has no effect in 1978-1981, so 06's fit stays exactly as it was and
# the two steps cannot fight each other.
#
# THE FIVE PARAMETERS
#   1. seal_rise_multiplier   how steeply seal mortality increases as seal
#                             effort rises above its 1978-1981 level
#   2. seal_saturation_effort the size of rise at which that increase is half
#                             saturated (seal mortality levelling off)
#   3. cod_extra_mortality    extra natural mortality on mature cod from 1990,
#                             per year. The seal gear takes no mature cod, and
#                             assessments have adult cod mortality rising then.
#   4. cod_moratorium_factor      what fraction of its catchability the bottom
#   5. haddock_moratorium_factor  trawl keeps on cod and on haddock from 1993.
#                                 The effort series is per gear, but the 1993
#                                 closure was per species, so trawl effort alone
#                                 cannot represent it.
#
# Each was chosen because an earlier hand-set test showed the model was wrong in
# exactly that respect (handover sections 11, 14 and 17). They are estimated
# here instead of picked.
#
# WHAT IS FITTED TO WHAT
# Survey biomass and fishery catch by species, 1978-2004, for the 16 species
# whose survey biomass carries no catchability multiplier. 2005-2010 is NOT
# used: the FishMIP regional protocol (Ortega-Cisneros et al. 2025, step 6)
# holds those years back, and they are reported at the end as the test of
# whether the estimates predict. Catch means the four fishing gears only; seals
# are a gear in this model and what they eat is not catch.
#
# The error is the squared difference of log10 values, averaged within each
# species and each variable first, so a species with many survey years does not
# outweigh one with few.
#
# Every run starts from the model state in 1978 taken from 07's projection. The
# parameters cannot change anything before then, so nothing is lost and one
# evaluation takes about a second.
#
# Run from the project root, after 06 and 07:
#   Rscript 10.Fit_dynamics_to_time_series.R

suppressMessages({
  library(mizer)
  library(therMizer)
  library(dplyr)
  library(tidyr)
  library(parallel)
  library(optimParallel)
})
source("Ess_Fishmip_data/Useful_functions.R")   # gaussian_sizeweight() for the Seals gear

region_area_m2 <- 146360555008
multiplied_forage_species <- c("Atlantic Herring", "Sand Lance", "Capelin")

first_fitted_year  <- 1978
last_fitted_year   <- 2004      # 2005-2010 are held out by the protocol
last_year          <- 2010
cod_mortality_starts <- 1990
moratorium_starts    <- 1993

fitted_values_file <- "Ess_Fishmip_data/fitted_dynamic_parameters.csv"
fitted_sim_file    <- "Ess_Fishmip_data/sims_GFDL_fitted.RDS"

quietly <- function(expression) suppressWarnings(suppressMessages(expression))

# ---- The model and its state in 1978 -----------------------------------------

baseline_sim <- readRDS("Ess_Fishmip_data/sims_GFDL_recalibrated.RDS")
params <- baseline_sim@params
species_names <- as.character(params@species_params$species)
fitted_species <- setdiff(species_names, multiplied_forage_species)

state_index_1978 <- which(dimnames(baseline_sim@n)[[1]] == as.character(first_fitted_year))
abundance_1978 <- baseline_sim@n[state_index_1978, , ]
resource_1978  <- baseline_sim@n_pp[state_index_1978, ]

effort_table <- read.csv("Ess_Fishmip_data/Effort_spin_1861_2010.csv",
                         check.names = FALSE)
effort_by_year <- as.matrix(effort_table[, setdiff(names(effort_table), "year")])
rownames(effort_by_year) <- as.character(effort_table$year)
effort_by_year <- effort_by_year[as.character(first_fitted_year:last_year),
                                 names(params@initial_effort)]

# Seal effort at the calibration state. The seal parameters only act on effort
# ABOVE this, which is what keeps 06's steady state untouched.
seal_effort_at_calibration <- mean(
  effort_table$Seals[effort_table$year %in% 1978:1981])

# ---- Observations ------------------------------------------------------------

observed_biomass <- readRDS("Ess_Fishmip_data/biomass_ts.RDS") |>
  transmute(Year = as.integer(Year), Species = as.character(Species),
            observed = Biomass) |>
  filter(observed > 0, Species %in% fitted_species,
         Year >= first_fitted_year, Year <= last_year)

observed_catch <- readRDS("Ess_Fishmip_data/new_NAFO_4VW_Catches.RDS") |>
  pivot_longer(-Year, names_to = "Species", values_to = "observed") |>
  filter(!is.na(observed), observed > 0, Species %in% fitted_species,
         Year >= first_fitted_year, Year <= last_year)

# ---- One model run for a given set of parameter values -----------------------

# Seal effort with a saturating rise. Up to the calibration level it is the
# original series. Above it, the rise is multiplied by seal_rise_multiplier and
# levels off on the scale of seal_saturation_effort.
reshape_seal_effort <- function(seal_effort, seal_rise_multiplier,
                                seal_saturation_effort) {
  rise <- pmax(seal_effort - seal_effort_at_calibration, 0)
  pmin(seal_effort, seal_effort_at_calibration) +
    seal_rise_multiplier * rise / (1 + rise / seal_saturation_effort)
}

run_model <- function(seal_rise_multiplier, seal_saturation_effort,
                      cod_extra_mortality, cod_moratorium_factor,
                      haddock_moratorium_factor) {
  effort <- effort_by_year
  effort[, "Seals"] <- reshape_seal_effort(effort[, "Seals"],
                                           seal_rise_multiplier,
                                           seal_saturation_effort)

  # Three stretches of years, because two of the parameters switch on part way.
  stretches <- list(c(first_fitted_year, cod_mortality_starts),
                    c(cod_mortality_starts, moratorium_starts),
                    c(moratorium_starts, last_year))
  abundance <- abundance_1978
  resource  <- resource_1978
  biomass_parts <- list()
  catch_parts   <- list()

  for (stretch in stretches) {
    stretch_params <- params
    if (stretch[1] >= cod_mortality_starts) {
      cod_is_mature <- stretch_params@w >=
        stretch_params@species_params["Atlantic Cod", "w_mat"]
      stretch_params@mu_b["Atlantic Cod", cod_is_mature] <-
        stretch_params@mu_b["Atlantic Cod", cod_is_mature] + cod_extra_mortality
    }
    if (stretch[1] >= moratorium_starts) {
      stretch_params@catchability["Bottom_Trawl", "Atlantic Cod"] <-
        stretch_params@catchability["Bottom_Trawl", "Atlantic Cod"] *
        cod_moratorium_factor
      stretch_params@catchability["Bottom_Trawl", "Haddock"] <-
        stretch_params@catchability["Bottom_Trawl", "Haddock"] *
        haddock_moratorium_factor
    }
    stretch_params@initial_n[]    <- abundance
    stretch_params@initial_n_pp[] <- resource
    sim <- quietly(project(
      stretch_params,
      effort = effort[as.character(stretch[1]:stretch[2]), , drop = FALSE],
      progress_bar = FALSE))

    yield_by_gear <- getYieldGear(sim)
    fishing_gears <- setdiff(dimnames(yield_by_gear)[[2]], "Seals")
    stretch_biomass <- getBiomass(sim)
    stretch_catch <- apply(yield_by_gear[, fishing_gears, , drop = FALSE],
                           c(1, 3), sum)
    # The last year of a stretch is the first of the next; keep it only once.
    is_final_stretch <- stretch[2] == last_year
    keep <- if (is_final_stretch) seq_len(nrow(stretch_biomass)) else
      seq_len(nrow(stretch_biomass) - 1)
    biomass_parts[[length(biomass_parts) + 1]] <- stretch_biomass[keep, , drop = FALSE]
    catch_parts[[length(catch_parts) + 1]]     <- stretch_catch[keep, , drop = FALSE]

    last_time <- dim(sim@n)[1]
    abundance <- sim@n[last_time, , ]
    resource  <- sim@n_pp[last_time, ]
  }
  list(biomass = do.call(rbind, biomass_parts),
       catch   = do.call(rbind, catch_parts))
}

# ---- Mismatch between a run and the observations -----------------------------

# Squared log10 error for each species, over the chosen years.
error_by_species <- function(modelled_matrix, observed_table, years) {
  observed_rows <- observed_table[observed_table$Year %in% years, ]
  observed_rows$modelled <- modelled_matrix[
    cbind(as.character(observed_rows$Year), observed_rows$Species)]
  observed_rows <- observed_rows[is.finite(observed_rows$modelled) &
                                   observed_rows$modelled > 0, ]
  observed_rows$squared_error <-
    (log10(observed_rows$modelled) - log10(observed_rows$observed))^2
  tapply(observed_rows$squared_error, observed_rows$Species, mean)
}

mismatch <- function(run, years) {
  mean(c(error_by_species(run$biomass, observed_biomass, years),
         error_by_species(run$catch,   observed_catch,   years)))
}

# The optimiser works on an unbounded-friendly scale: logs for the two seal
# parameters, which span orders of magnitude, and plain values for the rest.
parameter_names <- c("log_seal_rise_multiplier", "log_seal_saturation_effort",
                     "cod_extra_mortality", "cod_moratorium_factor",
                     "haddock_moratorium_factor")
lower_bounds <- c(log(0.01), log(0.005), 0,   0.001, 0.001)
upper_bounds <- c(log(5),    log(1000),  1.5, 1,     1)

run_from_vector <- function(values) {
  run_model(seal_rise_multiplier      = exp(values[1]),
            seal_saturation_effort    = exp(values[2]),
            cod_extra_mortality       = values[3],
            cod_moratorium_factor     = values[4],
            haddock_moratorium_factor = values[5])
}

fitted_years <- first_fitted_year:last_fitted_year
objective <- function(values) {
  run <- tryCatch(run_from_vector(values), error = function(e) NULL)
  if (is.null(run)) return(1e6)
  value <- mismatch(run, fitted_years)
  if (is.finite(value)) value else 1e6
}

# "No change": seal effort as in the file, no extra cod mortality, no
# moratorium effect. A saturation scale of 1000 is about a thousand times the
# largest rise in the effort series, so it is no saturation to within 0.1%.
no_change <- c(log(1), log(1000), 0, 1, 1)

# Check that the "no change" run is the 07 run, so any improvement below comes
# from the parameters and not from a difference in how the run is assembled.
baseline_run <- run_from_vector(no_change)
# Compared over the species that are fitted. Capelin and sculpin fall towards
# zero in this run, and a species at one part in a billion turns a 0.1% change
# in mortality into a large relative difference that means nothing.
baseline_biomass_from_07 <- getBiomass(baseline_sim)[rownames(baseline_run$biomass),
                                                     fitted_species]
restart_difference <- max(abs(log10(baseline_run$biomass[, fitted_species] /
                                      baseline_biomass_from_07)))
cat(sprintf("'no change' run against 07's projection: largest difference %.3f log10 units\n",
            restart_difference))
stopifnot(restart_difference < 0.02)

# ---- Fit ---------------------------------------------------------------------

# Several starting points, because a five-parameter fit to noisy series can have
# more than one basin. The best result is kept and the others are reported, so
# it is visible whether they agree.
starting_points <- list(
  no_change            = no_change,
  seals_saturate_early = c(log(1),   log(0.1), 0.3, 0.1, 0.1),
  weak_seals_high_codM = c(log(0.3), log(1),   0.6, 0.05, 0.05))

cluster <- makeCluster(min(length(lower_bounds) * 2 + 1, detectCores() - 1),
                       type = "FORK")
setDefaultCluster(cluster)
fits <- lapply(names(starting_points), function(start_name) {
  fit <- optimParallel(par = starting_points[[start_name]], fn = objective,
                       lower = lower_bounds, upper = upper_bounds,
                       control = list(maxit = 60, factr = 1e9))
  cat(sprintf("start '%s': mismatch %.4f -> %.4f after %d evaluations\n",
              start_name, objective(starting_points[[start_name]]), fit$value,
              fit$counts[["function"]]))
  fit
})
stopCluster(cluster)
names(fits) <- names(starting_points)

best_fit <- fits[[which.min(sapply(fits, function(fit) fit$value))]]
natural_scale <- function(values) c(
  seal_rise_multiplier      = exp(values[1]),
  seal_saturation_effort    = exp(values[2]),
  cod_extra_mortality       = values[3],
  cod_moratorium_factor     = values[4],
  haddock_moratorium_factor = values[5])

cat("\nfitted values from each starting point:\n")
print(round(t(sapply(fits, function(fit) c(natural_scale(fit$par),
                                            mismatch = fit$value))), 4))

fitted_values <- natural_scale(best_fit$par)
fitted_run <- run_from_vector(best_fit$par)

# ---- How much each parameter matters -----------------------------------------

# Put each fitted parameter back to its "no change" value in turn and see how
# much of the improvement is lost. A parameter whose removal costs almost
# nothing is not supported by the data, whatever value the optimiser gave it.
cat("\nmismatch over the fitted years (1978-2004); lower is better\n")
cat(sprintf("  %-34s %.4f\n", "no change (the 07 run)", objective(no_change)))
cat(sprintf("  %-34s %.4f\n", "all five fitted", best_fit$value))
parameter_groups <- list("without seal saturation (1 and 2)" = 1:2,
                         "without extra cod mortality (3)"   = 3,
                         "without cod moratorium (4)"        = 4,
                         "without haddock moratorium (5)"    = 5)
for (group_name in names(parameter_groups)) {
  values <- best_fit$par
  values[parameter_groups[[group_name]]] <- no_change[parameter_groups[[group_name]]]
  cat(sprintf("  %-34s %.4f\n", group_name, objective(values)))
}

# ---- Fit in the fitted years and in the held-out years -----------------------

skill_line <- function(modelled_matrix, observed_table, years) {
  observed_rows <- observed_table[observed_table$Year %in% years, ]
  observed_rows$modelled <- modelled_matrix[
    cbind(as.character(observed_rows$Year), observed_rows$Species)]
  per_species <- observed_rows |>
    filter(is.finite(modelled), modelled > 0) |>
    group_by(Species) |>
    filter(n() >= 4) |>
    summarise(bias = median(log10(modelled / observed)),
              correlation = suppressWarnings(cor(log(modelled), log(observed))),
              rmse = sqrt(mean((log10(modelled) - log10(observed))^2)),
              .groups = "drop")
  sprintf("bias x%.2f | within 2x %2d/%2d | median r %+.2f | median RMSE %.2f",
          10^median(per_species$bias), sum(abs(per_species$bias) < log10(2)),
          nrow(per_species), median(per_species$correlation, na.rm = TRUE),
          median(per_species$rmse))
}

held_out_years <- (last_fitted_year + 1):last_year
cat("\n16 species without survey multipliers; bias is modelled/observed\n")
for (run_name in c("no change", "fitted")) {
  run <- if (run_name == "no change") baseline_run else fitted_run
  cat(sprintf("%-10s biomass %d-%d: %s\n", run_name, first_fitted_year,
              last_fitted_year, skill_line(run$biomass, observed_biomass, fitted_years)))
  cat(sprintf("%-10s biomass %d-%d (held out): %s\n", run_name, last_fitted_year + 1,
              last_year, skill_line(run$biomass, observed_biomass, held_out_years)))
  cat(sprintf("%-10s catch   %d-%d: %s\n", run_name, first_fitted_year,
              last_fitted_year, skill_line(run$catch, observed_catch, fitted_years)))
  cat(sprintf("%-10s catch   %d-%d (held out): %s\n", run_name, last_fitted_year + 1,
              last_year, skill_line(run$catch, observed_catch, held_out_years)))
}

relative_to_calibration <- function(biomass_matrix, species) {
  mean(biomass_matrix[as.character(1994:2004), species]) /
    mean(biomass_matrix[as.character(1978:1981), species])
}
observed_relative <- function(species) {
  rows <- observed_biomass[observed_biomass$Species == species, ]
  mean(rows$observed[rows$Year %in% 1994:2004]) /
    mean(rows$observed[rows$Year %in% 1978:1981])
}
cat("\nbiomass in 1994-2004 relative to 1978-1981\n")
cat(sprintf("  %-20s %10s %8s %9s\n", "species", "no change", "fitted", "observed"))
for (species in c("Atlantic Cod", "Haddock", "American Plaice", "Witch Flounder",
                  "Yellowtail Flounder", "Winter Skate", "Thorny Skate", "Monkfish")) {
  cat(sprintf("  %-20s %10.2f %8.2f %9.2f\n", species,
              relative_to_calibration(baseline_run$biomass, species),
              relative_to_calibration(fitted_run$biomass, species),
              observed_relative(species)))
}
total_fish <- function(run, year) {
  sum(run$biomass[as.character(year), fitted_species]) / region_area_m2
}
cat(sprintf("\ntotal biomass of the 16 species in 2010: no change %.2f, fitted %.2f g/m2\n",
            total_fish(baseline_run, 2010), total_fish(fitted_run, 2010)))

# ---- Save --------------------------------------------------------------------

write.csv(data.frame(parameter = names(fitted_values),
                     fitted_value = signif(unname(fitted_values), 5),
                     no_change_value = c(1, Inf, 0, 1, 1)),
          fitted_values_file, row.names = FALSE)
saveRDS(list(biomass = fitted_run$biomass, catch = fitted_run$catch,
             baseline_biomass = baseline_run$biomass,
             baseline_catch = baseline_run$catch,
             fitted_values = fitted_values,
             all_fits = lapply(fits, function(fit) list(
               values = natural_scale(fit$par), mismatch = fit$value))),
        fitted_sim_file)
cat(sprintf("\nwrote %s and %s\n", fitted_values_file, fitted_sim_file))
