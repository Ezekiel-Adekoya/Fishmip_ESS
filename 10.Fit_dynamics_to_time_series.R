# Fit a small number of dynamic parameters to the catch and biomass time series.
#
# WHAT THIS IS, AND WHAT IT IS NOT
# 06 calibrates the model at a steady state (1978-1981 averages). Nothing in 06
# is tuned against how catch and biomass CHANGE through time. This script does
# that next step: it estimates five parameters that only act after the
# calibration period, plus a starting abundance for each species, by minimising
# the mismatch with the observed series.
#
# The calibrated parameters are not touched. The five dynamic parameters below
# are defined so that they have no effect in 1978-1981, and the starting
# abundances change where the model begins, not how it behaves.
#
# THE FIVE DYNAMIC PARAMETERS
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
# THE STARTING MULTIPLIERS
# One multiplier per fitted species, applied to its abundance at the start of
# 1961. Why they are needed: the model is calibrated to 1978-1981 but run from
# 1861, and effort before 1978 differs from the 1978-1981 average, so by the
# time the data begin some species have drifted far from the state they were
# calibrated at. Spiny dogfish arrives in 1978 at 1% of its calibrated biomass
# and thorny skate at nearly three times it, and the slow species never recover
# from that within the data period. A multiplier lets the data say where each
# species should start.
#
# They are applied in 1961, not 1978, because FishMIP output starts in 1961 and
# the correction should be in place for the whole submitted period. Species
# that turn over quickly forget their start within a decade; their multipliers
# are not pinned down by data from 1970 on and should stay near 1. The script
# reports which ones moved and how much the fit depends on them.
#
# THE PENALTY ON THE STARTING MULTIPLIERS
# Left free, the optimiser gave some species multipliers in the hundreds to
# chase one or two early survey points, for almost no gain in fit. So each
# multiplier now has to pay for itself: the objective is the mismatch plus
#     start_penalty_weight * mean(log10(multiplier)^2).
# With 16 fitted species, changing one species' start tenfold adds
# start_penalty_weight / 16 to the objective, so it only happens if the
# mismatch falls by more than that. Every species is treated the same way and
# no species list is chosen by hand. Multipliers the data support stay; the
# rest settle at 1.
#
# WHAT IS FITTED TO WHAT
# Survey biomass and fishery catch by species, 1970-2004, for the 16 species
# whose survey biomass carries no catchability multiplier. Catch means the four
# fishing gears only; seals are a gear in this model and what they eat is not
# catch.
#
# 2005-2010 IS HELD OUT, as the FishMIP regional protocol asks (Ortega-Cisneros
# et al. 2025, step 6: those years of catch are kept back to test predictive
# skill). Scores for 2005-2010 are printed separately at the end and are a
# genuine out-of-sample test.
#
# Fitting all years instead (last_fitted_year <- 2010) was tried once for
# comparison on 2026-11-27. It gave similar parameter values and a much better
# 2005-2010 fit, as it must, but leaves nothing to test the model against. The
# numbers are in handover section 22.
#
# The error is the squared difference of log10 values, averaged within each
# species and each variable first, so a species with many survey years does not
# outweigh one with few.
#
# Every run starts from the model state in 1961 taken from 07's projection, so
# the long spin-up is not repeated and one evaluation takes a second or two.
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

first_run_year     <- 1961      # start of the real forcing and of FishMIP output
first_fitted_year  <- 1970      # the survey and catch series begin here
last_fitted_year   <- 2004      # 2005-2010 are held out, as the FishMIP protocol asks
last_year          <- 2010
cod_mortality_starts <- 1990
moratorium_starts    <- 1993

# How strongly the starting multipliers are pulled towards 1 (see the header).
# 0 switches the penalty off. It can be set from outside to test other values,
# for example START_PENALTY_WEIGHT=0.01 Rscript 10.Fit_dynamics_to_time_series.R
start_penalty_weight <- as.numeric(Sys.getenv("START_PENALTY_WEIGHT", "0.1"))

fitted_values_file <- "Ess_Fishmip_data/fitted_dynamic_parameters.csv"
fitted_sim_file    <- "Ess_Fishmip_data/sims_GFDL_fitted.RDS"

quietly <- function(expression) suppressWarnings(suppressMessages(expression))

# ---- The model and its state in 1978 -----------------------------------------

baseline_sim <- readRDS("Ess_Fishmip_data/sims_GFDL_recalibrated.RDS")
params <- baseline_sim@params
species_names <- as.character(params@species_params$species)
fitted_species <- setdiff(species_names, multiplied_forage_species)

start_index <- which(dimnames(baseline_sim@n)[[1]] == as.character(first_run_year))
abundance_at_start <- baseline_sim@n[start_index, , ]
resource_at_start  <- baseline_sim@n_pp[start_index, ]

effort_table <- read.csv("Ess_Fishmip_data/Effort_spin_1861_2010.csv",
                         check.names = FALSE)
effort_by_year <- as.matrix(effort_table[, setdiff(names(effort_table), "year")])
rownames(effort_by_year) <- as.character(effort_table$year)
effort_by_year <- effort_by_year[as.character(first_run_year:last_year),
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
                      haddock_moratorium_factor, start_multipliers) {
  effort <- effort_by_year
  effort[, "Seals"] <- reshape_seal_effort(effort[, "Seals"],
                                           seal_rise_multiplier,
                                           seal_saturation_effort)

  # Three stretches of years, because two of the parameters switch on part way.
  stretches <- list(c(first_run_year, cod_mortality_starts),
                    c(cod_mortality_starts, moratorium_starts),
                    c(moratorium_starts, last_year))
  # Each fitted species starts 1961 at its spun-up abundance times its
  # multiplier, the same factor at every size.
  abundance <- abundance_at_start
  abundance[names(start_multipliers), ] <-
    abundance[names(start_multipliers), ] * start_multipliers
  resource  <- resource_at_start
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

# The optimiser works on scales where a unit step means something similar for
# every parameter: natural logs for the two seal parameters, plain values for
# the cod and moratorium ones, and log10 for the starting multipliers (so 1
# means ten times, -1 means a tenth).
dynamic_parameter_names <- c("log_seal_rise_multiplier",
                             "log_seal_saturation_effort",
                             "cod_extra_mortality", "cod_moratorium_factor",
                             "haddock_moratorium_factor")
number_of_dynamic <- length(dynamic_parameter_names)
multiplier_positions <- number_of_dynamic + seq_along(fitted_species)
lower_bounds <- c(log(0.01), log(0.005), 0,   0.001, 0.001,
                  rep(-3, length(fitted_species)))
upper_bounds <- c(log(5),    log(1000),  1.5, 1,     1,
                  rep(3, length(fitted_species)))

run_from_vector <- function(values) {
  run_model(seal_rise_multiplier      = exp(values[1]),
            seal_saturation_effort    = exp(values[2]),
            cod_extra_mortality       = values[3],
            cod_moratorium_factor     = values[4],
            haddock_moratorium_factor = values[5],
            start_multipliers = setNames(10^values[multiplier_positions],
                                         fitted_species))
}

fitted_years <- first_fitted_year:last_fitted_year

# The mismatch with the data alone. This is what the reports below quote, so
# the numbers can be compared between fits with different penalty weights.
data_mismatch <- function(values) {
  run <- tryCatch(run_from_vector(values), error = function(e) NULL)
  if (is.null(run)) return(1e6)
  value <- mismatch(run, fitted_years)
  if (is.finite(value)) value else 1e6
}

# What the starting multipliers cost: zero when every multiplier is 1.
start_penalty <- function(values) {
  start_penalty_weight * mean(values[multiplier_positions]^2)
}

# What the optimiser minimises.
objective <- function(values) data_mismatch(values) + start_penalty(values)

# "No change": seal effort as in the file, no extra cod mortality, no
# moratorium effect. A saturation scale of 1000 is about a thousand times the
# largest rise in the effort series, so it is no saturation to within 0.1%.
no_change <- c(log(1), log(1000), 0, 1, 1, rep(0, length(fitted_species)))

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

# Several starting points, because a fit with this many parameters to noisy
# series can have more than one basin. The best result is kept and the others
# are reported, so it is visible whether they agree. The third start gives the
# species known to be displaced a head start in the right direction.
informed_multipliers <- setNames(rep(0, length(fitted_species)), fitted_species)
informed_multipliers[c("Spiny Dogfish", "Thorny Skate", "Longhorn Sculpin",
                       "Winter Skate")] <- c(1.5, -0.4, -0.3, 0.3)
starting_points <- list(
  no_change            = no_change,
  seals_saturate_early = c(log(1), log(0.1), 0.3, 0.1, 0.1,
                           rep(0, length(fitted_species))),
  informed_starts      = c(log(0.5), log(1000), 0.6, 0.3, 0.2,
                           unname(informed_multipliers)))

cluster <- makeCluster(min(length(lower_bounds) * 2 + 1, detectCores() - 1),
                       type = "FORK")
setDefaultCluster(cluster)
fits <- lapply(names(starting_points), function(start_name) {
  fit <- optimParallel(par = starting_points[[start_name]], fn = objective,
                       lower = lower_bounds, upper = upper_bounds,
                       control = list(maxit = 150, factr = 1e8))
  cat(sprintf("start '%s': mismatch plus penalty %.4f -> %.4f after %d evaluations\n",
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
  haddock_moratorium_factor = values[5],
  setNames(10^values[multiplier_positions], paste0("start_x_", fitted_species)))

cat("\nfitted dynamic parameters from each starting point:\n")
print(round(t(sapply(fits, function(fit) c(
  natural_scale(fit$par)[seq_len(number_of_dynamic)],
  mismatch = data_mismatch(fit$par), penalty = start_penalty(fit$par)))), 4))
cat("\nfitted starting multipliers from each starting point:\n")
multiplier_table <- sapply(fits, function(fit) 10^fit$par[multiplier_positions])
rownames(multiplier_table) <- fitted_species
print(signif(multiplier_table, 3))

fitted_values <- natural_scale(best_fit$par)
fitted_run <- run_from_vector(best_fit$par)

# ---- How much each parameter matters -----------------------------------------

# Put each fitted parameter back to its "no change" value in turn and see how
# much of the improvement is lost. A parameter whose removal costs almost
# nothing is not supported by the data, whatever value the optimiser gave it.
cat(sprintf("\nmismatch over the fitted years (%d-%d); lower is better\n",
            first_fitted_year, last_fitted_year))
cat(sprintf("  (starting multiplier penalty weight %g; the penalty is not included below)\n",
            start_penalty_weight))
best_mismatch <- data_mismatch(best_fit$par)
cat(sprintf("  %-34s %.4f\n", "no change (the 07 run)", data_mismatch(no_change)))
cat(sprintf("  %-34s %.4f\n", "everything fitted", best_mismatch))
parameter_groups <- list("without starting multipliers"      = multiplier_positions,
                         "without the five dynamic ones"     = seq_len(number_of_dynamic),
                         "without seal saturation (1 and 2)" = 1:2,
                         "without extra cod mortality (3)"   = 3,
                         "without cod moratorium (4)"        = 4,
                         "without haddock moratorium (5)"    = 5)
for (group_name in names(parameter_groups)) {
  values <- best_fit$par
  values[parameter_groups[[group_name]]] <- no_change[parameter_groups[[group_name]]]
  cat(sprintf("  %-34s %.4f\n", group_name, data_mismatch(values)))
}

# How much each starting multiplier matters on its own: reset one species at a
# time to a multiplier of 1 and see what the fit loses. Those that cost nothing
# are not identified by the data and their fitted values should not be read as
# meaningful.
cat("\nstarting multipliers: fitted value, and the mismatch lost if reset to 1\n")
multiplier_cost <- sapply(seq_along(fitted_species), function(i) {
  values <- best_fit$par
  values[multiplier_positions[i]] <- 0
  data_mismatch(values) - best_mismatch
})
multiplier_report <- data.frame(
  species = fitted_species,
  start_multiplier = signif(10^best_fit$par[multiplier_positions], 3),
  mismatch_lost_if_reset = round(multiplier_cost, 4))
print(multiplier_report[order(-multiplier_report$mismatch_lost_if_reset), ],
      row.names = FALSE)

# ---- Fit in the earlier and the later years ----------------------------------

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

# The two periods are always reported separately so results stay comparable
# whichever years were fitted. The label says whether the later one was used.
earlier_years <- first_fitted_year:2004
later_years   <- 2005:last_year
later_label <- if (last_fitted_year >= last_year) "also fitted" else "held out"
cat("\n16 species without survey multipliers; bias is modelled/observed\n")
for (run_name in c("no change", "fitted")) {
  run <- if (run_name == "no change") baseline_run else fitted_run
  cat(sprintf("%-10s biomass %d-2004: %s\n", run_name, first_fitted_year,
              skill_line(run$biomass, observed_biomass, earlier_years)))
  cat(sprintf("%-10s biomass 2005-%d (%s): %s\n", run_name, last_year, later_label,
              skill_line(run$biomass, observed_biomass, later_years)))
  cat(sprintf("%-10s catch   %d-2004: %s\n", run_name, first_fitted_year,
              skill_line(run$catch, observed_catch, earlier_years)))
  cat(sprintf("%-10s catch   2005-%d (%s): %s\n", run_name, last_year, later_label,
              skill_line(run$catch, observed_catch, later_years)))
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
                  "Yellowtail Flounder", "Winter Skate", "Thorny Skate", "Monkfish",
                  "Spiny Dogfish", "Pollock")) {
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
                     no_change_value = c(1, Inf, 0, 1, 1,
                                         rep(1, length(fitted_species)))),
          fitted_values_file, row.names = FALSE)
saveRDS(list(biomass = fitted_run$biomass, catch = fitted_run$catch,
             baseline_biomass = baseline_run$biomass,
             baseline_catch = baseline_run$catch,
             fitted_values = fitted_values,
             start_penalty_weight = start_penalty_weight,
             all_fits = lapply(fits, function(fit) list(
               values = natural_scale(fit$par),
               mismatch = data_mismatch(fit$par)))),
        fitted_sim_file)
cat(sprintf("\nwrote %s and %s\n", fitted_values_file, fitted_sim_file))
