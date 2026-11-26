# 9. Project the model with the dynamic seal, 1861-2010
#
# Takes the model calibrated in script 8 and the seal growth and harvest rates
# fitted in script 5, and runs them through the same years, effort, temperature
# and plankton as script 07.
#
#   1861-1961  spin-up, with seals held near their 1960 number
#   1961-1967  seals grow at their fitted rate
#   1967-1990  the same, less the fitted harvest
#   1991-2010  no harvest
#
# Nothing from script 10 is used: no starting multipliers and no change in
# catchability after the moratorium. As it stands the run shows what the seal
# does on its own account.
#
# With FIT_COD_MORTALITY=1 an extra natural mortality on mature cod is added
# as an episode: it starts in 1990, as in script 10, and lasts a fitted number
# of years, after which it is switched off. Its level and its length are fitted
# to the cod survey biomass for 1970-2004. 2005-2010 stay held out.
#
# Writes dynamic_seal_run.rds to the SC folder and prints the comparison with
# the seal counts, the Bowen & Harrison diet and the survey and catch series.

suppressPackageStartupMessages({
  library(mizer); library(therMizer); library(dplyr); library(tidyr)
})
source("Ess_Fishmip_data/Useful_functions.R")
quietly <- function(expression) suppressWarnings(suppressMessages(expression))
region_area_m2 <- 146360555008
work_folder <- Sys.getenv("SC")
multiplied_forage_species <- c("Atlantic Herring", "Sand Lance", "Capelin")
harvest_years <- 1967:1990
last_fitted_year <- 2004           # 2005-2010 are held out, as in script 10
cod_mortality_starts <- 1990       # the same starting year as script 10

first_build <- readRDS(file.path(work_folder, "seal_calibrated.rds"))   # forcing arrays
calibrated  <- readRDS(file.path(work_folder, "seal_on_bowen_diet.rds"))
demography  <- readRDS(file.path(work_folder, "seal_demography.rds"))$with_harvest
params <- calibrated$params
all_species <- as.character(params@species_params$species)
fish <- setdiff(all_species, "Grey Seal")
fitted_species <- setdiff(fish, multiplied_forage_species)
is_seal <- all_species == "Grey Seal"
fish_size_bins <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)

seal_counts_table <- readRDS(
  "~/Desktop/Ez June 2026/ESS_Restructured/Ess_data/seal_abundance_digitised.RDS")
seal_counts <- setNames(seal_counts_table$scotian_mean, seal_counts_table$year)

effort_table <- read.csv("Ess_Fishmip_data/Effort_spin_1861_2010.csv", check.names = FALSE)
effort_by_year <- as.matrix(effort_table[, names(params@initial_effort)])
rownames(effort_by_year) <- effort_table$year

seals_in_thousands <- function(abundance) sum(abundance["Grey Seal", ] * params@dw) / 1000

# ---- Time-varying temperature and plankton, as in script 07 ------------------

params_dynamic <- quietly(upgradeTherParams(
  params, ocean_temp = first_build$tmat, n_pp_array = first_build$u217,
  vertical_migration_array = first_build$vm, exposure_array = first_build$ex,
  aerobic_effect = TRUE, metabolism_effect = TRUE))
stopifnot(!anyNA(scaled_temp_effect(params_dynamic, 1861)),
          isTRUE(all.equal(params_dynamic@search_vol, params@search_vol)))

# ---- Spin-up with few seals --------------------------------------------------

calibration_seals <- seals_in_thousands(params@initial_n)
spin_up_params <- params_dynamic
spin_up_params@species_params$R_max[is_seal] <-
  params@species_params$R_max[is_seal] * seal_counts[["1960"]] / calibration_seals
spin_up_params@initial_n["Grey Seal", ] <-
  params@initial_n["Grey Seal", ] * seal_counts[["1960"]] / calibration_seals
spin_up <- quietly(project(spin_up_params,
                           effort = effort_by_year[as.character(1861:1961), ],
                           progress_bar = FALSE))
abundance_in_1961 <- spin_up@n[dim(spin_up@n)[1], , ]
resource_in_1961  <- spin_up@n_pp[dim(spin_up@n_pp)[1], ]
abundance_in_1961["Grey Seal", ] <- abundance_in_1961["Grey Seal", ] *
  seal_counts[["1961"]] / seals_in_thousands(abundance_in_1961)

# ---- 1961-2010 with seal growth and harvest ----------------------------------

# The extra cod mortality in each year of the run. It is at its full level for
# every whole year of the episode, a fraction of it in the year the episode
# ends part-way through, and zero before and after.
cod_extra_mortality_by_year <- function(run_years, cod_extra_mortality, cod_episode_years) {
  years_into_episode <- run_years - cod_mortality_starts
  share_of_year_in_episode <- pmin(pmax(cod_episode_years - years_into_episode, 0), 1)
  ifelse(years_into_episode < 0, 0, cod_extra_mortality * share_of_year_in_episode)
}

# One run from the 1961 state for a given seal reproductive efficiency and
# harvest rate, and optionally an episode of extra mortality on mature cod.
# The years are grouped into stretches over which nothing changes, and the
# function returns one simulation per stretch.
run_from_1961 <- function(seal_erepro, seal_harvest_rate,
                          cod_extra_mortality = 0, cod_episode_years = 0) {
  growing_params <- params_dynamic
  growing_params@species_params$erepro[is_seal] <- seal_erepro
  growing_params@species_params$R_max[is_seal]  <- demography[["R_max"]]
  cod_is_mature <- growing_params@w >=
    growing_params@species_params$w_mat[all_species == "Atlantic Cod"]

  # The settings for each year that a projection step starts in.
  run_years <- 1961:2009
  seal_harvest_in_year <- ifelse(run_years %in% harvest_years, seal_harvest_rate, 0)
  cod_extra_in_year <- cod_extra_mortality_by_year(run_years, cod_extra_mortality,
                                                   cod_episode_years)
  # A new stretch starts whenever either setting changes.
  setting_changes <- c(TRUE, diff(seal_harvest_in_year) != 0 | diff(cod_extra_in_year) != 0)
  stretch_first_years <- run_years[setting_changes]
  stretch_last_years  <- c(stretch_first_years[-1], 2010)

  abundance <- abundance_in_1961
  resource  <- resource_in_1961
  stretch_simulations <- list()
  for (stretch_number in seq_along(stretch_first_years)) {
    first_year <- stretch_first_years[stretch_number]
    last_year  <- stretch_last_years[stretch_number]
    stretch_params <- growing_params
    stretch_params@mu_b["Grey Seal", ] <- growing_params@mu_b["Grey Seal", ] +
      seal_harvest_in_year[run_years == first_year]
    stretch_params@mu_b["Atlantic Cod", cod_is_mature] <-
      growing_params@mu_b["Atlantic Cod", cod_is_mature] +
      cod_extra_in_year[run_years == first_year]
    stretch_params@initial_n[]    <- abundance
    stretch_params@initial_n_pp[] <- resource
    simulation <- quietly(project(
      stretch_params, effort = effort_by_year[as.character(first_year:last_year), ],
      progress_bar = FALSE))
    stretch_simulations[[length(stretch_simulations) + 1]] <- simulation
    abundance <- simulation@n[dim(simulation@n)[1], , ]
    resource  <- simulation@n_pp[dim(simulation@n_pp)[1], ]
  }
  stretch_simulations
}
seal_numbers_in <- function(stretch_simulations) {
  numbers <- unlist(lapply(stretch_simulations, function(simulation) {
    n <- getN(simulation)[, "Grey Seal"] / 1000
    n[-length(n)]
  }))
  last <- stretch_simulations[[length(stretch_simulations)]]
  c(numbers, tail(getN(last)[, "Grey Seal"] / 1000, 1))
}

# The growth and harvest rates from script 5 were fitted with the seal on its
# own, fed at a fixed level. In the full model its food changes from year to
# year, so with those rates it does not follow the counts. With
# FIT_SEALS_IN_FULL_MODEL=1 the two rates are fitted again here, in the full
# model, to the counted numbers for 1961-2010.
seal_erepro <- demography[["erepro"]]
seal_harvest_rate <- demography[["harvest_rate"]]
# Seal rates from an earlier fit can be given from outside, to repeat a run
# or to start the fit close to its answer.
if (nzchar(Sys.getenv("SEAL_EREPRO"))) seal_erepro <- as.numeric(Sys.getenv("SEAL_EREPRO"))
if (nzchar(Sys.getenv("SEAL_HARVEST_RATE"))) {
  seal_harvest_rate <- as.numeric(Sys.getenv("SEAL_HARVEST_RATE"))
}
cod_extra_mortality <- 0
cod_episode_years <- 0
fit_seals <- Sys.getenv("FIT_SEALS_IN_FULL_MODEL", "0") == "1"
fit_cod_mortality <- Sys.getenv("FIT_COD_MORTALITY", "0") == "1"

# Fits the two seal rates to the counts, with the cod mortality as it stands.
fit_seal_rates_to_counts <- function() {
  count_years <- as.character(1961:2010)
  mismatch_with_counts <- function(values) {
    numbers <- tryCatch(
      seal_numbers_in(run_from_1961(exp(values[1]), max(values[2], 0),
                                    cod_extra_mortality, cod_episode_years)),
      error = function(e) NULL)
    if (is.null(numbers) || any(!is.finite(numbers)) || any(numbers <= 0)) return(1e6)
    sum((log(numbers[count_years]) - log(seal_counts[count_years]))^2) +
      if (values[2] < 0) 1e3 * values[2]^2 else 0
  }
  fit <- optim(c(log(seal_erepro), seal_harvest_rate), mismatch_with_counts,
               control = list(maxit = 120, reltol = 1e-4))
  cat(sprintf("seal rates fitted in the full model: erepro %.4f -> %.4f | harvest rate %.3f -> %.3f per year | log RMSE against counts %.3f\n",
              seal_erepro, exp(fit$par[1]), seal_harvest_rate, max(fit$par[2], 0),
              sqrt(fit$value / length(count_years))))
  seal_erepro <<- exp(fit$par[1])
  seal_harvest_rate <<- max(fit$par[2], 0)
}

# Cod survey biomass over the fitted years, and the model's mismatch with it
# (mean squared difference in log10 biomass).
cod_survey <- readRDS("Ess_Fishmip_data/biomass_ts.RDS") |>
  filter(Species == "Atlantic Cod", Biomass > 0,
         Year >= 1970, Year <= last_fitted_year)
cod_biomass_in <- function(stretch_simulations) {
  unlist(lapply(stretch_simulations, function(simulation) {
    cod <- getBiomass(simulation)[, "Atlantic Cod"]
    cod[setdiff(names(cod), "2010")]
  }))
}
cod_mismatch_with_survey <- function(level, episode_years) {
  modelled <- tryCatch(
    cod_biomass_in(run_from_1961(seal_erepro, seal_harvest_rate, level, episode_years)),
    error = function(e) NULL)
  if (is.null(modelled)) return(1e6)
  modelled <- modelled[as.character(cod_survey$Year)]
  if (any(!is.finite(modelled)) || any(modelled <= 0)) return(1e6)
  mean((log10(modelled) - log10(cod_survey$Biomass))^2)
}
# The episode can last between 1 year and the whole of 1990-2010 (21 years).
# Values outside the allowed ranges are pulled back and penalised.
longest_episode <- 2011 - cod_mortality_starts
fit_cod_mortality_to_survey <- function() {
  objective <- function(values) {
    level <- max(values[1], 0)
    episode_years <- min(max(values[2], 1), longest_episode)
    cod_mismatch_with_survey(level, episode_years) +
      (values[1] - level)^2 + 0.01 * (values[2] - episode_years)^2
  }
  starting_points <- list(c(0.6, 6), c(0.6, 14), c(1.5, 4))
  fits <- lapply(starting_points, function(start)
    optim(start, objective, control = list(maxit = 150, reltol = 1e-5)))
  for (fit in fits) {
    cat(sprintf("  cod episode fit: level %.3f per year for %.1f years | cod mismatch %.4f\n",
                max(fit$par[1], 0), min(max(fit$par[2], 1), longest_episode), fit$value))
  }
  best <- fits[[which.min(sapply(fits, function(fit) fit$value))]]
  cod_extra_mortality <<- max(best$par[1], 0)
  cod_episode_years   <<- min(max(best$par[2], 1), longest_episode)
}

if (fit_seals) fit_seal_rates_to_counts()
if (fit_cod_mortality) {
  cat(sprintf("cod mismatch with the survey, 1970-%d, without extra mortality: %.4f\n",
              last_fitted_year, cod_mismatch_with_survey(0, 0)))
  fit_cod_mortality_to_survey()
  cat(sprintf("extra mortality on mature cod: %.3f per year from %d, ending in %.1f\n",
              cod_extra_mortality, cod_mortality_starts,
              cod_mortality_starts + cod_episode_years))
  cat("  the same, with the episode running to the end of 2010: cod mismatch ",
      sprintf("%.4f\n", cod_mismatch_with_survey(cod_extra_mortality, longest_episode)))
  # Fewer cod change the seals' food a little, so their rates are fitted once
  # more with the cod episode in place.
  if (fit_seals) fit_seal_rates_to_counts()
  cat(sprintf("cod mismatch with the survey after the last seal fit: %.4f\n",
              cod_mismatch_with_survey(cod_extra_mortality, cod_episode_years)))
}
simulations <- c(list(spin_up), run_from_1961(seal_erepro, seal_harvest_rate,
                                              cod_extra_mortality, cod_episode_years))
growing_params <- params_dynamic
growing_params@species_params$erepro[is_seal] <- seal_erepro
growing_params@species_params$R_max[is_seal]  <- demography[["R_max"]]

# Join the stretches. The last year of one is the first of the next.
join_stretches <- function(per_simulation) {
  parts <- lapply(seq_along(simulations), function(i) {
    values <- per_simulation(simulations[[i]])
    if (i < length(simulations)) values[-nrow(values), , drop = FALSE] else values
  })
  do.call(rbind, parts)
}
biomass <- join_stretches(getBiomass)
catch   <- join_stretches(function(s) apply(getYieldGear(s), c(1, 3), sum))[, fish]
years <- as.integer(rownames(biomass))
state_in <- function(year) {
  for (simulation in rev(simulations)) {
    if (as.character(year) %in% dimnames(simulation@n)[[1]]) {
      return(list(n = simulation@n[as.character(year), , ],
                  n_pp = simulation@n_pp[as.character(year), ]))
    }
  }
}
seal_numbers <- sapply(years, function(year) seals_in_thousands(state_in(year)$n))
names(seal_numbers) <- years

# What the seals eat each year, kt by species, and the mortality they cause.
seal_pred_rate_in <- function(year) {
  state <- state_in(year)
  quietly(getPredRate(growing_params, n = state$n, n_pp = state$n_pp, t = year))[
    "Grey Seal", fish_size_bins]
}
eaten <- t(sapply(years, function(year) {
  state <- state_in(year); pred_rate <- seal_pred_rate_in(year)
  sapply(fish, function(species)
    sum(params@interaction["Grey Seal", species] * pred_rate *
          state$n[species, ] * params@w * params@dw)) / 1e9
}))
rownames(eaten) <- years
seal_mortality_on <- function(species, year, sizes) {
  state <- state_in(year); pred_rate <- seal_pred_rate_in(year)
  weights <- state$n[species, sizes] * params@w[sizes] * params@dw[sizes]
  sum(params@interaction["Grey Seal", species] * pred_rate[sizes] * weights) / sum(weights)
}

saveRDS(list(biomass = biomass, catch = catch, seal_numbers = seal_numbers,
             eaten = eaten, demography = demography,
             seal_erepro = seal_erepro, seal_harvest_rate = seal_harvest_rate,
             cod_extra_mortality = cod_extra_mortality,
             cod_mortality_starts = cod_mortality_starts,
             cod_episode_years = cod_episode_years),
        file.path(work_folder, Sys.getenv("RUN_FILE", "dynamic_seal_run.rds")))

# ---- Report ------------------------------------------------------------------

cod_juvenile_sizes <- params@w >= 20 & params@w <= 1000
all_sizes <- params@w > 0
cat("\nyear | seals model / count (thousand) | fish eaten kt | t per seal | sand lance %, cod % | sand lance kt | seal mortality: sand lance, cod 20-1000 g | 16 species g/m2\n")
for (year in c(1961, 1970, 1980, 1990, 1995, 2000, 2005, 2010)) {
  y <- as.character(year)
  cat(sprintf("%d | %5.0f / %3.0f | %6.1f | %.2f | %4.1f, %4.1f | %6.0f | %.3f, %.3f | %.2f\n",
              year, seal_numbers[y], seal_counts[y], sum(eaten[y, ]),
              sum(eaten[y, ]) / seal_numbers[y],
              100 * eaten[y, "Sand Lance"] / sum(eaten[y, ]),
              100 * eaten[y, "Atlantic Cod"] / sum(eaten[y, ]),
              biomass[y, "Sand Lance"] / 1e9,
              seal_mortality_on("Sand Lance", year, all_sizes),
              seal_mortality_on("Atlantic Cod", year, cod_juvenile_sizes),
              sum(biomass[y, fitted_species]) / region_area_m2))
}

observed_biomass <- readRDS("Ess_Fishmip_data/biomass_ts.RDS") |>
  transmute(Year = as.integer(Year), Species = as.character(Species), observed = Biomass) |>
  filter(observed > 0, Species %in% fitted_species)
observed_catch <- readRDS("Ess_Fishmip_data/new_NAFO_4VW_Catches.RDS") |>
  pivot_longer(-Year, names_to = "Species", values_to = "observed") |>
  filter(!is.na(observed), observed > 0, Species %in% fitted_species)
skill_line <- function(modelled_matrix, observed_table, period) {
  rows <- observed_table[observed_table$Year %in% period, ]
  rows$modelled <- modelled_matrix[cbind(as.character(rows$Year), rows$Species)]
  per_species <- rows |> filter(is.finite(modelled), modelled > 0) |>
    group_by(Species) |> filter(n() >= 4) |>
    summarise(bias = median(log10(modelled / observed)),
              correlation = suppressWarnings(cor(log(modelled), log(observed))),
              rmse = sqrt(mean((log10(modelled) - log10(observed))^2)), .groups = "drop")
  sprintf("bias x%.2f | within 2x %2d/%2d | median r %+.2f | median RMSE %.2f",
          10^median(per_species$bias), sum(abs(per_species$bias) < log10(2)),
          nrow(per_species), median(per_species$correlation, na.rm = TRUE),
          median(per_species$rmse))
}
cat("\ndynamic seal run against data, 16 species; bias is modelled/observed\n")
for (period in list(1970:2004, 2005:2010)) {
  cat(sprintf("  biomass %d-%d: %s\n", min(period), max(period),
              skill_line(biomass, observed_biomass, period)))
  cat(sprintf("  catch   %d-%d: %s\n", min(period), max(period),
              skill_line(catch, observed_catch, period)))
}
for (period in list(1990:1993, 1994:2004, 2005:2010)) {
  cod_rows <- readRDS("Ess_Fishmip_data/biomass_ts.RDS") |> filter(Species == "Atlantic Cod")
  cat(sprintf("cod in %d-%d relative to 1978-1981: model %.2f, observed %.2f\n",
              min(period), max(period),
              mean(biomass[as.character(period), "Atlantic Cod"]) /
                mean(biomass[as.character(1978:1981), "Atlantic Cod"]),
              mean(cod_rows$Biomass[cod_rows$Year %in% period]) /
                mean(cod_rows$Biomass[cod_rows$Year %in% 1978:1981])))
}
cat("\nbiomass in 1994-2004 relative to 1978-1981: model (observed)\n")
all_observed <- readRDS("Ess_Fishmip_data/biomass_ts.RDS")
for (species in c("Atlantic Cod", "Haddock", "American Plaice", "Silver Hake",
                  "Redfish Unseparated", "Winter Skate", "Thorny Skate", "Sand Lance")) {
  rows <- all_observed[all_observed$Species == species, ]
  cat(sprintf("  %-20s %5.2f  (%.2f)\n", species,
              mean(biomass[as.character(1994:2004), species]) /
                mean(biomass[as.character(1978:1981), species]),
              mean(rows$Biomass[rows$Year %in% 1994:2004]) /
                mean(rows$Biomass[rows$Year %in% 1978:1981])))
}
