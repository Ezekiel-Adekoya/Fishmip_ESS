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
# Nothing from script 10 is used: no starting multipliers, no extra cod
# mortality, no change in catchability after the moratorium. This run shows
# what the seal does on its own account.
#
# Prey switching (seal_prey_switching.R, after Smout et al. 2014) is on unless
# PREY_SWITCHING=0: the seal's preference for a prey falls as that prey becomes
# scarce and rises as it becomes common. It has no effect in the calibrated
# state.
#
# The seal's diet was measured in 1991-1993 and in 1993-1998
# (seal_diet_data.R), so the script reports the model's diet in those years
# against the measurements. Script 10 uses that to adjust the calibration.
# NO_SWITCHING_FOR names prey (separated by commas) the seal does not switch
# on, for example NO_SWITCHING_FOR="Atlantic Cod".
#
# SEAL_RATES_FROM names an earlier run file to take the seal's three fitted
# values from, instead of fitting them again.
#
# Writes dynamic_seal_run.rds to the SC folder and prints the comparison with
# the seal counts, the Bowen & Harrison diet and the survey and catch series.

suppressPackageStartupMessages({
  library(mizer); library(therMizer); library(dplyr); library(tidyr)
})
source("Ess_Fishmip_data/Useful_functions.R")
source("dynamic_seal_trial/seal_prey_switching.R")
source("dynamic_seal_trial/seal_diet_data.R")
quietly <- function(expression) suppressWarnings(suppressMessages(expression))
region_area_m2 <- 146360555008
work_folder <- Sys.getenv("SC")
multiplied_forage_species <- c("Atlantic Herring", "Sand Lance", "Capelin")
harvest_years <- 1967:1990
prey_without_switching <- trimws(strsplit(Sys.getenv("NO_SWITCHING_FOR", ""), ",")[[1]])
prey_switching_on <- Sys.getenv("PREY_SWITCHING", "1") != "0"

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

# The seal's preference shifts with prey availability, measured against the
# calibrated state, where the shift is nil.
if (prey_switching_on) {
  # SWITCHING_EXPONENT tries a value other than Smout et al.'s 1.5.
  params_dynamic <- add_seal_prey_switching(
    params_dynamic, reference_n = params@initial_n,
    exponent = as.numeric(Sys.getenv("SWITCHING_EXPONENT", smout_switching_exponent)),
    no_switching_for = prey_without_switching)
}
# The seal's preference for each species given the abundances at one moment.
seal_preference_given <- function(abundance) {
  if (prey_switching_on) seal_preference_now(params_dynamic, abundance)["Grey Seal", ]
  else params@interaction["Grey Seal", ]
}

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

# One run from the 1961 state for a given seal reproductive efficiency and
# harvest rate. Returns the three stretches of simulation.
run_from_1961 <- function(seal_erepro, seal_harvest_rate, seal_pup_limit) {
  growing_params <- params_dynamic
  growing_params@species_params$erepro[is_seal] <- seal_erepro
  growing_params@species_params$R_max[is_seal]  <- seal_pup_limit
  stretches <- list(list(from = 1961, to = 1967, harvest = 0),
                    list(from = 1967, to = 1991, harvest = seal_harvest_rate),
                    list(from = 1991, to = 2010, harvest = 0))
  abundance <- abundance_in_1961
  resource  <- resource_in_1961
  stretch_simulations <- list()
  for (stretch in stretches) {
    stretch_params <- growing_params
    stretch_params@mu_b["Grey Seal", ] <- growing_params@mu_b["Grey Seal", ] + stretch$harvest
    stretch_params@initial_n[]    <- abundance
    stretch_params@initial_n_pp[] <- resource
    simulation <- quietly(project(
      stretch_params, effort = effort_by_year[as.character(stretch$from:stretch$to), ],
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
#
# The counts level off after about 2000 while food per seal is still falling
# only slowly, so food alone does not stop the model seals. The third fitted
# value is the limit on the number of pups recruited each year (the seal's
# R_max): as the colony fills, each extra adult adds fewer surviving pups.
#
# The fit is to the counts themselves (thousands of seals), not their logs.
# On the log scale the many early years of few seals dominate and the
# levelling-off after 2000 carries almost no weight. SEAL_FIT_SCALE=log gives
# the log fit back.
seal_erepro <- demography[["erepro"]]
seal_harvest_rate <- demography[["harvest_rate"]]
seal_pup_limit <- demography[["R_max"]]
seal_fit_scale <- Sys.getenv("SEAL_FIT_SCALE", "natural")
if (nzchar(Sys.getenv("SEAL_RATES_FROM"))) {
  earlier_run <- readRDS(file.path(work_folder, Sys.getenv("SEAL_RATES_FROM")))
  seal_erepro <- earlier_run$seal_erepro
  seal_harvest_rate <- earlier_run$seal_harvest_rate
  seal_pup_limit <- earlier_run$seal_pup_limit
}
if (Sys.getenv("FIT_SEALS_IN_FULL_MODEL", "0") == "1") {
  count_years <- as.character(1961:2010)
  mismatch_with_counts <- function(values) {
    numbers <- tryCatch(seal_numbers_in(run_from_1961(exp(values[1]), max(values[2], 0),
                                                      exp(values[3]))),
                        error = function(e) NULL)
    if (is.null(numbers) || any(!is.finite(numbers)) || any(numbers <= 0)) return(1e6)
    mismatch <- if (seal_fit_scale == "log") {
      sum((log(numbers[count_years]) - log(seal_counts[count_years]))^2)
    } else {
      # divided by the mean count squared, so the size is comparable to the log fit
      sum((numbers[count_years] - seal_counts[count_years])^2) /
        mean(seal_counts[count_years])^2
    }
    mismatch + if (values[2] < 0) 1e3 * values[2]^2 else 0
  }
  # Two starting points for the pup limit, since a limit set far too high
  # has no effect and gives the fit nothing to follow.
  fits <- lapply(c(1, 0.1), function(pup_limit_start_factor)
    optim(c(log(seal_erepro), seal_harvest_rate, log(seal_pup_limit * pup_limit_start_factor)),
          mismatch_with_counts, control = list(maxit = 300, reltol = 1e-5)))
  fit <- fits[[which.min(sapply(fits, function(one_fit) one_fit$value))]]
  cat(sprintf("seal rates fitted in the full model (%s scale): erepro %.4f -> %.4f | harvest rate %.3f -> %.3f per year | pup limit %.3g -> %.3g | mismatch %.4f (other start %.4f)\n",
              seal_fit_scale, seal_erepro, exp(fit$par[1]), seal_harvest_rate,
              max(fit$par[2], 0), seal_pup_limit, exp(fit$par[3]), fit$value,
              max(sapply(fits, function(one_fit) one_fit$value))))
  seal_erepro <- exp(fit$par[1])
  seal_harvest_rate <- max(fit$par[2], 0)
  seal_pup_limit <- exp(fit$par[3])
}
simulations <- c(list(spin_up), run_from_1961(seal_erepro, seal_harvest_rate, seal_pup_limit))
growing_params <- params_dynamic
growing_params@species_params$erepro[is_seal] <- seal_erepro
growing_params@species_params$R_max[is_seal]  <- seal_pup_limit

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
count_years <- as.character(1961:2010)
cat(sprintf("seals against counts, 1961-2010: log RMSE %.3f | RMSE %.1f thousand\n",
            sqrt(mean((log(seal_numbers[count_years]) - log(seal_counts[count_years]))^2)),
            sqrt(mean((seal_numbers[count_years] - seal_counts[count_years])^2))))

# What the seals eat each year, kt by species, and the mortality they cause.
seal_pred_rate_in <- function(year) {
  state <- state_in(year)
  quietly(getPredRate(growing_params, n = state$n, n_pp = state$n_pp, t = year))[
    "Grey Seal", fish_size_bins]
}
eaten <- t(sapply(years, function(year) {
  state <- state_in(year); pred_rate <- seal_pred_rate_in(year)
  preference <- seal_preference_given(state$n)
  sapply(fish, function(species)
    sum(preference[species] * pred_rate *
          state$n[species, ] * params@w * params@dw)) / 1e9
}))
rownames(eaten) <- years
seal_mortality_on <- function(species, year, sizes) {
  state <- state_in(year); pred_rate <- seal_pred_rate_in(year)
  weights <- state$n[species, sizes] * params@w[sizes] * params@dw[sizes]
  sum(seal_preference_given(state$n)[species] * pred_rate[sizes] * weights) / sum(weights)
}
# The seal's diet in each period the diet was measured in.
model_diet_share <- t(sapply(seal_diet_years, function(period_years) {
  eaten_in_period <- colSums(eaten[as.character(period_years), ])
  eaten_in_period / sum(eaten_in_period)
}))

saveRDS(list(biomass = biomass, catch = catch, seal_numbers = seal_numbers,
             eaten = eaten, demography = demography,
             seal_erepro = seal_erepro, seal_harvest_rate = seal_harvest_rate,
             seal_pup_limit = seal_pup_limit, prey_switching_on = prey_switching_on,
             prey_without_switching = prey_without_switching,
             model_diet_share = model_diet_share, measured_diet_share = seal_diet_share),
        file.path(work_folder, "dynamic_seal_run.rds"))

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

cat(sprintf("\nprey switching: %s%s\n", if (prey_switching_on) "on" else "off",
            if (prey_switching_on && length(prey_without_switching) > 0)
              paste0(", but not for ", paste(prey_without_switching, collapse = ", ")) else ""))
cat("seal diet, percent of fish eaten: model (measured)\n")
cat(sprintf("  %-20s %-16s %-16s\n", "", rownames(seal_diet_share)[1], rownames(seal_diet_share)[2]))
for (species in colnames(seal_diet_share)) {
  cat(sprintf("  %-20s %5.1f (%4.1f)     %5.1f (%4.1f)\n", species,
              100 * model_diet_share[1, species], 100 * seal_diet_share[1, species],
              100 * model_diet_share[2, species], 100 * seal_diet_share[2, species]))
}
cat(sprintf("  mean absolute gap over the named species: %.1f and %.1f percent of the diet\n",
            100 * mean(abs(model_diet_share[1, colnames(seal_diet_share)] - seal_diet_share[1, ]), na.rm = TRUE),
            100 * mean(abs(model_diet_share[2, colnames(seal_diet_share)] - seal_diet_share[2, ]), na.rm = TRUE)))
cat("seal preference for cod, relative to the calibrated state:",
    paste(sprintf("%d %.2f", c(1980, 1990, 2000, 2010),
                  sapply(c(1980, 1990, 2000, 2010), function(year)
                    seal_preference_given(state_in(year)$n)["Atlantic Cod"] /
                      params@interaction["Grey Seal", "Atlantic Cod"])), collapse = " | "), "\n")

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
