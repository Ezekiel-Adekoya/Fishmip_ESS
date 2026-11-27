# 8. Grey seal as a model species, eating what Bowen & Harrison measured
#
# Takes the 20-species model from scripts 1 and 2 (run on the raised sand lance
# base from script 7) and calibrates the seal's feeding so that, at the
# 1978-1981 state:
#   - each seal eats 1.73 t of fish a year (12 x 170 kg body mass, 85% fish);
#   - the fish it eats are in the Bowen & Harrison (1994) proportions by
#     weight: 80.9% sand lance, 11.0% cod, 2.6% silver hake, 1.4% plaice,
#     1.3% redfish, 1.1% yellowtail, 1.7% everything else. The species
#     Bowen & Harrison do not name share one preference, so the seal takes
#     them in proportion to how much of each is there. (Giving each of them
#     an equal share of the diet forced a very high preference for rare
#     species: halibut ended at the highest preference of all and was eaten
#     out in the projection.)
#   - the other 15% of its food comes from outside the model.
# The fish are re-matched to their observed biomass each time the seal changes.
#
# The diets were measured in 1991-1993 and 1993-1998 (seal_diet_data.R), not
# in 1978-1981. The shares the seal is calibrated to here can therefore be
# given from outside (REFERENCE_DIET_FILE, written by script 10), chosen so
# that the diet comes out close to the measured proportions in those years of
# the projection. Without that file the average of the two measured diets is
# used.
#
# Writes the calibrated model and the seal's feeding level to the SC folder.
# Script 5 then fits the seal's growth and harvest at that feeding level, and
# script 9 projects.

suppressPackageStartupMessages({library(mizer); library(therMizer)})
source("Ess_Fishmip_data/Useful_functions.R")
source("dynamic_seal_trial/seal_diet_data.R")
quietly <- function(expression) suppressWarnings(suppressMessages(expression))
region_area_m2 <- 146360555008
work_folder <- Sys.getenv("SC")

fish_per_seal_tonnes <- 12 * 170 * 0.85 / 1000      # 1.73 t of fish per seal a year
outside_food_per_fish <- 0.15 / 0.85                 # non-fish food, relative to fish
bowen_share <- starting_seal_diet_share
measured_share <- bowen_share
reference_diet_file <- Sys.getenv("REFERENCE_DIET_FILE")
if (nzchar(reference_diet_file) && file.exists(file.path(work_folder, reference_diet_file))) {
  reference_diet <- readRDS(file.path(work_folder, reference_diet_file))
  bowen_share <- reference_diet[names(measured_share)]
  cat("calibrating to the shares in", reference_diet_file, ":\n")
  print(round(100 * bowen_share, 2))
}
# From here on bowen_share is the diet to reach in the calibrated state.
# Preferred prey about 20 g with a wide spread: fitted in handover section 21
# to the mean prey weights Bowen & Harrison measured.
seal_weight_over_prey_weight <- 225000 / 20
seal_prey_size_spread <- 1.5

calibrated <- readRDS(file.path(work_folder, "seal_calibrated.rds"))
params <- calibrated$params
seal_target_thousand <- calibrated$target_seals
all_species <- as.character(params@species_params$species)
fish <- setdiff(all_species, "Grey Seal")
is_seal <- all_species == "Grey Seal"
fish_size_bins <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
adult_seal_sizes <- params@w >= 16000

# The share each named species should make of the fish eaten, and the 1.7%
# that the species Bowen & Harrison do not name should make together.
unnamed_species <- setdiff(fish, names(bowen_share))
unnamed_share_together <- 1 - sum(bowen_share)
target_share <- c(bowen_share, "unnamed species together" = unnamed_share_together)
target_total_kt <- fish_per_seal_tonnes * seal_target_thousand

seals_in_thousands <- function(p) sum(p@initial_n["Grey Seal", ] * p@dw) / 1000
# Weight of each fish species eaten by seals, kt a year.
eaten_by_seals <- function(p) {
  seal_pred_rate <- quietly(getPredRate(p))["Grey Seal", fish_size_bins]
  sapply(fish, function(species)
    sum(p@interaction["Grey Seal", species] * seal_pred_rate *
          p@initial_n[species, ] * p@w * p@dw)) / 1e9
}
solve_steady <- function(p) quietly(steady(p, t_max = 300, progress_bar = FALSE,
                                           preserve = "erepro"))

# ---- The seal eats fish and outside food only, at the measured prey sizes ----

given <- params@given_species_params
given$temp_min <- params@species_params$temp_min   # kept here so they survive
given$temp_max <- params@species_params$temp_max   # therMizer rebuilding the table
params@species_params$interaction_resource[is_seal] <- 0
params@species_params$beta[is_seal]  <- seal_weight_over_prey_weight
params@species_params$sigma[is_seal] <- seal_prey_size_spread
given$interaction_resource[given$species == "Grey Seal"] <- 0
given$beta[given$species == "Grey Seal"]  <- seal_weight_over_prey_weight
given$sigma[given$species == "Grey Seal"] <- seal_prey_size_spread
params@given_species_params <- given
params <- quietly(setPredKernel(params))

outside_food <- ext_encounter(params)
outside_food["Grey Seal", ] <- 0
ext_encounter(params) <- outside_food
interaction_matrix <- params@interaction
interaction_matrix["Grey Seal", names(bowen_share)] <- bowen_share / max(bowen_share)
interaction_matrix["Grey Seal", unnamed_species] <-
  unnamed_share_together / length(unnamed_species) / max(bowen_share)
params <- quietly(setInteraction(params, interaction = interaction_matrix))

# ---- Adjust preference and search rate until amount and proportions match ----

for (pass in 1:40) {
  eaten <- eaten_by_seals(params)
  share_now <- eaten / sum(eaten)
  amount_error <- sum(eaten) / target_total_kt
  share_error <- max(abs(share_now[names(bowen_share)] - bowen_share))
  if (pass > 1 && abs(amount_error - 1) < 0.03 && share_error < 0.01) break

  # Preference for each species moves towards its target share. The largest
  # preference is kept at 1; the search rate carries the overall amount.
  # The unnamed species all move by the same factor, so they keep one
  # preference between them.
  preference <- params@interaction["Grey Seal", fish]
  preference[names(bowen_share)] <- preference[names(bowen_share)] *
    (bowen_share / pmax(share_now[names(bowen_share)], 1e-12))^0.7
  preference[unnamed_species] <- preference[unnamed_species] *
    (unnamed_share_together / max(sum(share_now[unnamed_species]), 1e-12))^0.7
  interaction_matrix <- params@interaction
  interaction_matrix["Grey Seal", fish] <- preference / max(preference)
  params <- quietly(setInteraction(params, interaction = interaction_matrix))

  search_step <- min(2.5, max(0.5, (1 / amount_error)^0.6))
  search_volume <- params@search_vol
  search_volume["Grey Seal", ] <- search_volume["Grey Seal", ] * search_step
  params <- setSearchVolume(params, search_vol = search_volume)

  # Outside food in fixed proportion to the fish the seal encounters.
  outside_food <- ext_encounter(params)
  fish_encounter <- getEncounter(params)["Grey Seal", ] - outside_food["Grey Seal", ]
  outside_food["Grey Seal", ] <- outside_food_per_fish * fish_encounter
  ext_encounter(params) <- outside_food

  params <- solve_steady(params)
  seal_scale <- seal_target_thousand / seals_in_thousands(params)
  params@species_params$R_max[is_seal] <- params@species_params$R_max[is_seal] * seal_scale
  params@initial_n["Grey Seal", ] <- params@initial_n["Grey Seal", ] * seal_scale
  params <- quietly(matchBiomasses(params, species = fish))
  eaten <- eaten_by_seals(params)
  cat(sprintf("pass %2d: seals eat %5.1f kt of fish (target %.1f) | sand lance %4.1f%% cod %4.1f%% | fish %.2f g/m2 | seal feeding level %.2f\n",
              pass, sum(eaten), target_total_kt, 100 * eaten["Sand Lance"] / sum(eaten),
              100 * eaten["Atlantic Cod"] / sum(eaten),
              sum(getBiomass(params)[fish]) / region_area_m2,
              mean(getFeedingLevel(params)["Grey Seal", adult_seal_sizes])))
}
# Settle with the seals at their target number and the fish at theirs. The
# seal's recruitment limit is scaled with its numbers, or the next settling
# would undo the change.
for (pass in 1:5) {
  params <- solve_steady(params)
  seal_scale <- seal_target_thousand / seals_in_thousands(params)
  params@species_params$R_max[is_seal] <- params@species_params$R_max[is_seal] * seal_scale
  params@initial_n["Grey Seal", ] <- params@initial_n["Grey Seal", ] * seal_scale
  params <- quietly(matchBiomasses(params, species = fish))
}
params <- solve_steady(params)
stopifnot(abs(seals_in_thousands(params) / seal_target_thousand - 1) < 0.05)

# ---- Put each fish back at its own reproduction level --------------------------

# While the seal was being adjusted it passed through states with far too many
# seals, and steady() answered by raising the reproductive efficiency of the
# fish they were eating. With preserve = "erepro" that never comes back down:
# dogfish ended with erepro 20. Each fish is therefore returned to the
# reproduction level it has in the base model. That recomputes erepro and R_max
# for the state the model is in now, and does not move that state.
base_params <- readRDS(Sys.getenv("BASE_PARAMS"))
base_reproduction_level <- getReproductionLevel(base_params)[fish]
biomass_before <- getBiomass(params)
params <- quietly(setBevertonHolt(params, reproduction_level = base_reproduction_level))
# Dogfish sits close to the limit of erepro = 1 and the seal's small take of
# its young pushes it over. Any fish that would need erepro above 0.95 is given
# a lower reproduction level instead (a higher recruitment ceiling), just far
# enough to bring erepro back to 0.95.
largest_allowed_erepro <- 0.95
erepro_with_no_limit <- params@species_params[fish, "erepro"] * (1 - base_reproduction_level)
adjusted_level <- pmin(base_reproduction_level,
                       1 - erepro_with_no_limit / largest_allowed_erepro)
if (any(adjusted_level < base_reproduction_level - 1e-12)) {
  changed <- fish[adjusted_level < base_reproduction_level - 1e-12]
  cat("reproduction level lowered to keep erepro at or below", largest_allowed_erepro, ":",
      paste(sprintf("%s %.2f -> %.2f", changed, base_reproduction_level[changed],
                    adjusted_level[changed]), collapse = "; "), "\n")
  params <- setBevertonHolt(params, reproduction_level = adjusted_level)
}
stopifnot(max(params@species_params[fish, "erepro"]) <= 1,
          max(abs(getBiomass(params) / biomass_before - 1)) < 1e-8)

# ---- Report -------------------------------------------------------------------

species_params <- params@species_params
biomass_ratio <- getBiomass(params)[fish] / species_params[fish, "biomass_observed"]
yield_modelled <- colSums(getYieldGear(params))[fish]
yield_observed <- species_params[fish, "yield_observed"]
has_yield <- !is.na(yield_observed) & yield_observed > 0
yield_ratio <- (yield_modelled / yield_observed)[has_yield]
eaten <- eaten_by_seals(params)
seal_feeding_level <- mean(getFeedingLevel(params)["Grey Seal", adult_seal_sizes])

cat(sprintf("\ncalibrated with %.1f thousand seals: fish %.2f g/m2 | biomass within 2x %d/19 | yield within 2x %d/%d | fish with erepro > 1: %s\n",
            seals_in_thousands(params), sum(getBiomass(params)[fish]) / region_area_m2,
            sum(biomass_ratio > 0.5 & biomass_ratio < 2),
            sum(yield_ratio > 0.5 & yield_ratio < 2), length(yield_ratio),
            if (any(species_params[fish, "erepro"] > 1))
              paste(fish[species_params[fish, "erepro"] > 1], collapse = ", ") else "none"))
cat(sprintf("seals eat %.1f kt of fish a year, %.2f t per seal (target %.2f); feeding level %.2f\n",
            sum(eaten), sum(eaten) / seals_in_thousands(params), fish_per_seal_tonnes,
            seal_feeding_level))
cat("share of the fish eaten, model against the shares calibrated to:\n")
for (species in names(bowen_share)) {
  cat(sprintf("  %-20s %5.1f%%  (%.1f%%)   preference %.3f\n", species,
              100 * eaten[species] / sum(eaten), 100 * bowen_share[species],
              params@interaction["Grey Seal", species]))
}
cat(sprintf("  %-20s %5.1f%%  (%.1f%%)   preference %.3f, the same for each\n",
            "unnamed species", 100 * sum(eaten[unnamed_species]) / sum(eaten),
            100 * unnamed_share_together, params@interaction["Grey Seal", unnamed_species[1]]))
for (species in unnamed_species) {
  cat(sprintf("    %-20s %5.2f%%\n", species, 100 * eaten[species] / sum(eaten)))
}
# Mean weight of the individual fish eaten, to compare with the measured ones.
seal_pred_rate <- quietly(getPredRate(params))["Grey Seal", fish_size_bins]
measured_prey_weight <- c("Atlantic Cod" = 201, "American Plaice" = 132,
                          "Yellowtail Flounder" = 89, "Silver Hake" = 66,
                          "Redfish Unseparated" = 32, "Sand Lance" = 15)
cat("mean weight of fish eaten, g, model (measured):\n")
for (species in names(measured_prey_weight)) {
  numbers_eaten <- seal_pred_rate * params@initial_n[species, ] * params@dw
  cat(sprintf("  %-20s %5.0f  (%d)\n", species,
              sum(numbers_eaten * params@w) / sum(numbers_eaten),
              measured_prey_weight[species]))
}
sand_lance_biomass <- params@initial_n["Sand Lance", ] * params@w * params@dw
cat(sprintf("seal predation mortality on sand lance: %.3f per year (biomass weighted)\n",
            sum(params@interaction["Grey Seal", "Sand Lance"] * seal_pred_rate *
                  sand_lance_biomass) / sum(sand_lance_biomass)))

saveRDS(list(params = params, seal_feeding_level = seal_feeding_level,
             target_share = target_share, calibrated_share = bowen_share),
        file.path(work_folder, "seal_on_bowen_diet.rds"))
