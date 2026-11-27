# 7. A fish base with realistic stocks of sand lance and herring
#
# The pipeline model is calibrated to the raw trawl survey figures for the
# forage fish: 8.2 kt of sand lance and 6.1 kt of herring. Trawls barely catch
# either.
#
#   Sand lance: Bowen & Harrison found it to be four fifths of what grey seals
#   eat. At 1.73 t of fish per seal a year, the 15 thousand seals of 1978-1981
#   need about 21 kt of sand lance a year and the 275 thousand of 2010 about
#   385 kt. A dynamic seal cannot be built on 8 kt.
#
#   Herring: reported 4VW catches were 230 kt in 1970 and 24-39 kt in 1980 and
#   1990, several times the 6 kt stock. The model can only take a fraction of
#   that, with a fishing mortality near 1 a year, and the stock cannot respond
#   when fishing and predation ease.
#
# This script raises the two targets in steps, recalibrating the whole
# community each time exactly as script 06 does. For herring it also lowers
# the catchability of the pelagic gear as the stock is raised, so the catch in the
# calibrated state ends at the observed one. It then loosens the limit on the
# recruitment of each, so the stocks can respond when their predators and
# fishing change.
#
# The target levels are working assumptions, not measurements: no reliable
# biomass estimate for either on the eastern Scotian Shelf was found. They are
# set with SAND_LANCE_TARGET_KT (default 400) and HERRING_TARGET_KT (default
# 700, about the level of the Paper 1 model; 0 leaves herring as it is).
#
# Reads Ess_Fishmip_data/params_GFDL_steady.RDS. Writes the new base to the
# folder named by the environment variable SC.

suppressPackageStartupMessages({library(mizer); library(therMizer)})
source("Ess_Fishmip_data/Useful_functions.R")
quietly <- function(expression) suppressWarnings(suppressMessages(expression))
region_area_m2 <- 146360555008
output_folder <- Sys.getenv("SC")
stopifnot(nzchar(output_folder), dir.exists(output_folder))

sand_lance_target <- 1e9 * as.numeric(Sys.getenv("SAND_LANCE_TARGET_KT", "400"))
herring_target    <- 1e9 * as.numeric(Sys.getenv("HERRING_TARGET_KT", "700"))
# Reproduction level after the recalibration. The pipeline model has both
# species at 0.99 or above, which pins recruitment at its ceiling: the stock is
# then very stable but can never grow. A lower level lets recruitment rise
# (0.5 lets it double) but also lets it fall when the spawning stock falls.
# At 0.5 sand lance collapsed during the spin-up, and herring sat at about
# 100 kt from 1961 to 1980 with a 1970 catch of 30 kt against 230 kt reported.
# At 0.9 herring holds about 420 kt before 1980, rises to 700-760 kt in the
# 1990s, and its catches come out near the reported ones. Both default to 0.9.
sand_lance_reproduction_level <- as.numeric(
  Sys.getenv("SAND_LANCE_REPRODUCTION_LEVEL", "0.9"))
herring_reproduction_level <- as.numeric(
  Sys.getenv("HERRING_REPRODUCTION_LEVEL", "0.9"))

params <- readRDS("Ess_Fishmip_data/params_GFDL_steady.RDS")
solve_steady <- function(p) quietly(steady(p, t_max = 300, progress_bar = FALSE,
                                           preserve = "erepro"))
fishing_gears <- setdiff(names(params@initial_effort), "Seals")
modelled_yield <- function(p, species) sum(getYieldGear(p)[fishing_gears, species])

# Sets the catchability of the gear that takes herring so that the catch in
# the calibrated state is the observed one.
match_herring_catch <- function(p) {
  gears <- gear_params(p)
  herring_row <- which(gears$species == "Atlantic Herring" & gears$gear == "Pelagic")
  observed_yield <- p@species_params["Atlantic Herring", "yield_observed"]
  for (pass in 1:5) {
    gears$catchability[herring_row] <- gears$catchability[herring_row] *
      observed_yield / modelled_yield(p, "Atlantic Herring")
    gear_params(p) <- gears
  }
  p
}

# Raises the calibration target of one species in steps of a half, with the
# whole community recalibrated at each step.
raise_target <- function(p, species, final_target, keep_catch = FALSE) {
  target_now <- p@species_params[species, "biomass_observed"]
  while (target_now < final_target) {
    previous_target <- target_now
    target_now <- min(final_target, target_now * 1.5)
    p@species_params[species, "biomass_observed"] <- target_now
    if (keep_catch) {
      # Lower the catchability by the size of the step before recalibrating,
      # so the catch asked of the stock does not jump with its biomass. (The
      # observed herring catch is several times the starting stock, and
      # without this the first steps do not settle.)
      gears <- gear_params(p)
      herring_row <- which(gears$species == species & gears$gear == "Pelagic")
      gears$catchability[herring_row] <- gears$catchability[herring_row] *
        previous_target / target_now
      gear_params(p) <- gears
    }
    given <- p@given_species_params
    given$biomass_observed[given$species == species] <- target_now
    p@given_species_params <- given
    for (pass in 1:3) p <- solve_steady(quietly(matchBiomasses(p)))
    cat(sprintf("%s target %4.0f kt: model %4.0f kt | catch %5.1f kt | all species %.2f g/m2 | worst erepro %.3f\n",
                species, target_now / 1e9, getBiomass(p)[species] / 1e9,
                modelled_yield(p, species) / 1e9,
                sum(getBiomass(p)) / region_area_m2, max(p@species_params$erepro)))
  }
  p
}

params <- raise_target(params, "Sand Lance", sand_lance_target)
if (herring_target > 0) {
  params <- raise_target(params, "Atlantic Herring", herring_target, keep_catch = TRUE)
}
# Bring the herring catch to the observed one, then settle. Changing the
# catchability moves the stock and re-matching the stock moves the catch, so
# the two are repeated until both hold.
for (pass in 1:8) {
  if (herring_target > 0) params <- solve_steady(match_herring_catch(params))
  params <- solve_steady(quietly(matchBiomasses(params)))
}

biomass_before <- getBiomass(params)
new_levels <- c("Sand Lance" = sand_lance_reproduction_level)
if (herring_target > 0) new_levels["Atlantic Herring"] <- herring_reproduction_level
params <- setBevertonHolt(params, reproduction_level = new_levels)
stopifnot(max(abs(getBiomass(params) / biomass_before - 1)) < 1e-8)

species_params <- params@species_params
biomass_ratio <- getBiomass(params) / species_params$biomass_observed
yield_modelled <- colSums(getYieldGear(params)[fishing_gears, ])
has_yield <- !is.na(species_params$yield_observed) & species_params$yield_observed > 0
yield_ratio <- (yield_modelled / species_params$yield_observed)[has_yield]
cat(sprintf("\nbase with raised forage fish: %.0f kt sand lance | %.0f kt herring | all species %.2f g/m2\n",
            getBiomass(params)["Sand Lance"] / 1e9, getBiomass(params)["Atlantic Herring"] / 1e9,
            sum(getBiomass(params)) / region_area_m2))
cat(sprintf("biomass within 2x: %d/19 | yield within 2x: %d/%d | species with erepro > 1: %d\n",
            sum(biomass_ratio > 0.5 & biomass_ratio < 2),
            sum(yield_ratio > 0.5 & yield_ratio < 2), length(yield_ratio),
            sum(species_params$erepro > 1)))
cat(sprintf("sand lance erepro %.3f | herring erepro %.3f, catch %.1f kt (observed %.1f), fishing mortality on adults %.2f a year\n",
            species_params["Sand Lance", "erepro"], species_params["Atlantic Herring", "erepro"],
            yield_modelled["Atlantic Herring"] / 1e9,
            species_params["Atlantic Herring", "yield_observed"] / 1e9,
            max(getFMort(params)["Atlantic Herring", ])))
# Left alone for 100 years the base should stay where it is.
left_alone <- quietly(project(params, t_max = 100, progress_bar = FALSE))
biomass_left_alone <- getBiomass(left_alone)
drift <- biomass_left_alone[nrow(biomass_left_alone), ] / biomass_left_alone[1, ]
cat(sprintf("largest change when left alone for 100 years: %s x%.3f\n",
            names(drift)[which.max(abs(log(drift)))], drift[which.max(abs(log(drift)))]))
saveRDS(params, file.path(output_folder, "base_params_raised_forage_fish.RDS"))
