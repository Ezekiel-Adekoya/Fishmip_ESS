# 7. A fish base with enough sand lance to feed the seals
#
# The pipeline model is calibrated to 8.2 kt of sand lance, the raw trawl
# survey figure. Trawls barely catch sand lance, and Bowen & Harrison found it
# to be four fifths of what grey seals eat. At 1.73 t of fish per seal a year,
# the 15 thousand seals of 1978-1981 need about 21 kt of sand lance a year and
# the 275 thousand of 2010 about 385 kt. A dynamic seal cannot be built on
# 8 kt.
#
# This script raises the sand lance target in steps, recalibrating the whole
# community each time exactly as script 06 does, and then loosens the limit on
# sand lance recruitment so the stock can respond when its predators change.
#
# The target level is a working assumption, not a measurement: no reliable
# biomass estimate for sand lance on the Scotian Shelf was found. It is set
# with SAND_LANCE_TARGET_KT (default 400).
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
# Reproduction level for sand lance after the recalibration. The pipeline model
# has it at 0.99, which pins recruitment at its ceiling: the stock is then very
# stable but can never grow. A lower level lets recruitment rise (0.5 lets it
# double) but also lets it fall when the spawning stock falls, and at 0.5 the
# stock collapsed during the spin-up. Set with SAND_LANCE_REPRODUCTION_LEVEL.
sand_lance_reproduction_level <- as.numeric(
  Sys.getenv("SAND_LANCE_REPRODUCTION_LEVEL", "0.9"))

params <- readRDS("Ess_Fishmip_data/params_GFDL_steady.RDS")
solve_steady <- function(p) quietly(steady(p, t_max = 300, progress_bar = FALSE,
                                           preserve = "erepro"))

target_now <- params@species_params["Sand Lance", "biomass_observed"]
while (target_now < sand_lance_target) {
  target_now <- min(sand_lance_target, target_now * 1.5)
  params@species_params["Sand Lance", "biomass_observed"] <- target_now
  given <- params@given_species_params
  given$biomass_observed[given$species == "Sand Lance"] <- target_now
  params@given_species_params <- given
  for (pass in 1:3) params <- solve_steady(quietly(matchBiomasses(params)))
  cat(sprintf("sand lance target %4.0f kt: model %4.0f kt | all species %.2f g/m2 | worst erepro %.3f\n",
              target_now / 1e9, getBiomass(params)["Sand Lance"] / 1e9,
              sum(getBiomass(params)) / region_area_m2, max(params@species_params$erepro)))
}
for (pass in 1:3) params <- solve_steady(quietly(matchBiomasses(params)))

biomass_before <- getBiomass(params)
params <- setBevertonHolt(params, reproduction_level =
                            c("Sand Lance" = sand_lance_reproduction_level))
stopifnot(max(abs(getBiomass(params) / biomass_before - 1)) < 1e-8)

species_params <- params@species_params
biomass_ratio <- getBiomass(params) / species_params$biomass_observed
yield_modelled <- colSums(getYieldGear(params)[setdiff(names(params@initial_effort), "Seals"), ])
has_yield <- !is.na(species_params$yield_observed) & species_params$yield_observed > 0
yield_ratio <- (yield_modelled / species_params$yield_observed)[has_yield]
cat(sprintf("\nbase with raised sand lance: %.0f kt sand lance | all species %.2f g/m2\n",
            getBiomass(params)["Sand Lance"] / 1e9, sum(getBiomass(params)) / region_area_m2))
cat(sprintf("biomass within 2x: %d/19 | yield within 2x: %d/%d | species with erepro > 1: %d | sand lance erepro %.3f\n",
            sum(biomass_ratio > 0.5 & biomass_ratio < 2),
            sum(yield_ratio > 0.5 & yield_ratio < 2), length(yield_ratio),
            sum(species_params$erepro > 1), species_params["Sand Lance", "erepro"]))
saveRDS(params, file.path(output_folder, "base_params_raised_sand_lance.RDS"))
