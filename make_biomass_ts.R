# Build Ess_Fishmip_data/biomass_ts.RDS -- the observed survey biomass time
# series that 05.Thermizer_with_Runs.Rmd reads at its "Survey data" step and
# then compares the FishMIP biomass outputs (tcb, tpb, tdb) against.
#
# WHY THIS SCRIPT EXISTS
# 05 builds the same survey table internally earlier on, but renames its value
# column to Biomass_Obs, while everything after the "Survey data" step sums a
# column called Biomass. So the readRDS is filling a real gap rather than being
# redundant. Copies of biomass_ts.RDS in the sibling ESS repos cannot be used:
# they hold 13 species under raw upper-case survey names (AMERICAN PLAICE,
# HERRING(ATLANTIC), ...) plus Greenland Halibut, which is not in the model, so
# they fail 05's own check that the three pelagic species are present.
#
# This script reproduces 05's construction step for step from the raw survey
# file, so the species names and the numbers agree with what 05 builds itself.
#
# Run from the project root:  Rscript make_biomass_ts.R

suppressMessages({
  library(dplyr)
  library(tidyr)
})

raw_survey_file  <- "Ess_Thermizer_data/LF_Biomass_Abund.csv"
output_file      <- "Ess_Fishmip_data/biomass_ts.RDS"

# The 4VsW strata and the year span 05 uses. Strata codes outside 443-466 are a
# different survey footprint and would not be comparable across years.
first_stratum    <- 443
last_stratum     <- 466
first_survey_year <- 1970
last_survey_year  <- 2024

# Raw DFO survey common names mapped to the model's species names. The raw file
# uses its own capitalised codes, so this mapping is what makes the output
# comparable with species_params(params)$species.
survey_name_to_model_name <- c(
  "COD(ATLANTIC)"             = "Atlantic Cod",
  "HADDOCK"                   = "Haddock",
  "WHITE HAKE"                = "White Hake",
  "SILVER HAKE"               = "Silver Hake",
  "POLLOCK"                   = "Pollock",
  "REDFISH UNSEPARATED"       = "Redfish Unseparated",
  "AMERICAN PLAICE"           = "American Plaice",
  "YELLOWTAIL FLOUNDER"       = "Yellowtail Flounder",
  "SPINY DOGFISH"             = "Spiny Dogfish",
  "HALIBUT(ATLANTIC)"         = "Atlantic Halibut",
  "SQUIRREL OR RED HAKE"      = "Red Hake",
  "MONKFISH,GOOSEFISH,ANGLER" = "Monkfish",
  "NORTHERN SAND LANCE"       = "Sand Lance",
  "WITCH FLOUNDER"            = "Witch Flounder",
  "HERRING(ATLANTIC)"         = "Atlantic Herring",
  "THORNY SKATE"              = "Thorny Skate",
  "WINTER SKATE"              = "Winter Skate",
  "LONGHORN SCULPIN"          = "Longhorn Sculpin",
  "CAPELIN"                   = "Capelin")

model_species_in_order <- c(
  "Atlantic Cod", "Haddock", "White Hake", "Red Hake", "Silver Hake",
  "Pollock", "Redfish Unseparated", "American Plaice", "Yellowtail Flounder",
  "Witch Flounder", "Atlantic Herring", "Capelin", "Thorny Skate",
  "Winter Skate", "Spiny Dogfish", "Monkfish", "Sand Lance",
  "Longhorn Sculpin", "Atlantic Halibut")

# Survey catchability multipliers for the three small pelagics, copied from 05
# so the two agree. The bottom-trawl survey misses these species at all sizes
# (they are off-bottom, or bury in sediment), so the raw survey biomass is not a
# usable estimate of standing stock for them.
#
# CAVEAT, recorded deliberately: these multipliers were tuned to let the seal
# predation terms balance, not derived as independent catchability estimates.
# Including them here means the FishMIP pelagic biomass output (tpb) is compared
# against survey values inflated by up to 3000x. That is the intended behaviour
# as of this writing -- it keeps the comparison on the same footing as the
# biomass the model was calibrated to match -- but it is a modelling choice, not
# an observation. Set apply_forage_multipliers to FALSE for raw survey values.
apply_forage_multipliers <- TRUE
forage_catchability_mult <- c(`Atlantic Herring` = 50,
                              `Sand Lance`       = 3000,
                              `Capelin`          = 200)

stopifnot(file.exists(raw_survey_file))
raw_survey <- read.csv(raw_survey_file)

# Sum biomass at length into an annual total per species, converting the raw
# kilograms to grams, then widen so each species is a column.
survey_biomass_wide <- raw_survey %>%
  mutate(stratum_number = suppressWarnings(as.numeric(as.character(STRAT)))) %>%
  filter(!is.na(stratum_number),
         stratum_number >= first_stratum, stratum_number <= last_stratum,
         YEAR >= first_survey_year, YEAR <= last_survey_year) %>%
  mutate(model_species_name = survey_name_to_model_name[COMMON_NAME]) %>%
  filter(!is.na(model_species_name)) %>%
  group_by(Year = YEAR, model_species_name) %>%
  summarise(biomass_g = sum(Biomass_Strata_KG, na.rm = TRUE) * 1e3,
            .groups = "drop") %>%
  pivot_wider(names_from = model_species_name, values_from = biomass_g,
              values_fill = 0) %>%
  select(Year, any_of(model_species_in_order))

# A species absent from a year's survey means zero caught, not unknown.
survey_biomass_wide[is.na(survey_biomass_wide)] <- 0

if (apply_forage_multipliers) {
  for (species_name in names(forage_catchability_mult)) {
    stopifnot(species_name %in% names(survey_biomass_wide))
    survey_biomass_wide[[species_name]] <-
      survey_biomass_wide[[species_name]] *
      forage_catchability_mult[[species_name]]
  }
}

# Long format with a column called Biomass, which is the name 05 sums over.
biomass_ts <- survey_biomass_wide %>%
  pivot_longer(cols = -Year, names_to = "Species", values_to = "Biomass") %>%
  arrange(Species, Year) %>%
  as.data.frame()

# The checks 05 itself makes on this object, run here so a bad file cannot
# reach the model run. The three pelagic names must be present verbatim,
# because 05 selects pelagic species by exact name match.
pelagic_species <- c("Atlantic Herring", "Capelin", "Sand Lance")
stopifnot(
  all(c("Year", "Species", "Biomass") %in% names(biomass_ts)),
  all(pelagic_species %in% unique(biomass_ts$Species)),
  all(unique(biomass_ts$Species) %in% model_species_in_order),
  !anyNA(biomass_ts$Biomass),
  all(biomass_ts$Biomass >= 0)
)

saveRDS(biomass_ts, output_file)

cat(sprintf("wrote %s\n  %d rows, %d species, years %d-%d\n",
            output_file, nrow(biomass_ts),
            length(unique(biomass_ts$Species)),
            min(biomass_ts$Year), max(biomass_ts$Year)))
cat(sprintf("  forage multipliers applied: %s\n",
            if (apply_forage_multipliers)
              paste(sprintf("%s x%g", names(forage_catchability_mult),
                            forage_catchability_mult), collapse = ", ")
            else "none"))
