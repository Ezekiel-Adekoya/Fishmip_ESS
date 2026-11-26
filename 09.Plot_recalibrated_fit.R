# Plot the recalibrated dynamic run (from 07) against the observations.
#
# Three figures are written to figures_recalibrated/:
#   biomass_by_species.png  modelled and survey biomass, one panel per species
#   catch_by_species.png    modelled fishery catch and NAFO landings
#   totals.png              community totals for the species without multipliers
#
# Catch is fishery catch only: seals are a gear in this model, so their
# consumption is left out. The grey band marks 2005-2010, the years the FishMIP
# protocol holds back from calibration. The dashed line marks 1978-1981, the
# window the model was calibrated on.
#
# Run from the project root after 07:  Rscript 09.Plot_recalibrated_fit.R

suppressMessages({
  library(mizer); library(therMizer); library(dplyr); library(tidyr)
  library(ggplot2)
})
source("Ess_Fishmip_data/Useful_functions.R")

region_area_m2 <- 146360555008
multiplied_forage_species <- c("Atlantic Herring", "Sand Lance", "Capelin")
output_folder <- "figures_recalibrated"
dir.create(output_folder, showWarnings = FALSE)

sim <- readRDS("Ess_Fishmip_data/sims_GFDL_recalibrated.RDS")

to_long <- function(year_by_species, value_name) {
  as.data.frame(as.table(year_by_species)) |>
    setNames(c("Year", "Species", value_name)) |>
    mutate(Year = as.integer(as.character(Year)), Species = as.character(Species))
}

fishing_gears <- setdiff(dimnames(getYieldGear(sim))[[2]], "Seals")
modelled_catch <- apply(getYieldGear(sim)[, fishing_gears, , drop = FALSE],
                        c(1, 3), sum)

modelled <- inner_join(to_long(getBiomass(sim), "biomass"),
                       to_long(modelled_catch, "catch"),
                       by = c("Year", "Species")) |>
  filter(Year >= 1961)

observed_biomass <- readRDS("Ess_Fishmip_data/biomass_ts.RDS") |>
  transmute(Year = as.integer(Year), Species = as.character(Species),
            biomass = Biomass) |>
  filter(Year <= 2010, biomass > 0)

observed_catch <- readRDS("Ess_Fishmip_data/new_NAFO_4VW_Catches.RDS") |>
  pivot_longer(-Year, names_to = "Species", values_to = "catch") |>
  filter(!is.na(catch), catch > 0, Year <= 2010)

# Mark the species whose survey biomass carries a catchability multiplier, so a
# reader does not take those panels as a comparison against raw data.
label_species <- function(species) {
  ifelse(species %in% multiplied_forage_species,
         paste0(species, " (survey x multiplier)"), species)
}

panel_plot <- function(modelled_values, observed_values, value_column, y_title,
                       plot_title) {
  modelled_values$panel <- label_species(modelled_values$Species)
  observed_values$panel <- label_species(observed_values$Species)
  ggplot() +
    annotate("rect", xmin = 2004.5, xmax = 2010.5, ymin = 0, ymax = Inf,
             fill = "grey80", alpha = 0.5) +
    geom_vline(xintercept = c(1978, 1981), linetype = "dashed",
               linewidth = 0.3, colour = "grey40") +
    geom_line(data = modelled_values,
              aes(Year, .data[[value_column]] / region_area_m2,
                  colour = "Model"), linewidth = 0.7) +
    geom_point(data = observed_values,
               aes(Year, .data[[value_column]] / region_area_m2,
                   colour = "Observed"), size = 0.9) +
    facet_wrap(~ panel, scales = "free_y", ncol = 4) +
    scale_y_log10() +
    scale_colour_manual(values = c(Model = "#1b6ca8", Observed = "black"),
                        name = NULL) +
    labs(title = plot_title, x = NULL, y = y_title,
         caption = paste("Dashed lines: 1978-1981 calibration window.",
                         "Grey band: 2005-2010, held out by the protocol.",
                         "Log scale.")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", strip.text = element_text(size = 8))
}

ggsave(file.path(output_folder, "biomass_by_species.png"),
       panel_plot(modelled, observed_biomass, "biomass",
                  expression("Biomass (g m"^-2*")"),
                  "Recalibrated model on the COBALT resource: biomass"),
       width = 13, height = 10, dpi = 150)

species_with_catch <- unique(observed_catch$Species)
ggsave(file.path(output_folder, "catch_by_species.png"),
       panel_plot(filter(modelled, Species %in% species_with_catch, catch > 0),
                  observed_catch, "catch",
                  expression("Catch (g m"^-2*" year"^-1*")"),
                  "Recalibrated model on the COBALT resource: fishery catch"),
       width = 13, height = 10, dpi = 150)

# Community totals over the species with no multiplier, summed only over the
# species-years that have an observation so the two lines cover the same things.
total_of <- function(observed_values, value_column) {
  inner_join(modelled, observed_values, by = c("Year", "Species"),
             suffix = c("_model", "_observed")) |>
    filter(!Species %in% multiplied_forage_species) |>
    group_by(Year) |>
    summarise(Model    = sum(.data[[paste0(value_column, "_model")]]),
              Observed = sum(.data[[paste0(value_column, "_observed")]]),
              .groups = "drop") |>
    pivot_longer(-Year, names_to = "source", values_to = "value") |>
    mutate(quantity = value_column)
}
totals <- bind_rows(total_of(observed_biomass, "biomass"),
                    total_of(observed_catch, "catch")) |>
  mutate(quantity = ifelse(quantity == "biomass",
                           "Total biomass (16 species)",
                           "Total fishery catch (16 species)"))

ggsave(file.path(output_folder, "totals.png"),
       ggplot(totals, aes(Year, value / region_area_m2, colour = source)) +
         annotate("rect", xmin = 2004.5, xmax = 2010.5, ymin = 0, ymax = Inf,
                  fill = "grey80", alpha = 0.5) +
         geom_vline(xintercept = c(1978, 1981), linetype = "dashed",
                    linewidth = 0.3, colour = "grey40") +
         geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
         facet_wrap(~ quantity, scales = "free_y", ncol = 1) +
         scale_colour_manual(values = c(Model = "#1b6ca8", Observed = "black"),
                             name = NULL) +
         labs(title = "Recalibrated model on the COBALT resource: community totals",
              x = NULL, y = expression("g m"^-2*"  (catch: per year)"),
              caption = paste("Species without survey multipliers, summed over",
                              "species-years that have an observation.")) +
         theme_bw(base_size = 11) + theme(legend.position = "top"),
       width = 9, height = 7.5, dpi = 150)

cat("wrote", paste(list.files(output_folder, full.names = TRUE), collapse = ", "), "\n")
