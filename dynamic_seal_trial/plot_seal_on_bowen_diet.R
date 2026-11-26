# Figures for the dynamic seal run from script 9, against the data and against
# the seal gear model as fitted in script 10.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
region_area_m2 <- 146360555008
run <- readRDS(file.path(Sys.getenv("SC"), "dynamic_seal_run.rds"))
gear_fit <- readRDS("Ess_Fishmip_data/sims_GFDL_fitted.RDS")
label <- Sys.getenv("RUN_LABEL", "dynamic seal")
plot_years <- 1961:2010
seal_counts_table <- readRDS(
  "~/Desktop/Ez June 2026/ESS_Restructured/Ess_data/seal_abundance_digitised.RDS")
observed <- readRDS("Ess_Fishmip_data/biomass_ts.RDS") |>
  transmute(Year = as.integer(Year), Species = as.character(Species),
            value = Biomass / region_area_m2) |>
  filter(value > 0, Year <= 2010)
run_colours <- c("dynamic seal" = "#b2182b", "seal gear, fitted (script 10)" = "#4a90d9")

long_biomass <- function(biomass_matrix, run_name) {
  as.data.frame(biomass_matrix / region_area_m2) |>
    mutate(Year = as.integer(rownames(biomass_matrix)), run = run_name) |>
    pivot_longer(-c(Year, run), names_to = "Species", values_to = "value") |>
    filter(Year %in% plot_years)
}
modelled <- bind_rows(long_biomass(run$biomass, "dynamic seal"),
                      long_biomass(gear_fit$biomass, "seal gear, fitted (script 10)")) |>
  filter(Species != "Grey Seal")

species_plot <- ggplot(modelled, aes(Year, value)) +
  geom_line(aes(colour = run), linewidth = 0.6) +
  geom_point(data = filter(observed, Species %in% unique(modelled$Species)), size = 0.5) +
  facet_wrap(~Species, scales = "free_y", ncol = 5) + scale_y_log10() +
  scale_colour_manual(values = run_colours) +
  labs(title = paste0("Biomass by species: ", label),
       subtitle = "Points: survey (herring, sand lance and capelin points are raw trawl figures, far below the true stock).",
       y = "g per m2 (log scale)", x = NULL, colour = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave("figures_recalibrated/dynseal_bowen_biomass_by_species.png", species_plot,
       width = 14, height = 9, dpi = 150)

years <- as.character(plot_years)
summary_panels <- bind_rows(
  data.frame(Year = plot_years, panel = "Seals (thousands)", series = "model",
             value = run$seal_numbers[years]),
  data.frame(Year = seal_counts_table$year, panel = "Seals (thousands)", series = "data",
             value = seal_counts_table$scotian_mean),
  data.frame(Year = plot_years, panel = "Sand lance stock (kt)", series = "model",
             value = run$biomass[years, "Sand Lance"] / 1e9),
  data.frame(Year = plot_years, panel = "Fish eaten per seal (t a year; target 1.73)",
             series = "model", value = rowSums(run$eaten[years, ]) / run$seal_numbers[years]),
  data.frame(Year = plot_years, panel = "Sand lance share of fish eaten (%; Bowen 81)",
             series = "model",
             value = 100 * run$eaten[years, "Sand Lance"] / rowSums(run$eaten[years, ])),
  data.frame(Year = plot_years, panel = "Cod share of fish eaten (%; Bowen 11)",
             series = "model",
             value = 100 * run$eaten[years, "Atlantic Cod"] / rowSums(run$eaten[years, ])),
  data.frame(Year = plot_years, panel = "Cod biomass (g per m2)", series = "model",
             value = run$biomass[years, "Atlantic Cod"] / region_area_m2),
  filter(observed, Species == "Atlantic Cod") |>
    transmute(Year, panel = "Cod biomass (g per m2)", series = "data", value)) |>
  filter(Year %in% plot_years)
summary_plot <- ggplot(summary_panels, aes(Year, value)) +
  geom_line(data = filter(summary_panels, series == "model"), colour = "#b2182b",
            linewidth = 0.8) +
  geom_point(data = filter(summary_panels, series == "data"), size = 1) +
  facet_wrap(~panel, scales = "free_y", ncol = 3) +
  labs(title = paste0("Seals, their food and cod: ", label),
       subtitle = "Line: model. Points: seal counts and cod survey.", y = NULL, x = NULL) +
  theme_bw(base_size = 10)
ggsave("figures_recalibrated/dynseal_bowen_seals_food_cod.png", summary_plot,
       width = 12, height = 7, dpi = 150)
