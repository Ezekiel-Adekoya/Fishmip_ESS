# Figures for the dynamic seal run from script 9, against the data and against
# the seal gear model as fitted in script 10.
#
# Set from outside:
#   RUN_FILE          the run to draw (default dynamic_seal_run.rds)
#   EARLIER_RUN_FILE  optional: a second dynamic seal run drawn in grey for
#                     comparison, for example the run without extra cod mortality
#   EARLIER_RUN_LABEL the legend name of that second run
#   FIGURE_TAG        added to the figure file names, so runs do not overwrite
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
region_area_m2 <- 146360555008
run <- readRDS(file.path(Sys.getenv("SC"), Sys.getenv("RUN_FILE", "dynamic_seal_run.rds")))
figure_tag <- Sys.getenv("FIGURE_TAG", "")
earlier_run_label <- Sys.getenv("EARLIER_RUN_LABEL", "dynamic seal, earlier run")
earlier_run <- if (nzchar(Sys.getenv("EARLIER_RUN_FILE"))) {
  readRDS(file.path(Sys.getenv("SC"), Sys.getenv("EARLIER_RUN_FILE")))
}
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
run_colours[earlier_run_label] <- "grey55"

long_biomass <- function(biomass_matrix, run_name) {
  as.data.frame(biomass_matrix / region_area_m2) |>
    mutate(Year = as.integer(rownames(biomass_matrix)), run = run_name) |>
    pivot_longer(-c(Year, run), names_to = "Species", values_to = "value") |>
    filter(Year %in% plot_years)
}
modelled <- bind_rows(long_biomass(run$biomass, "dynamic seal"),
                      long_biomass(gear_fit$biomass, "seal gear, fitted (script 10)")) |>
  filter(Species != "Grey Seal")
if (!is.null(earlier_run)) {
  modelled <- bind_rows(modelled, long_biomass(earlier_run$biomass, earlier_run_label)) |>
    filter(Species != "Grey Seal")
}

species_plot <- ggplot(modelled, aes(Year, value)) +
  geom_line(aes(colour = run), linewidth = 0.6) +
  geom_point(data = filter(observed, Species %in% unique(modelled$Species)), size = 0.5) +
  facet_wrap(~Species, scales = "free_y", ncol = 5) + scale_y_log10() +
  scale_colour_manual(values = run_colours) +
  labs(title = paste0("Biomass by species: ", label),
       subtitle = "Points: survey (herring, sand lance and capelin points are raw trawl figures, far below the true stock).",
       y = "g per m2 (log scale)", x = NULL, colour = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(paste0("figures_recalibrated/dynseal_bowen_biomass_by_species", figure_tag, ".png"),
       species_plot,
       width = 14, height = 9, dpi = 150)

years <- as.character(plot_years)
# The model lines of the summary figure for one run.
model_panels <- function(one_run, series_name) {
  bind_rows(
    data.frame(Year = plot_years, panel = "Seals (thousands)",
               value = one_run$seal_numbers[years]),
    data.frame(Year = plot_years, panel = "Sand lance stock (kt)",
               value = one_run$biomass[years, "Sand Lance"] / 1e9),
    data.frame(Year = plot_years, panel = "Fish eaten per seal (t a year; target 1.73)",
               value = rowSums(one_run$eaten[years, ]) / one_run$seal_numbers[years]),
    data.frame(Year = plot_years, panel = "Sand lance share of fish eaten (%; Bowen 81)",
               value = 100 * one_run$eaten[years, "Sand Lance"] / rowSums(one_run$eaten[years, ])),
    data.frame(Year = plot_years, panel = "Cod share of fish eaten (%; Bowen 11)",
               value = 100 * one_run$eaten[years, "Atlantic Cod"] / rowSums(one_run$eaten[years, ])),
    data.frame(Year = plot_years, panel = "Cod biomass (g per m2)",
               value = one_run$biomass[years, "Atlantic Cod"] / region_area_m2)) |>
    mutate(series = series_name)
}
data_panels <- bind_rows(
  data.frame(Year = seal_counts_table$year, panel = "Seals (thousands)",
             value = seal_counts_table$scotian_mean),
  filter(observed, Species == "Atlantic Cod") |>
    transmute(Year, panel = "Cod biomass (g per m2)", value)) |>
  filter(Year %in% plot_years)
model_lines <- model_panels(run, "dynamic seal")
if (!is.null(earlier_run)) {
  model_lines <- bind_rows(model_panels(earlier_run, earlier_run_label), model_lines)
}
summary_plot <- ggplot(model_lines, aes(Year, value)) +
  geom_line(aes(colour = series), linewidth = 0.8) +
  geom_point(data = data_panels, size = 1) +
  scale_colour_manual(values = run_colours) +
  facet_wrap(~panel, scales = "free_y", ncol = 3) +
  labs(title = paste0("Seals, their food and cod: ", label),
       subtitle = "Lines: model. Points: seal counts and cod survey.", y = NULL, x = NULL,
       colour = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(paste0("figures_recalibrated/dynseal_bowen_seals_food_cod", figure_tag, ".png"),
       summary_plot, width = 12, height = 7, dpi = 150)
