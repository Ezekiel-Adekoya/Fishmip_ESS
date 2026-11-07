# Fishmip outputs as requested

A five-stage R Markdown pipeline that builds temperature, plankton (resource) and depth-use forcing for a mizer/therMizer model of the ESS (19 fish species plus seals), runs it **with and without fishing**, then turns the saved simulations into FishMIP output variables (1961-2010, GFDL-MOM6-COBALT2 obsclim). Author: Ezekiel (resource script adapted from Phoebe Woodworth-Jefcoats, NOAA).

## Pipeline overview

The first three scripts are independent of each other; the resource script needs the plankton script's output; stage E needs everything. Stage E is run as two documents: the dynamics/ensemble Rmd (E1) and the FishMIP export Rmd (E2).

| Step | Script | What it does | Key output(s) |
|---|---|---|---|
| A | Bias correction of GFDL thetao | Bias-corrects temperature, averages to two layers, adds spin-up | `Ess_Fishmip_data/Temp_by_depth_spin.csv` |
| B | GFDL obsclim plankton processing | Annual regional plankton totals, with spin-up | `Ess_Fishmip_data/GFDL_obsclim_{diat,diaz,pico,meso,micro}_1861_2010.csv` |
| C | Realms set-up | Depth-use (time-in-layer) and exposure arrays | `GFDL_ESS_vertical_migration_array.Rdata`, `GFDL_ESS_exposure_array.csv` |
| D | Prep plankton for therMizer | Builds the resource spectrum per year from B | `GFDL_resource_spectra.dat` (plus `_slope`, `_intercept`) |
| E1 | Exploring thermal effects using therMizer (dynamics) | Loads all forcing, builds `params_dyn`, runs the randomised-web ensemble with and without fishing | `sims_with_fishing`, `sims_with_no_fishing` (cached RDS) |
| E2 | FishMIP outputs | Turns both ensembles into FishMIP variables, plots them against observations | CSVs in `Fishmip_output/with_fishing/` and `Fishmip_output/without_fishing/` |

All scripts use the same base model, `params_stage1_steady.RDS` (Stage 2 is not forward-stable), and must agree on **two realms in a fixed order: `r0to50`, `r50to150`**. therMizer matches temperature columns to realms by position, so a mismatch will not error, it will silently pair the wrong temperature with the wrong layer.

---

## A. Bias correction of GFDL-MOM6-COBALT2 temperature

**Method:** (1) load per-depth GFDL and WOA monthly climatologies (1981-2010) and take cos-latitude-weighted regional means; (2) bias = WOA minus GFDL for each of 9 depths x 12 months; (3) read the raw monthly GFDL series (1961-2010) from the zarr store, matching depths to WOA (within 0.5 m); (4) add the bias (same offset every year, so trend is untouched); (5) annual means, requiring all 12 months; (6) average to `r0to50` (2.5-32.5 m) and `r50to150` (51.25-156.25 m); (7) spin-up: repeat the first `n_cycle` years back to 1861, giving 150 rows.

**Inputs:** `Gfdl_temp forcing/isimip_3a/<depth>_.../` (one GFDL and one WOA parquet per depth), and the GFDL thetao zarr store.
**Outputs:** `Gfdl_temp forcing/gfdl_temp_Annual_temp_bias-corrected_all-depths_1961-2010.csv`, `Temp_by depth.csv` (project root), `Ess_Fishmip_data/Temp_by_depth_spin.csv` (used downstream), and diagnostic plots.

**Known limitation: level, not trend.** The WOA reference is a climatology, so GFDL's own warming rate is inherited unchanged (about +0.26 C/decade at 0-50 m and +0.21 at 50-150 m, 1961-2024, with the surface-deep gap widening). Observations (Zhai et al. 2025) show stronger bottom than surface warming and a regime shift around 2012. Because most species are weighted to the deep layer, temperature effects are likely biased downward. Fixing this needs a time-resolved reference such as GLORYS12V1, not different WOA offsets. July warming confined to the top 10 m is also diluted by the 0-50 m average.

**Check:**
- `library(ggplot2)` is commented out, so the first plot chunk fails until a later chunk loads it.
- The comment says the spin-up repeats 1961-1970, but `n_cycle <- 5` repeats 1961-1965.
- Layer means are simple (not thickness-weighted) averages.
- `Temp_by depth.csv` (with a space) goes to the project root while the spin-up file goes to `Ess_Fishmip_data/`.
- The Overview mentions `Ocean_temp_by_Depth_1961_2100.csv` and an ISIMIP3b splice, neither of which is made here. `jsonlite` is used without being loaded.

---

## B. GFDL obsclim plankton processing

Converts monthly, vertically integrated plankton carbon (mol C m-2, ISIMIP `-vint` files) for five groups (diatoms, diazotrophs, picophytoplankton, mesozooplankton, microzooplankton) into an **annual regional total in g wet weight**: unweighted mean across grid cells, then annual mean, then x 12.011 x 10 x region area (146,360,555,008 m2). Then builds a 1861-1960 spin-up by repeating the first `n_cycle` years (aligned backward so 1961 restarts the cycle), plots it, and writes one `V1` (date) / `V2` (value) CSV per group for script D.

**Inputs:** `Gfdl_Plankton Forcing/Plankton_data/gfdl-mom6-cobalt2_obsclim_<var>-vint_15arcmin_eastern-scotian-shelf_monthly_1961_2010.csv`
**Outputs:** annual CSVs in `Gfdl_Plankton Forcing/plankton/`, `GFDL_obsclim_<group>_1861_2010.csv`, and `plankton_spinup_1861_2010.png` in `Ess_Fishmip_data/`.

**Check:**
- The spin-up chunk needs a `plankton_all` data frame (Year, val, variable) that is not built anywhere in this file. The group labels there (`Diatoms`, ...) also differ from the labels written above (`Diat`, ...). A combining step appears to be missing.
- The 10-year spin-up cycle here differs from the 5-year cycle in script A, so temperature and plankton spin-ups are out of phase.
- The x10 carbon-to-wet-weight factor was set for phytoplankton; script D divides it back out and applies group-specific factors.

---

## C. Realms set-up (vertical migration and exposure arrays)

`vertical_migration_array[realm, species, size]` holds the **proportion of time** each species spends in each layer; it must sum to 1 across realms, because therMizer evaluates the thermal curve per layer and sums the weighted results. Weights are size-constant. Examples: most gadoids and Spiny Dogfish 0.2/0.8, flatfish, skates, redfish and monkfish 0/1 (deep), Herring and Capelin 1/0 (shallow), Sand Lance 0.9/0.1. `exposure_array[realm, species]` is 1 wherever the weight is above zero. The script asserts every species is covered and every row sums to 1.

**Output:** `GFDL_ESS_vertical_migration_array.Rdata`, `GFDL_ESS_exposure_array.csv`.
**Note:** the old 18-realm (Hawaii longline) scheme is kept commented out for reference only. The header text still describes the older ISIMIP3b workflow and lists file names without the `GFDL_` prefix.

---

## D. Prep plankton for therMizer (resource spectra)

Builds, for each year 1861-2010, the background resource on the model's `w_full` size grid.
1. Convert each group's wet weight to carbon-corrected wet weight: divide out script B's uniform x10, then apply group-specific factors (pico 4.74, large 26.5, zmicro 10, zmeso 10; Menden-Deuer & Lessard 2000 for phytoplankton). Size classes (ESD): pico 0.2-10 um, large 10-200 um, zmicro 2-200 um, zmeso 200-20,000 um.
2. Spread each group's biomass over its own size range as a power law `c_i w^b` whose integral equals the group biomass exactly, then **sum the groups** (overlapping classes add, no mid-point size needed).
3. Find the within-class exponent `b` self-consistently (build, fit aggregate slope, repeat; about 5 passes, b near -1.99).
4. Continue the power law above the 4.19 g mesozooplankton limit at the same slope. This stands in for unmodelled prey (benthos, small non-modelled fish) that mizer's Stage-1 resource implicitly represented; therMizer's `w_pp_cutoff` stops it at 100 g. This deliberately breaks biomass conservation above 4.19 g.
5. Write `log10(numbers per bin)`, which is what `plankton_forcing()` expects (empty bins are `-Inf`, meaning zero resource). No slope or intercept scaling is applied.

The chunk ends with verification: slope vs the calibrated resource, level ratio, per-group biomass conservation (< 1e-8), and the size of the extended band.

**Inputs:** `GFDL_obsclim_<group>_1861_2010.csv` (from B), `params_stage1_steady.RDS`.
**Outputs:** `GFDL_resource_spectra.dat` (150 years x 211 bins, row names are years), `GFDL_resource_slope.dat`, `GFDL_resource_intercept.dat`.

**Check:**
- `dplyr` is used (`%>%`, `mutate`) but not loaded in this file; it depends on the session.
- Comments say the committed `.dat` has 197 columns (`stage2_factual`), but the code uses Stage 1 (211 columns), which is what script E asserts. Regenerate before running E.
- `save_dir` is built with `paste0(here(), "Ess_Fishmip_data/")` and has no `/` separator (unused, but wrong). One `stopifnot` guard is duplicated.
- The mizer-vs-ESM level is not expected to agree; the slope comparison is the informative one.

---

## E1. Dynamics and ensemble (exploring thermal effects using therMizer)

**Set-up:**
- Loads `Useful_functions.R` first, before anything that triggers `setParams()`, because `gaussian_sizeweight` is the `sel_func` for the Seals gear rows. The Rmd defines a fallback only if the file does not supply it.
- Loads `params_stage1_steady.RDS`; attaches `temp_min` / `temp_max` from `species_params_with_temp.csv` (FishBase, with an order assertion); overrides Spiny Dogfish (0-12 C) and Thorny Skate (-1.3-14 C), because FishBase/AquaMaps ranges are occurrence envelopes rather than physiological limits (Capelin and Red Hake overrides are commented out).
- Loads the resource (D) and temperature (A) as **full 1861-2010 series**, plus the arrays from C, and builds `params_dyn` with `upgradeTherParams()`.
- `reference_years <- 1978:1981` is used only for sanity checks (resource level against Paper 1, and that every species has at least one realm inside its thermal range at the reference temperature). It no longer drives the dynamic runs. No steady state is solved here, because one does not exist under time-varying forcing; calibration belongs in the separate calibration script.

**Scenarios:**

| Scenario | Effort | Object |
|---|---|---|
| With fishing (factual) | `effort_factual` from `Effort_spin_1861_2010.csv` | `sims_with_fishing` |
| Without fishing (counterfactual) | `effort_no_fishing`: same matrix with every gear, including seals, set to 0 | `sims_with_no_fishing` |

**Ensemble:** `run_paired_sims()` (in `Useful_functions.R`) randomises the interaction matrix and seal catchability once per member (lognormal noise, SD `sd_web`; zeros stay zero) and projects each web under both efforts, so index `i` is the same web in both lists. It takes `effort_factual` and `effort_counter` and returns `list(factual, counter)`. Failed members are stored as `NULL` with `x[i] <- list(...)`, so slots are never deleted and the pairing is preserved. Members that failed in either scenario are dropped together.
- `n_sims` is currently 1 (set to 1000 for the full ensemble); `sd_web = 0.1` (about 10% CV per link), 0.2 for a stress test; `seed = 42`.
- Results are cached in RDS files whose names encode `n_sims` and `sd_web` (for example `sims_with_fishing_n1_sd0.10.RDS`), so a small test run cannot be loaded as the full ensemble. The cache cannot detect changed params, effort or forcing: set `force_resim <- TRUE` or delete the RDS after any such change.
- The simulations cover 1861-2010 only, so plots cannot extend past 2010 unless the effort, temperature and plankton series are extended first.

---

## E2. FishMIP output variables

Built from `sims_with_fishing` and `sims_with_no_fishing`. The script uses no helper functions: one loop builds annual per-species biomass and yield tables for every simulation, and each output is a short explicit pipeline on those tables (sum species within each sim, then median and min/max across sims, converted from g to g m-2). Every output is computed, exported and plotted for both scenarios.

**Output files**, one set per scenario, in separate folders:

```
Fishmip_output/with_fishing/mizer_gfdl-mom6-cobalt2_obsclim_histsoc_15arcmin_<var>_ESS_Yearly_1961_2010.csv
Fishmip_output/without_fishing/mizer_gfdl-mom6-cobalt2_obsclim_nat_15arcmin_<var>_ESS_Yearly_1961_2010.csv
```

The tags (`histsoc`, `nat`) are set in `scen_tag` at the top; `nat` is a placeholder, so change it to whatever your protocol calls the no-fishing run.

| Variable | Description |
|---|---|
| `tcb` | Total consumer biomass |
| `tcblog10` | Biomass in log10 weight bins (monthly, first sim only) |
| `tpb`, `tdb` | Pelagic (Herring, Capelin, Sand Lance) and demersal biomass |
| `tc` | Total catch (vs NAFO 4VW) |
| `tclog10` | Catch in log10 weight bins (no header, no Year column; first sim only) |
| `tpc`, `tdc` | Pelagic and demersal catch |
| `bp30cm` ... `bd90cm`, `cp30cm` ... `cd90cm` | Biomass and catch by size class (L_inf <30, 30-90, >=90 cm) x pelagic/demersal |

**Catch accounting:**
- Seals act in the model but are not fishery catch, so every catch output sums only the gears whose name does not contain "seal". The script prints `gear_names` and `fishing_gears` at the top; check them.
- `tclog10` uses `getFMortGear()` summed over those gears (it replaces `getFMort()`, which includes seals). It holds a time x gear x species x weight array for one sim at a time.
- **Without fishing every gear has zero effort, so all catch outputs for that scenario are zero by construction.** They are still exported and plotted for completeness.
- Scale checks (`log10(median obs/sim)`; about 0 is consistent, about 3 or 6 flags a unit problem) use the with-fishing run only.

**Required inputs:** `Useful_functions.R`, `params_stage1_steady.RDS`, `species_params_with_temp.csv`, `GFDL_resource_spectra.dat`, `Temp_by_depth_spin.csv`, the two realm files from C, `Effort_spin_1861_2010.csv`, `biomass_ts.RDS`, `yield_ts.RDS`, `new_NAFO_4VW_Catches.RDS`, and both simulation ensembles in memory (`sims_with_fishing`, `sims_with_no_fishing`).

**Check:**
- `base_dir` in the export Rmd is still an absolute Dropbox path. The dynamics Rmd uses `here()`.
- `catch_to_g <- 1`, but NAFO catches are usually tonnes (likely 1e6); verify with the scale check.
- Sand Lance survey biomass looks implausibly high.
- `tcblog10` and `tclog10` use the first simulation only.
- The tag `nat` is a guess (see above).

---

## Diagnostic figures

Ensemble median and min/max ribbons for each scenario, with observations as black points, using shared colours (blue = without fishing, red = with fishing).
- **Biomass by species:** 19 panels, free scales. Simulations are shown over the span of the observations.
- **Total biomass:** species are summed within each ensemble member before taking the median and range (summing per-species medians would not give the median of the total). The observed total is only comparable across years if the same species are present each year; the script prints how many species each year has.
- **Catch by species:** the 16 species with catches (Sand Lance, Capelin and Longhorn Sculpin excluded). Seal gears are excluded. The without-fishing line is zero by construction.

Plot fixes applied: the panel border needs `fill = NA` or it paints over the data; y-axis labels use automatic accuracy (a fixed `accuracy = 1` rounds small species to "0"); colours are a named vector so they cannot drift from the legend labels; observations are drawn once, not once per scenario.

---

## Open issue: modelled catches exceed observed

In the early years (1980s) modelled catches are well above observed for Atlantic Cod (about 400 kt against about 100 kt), Pollock, Redfish, Atlantic Herring and Atlantic Halibut, and Silver Hake shows spikes near 4,000 kt around 1982 and 1990 against negligible observations. This is not a plotting issue. Possible causes: biomass drift over the 117 years from the calibrated state to 1978, the resource level (the GFDL-derived resource is about 0.74x Paper 1's), or effort/catchability. Compare the biomass panels for the same species: if biomass is also too high, suspect drift or resource level; if only catches are, suspect effort or catchability.

---

## Cross-script issues to resolve first

1. Align spin-up cycle lengths: A uses 5 years, B uses 10.
2. Add the missing `plankton_all` construction to B (or document where it comes from).
3. Regenerate `GFDL_resource_spectra.dat` on Stage 1 (211 bins) before running E1.
4. Replace hardcoded Dropbox paths with `here()` (the FishMIP export Rmd still sets `base_dir` to an absolute path; several upstream scripts do too).
5. Resolve the early-year catch overshoot above before using the catch outputs.
6. Confirm the no-fishing file tag with the FishMIP protocol.

*Resolved since the previous version:* the duplicate `setup` chunk, the reference period (1978-81 for sanity checks only), the ensemble object names and `effort_counter` handling, and the `sim0` comparison run, which now uses `params_dyn`.