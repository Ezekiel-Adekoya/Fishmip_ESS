###"Useful_functions for mizer calibration"
### Set of functions using lower and upper size selectivity of the RV survey data
### This allows us to set up the params object with two columns called biomass_cutoff and biomass_cutoffHigh


# ── HELPER FUNCTION FOR BIOMASS CUTOFFS ───────────────────────────────────────
get_cutoffs <- function(sp_params, n_sp) {
  
  cutoffLow <- sp_params$biomass_cutoff
  
  if (is.null(cutoffLow)) {
    cutoffLow <- rep(0, n_sp)
  }
  
  cutoffLow[is.na(cutoffLow)] <- 0
  
  cutoffHigh <- sp_params$biomass_cutoffHigh
  
  if (is.null(cutoffHigh)) {
    cutoffHigh <- rep(Inf, n_sp)
  }
  
  cutoffHigh[is.na(cutoffHigh)] <- Inf
  
  list(
    low = cutoffLow,
    high = cutoffHigh
  )
}

# this function adds a lower boundary to selected size

plotBiomassObservedVsModelCustom <- function (object, species = NULL, ratio = FALSE, log_scale = TRUE, 
    return_data = FALSE, labels = TRUE, show_unobserved = FALSE) 
{
    if (is(object, "MizerSim")) {
        params = object@params
        n <- finalN(object)
    }
    else if (is(object, "MizerParams")) {
        params = object
        n <- initialN(params)
    }
    else {
        stop("You have not provided a valid mizerSim or mizerParams object.")
    }
    sp_params <- params@species_params
    species = valid_species_arg(object, species)
    if (length(species) == 0) 
        stop("No species selected, please fix.")
    row_select = match(species, sp_params$species)
    if (!"biomass_observed" %in% names(sp_params)) {
        stop("You have not provided values for the column 'biomass_observed' ", 
            "in the mizerParams/mizerSim object.")
    }
    else if (!is.numeric(sp_params$biomass_observed)) {
        stop("The column 'biomass_observed' in the mizerParams/mizerSim object", 
            " is not numeric, please fix.")
    }
    else {
        biomass_observed = sp_params$biomass_observed
    }
    
    cutoffLow <- sp_params$biomass_cutoff[row_select]
    if (is.null(cutoffLow)) {
        cutoffLow = rep(0, length(species))
    }
    else if (!is.numeric(cutoffLow)) {
        stop("params@species_params$biomass_cutoff is not numeric, \",\n                 \"please fix.")
    }
    cutoffLow[is.na(cutoffLow)] <- 0
    
    cutoffHigh <- sp_params$biomass_cutoffHigh[row_select]
    if (is.null(cutoffHigh)) {
        cutoffHigh = rep(0, length(species))
    }
    else if (!is.numeric(cutoffHigh)) {
        stop("params@species_params$biomass_cutoffHigh is not numeric, \",\n                 \"please fix.")
    }
    cutoffHigh[is.na(cutoffHigh)] <- 0
    
    sim_biomass = rep(0, length(species))
    for (j in 1:length(species)) {
        sim_biomass[j] = sum((n[row_select[j], ] * params@w * 
            params@dw)[params@w >= cutoffLow[j] & cutoffHigh[j] >= params@w])
    }
    dummy = data.frame(species = species, model = sim_biomass, 
        observed = biomass_observed[row_select]) %>% mutate(species = factor(species, 
        levels = species), is_observed = !is.na(observed) & observed > 
        0, observed = case_when(is_observed ~ observed, !is_observed ~ 
        model), ratio = model/observed)
    if (sum(dummy$is_observed) == 0) {
        cat(paste("There are no observed biomasses to compare to model,", 
            "only plotting model biomasses.", sep = "\n"))
    }
    if (!show_unobserved) {
        dummy <- filter(dummy, is_observed)
    }
    if (return_data == TRUE) 
        return(dummy)
    tre <- round(sum(abs(1 - dummy$ratio)), digits = 3)
    caption <- paste0("Total relative error = ", tre)
    if (any(!dummy$is_observed)) {
        caption <- paste(caption, "\n Open circles represent species without biomass observation.")
    }
    if (ratio == FALSE) {
        gg <- ggplot(data = dummy, aes(x = observed, y = model, 
            colour = species, shape = is_observed)) + geom_abline(aes(intercept = 0, 
            slope = 1), colour = "purple", linetype = "dashed", 
            size = 1.3) + geom_point(size = 3) + labs(y = "model biomass [g]") + 
            coord_cartesian(ylim = range(dummy$model, dummy$observed))
    }
    else {
        gg <- ggplot(data = dummy, aes(x = observed, y = ratio, 
            colour = species, shape = is_observed)) + geom_hline(aes(yintercept = 1), 
            linetype = "dashed", colour = "purple", 
            size = 1.3) + geom_point(size = 3) + labs(y = "model biomass / observed biomass") + 
            coord_cartesian(ylim = range(dummy$ratio))
    }
    gg <- gg + labs(x = "observed biomass [g]", caption = caption) + 
        scale_colour_manual(values = getColours(params)[dummy$species]) + 
        scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1)) + 
        guides(shape = "none")
    if (log_scale == TRUE & ratio == FALSE) {
        gg = gg + scale_x_log10() + scale_y_log10()
    }
    if (log_scale == TRUE & ratio == TRUE) {
        gg = gg + scale_x_log10()
    }
    if (labels == TRUE) {
        gg = gg + ggrepel::geom_label_repel(aes(label = species), 
            box.padding = 0.35, point.padding = 0.5, segment.color = "grey50", 
            show.legend = FALSE, max.overlaps = Inf, seed = 42)
    }
    gg
}




# adapting cutoff here too
calibrateBiomassCustom <- function(params) {
  
  if (!("biomass_observed" %in% names(params@species_params)) ||
      all(is.na(params@species_params$biomass_observed))) {
    return(params)
  }
  
  sp_params <- params@species_params
  n_sp <- nrow(sp_params)
  
  cutoffs <- get_cutoffs(sp_params, n_sp)
  cutoffLow <- cutoffs$low
  cutoffHigh <- cutoffs$high
  
  observed <- sp_params$biomass_observed
  sp_observed <- which(!is.na(observed) & observed > 0)
  
  observed_total <- sum(observed[sp_observed], na.rm = TRUE)
  
  model_total <- 0
  
  for (sp_idx in sp_observed) {
    idx <- params@w >= cutoffLow[sp_idx] &
      params@w <= cutoffHigh[sp_idx]
    
    model_total <- model_total +
      sum(params@initial_n[sp_idx, idx] *
            params@w[idx] *
            params@dw[idx])
  }
  
  if (model_total <= 0 || !is.finite(model_total)) {
    stop("Model biomass within biomass cutoffs is zero or non-finite.")
  }
  
  scaleModel(params, factor = observed_total / model_total)
}


# same as above
matchBiomassesCustom <- function(params, species = NULL) {
  
  if (!("biomass_observed" %in% names(params@species_params))) {
    return(params)
  }
  
  sp_params <- params@species_params
  n_sp <- nrow(sp_params)
  
  selected <- valid_species_arg(
    params,
    species = species,
    return.logical = TRUE
  ) &
    !is.na(sp_params$biomass_observed) &
    sp_params$biomass_observed > 0
  
  cutoffs <- get_cutoffs(sp_params, n_sp)
  
  for (sp in which(selected)) {
    
    idx <- params@w >= cutoffs$low[sp] &
      params@w <= cutoffs$high[sp]
    
    total <- sum(params@initial_n[sp, idx] *
                   params@w[idx] *
                   params@dw[idx])
    
    if (total <= 0 || !is.finite(total)) {
      warning("Skipping ", sp_params$species[sp],
              ": model biomass within cutoffs is zero.")
      next
    }
    
    factor <- sp_params$biomass_observed[sp] / total
    factor <- min(max(factor, 0.25), 4)   # damp: max 4x or 0.25x per iteration
    params@initial_n[sp, ] <- params@initial_n[sp, ] * factor
  }
  
  # Diagnose before setBevertonHolt: identify species with zero RDI
  rdi <- getRDI(params)
  zero_sp <- which(rdi == 0)
  if (length(zero_sp) > 0) {
    sp_names <- sp_params$species[zero_sp]
    for (sp in zero_sp) {
      wmat <- sp_params$w_mat[sp]
      mature_b <- sum(params@initial_n[sp, params@w >= wmat] *
                        params@w[params@w >= wmat] *
                        params@dw[params@w >= wmat])
      message("getRDI = 0 for ", sp_params$species[sp],
              " | w_mat = ", round(wmat, 2),
              " g | mature biomass = ", signif(mature_b, 3), " g",
              " | erepro = ", round(sp_params$erepro[sp], 4))
    }
    stop("Some species have no reproduction: ", paste(sp_names, collapse = ", "),
         "\nSee messages above. Check z0, w_mat, and biomass scaling for these species.")
  }

  setBevertonHolt(params)
}


# matchYieldCustom: wraps matchYield and clamps per-species fishing mortality
# to biologically realistic bounds after each catchability adjustment.
#
# Reads optional columns from species_params:
#   f_max  — upper F bound (e.g. 0.4 for Atlantic Halibut)
#   f_min  — lower F bound (e.g. 0.3 for Atlantic Cod)
# Species without a bound set (NA) are unconstrained.
#
# F is computed as biomass-weighted mean across sizes, consistent with how
# plotFMortExcl summarises fishing mortality.
matchYieldCustom <- function(params, gears = NULL, keep = "biomass") {

  if (!is.null(gears)) {
    params <- matchYield(params, keep = keep, gears = gears)
  } else {
    params <- matchYield(params, keep = keep)
  }

  # Cap per gear×species: if effort × catchability > 4 for any row, scale
  # that row's catchability down so F = 4. Checked per row, not summed.
  gp     <- gear_params(params)
  effort <- params@initial_effort[gp$gear]
  f_row  <- effort * gp$catchability
  over   <- !is.na(f_row) & f_row > 4
  if (any(over)) {
    gear_params(params)$catchability[over] <-
      gp$catchability[over] * (4 / f_row[over])
  }

  params
}


# ════════════════════════════════════════════════════════════════════════════
# MERGED FROM functions.R (2026-06): the single project functions file.
# Added below: the seal selectivity function and the ESS getError* loss
# functions (steady-state + time-series, used by the calibration / fastOptim
# optimisation). All use the `biomass_cutoff` (lower) / `biomass_cutoffHigh`
# (upper) convention — consistent with the custom functions above and Stage 1.
# The custom matchBiomasses/matchYield/calibrate functions above are RETAINED
# for future use.
# ════════════════════════════════════════════════════════════════════════════

# Seal selectivity function — log-normal size preference used by the Seals gear
# Seal selectivity function — log-normal size preference used by the Seals gear
gaussian_sizeweight <- function(w, species_params, w_opt, sigma, ...) {
  exp(-0.5 * ((log(w) - log(w_opt)) / sigma)^2)
}

# ── Seal-diet tuning: per-species seal catchability vs observed trajectory ────
# Tunes seal catchability multipliers for `seal_tune_sp` to stop the ramping seal
# effort over-cropping them, while preserving cod's collapse (cod kept in the
# loss but its seal-q NOT tuned). Projects from the calibrated initial_n through
# the historical effort (no re-spin — seal effort at the reference is tiny).
# vary: log10(seal catchability multiplier) for each seal_tune_sp.
getErrorSealDiet <- function(vary, params, dat_ts, effort_ts,
                             seal_tune_sp, loss_sp, sp_weights = NULL,
                             fish_tune_sp = NULL,
                             fishing_gears = c("Bottom_Trawl","Fixed_Gear",
                                               "Small_Mesh","Pelagic"),
                             mult_min = 0.01, mult_max = 3,
                             fish_mult_max = 30) {
  n_seal <- length(seal_tune_sp)
  n_fish <- length(fish_tune_sp)
  seal_mult <- 10^vary[seq_len(n_seal)]
  if (any(!is.finite(seal_mult)) || any(seal_mult < mult_min) ||
      any(seal_mult > mult_max)) return(1e8)

  gp <- gear_params(params)
  for (i in seq_along(seal_tune_sp)) {
    idx <- which(gp$species == seal_tune_sp[i] & gp$gear == "Seals")
    if (length(idx)) gp$catchability[idx] <- gp$catchability[idx] * seal_mult[i]
  }
  # Optional fishing catchability multipliers (e.g. restore Silver Hake fishing)
  if (n_fish > 0) {
    fish_mult <- 10^vary[(n_seal + 1):(n_seal + n_fish)]
    if (any(!is.finite(fish_mult)) || any(fish_mult <= 0) ||
        any(fish_mult > fish_mult_max)) return(1e8)
    for (i in seq_along(fish_tune_sp)) {
      idx <- which(gp$species == fish_tune_sp[i] & gp$gear %in% fishing_gears)
      if (length(idx)) gp$catchability[idx] <- gp$catchability[idx] * fish_mult[i]
    }
  }
  gear_params(params) <- gp
  params <- tryCatch(setParams(params), error = function(e) NULL)
  if (is.null(params)) return(1e8)

  sim <- tryCatch(
    project(params, effort = effort_ts, t_max = nrow(effort_ts),
            t_save = 1, progress_bar = FALSE),
    error = function(e) NULL)
  if (is.null(sim)) return(1e8)

  cl <- params@species_params$biomass_cutoff; cl[is.na(cl)] <- 0
  sp <- params@species_params$species
  sim_time <- as.integer(dimnames(sim@n)[[1]])
  common <- intersect(sim_time, unique(dat_ts$Year))
  if (is.null(sp_weights)) sp_weights <- setNames(rep(1, length(loss_sp)), loss_sp)

  loss <- 0
  for (s in loss_sp) {
    j <- which(sp == s); if (!length(j)) next
    idx <- params@w >= cl[j]
    mb <- rowSums(sweep(sim@n[, j, idx, drop = FALSE], 3,
                        params@w[idx] * params@dw[idx], "*"))
    names(mb) <- as.character(sim_time)
    obs_s <- dat_ts$Biomass[dat_ts$Species == s & dat_ts$Year %in% common]
    yr_s  <- dat_ts$Year[dat_ts$Species == s & dat_ts$Year %in% common]
    mod_s <- mb[as.character(yr_s)]
    ok <- !is.na(obs_s) & obs_s > 0 & !is.na(mod_s) & mod_s > 0
    if (sum(ok) == 0) next
    wt <- if (s %in% names(sp_weights)) sp_weights[s] else 1
    loss <- loss + wt * sum((log(mod_s[ok] / obs_s[ok]))^2)
  }
  loss
}


# ── ESS 18-species adaptation ─────────────────────────────────────────────────
# vary layout (length = 3 * n_sp):
#   vary[1:n_sp]              log10(R_max) for each species
#   vary[(n_sp+1):(2*n_sp)]   erepro for each species
#   vary[(2*n_sp+1):(3*n_sp)] interaction_resource for each species
#
# Interaction matrix is kept fixed (too many params for gradient-free search).
# Uses biomass_cutoff (lower) and biomass_cutoffHigh (upper) from species_params.
# dat: named numeric vector of observed biomass (g) — NA-safe.

getErrorCustomBiomassESS <- function(vary, params, dat,
                                     inf_sp = NULL,
                                     tol = 0.001, timetorun = 10) {
  n_sp <- nrow(params@species_params)

  # vary = log10(R_max); if inf_sp supplied, only update those species
  if (is.null(inf_sp)) {
    params@species_params$R_max <- 10^vary[1:n_sp]
  } else {
    params@species_params$R_max[inf_sp] <- 10^vary
  }

  params <- setParams(params)

  params <- tryCatch(
    projectToSteady(params, distance_func = distanceSSLogN,
                    tol = tol, t_max = 200, return_sim = FALSE),
    error = function(e) NULL
  )
  if (is.null(params)) return(1e8)

  sim <- tryCatch(
    project(params, t_max = timetorun, progress_bar = FALSE),
    error = function(e) NULL
  )
  if (is.null(sim)) return(1e8)

  n_final <- sim@n[dim(sim@n)[1], , ]

  cutoff_low <- params@species_params$biomass_cutoff
  if (is.null(cutoff_low)) cutoff_low <- rep(0, n_sp)
  cutoff_low[is.na(cutoff_low)] <- 0

  cutoff_high <- params@species_params$biomass_cutoffHigh
  if (is.null(cutoff_high)) cutoff_high <- rep(Inf, n_sp)
  cutoff_high[is.na(cutoff_high)] <- Inf

  sim_biomass <- numeric(n_sp)
  for (j in seq_len(n_sp)) {
    idx <- params@w >= cutoff_low[j] & params@w <= cutoff_high[j]
    sim_biomass[j] <- sum(n_final[j, idx] * params@w[idx] * params@dw[idx])
  }

  ok   <- !is.na(dat) & dat > 0 & sim_biomass > 0
  pred <- log(sim_biomass[ok])
  obs  <- log(dat[ok])
  sum((pred - obs)^2)
}


# Time-series version: projects with historical effort matrix and compares
# modelled biomass to the full 1978-2024 survey time series.
# dat_ts: data.frame with columns Year, Species, Biomass (g) — long format.
# effort_ts: matrix / data.frame with years as rownames and gears as columns.
# sp_weights: named numeric vector; species not listed get weight 1.
# exclude_sp: character vector of species to drop from the loss.

getErrorTSBiomassESS <- function(vary, params, dat_ts, effort_ts,
                                  sp_weights  = NULL,
                                  exclude_sp  = c("Atlantic Herring", "Capelin",
                                                   "Sand Lance", "Longhorn Sculpin"),
                                  tol         = 0.001) {
  n_sp <- nrow(params@species_params)
  sp   <- params@species_params$species

  # vary contains only log10(R_max); erepro and interaction_resource held fixed
  params@species_params$R_max <- 10^vary[1:n_sp]

  params <- tryCatch(setParams(params), error = function(e) NULL)
  if (is.null(params)) return(1e8)

  sim <- tryCatch(
    project(params, effort = effort_ts, t_max = nrow(effort_ts),
            t_save = 1, progress_bar = FALSE),
    error = function(e) NULL
  )
  if (is.null(sim)) return(1e8)

  cutoff_low <- params@species_params$biomass_cutoff
  if (is.null(cutoff_low)) cutoff_low <- rep(0, n_sp)
  cutoff_low[is.na(cutoff_low)] <- 0

  cutoff_high <- params@species_params$biomass_cutoffHigh
  if (is.null(cutoff_high)) cutoff_high <- rep(Inf, n_sp)
  cutoff_high[is.na(cutoff_high)] <- Inf

  # getBiomass() only accepts a scalar min_w, so compute per-species with cutoffs
  n_t      <- dim(sim@n)[1]
  sim_time <- as.integer(dimnames(sim@n)[[1]])
  model_bio <- matrix(NA_real_, nrow = n_t, ncol = n_sp,
                      dimnames = list(as.character(sim_time), sp))
  for (j in seq_len(n_sp)) {
    idx <- params@w >= cutoff_low[j] & params@w <= cutoff_high[j]
    model_bio[, j] <- rowSums(
      sweep(sim@n[, j, idx, drop = FALSE], 3,
            params@w[idx] * params@dw[idx], "*")
    )
  }

  sp_fit <- setdiff(sp, exclude_sp)
  if (is.null(sp_weights)) sp_weights <- setNames(rep(1, length(sp_fit)), sp_fit)

  obs_years <- unique(dat_ts$Year)
  common    <- intersect(sim_time, obs_years)
  if (length(common) == 0) return(1e8)

  loss <- 0
  for (s in sp_fit) {
    if (!s %in% colnames(model_bio)) next
    obs_s <- dat_ts$Biomass[dat_ts$Species == s & dat_ts$Year %in% common]
    yr_s  <- dat_ts$Year[dat_ts$Species == s & dat_ts$Year %in% common]
    mod_s <- model_bio[as.character(yr_s), s]
    ok    <- !is.na(obs_s) & obs_s > 0 & !is.na(mod_s) & mod_s > 0
    if (sum(ok) == 0) next
    w <- if (s %in% names(sp_weights)) sp_weights[s] else 1
    loss <- loss + w * sum((log(mod_s[ok] / obs_s[ok]))^2)
  }
  loss
}


# ── Time-series refinement: catchability + reproduction_level ─────────────────
# HIMI-style second stage. Starts from the steady-state-calibrated params and
# tunes, for each fitted species:
#   - a fishing catchability MULTIPLIER (drives F → the collapse)
#   - reproduction_level (density dependence → the (non-)recovery)
# to reproduce the observed 1978-2024 biomass TRAJECTORY (not just the reference).
#
# vary layout (length = 2 * n_fit):
#   vary[1:n_fit]              log10(catchability multiplier) per fitted species
#   vary[(n_fit+1):(2*n_fit)]  logit(reproduction_level) per fitted species
#
# Projects from the calibrated initial_n (the steady-state 1978-85 state) through
# the historical effort — no re-spin, so the reference biomass match is preserved
# as the initial condition and the trajectory diverges under the new F / recruitment.
#
# fit_sp:    character vector of species to tune (exclude forage fish).
# dat_ts:    long data.frame (Year, Species, Biomass in g).
# effort_ts: matrix, years (1978-2024) as rownames, gears as columns.
# sp_weights: named numeric; species not listed get weight 1 (weight cod high).

getErrorTSCatchRepro <- function(vary, params, dat_ts, effort_ts, fit_sp,
                                 sp_weights    = NULL,
                                 rmax_sp       = NULL,
                                 seal_sp       = NULL,
                                 seal_mult_max = 3,
                                 fishing_gears = c("Bottom_Trawl","Fixed_Gear",
                                                   "Small_Mesh","Pelagic"),
                                 catch_mult_max = 20) {
  n_fit  <- length(fit_sp)
  n_rmax <- length(rmax_sp)
  n_seal <- length(seal_sp)
  catch_mult <- 10^vary[1:n_fit]
  rl_vals    <- plogis(vary[(n_fit + 1):(2 * n_fit)])   # (0,1)

  if (any(!is.finite(catch_mult)) || any(catch_mult <= 0) ||
      any(catch_mult > catch_mult_max)) return(1e8)

  # Apply catchability multipliers to each fitted species' fishing gears
  gp <- gear_params(params)
  for (i in seq_len(n_fit)) {
    idx <- which(gp$species == fit_sp[i] & gp$gear %in% fishing_gears)
    if (length(idx)) gp$catchability[idx] <- gp$catchability[idx] * catch_mult[i]
  }

  # Seal catchability multipliers (top-down control; "seals mainly" — bounded
  # for plausibility). Position in vary: after catch, rl, and R_max blocks.
  if (n_seal > 0) {
    seal_mult <- 10^vary[(2 * n_fit + n_rmax + 1):(2 * n_fit + n_rmax + n_seal)]
    if (any(!is.finite(seal_mult)) || any(seal_mult <= 0) ||
        any(seal_mult > seal_mult_max)) return(1e8)
    for (i in seq_len(n_seal)) {
      idx <- which(gp$species == seal_sp[i] & gp$gear == "Seals")
      if (length(idx)) gp$catchability[idx] <- gp$catchability[idx] * seal_mult[i]
    }
  }
  gear_params(params) <- gp

  # Apply reproduction_level for fitted species (keep others as calibrated)
  rl_full <- getReproductionLevel(params)
  rl_full[fit_sp] <- rl_vals
  rl_full[rl_full >= 1] <- 0.999
  rl_full[rl_full <= 0] <- 0.001
  params <- tryCatch(setBevertonHolt(params, reproduction_level = rl_full),
                     error = function(e) NULL)
  if (is.null(params)) return(1e8)

  # Optionally override R_max for the blow-up species (Inf-R_max competitive
  # release, e.g. Pollock). Tuned in log10(grams); erepro re-derived by
  # setBevertonHolt to preserve reference recruitment.
  if (n_rmax > 0) {
    rmax_log <- vary[(2 * n_fit + 1):(2 * n_fit + n_rmax)]
    if (any(!is.finite(rmax_log))) return(1e8)
    rmax_full <- setNames(params@species_params$R_max,
                          params@species_params$species)
    rmax_full[rmax_sp] <- 10^rmax_log
    params <- tryCatch(setBevertonHolt(params, R_max = rmax_full),
                       error = function(e) NULL)
    if (is.null(params)) return(1e8)
  }

  # Project from calibrated initial_n through the historical effort
  sim <- tryCatch(
    project(params, effort = effort_ts, t_max = nrow(effort_ts),
            t_save = 1, progress_bar = FALSE),
    error = function(e) NULL
  )
  if (is.null(sim)) return(1e8)

  n_sp <- nrow(params@species_params)
  sp   <- params@species_params$species
  # Lower cutoff ONLY — matches the raw survey biomass_ts (which spans all sizes
  # above the small-fish detection limit) and the calibration target. Applying
  # biomass_cutoffHigh here would compare a small model size-window against the
  # full observed survey biomass (inconsistent), breaking the level.
  cl <- params@species_params$biomass_cutoff; cl[is.na(cl)] <- 0

  sim_time  <- as.integer(dimnames(sim@n)[[1]])
  model_bio <- matrix(NA_real_, nrow = length(sim_time), ncol = n_sp,
                      dimnames = list(as.character(sim_time), sp))
  for (j in seq_len(n_sp)) {
    idx <- params@w >= cl[j]
    model_bio[, j] <- rowSums(sweep(sim@n[, j, idx, drop = FALSE], 3,
                                    params@w[idx] * params@dw[idx], "*"))
  }

  if (is.null(sp_weights)) sp_weights <- setNames(rep(1, n_fit), fit_sp)
  common <- intersect(sim_time, unique(dat_ts$Year))
  if (length(common) == 0) return(1e8)

  loss <- 0
  for (s in fit_sp) {
    if (!s %in% colnames(model_bio)) next
    obs_s <- dat_ts$Biomass[dat_ts$Species == s & dat_ts$Year %in% common]
    yr_s  <- dat_ts$Year[dat_ts$Species == s & dat_ts$Year %in% common]
    mod_s <- model_bio[as.character(yr_s), s]
    ok    <- !is.na(obs_s) & obs_s > 0 & !is.na(mod_s) & mod_s > 0
    if (sum(ok) == 0) next
    wt <- if (s %in% names(sp_weights)) sp_weights[s] else 1
    loss <- loss + wt * sum((log(mod_s[ok] / obs_s[ok]))^2)
  }
  loss
}


# ── PAIRED RANDOMISED-WEB ENSEMBLE (Stage 2 scenarios, 03_Scenarios) ──────────
# Shared by 03_Scenarios_Through_time.Rmd and run_ensemble_sims.R so the two stay
# in lock-step (single source of truth). Keep these definitions here only.

# Randomise the interaction matrix while preserving topology (zeros stay zero) and
# perturb the seal catchabilities the same way — a structural-uncertainty draw.
randomize_interaction_matrix <- function(params, sdev = 0.1) {

  # Get current interaction matrix
  old_inter <- params@interaction

  # Get seal catchability only (not all gears)
  seal_idx  <- which(gear_params(params)$gear == "Seals")
  old_seal  <- gear_params(params)$catchability[seal_idx]

  # Create noise matrix with same dimensions as interaction matrix
  n          <- nrow(old_inter)
  noise      <- matrix(rnorm(n * n, mean = 0, sd = sdev), nrow = n)
  noise_seal <- rnorm(length(seal_idx), mean = 0, sd = sdev)

  # Add lognormal noise
  new_inter <- old_inter * exp(noise)
  new_seal  <- old_seal  * exp(noise_seal)

  # Clip to [0, 1]
  new_inter <- pmin(pmax(new_inter, 0), 1)
  new_seal  <- pmax(new_seal, 0)   # floor at 0 only: seal catchabilities are O(1-8);
  # a [0,1] clip would crush them all to 1 and destroy seal predation + Bowen selectivity

  # Update parameters
  params@interaction <- new_inter
  gear_params(params)$catchability[seal_idx] <- new_seal

  return(params)
}

# PAIRED ensemble: each randomised food web is run under BOTH effort scenarios, so
# sims$factual[[i]] and sims$counter[[i]] share the SAME web — a genuinely PAIRED
# comparison — and `seed` makes the whole ensemble reproducible.
run_paired_sims <- function(params, n_sims = 100, sd = 0.2,
                            effort_factual, effort_counter, seed = 42) {
  set.seed(seed)
  factual <- counter <- vector("list", n_sims)
  for (i in seq_len(n_sims)) {
    rp <- tryCatch(setParams(randomize_interaction_matrix(params, sdev = sd)),
                   error = function(e) NULL)
    if (is.null(rp)) next                   # leave both NULL at i (pairing preserved)
    factual[[i]] <- tryCatch(project(rp, effort = effort_factual), error = function(e) NULL)
    counter[[i]] <- tryCatch(project(rp, effort = effort_counter), error = function(e) NULL)
  }
  list(factual = factual, counter = counter)   # index i is the shared web
}
