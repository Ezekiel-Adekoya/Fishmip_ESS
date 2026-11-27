# Prey switching for the grey seal
#
# In plain mizer a predator's preference for each prey species is fixed, so
# the seal takes each prey in direct proportion to how much of it is there.
# Smout et al. (2014, ICES J. Mar. Sci. 71: 81-89) fitted the diets of UK grey
# seals and found a switching response clearly favoured: the rate at which
# seals take a prey rises with that prey's availability to the power m = 1.5
# (m = 1 is no switching). Seals therefore take disproportionately little of a
# prey that has become scarce and disproportionately much of one that has
# become common.
#
# Here that is done by multiplying the seal's fixed preference for each fish
# species, at every time step, by
#
#     (availability now / availability in the calibrated state) ^ (m - 1)
#
# so the multiplier is 1 in the calibrated state and the calibration is not
# disturbed. Availability is the biomass of the species weighted by the sizes
# a typical adult seal eats. m is taken from Smout et al. and is not fitted.
#
# The exponent can be set to other values to try other responses:
#   m = 1.5  Smout et al.: seals turn away from a prey as it becomes scarce.
#   m = 1    plain mizer: seals take each prey in proportion to what is there.
#   m = 0    fixed shares: each seal eats the same amount of a prey however
#            much of it there is, so the diet keeps its proportions and the
#            mortality seals cause on a prey RISES as the prey becomes scarce.
#            The two Sable Island diets (1991-93 and 1993-98) look more like
#            this for cod and flatfish, whose shares held or rose while their
#            stocks fell. Two diet samples are thin support, so treat it as a
#            hypothesis. Because the mortality grows without limit as a prey
#            vanishes, `largest_multiplier` caps how far the preference can
#            rise above its calibrated value.
#
# A Type III response for each prey (`half_saturation` given):
#   The seal's intake of a prey follows  N^m / (K^m + N^m), with N the prey's
#   availability and K the availability at which intake is half its maximum.
#   With m = 1.5 this has both halves of what the data seem to need. While a
#   prey is well above K the intake is nearly flat, so seals take a rising
#   fraction of a prey that is becoming scarce (the squeeze). Below about
#   0.6 K the intake falls away faster than the prey does, so seals take a
#   falling fraction of a prey that is nearly gone (the relief). The fraction
#   taken is highest near 0.6 K: the "predator pit".
#   `half_saturation` is K as a fraction of the prey's availability in the
#   calibrated state, the same for every prey. It is not measured for these
#   prey. In Smout et al.'s model the saturation is shared across all prey;
#   saturation prey by prey is a departure from it, suggested by the Sable
#   Island diets.
#
# Use:  params <- add_seal_prey_switching(params)
# The model is given the extra class "sealPreySwitching". mizer then calls the
# two functions at the bottom of this file whenever it works out what
# predators encounter and the mortality they cause, and those hand on to
# therMizer's and mizer's own functions with the seal's preference replaced.
# This is the same way therMizer adds its temperature effect. The two
# functions must stay defined while the model is run, so every script that
# projects with switching has to source this file.

smout_switching_exponent <- 1.5

# The seal's preference for every species at one moment, given the abundances.
seal_preference_now <- function(params, n) {
  switching <- params@other_params$other$seal_prey_switching
  availability_now <- as.vector(n[switching$prey_rows, , drop = FALSE] %*%
                                  switching$weight_by_size)
  relative_availability <- availability_now / switching$reference_availability
  if (is.null(switching$half_saturation)) {
    preference_multiplier <- pmin(relative_availability^(switching$exponent - 1),
                                  switching$largest_multiplier)
  } else {
    # Type III for each prey, scaled to be 1 in the calibrated state.
    m <- switching$exponent
    k <- switching$half_saturation
    preference_multiplier <- relative_availability^(m - 1) * (1 + k^m) /
      (relative_availability^m + k^m)
  }
  preference <- params@interaction
  preference[switching$seal_row, switching$prey_rows] <-
    preference[switching$seal_row, switching$prey_rows] * preference_multiplier
  preference
}

# Stores what the switching needs in the model and tells mizer to use the two
# functions below. `reference_n` is the abundance at which the multiplier is 1.
# `no_switching_for` names prey the seal does not switch on: its preference
# for them stays fixed.
add_seal_prey_switching <- function(params, reference_n = params@initial_n,
                                    exponent = smout_switching_exponent,
                                    no_switching_for = character(0),
                                    largest_multiplier = Inf,
                                    half_saturation = NULL,
                                    typical_seal_weight = 170000) {
  species <- as.character(params@species_params$species)
  seal_row <- which(species == "Grey Seal")
  prey_rows <- which(species != "Grey Seal")
  stopifnot(all(no_switching_for %in% species))
  exponent_by_prey <- ifelse(species[prey_rows] %in% no_switching_for, 1, exponent)
  # How strongly a seal of typical weight selects each prey size, times the
  # weight and width of each size class: abundance times this is the biomass
  # in the sizes seals eat.
  fish_size_bins <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
  seal_size_bin <- which.min(abs(params@w - typical_seal_weight))
  size_selection <- getPredKernel(params)[seal_row, seal_size_bin, fish_size_bins]
  weight_by_size <- size_selection * params@w * params@dw
  params@other_params$other$seal_prey_switching <- list(
    seal_row = seal_row, prey_rows = prey_rows, exponent = exponent_by_prey,
    largest_multiplier = largest_multiplier, half_saturation = half_saturation,
    weight_by_size = weight_by_size,
    reference_availability = as.vector(reference_n[prey_rows, , drop = FALSE] %*% weight_by_size))
  # Registered as a mizer extension, so that mizer keeps the class when it
  # checks the model at the start of a projection.
  if (!inherits(params, "sealPreySwitching")) {
    params <- recordExtension(params, "sealPreySwitching")
  }
  params
}

# The same encounter and predation mortality as before, with the seal's
# preference replaced by its preference at this moment.
projectEncounter.sealPreySwitching <- function(params, n, n_pp, n_other, t = 0, ...) {
  params@interaction <- seal_preference_now(params, n)
  NextMethod()
}
projectPredMort.sealPreySwitching <- function(params, n, n_pp, n_other, t = 0, pred_rate, ...) {
  params@interaction <- seal_preference_now(params, n)
  NextMethod()
}
