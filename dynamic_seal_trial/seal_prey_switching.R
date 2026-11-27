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
  preference_multiplier <-
    (availability_now / switching$reference_availability)^(switching$exponent - 1)
  preference <- params@interaction
  preference[switching$seal_row, switching$prey_rows] <-
    preference[switching$seal_row, switching$prey_rows] * preference_multiplier
  preference
}

# Stores what the switching needs in the model and tells mizer to use the two
# functions below. `reference_n` is the abundance at which the multiplier is 1.
add_seal_prey_switching <- function(params, reference_n = params@initial_n,
                                    exponent = smout_switching_exponent,
                                    typical_seal_weight = 170000) {
  species <- as.character(params@species_params$species)
  seal_row <- which(species == "Grey Seal")
  prey_rows <- which(species != "Grey Seal")
  # How strongly a seal of typical weight selects each prey size, times the
  # weight and width of each size class: abundance times this is the biomass
  # in the sizes seals eat.
  fish_size_bins <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
  seal_size_bin <- which.min(abs(params@w - typical_seal_weight))
  size_selection <- getPredKernel(params)[seal_row, seal_size_bin, fish_size_bins]
  weight_by_size <- size_selection * params@w * params@dw
  params@other_params$other$seal_prey_switching <- list(
    seal_row = seal_row, prey_rows = prey_rows, exponent = exponent,
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
