# 10. Make the seal's diet match the measurements in the years they were made
#
# The seal's diet was measured in 1991-1993 and in 1993-1998
# (seal_diet_data.R). Script 8 calibrates the seal in the 1978-1981 state, and
# by the 1990s the fish community has changed. This script repeats
#
#   script 8 (calibrate the seal) -> script 5 (seal growth) -> script 9 (project)
#
# and after each round moves the shares that script 8 calibrates to, so that
# the projected diet comes closer to the measured one in both periods. A
# species that is too small a share of the projected diet gets a larger share
# in the calibrated state, and the other way round.
#
# The seal has one preference for each prey, so the two periods cannot both be
# matched exactly: each species is moved by the average of what the two
# periods ask for, and the result is a compromise between them.
#
# The seal's fitted growth, harvest and pup limit are taken from an earlier
# run (SEAL_RATES_FROM) during the rounds, to save time. Run script 9 once
# more with FIT_SEALS_IN_FULL_MODEL=1 afterwards.
#
# Needs the same SC, BASE_PARAMS and SEAL_SOURCE settings as scripts 8 and 9.

source("dynamic_seal_trial/seal_diet_data.R")
work_folder <- Sys.getenv("SC")
number_of_rounds <- as.integer(Sys.getenv("DIET_ROUNDS", "6"))
# No calibrated share may move further than this factor from the average of
# the measured ones, and each round moves only part of the way. Without these
# limits the rounds run away: a larger share of cod in the calibrated diet
# makes the seals reduce cod, cod then becomes a smaller share of the later
# diet, and the next round asks for more cod again.
largest_adjustment <- 2
step_towards_measured <- 0.4
reference_diet_file <- "reference_diet.rds"
Sys.setenv(REFERENCE_DIET_FILE = reference_diet_file)
named_species <- colnames(seal_diet_share)

run_script <- function(script, log_name) {
  status <- system2("Rscript", file.path("dynamic_seal_trial", script),
                    stdout = file.path(work_folder, log_name), stderr = file.path(work_folder, log_name))
  if (status != 0) stop(script, " failed; see ", file.path(work_folder, log_name))
}
# How far the projected diet is from the measured one: the mean, over species
# and periods, of the squared difference in log share.
diet_mismatch <- function(model_share) {
  mean((log(pmax(model_share[, named_species], 1e-4)) -
          log(pmax(seal_diet_share, 1e-4)))^2, na.rm = TRUE)
}

best_mismatch <- Inf
for (round_number in 1:number_of_rounds) {
  run_script("8_seal_on_bowen_diet.R", "diet_round_8.log")
  run_script("5_fit_seal_growth_and_harvest.R", "diet_round_5.log")
  run_script("9_project_seal_on_bowen_diet.R", "diet_round_9.log")
  run <- readRDS(file.path(work_folder, "dynamic_seal_run.rds"))
  model_share <- run$model_diet_share
  mismatch <- diet_mismatch(model_share)
  cat(sprintf("round %d: diet mismatch %.3f | sand lance %.0f/%.0f, cod %.0f/%.0f, plaice %.1f/%.1f, yellowtail %.1f/%.1f, silver hake %.1f/%.1f (percent in 1991-93/1993-98)\n",
              round_number, mismatch,
              100 * model_share[1, "Sand Lance"], 100 * model_share[2, "Sand Lance"],
              100 * model_share[1, "Atlantic Cod"], 100 * model_share[2, "Atlantic Cod"],
              100 * model_share[1, "American Plaice"], 100 * model_share[2, "American Plaice"],
              100 * model_share[1, "Yellowtail Flounder"], 100 * model_share[2, "Yellowtail Flounder"],
              100 * model_share[1, "Silver Hake"], 100 * model_share[2, "Silver Hake"]))

  calibrated_share <- if (file.exists(file.path(work_folder, reference_diet_file))) {
    readRDS(file.path(work_folder, reference_diet_file))[named_species]
  } else starting_seal_diet_share[named_species]
  # Keep the calibration that gave the closest diet so far.
  if (mismatch < best_mismatch) {
    best_mismatch <- mismatch
    saveRDS(calibrated_share, file.path(work_folder, "best_reference_diet.rds"))
  }
  if (round_number == number_of_rounds) break

  # What each period asks for, as a log ratio, averaged over the periods that
  # have a figure for the species.
  log_ratio <- log(pmax(seal_diet_share, 1e-4)) -
    log(pmax(model_share[, named_species], 1e-4))
  average_log_ratio <- colMeans(log_ratio, na.rm = TRUE)
  new_share <- calibrated_share * exp(step_towards_measured * average_log_ratio)
  new_share <- pmin(pmax(new_share, starting_seal_diet_share / largest_adjustment),
                    starting_seal_diet_share * largest_adjustment)
  new_share <- new_share * (1 - share_of_unnamed_species_together) / sum(new_share)
  saveRDS(new_share, file.path(work_folder, reference_diet_file))
}

# Finish on the best calibration found, not simply the last one.
file.copy(file.path(work_folder, "best_reference_diet.rds"),
          file.path(work_folder, reference_diet_file), overwrite = TRUE)
run_script("8_seal_on_bowen_diet.R", "diet_round_8.log")
run_script("5_fit_seal_growth_and_harvest.R", "diet_round_5.log")
run_script("9_project_seal_on_bowen_diet.R", "diet_round_9.log")
cat(sprintf("finished on the best calibration: diet mismatch %.3f\n",
            diet_mismatch(readRDS(file.path(work_folder, "dynamic_seal_run.rds"))$model_diet_share)))
cat("shares calibrated to in the 1978-1981 state (%):\n")
print(round(100 * readRDS(file.path(work_folder, reference_diet_file)), 2))
