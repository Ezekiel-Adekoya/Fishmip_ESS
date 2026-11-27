# 10. Make the seal's diet match Bowen & Harrison in the years they sampled
#
# Bowen & Harrison (1994) collected their samples from July 1991 to January
# 1993. Script 8 calibrates the seal in the 1978-1981 state, and by 1991-1993
# the fish community has changed, so the diet there is no longer what they
# measured. This script repeats
#
#   script 8 (calibrate the seal) -> script 5 (seal growth) -> script 9 (project)
#
# and after each round moves the shares that script 8 calibrates to, until the
# projected diet in 1991-1993 is the measured one. A species that is too small
# a share of the 1991-1993 diet gets a larger share in the calibrated state,
# and the other way round.
#
# The seal's fitted growth, harvest and pup limit are taken from an earlier
# run (SEAL_RATES_FROM) during the rounds, to save time. Run script 9 once
# more with FIT_SEALS_IN_FULL_MODEL=1 afterwards.
#
# Needs the same SC, BASE_PARAMS and SEAL_SOURCE settings as scripts 8 and 9.

work_folder <- Sys.getenv("SC")
largest_number_of_rounds <- as.integer(Sys.getenv("DIET_ROUNDS", "8"))
close_enough <- 0.01        # every named species within 1% of the diet
# No calibrated share may move further than this factor from the measured one.
# It stops a species that nearly vanishes from the model in 1991-1993 from
# taking over the calibrated diet.
largest_adjustment <- 4
reference_diet_file <- "reference_diet.rds"
Sys.setenv(REFERENCE_DIET_FILE = reference_diet_file)

run_script <- function(script, log_name) {
  status <- system2("Rscript", file.path("dynamic_seal_trial", script),
                    stdout = file.path(work_folder, log_name), stderr = file.path(work_folder, log_name))
  if (status != 0) stop(script, " failed; see ", file.path(work_folder, log_name))
}

measured_share <- NULL
for (round_number in 1:largest_number_of_rounds) {
  run_script("8_seal_on_bowen_diet.R", "diet_round_8.log")
  run_script("5_fit_seal_growth_and_harvest.R", "diet_round_5.log")
  run_script("9_project_seal_on_bowen_diet.R", "diet_round_9.log")
  run <- readRDS(file.path(work_folder, "dynamic_seal_run.rds"))
  measured_share <- run$bowen_share
  projected_share <- run$diet_in_sampling_years[names(measured_share)]
  largest_gap <- max(abs(projected_share - measured_share))
  cat(sprintf("round %d: diet in 1991-93, model (measured): %s | largest gap %.3f\n",
              round_number,
              paste(sprintf("%s %.1f (%.1f)", names(measured_share), 100 * projected_share,
                            100 * measured_share), collapse = ", "), largest_gap))
  if (largest_gap < close_enough) break

  calibrated_share <- if (file.exists(file.path(work_folder, reference_diet_file))) {
    readRDS(file.path(work_folder, reference_diet_file))[names(measured_share)]
  } else measured_share
  new_share <- calibrated_share * (measured_share / pmax(projected_share, 1e-6))^0.8
  new_share <- pmin(pmax(new_share, measured_share / largest_adjustment),
                    measured_share * largest_adjustment)
  # The named species keep the same total as measured, so the unnamed species
  # keep their 1.7% between them.
  new_share <- new_share * sum(measured_share) / sum(new_share)
  saveRDS(new_share, file.path(work_folder, reference_diet_file))
}
cat("shares calibrated to in the 1978-1981 state (%):\n")
print(round(100 * readRDS(file.path(work_folder, reference_diet_file)), 2))
