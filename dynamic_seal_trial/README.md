# Dynamic seal trial (working scripts, not part of the pipeline)

These three scripts replaced the Seals gear in the FishMIP recalibrated model
(`Ess_Fishmip_data/params_GFDL_steady.RDS`) with Grey Seal as a dynamic species,
following Paper 2's recipe. They are kept so the trial can be rerun or picked
up; they are **not** cleaned up to the standard of the numbered scripts, and
they were written to be run once, in order.

What needs changing before they will run anywhere else:

- they read and write intermediate files through the environment variable
  `SC` (a scratch folder); set it to any writable folder first;
- they read `paper2_base_LOCKED.rds`, `seal_abundance_digitised.RDS` from
  absolute paths under `~/Desktop/Ez June 2026/`;
- `3_project_and_score.R` reads the fish-eating calibration from
  `seal_qB.rds`; the line that creates it
  (`qB <- make_fish_eating(C$params)`) was swapped for a `readRDS` after the
  first successful run and needs swapping back for a fresh run.

Results and the reasons the trial does not fit are in handover section 20.
Three things learned the hard way are recorded there too, because they will
bite anyone adding a species to a therMizer model.

## Later scripts (handover section 21)

- `4_revised_seal_prescribed_numbers.R`: seal eats fish and outside food only,
  preferred prey about 200 g, numbers set from the counts each year. Fits the
  fish well but for the wrong reason (prey four times too large).
- `5_fit_seal_growth_and_harvest.R`: seal-only model; fits an unharvested
  growth rate, a harvest rate for 1967-1990 and a recruitment ceiling to the
  digitised counts.
- `6_data_based_seal_exponential_harvest.R`: prey size fitted to the measured
  prey weights, consumption target from the Bowen & Harrison diet shares, seal
  numbers from the fitted demography. Needs the outputs of 2 and 5.
