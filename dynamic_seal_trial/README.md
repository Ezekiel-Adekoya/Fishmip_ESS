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

## Seal on the measured diets, on raised sand lance and herring stocks (handover sections 27-32)

Run in this order from the repository root. `SC` is any writable scratch
folder; the other variables are read by the scripts named.

```
export SC=/some/scratch/folder
export SAND_LANCE_TARGET_KT=1500
export BASE_PARAMS="$SC/base_params_raised_forage_fish.RDS"
export SEAL_SOURCE=seal_on_bowen_diet.rds
Rscript dynamic_seal_trial/7_fish_base_with_raised_forage_fish.R
Rscript dynamic_seal_trial/1_build_structure.R
Rscript dynamic_seal_trial/2_calibrate_with_seals.R
DIET_ROUNDS=4 Rscript dynamic_seal_trial/10_match_seal_diet_to_both_periods.R
FIT_SEALS_IN_FULL_MODEL=1 Rscript dynamic_seal_trial/9_project_seal_on_bowen_diet.R
Rscript dynamic_seal_trial/plot_seal_on_bowen_diet.R
```

- `7_fish_base_with_raised_forage_fish.R`: raises the sand lance and herring
  targets in steps, recalibrating the fish each time, keeps the herring catch
  at the observed one, sets their reproduction levels (0.9), and checks that
  the base stays put when left alone for 100 years. `HERRING_TARGET_KT`
  (default 700; 0 leaves herring alone), `SAND_LANCE_REPRODUCTION_LEVEL` and
  `HERRING_REPRODUCTION_LEVEL` change the settings.
- `seal_diet_data.R`: the two measured diets (Bowen & Harrison 1994 for
  1991-1993, and 2007 for 1993-1998) and the seven prey the seal is calibrated
  on one by one.
- `8_seal_on_bowen_diet.R`: calibrates the seal to 1.73 t of fish per seal a
  year in the proportions it is given (the average of the two measured diets,
  or the shares script 10 has written).
- `9_project_seal_on_bowen_diet.R`: projects 1861-2010 and reports the diet in
  both measured periods. With `FIT_SEALS_IN_FULL_MODEL=1` it first fits the
  seal's reproductive efficiency and harvest rate, and the limit on pups
  recruited a year, to the counts inside the full model. `SEAL_RATES_FROM`
  takes those three from an earlier run file instead.
- `seal_prey_switching.R`: prey switching for the seal, after Smout et al.
  (2014, ICES J. Mar. Sci. 71: 81-89). The seal's preference for a prey is
  multiplied by (availability now / availability in the calibrated state)^0.5,
  their exponent of 1.5 less 1. Sourced by script 9; `PREY_SWITCHING=0` turns
  it off, `NO_SWITCHING_FOR="Atlantic Cod"` turns it off for named prey, and
  `SWITCHING_EXPONENT` tries another exponent.
- `10_match_seal_diet_to_both_periods.R`: repeats scripts 8, 5 and 9, moving
  the diet script 8 calibrates to in small bounded steps so the projected diet
  comes closer to the measured one in both periods, and finishes on the best
  round. It runs scripts 8 and 5 itself.
- A check on any of these models must load therMizer as well as mizer.
  Without it the temperature effect is silently left out and a steady model
  looks as if it drifts.
- Scripts 1, 2 and 5 are the earlier ones, now reading their starting model
  from `BASE_PARAMS` and `SEAL_SOURCE`. They and scripts 8 and 9 still read
  Paper 2's seal parameters and the digitised seal counts from
  `~/Desktop/Ez June 2026/`.
