# Results

Running `run_two_stage_ga(seed)` writes two files to this directory:

- `latest_run.xlsx`: decoded UAV parameters and coverage durations.
- `latest_run.png`: the global and per-object coverage timelines.

Generated results are ignored by Git by default because they depend on the random seed and can be regenerated locally.

## Current baseline

The previously saved two-stage exploratory run reported `0.00 s` for the strict global objective, which requires all 33 lines of sight to be blocked simultaneously within the 67-second horizon.

This baseline is reported without embellishment. A zero result does not prove that no feasible positive solution exists; it only shows that the recorded run did not find one.
