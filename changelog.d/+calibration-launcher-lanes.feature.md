The full calibration launcher now submits one job per suite, sized from the checkout's
`calibration_plan.tsv` (lanes x CPUs, memory and time per suite), instead of a throttled 4-CPU
array. It validates the plan first, refuses a launch whose CPUs would pass the shared-cluster cap
or whose longest suite would run into the nightly, and cancels every job already submitted when a
submission fails. It prints `jobs=<ids> rows=<id> log=<dir> cpus=<total>`.
