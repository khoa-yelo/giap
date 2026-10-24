# Exact Wallenius reference (M8-12b; ADR-0044 decision 1; GIAP_SPEC.md §7.7.5, §11.4)

- Command: `sbatch scripts/slurm/golden_wallenius.sbatch` (runs `scripts/golden/make_wallenius_exact.py --workers 4` twice).
- Slurm job: 1879871 on wheat-01, 2026-10-24 (03:12 UTC).
- Reason: M8-12b: the exact Wallenius reference (ADR-0044 decision 1)
- Scripts: git commit `352446e8a932c6cf5a6dcf9802ad6f23dcdac550`.
- Environment: Python 3.12.3, numpy 2.5.4, np.longdouble nmant 63, bits 128, machine x86_64.
- Checks (M8-12b §3.4):
  - R1 80-bit vs float64, all tails: max rel 7.021e-14 (limit 1e-13);
  - R1 80-bit vs R2 on the R2 tails (CAT_003 over, CAT_004 over, CAT_010 over, CAT_011 over, CAT_012 over, CAT_013 under): max rel 0.000e+00 (limit 1e-13);
  - R2 spot checks on 10 tails >= 1e-3 (limit 1e-13): CAT_008 over rel 0.000e+00; CAT_009 under rel 0.000e+00; CAT_019 under rel 0.000e+00; CAT_047 under rel 0.000e+00; CAT_050 under rel 0.000e+00; CAT_073 under rel 0.000e+00; CAT_079 under rel 0.000e+00; CAT_082 under rel 0.000e+00; CAT_123 under rel 0.000e+00; CAT_131 under rel 0.000e+00;
  - Σf = 1 by R2 on 3 categories (limit 1e-50): CAT_019 |Σ - 1| 0.000e+00; CAT_031 |Σ - 1| 0.000e+00; CAT_143 |Σ - 1| 0.000e+00;
  - vacuity floor (ADR-0044 decision 4): {'1e-20_to_1e-8': 2, 'below_1e-20': 2};
  - ω max ulp difference against the numpy mean ratio: 3.0.
- Inputs: `tests/data/golden/r/inputs/ora_genes.tsv.gz` sha256 `e94df51857577faf669c1e59e16454eb911ffdc7149d0c38b8956683254b1023`;
  `tests/data/golden/r/inputs/ora_categories.tsv.gz` sha256 `114fcef698799d06250a8561e14022d9ea83af0cc6602193b48c69b35f9fef34`.
- Method: R1 (urn recursion in 80-bit long double) for tails >= 1e-3; R2 (the exact rational
  integral, fixed point with >= 200 surviving bits) for tails < 1e-3; the larger tail is
  1 - small + pmf(k); ω is the exact rational mean ratio, correctly rounded once to float64.
  The reference is ADR-0044 decision 1. No goseq, BiasedUrn or scipy code is used.
- Reproducibility: two runs in this job were byte-identical.
- Regeneration only through this job and only with a reason; a change to the ORA inputs requires
  it. CI fails a change to `tests/data/golden/wallenius/` that does not update this file.
