# M8-12b: Median-exact golden GSEA inputs, the R rerun, and the exact Wallenius reference

Milestone M8 · track G (golden data) · **pass 2** · depends on: M8-12 (#192, merged), ADR-0038,
ADR-0044 (accepted 2026-10-23) · blocks: M8-13 G1 exact-score equality, example 9, G3a, the G3b
golden comparisons, and A1 (Wallenius) · estimate ≈ 270 code lines (Python ≈ 200, shell ≈ 70)
(+ ≈ 200 test lines) · **statistical** (the exact reference) · **Slurm** (two short jobs plus one
test job)

## Spec sections implemented

- §11.5: "`scripts/golden/make_golden_r.R` … writes: fgsea multilevel results … for frozen ranked
  TU scores and gene sets; goseq Wallenius p-values …". M8-12 items 1–4 are kept, with one
  changed property of the GSEA inputs: their median is **exactly 0.0**.
- §7.7.5 criterion 1, as overridden by **ADR-0044** decision 1: the reference becomes "the exact
  Wallenius distribution", committed as `tests/data/golden/wallenius/exact_tails.tsv.gz` by
  `scripts/golden/make_wallenius_exact.py`. This uses pure Python, runs once on Slurm (≤ 4 CPUs,
  ≤ 1 h) and has a PROVENANCE entry.
- §11.4 (the recorded-output provenance template), R3 (`tolerances.yaml` is protected; this
  issue changes no key), R11 (≤ 5 MB per file, gzip text fixtures), R14 (job sizing).
- M8-13 lead decision 1: "regenerate the golden GSEA inputs with an exactly-zero median, ties and
  zeros kept, ORA files byte-identical, A1 rechecked."

**Out of scope.** The `golden_r.wallenius_goseq_abs_tol: 1.0e-5` key (ADR-0044 decision 2) goes
into `tests/calibration/tolerances.yaml`. That file is protected, so **the lead** commits the key
with the trailer `Gate-Change-Approved: ADR-0044`, in the M8-13 test PR (ADR-0044 "Consequences").
No file of this issue reads that key. The M8-13 tests themselves are not part of this issue either.

## Design in one paragraph

The GSEA scores have 6 exact zeros. So the median of 2,000 values is exactly 0.0 **if and only
if** exactly 997 scores are negative and 997 are positive. The generator gets there by negating
the |d| non-special scores closest to zero on the majority side, where 2d = n_neg − n_pos. At
`GOLDEN_SEED` d = 8, so 8 rows of `gsea_stats` change, each by ≤ 0.0175. No RNG draw is added,
and the ORA stream is a separate `SeedSequence` child, so the ORA files keep their bytes.
`golden_r.sbatch` then reruns once (about 1 min of R) and gains a guard: it fails if the ORA
inputs or `goseq.tsv.gz` change. The exact Wallenius reference is a new pure-Python generator:
- the urn recursion (R1) in 80-bit long double for every category;
- the exact rational integral (R2), in fixed-point integers with ≥ 200 surviving bits, for the
  6 tails below 1e-3;
- the complementary tail from the identity Σf = 1.

It takes about 25 min serial, or about 8 min on 4 workers. Each golden directory then has its
own PROVENANCE, enforced in CI.

## 1. Making the median exactly 0.0

### 1.1 What is true now

`make_golden_inputs.py` works in four steps:
1. it draws N(0, 1) and subtracts its median;
2. it rounds to integer micro-units;
3. it overwrites 106 random TUs (the "special" ones, `frozen`): 50 at +0.612345, 50 at −0.734567
   and 6 at 0;
4. `_make_unique` nudges every other value up by 1 µ until it is distinct, non-zero and not a tie
   value.

The overwrite is not sign-balanced, so the committed file has **1,005 negative, 6 zero and 989
positive** scores. Its median is (−0.005775 + −0.005635)/2 = −0.005705 (M8-13 F1).

### 1.2 The condition, and why "±x" cannot work here

Let N = 2000 (even) and let s₀ ≤ … ≤ s₁₉₉₉ be the sorted scores. `numpy.median` returns
(s₉₉₉ + s₁₀₀₀)/2, computed as the mean of the two middle elements.

- **±x is impossible while zeros exist.** Suppose s₉₉₉ = −x and s₁₀₀₀ = +x with x > 0. Every zero
  lies strictly between them in sorted order, but they are adjacent. So with Z = 6 zeros, a
  median of exactly 0.0 needs **both middle order statistics to be 0**.
- **Both middle values are 0 ⇔ n_neg ≤ 999 and n_neg + Z ≥ 1001,** that is n_neg ∈ [995, 999].
  The generator takes the symmetric point, **n_neg = n_pos = (N − Z)/2 = 997**. The zeros then
  sit at sorted positions 997–1002, so the median is 0.0 even if one middle value moves by up to
  2 positions under a future edit.
- **Exact in float64:**
  - every score is parsed from a `d.dddddd` string, and `micro_str(0)` writes `0.000000`, so every
    zero is +0.0 (never −0.0);
  - s₉₉₉ = s₁₀₀₀ = +0.0, and (+0.0 + +0.0) = +0.0 exactly, as is +0.0 / 2;
  - every median convention gives +0.0: lower, upper, midpoint, `statistics.median`,
    `np.quantile(·, 0.5)` with linear interpolation (a + (b − a)·t with a = b = 0), and R's
    `median`;
  - `giap.function.gsea.tu_ranking` computes `means − np.median(means)` = `means − 0.0`, and
    x − 0.0 == x bitwise for every finite x, including +0.0. With one gene per TU (the M8-13
    fixture), the TU mean is the score itself (a mean over one float is that float), so GIAP sends
    R **exactly** the golden scores.
- **General rule** (used in the toy examples): with N even you need Z ≥ 2 zeros and
  |n_neg − n_pos| ≤ Z − 2. With N odd you need Z ≥ 1 and |n_neg − n_pos| ≤ Z − 1.

### 1.3 The generator change (`scripts/golden/make_golden_inputs.py`, ≈ 25 lines)

In `gsea_inputs`, after `_make_unique` and before `values = micro / SCALE`, add a new pure
function `_balance_signs(micro, frozen, taken) -> np.ndarray`. `_make_unique` itself is unchanged.

1. n_neg = #(micro < 0) and n_pos = #(micro > 0) over all 2,000 TUs; the frozen tie values add
   50 each. n_neg + n_pos = 1994, which is even, so d = (n_neg − n_pos) / 2 is an integer.
2. If d ≠ 0, the candidates are the **non-frozen** TUs on the majority side (negative when d > 0).
   Order them by (|micro|, index) ascending, stably, and take the first |d|.
3. Process those TUs in that order. Set new = −micro. While new is a value currently held by
   another TU (any TU, frozen or not), step new **away from zero** by 1 µ (`new += sign(new)`).
   The step never reaches 0 or crosses it, so the sign is kept. Then assign it.
4. **Post-conditions.** The script exits non-zero (`SystemExit` with the failing condition) unless
   all of these hold:
   - n_neg == n_pos == 997 and #(micro == 0) == 6;
   - 50 TUs equal 612345 and 50 equal −734567;
   - the non-frozen values are distinct, non-zero and not tie values;
   - `float(np.median(micro / SCALE)) == 0.0`, an **exact** `==`. This replaces the old
     `abs(median) > 0.05` check.
5. Everything after that is unchanged. `values = micro / SCALE` feeds the set planting (`sd`, the
   tilted weights) as before, and the docstring says "median exactly 0.0 (M8-12b)".

**Effect at GOLDEN_SEED = 50000** (checked by the spec-writer on a scratch copy; this is not
committed code). d = 8, and the 8 negatives closest to zero are negated:

| TU | before | after |
|---|---|---|
| TU_0807 | −0.000554 | **0.000555** (+0.000554 is taken: one step away from zero) |
| TU_1362 | −0.001995 | 0.001995 |
| TU_0809 | −0.005135 | 0.005135 |
| TU_1686 | −0.005627 | 0.005627 |
| TU_1539 | −0.005635 | 0.005635 |
| TU_0883 | −0.005775 | 0.005775 |
| TU_1339 | −0.007387 | 0.007387 |
| TU_0108 | −0.008749 | 0.008749 |

- The sorted scores 995–1004 are −0.011526, −0.011246, 0 × 6, 0.000554, 0.000555.
- The tie blocks and the zeros are the **same TUs** as before, because the special draw comes
  first and is unchanged.
- `ora_genes.tsv.gz` and `ora_categories.tsv.gz` are byte-identical.
- `gsea_sets.tsv.gz` changes: about 110 of the 120 sets get different members. `rng.choice(p = w)`
  consumes a different number of draws once w moves. SET_026 (30 names, 25 matching) and SET_027
  (12 names, 8 matching, so dropped by minSize) are structural and stay the same.
- The planted |z| still span about 2.3 to 13; 24 sets have |z| ≥ 2.

### 1.4 Determinism

- No RNG draw is added or moved. `_balance_signs` is a deterministic function of `micro` and
  `frozen`, and its sort is stable.
- The ORA inputs come from the second `SeedSequence(seed).spawn(2)` child, so their bytes do not
  depend on any GSEA change. That is why they stay identical (sha256 below).
- Two runs give identical bytes: gzip mtime 0, no file name, as in M8-12 example 3.

## 2. The R rerun (`scripts/slurm/golden_r.sbatch`)

### 2.1 Changes to the job (≈ 20 lines of shell)

1. **ORA guard.** Before `rm -rf "$DEST/inputs" "$DEST/outputs"`, the job compares the new
   `ora_genes.tsv.gz`, `ora_categories.tsv.gz` and `outputs/goseq.tsv.gz` with the committed ones
   (`cmp`). If they differ it exits 1, unless `GOLDEN_R_ALLOW_ORA_CHANGE=1` is set. A first run with
   no committed files passes. This makes "ORA byte-identical" a property of the job and not only of
   a test.
2. **Reason and lineage in PROVENANCE.md.** The heredoc gains two lines:
   - `- Reason: $GOLDEN_R_REASON`, which is required. The job fails if it is unset;
   - `- Supersedes: Slurm job <id>`, where the id is read from the old PROVENANCE.md before it is
     overwritten, or `none`.

   Every other line keeps its current form, so `test_provenance_records_the_run` still passes.

The job keeps building nothing new. `giap-golden-r` exists and its `.giap-built` marker matches
the yml, so there is no rebuild. The lock re-export must come out byte-identical
(sha256 `eb4448ff…0b7c`).

### 2.2 Running it

- **Order:**
  1. commit and push the generator and sbatch changes first. PROVENANCE records
     `git rev-parse HEAD` and marks "(scripts modified after this commit)" if the scripts are
     dirty, which must not happen;
  2. sum the CPUs of running and pending `giap-*` jobs (R14, ≤ 128);
  3. submit.
- **Command:**
  `GOLDEN_R_REASON="M8-12b: median-exact GSEA inputs (M8-13 lead decision 1)" sbatch --export=ALL --cpus-per-task=2 --mem=8G --time=00:20:00 scripts/slurm/golden_r.sbatch /home/classes/bios/270/khoa/giap/local/worktrees/m8-12b`
- **Sizing.** M8-12's run (job 1868639) took 42 s, single-threaded R with `nproc = 1`, two R runs
  plus two input generations. So 2 CPUs, 8 GB and 20 min leave ≥ 10× headroom. The script header
  (4 CPUs, 16 GB, 2 h) is not edited; the command-line flags override it.
- Record the job ID in PROGRESS.md.

### 2.3 Expected outputs

| File | Expected |
|---|---|
| `inputs/gsea_stats.tsv.gz` | changed in exactly 8 rows (the table in §1.3) |
| `inputs/gsea_sets.tsv.gz` | changed (new memberships, same 120 set IDs and sizes rule) |
| `inputs/ora_genes.tsv.gz` | byte-identical, sha256 `e94df51857577faf669c1e59e16454eb911ffdc7149d0c38b8956683254b1023` |
| `inputs/ora_categories.tsv.gz` | byte-identical, sha256 `114fcef698799d06250a8561e14022d9ea83af0cc6602193b48c69b35f9fef34` |
| `outputs/fgsea.tsv.gz` | regenerated, 119 rows (SET_027 dropped) |
| `outputs/calc_gsea_stat.tsv.gz` | regenerated, 120 rows |
| `outputs/goseq.tsv.gz` | byte-identical, sha256 `293f78ac37c55b32849913c88f87fa45df2aa1c2d0e18bdb536c2ac5e906b879` (same inputs, lock and `set.seed`) |
| `VERSIONS.txt` | expected byte-identical (`a07134…2e4e`). A difference only in the warnings section is acceptable and is noted in the PR; any other difference stops the work (env drift) |
| `PROVENANCE.md` | new job ID, date, commit, Reason and Supersedes: 1868639 |
| `scripts/golden/golden_r.explicit.txt` | unchanged |

## 3. The exact Wallenius reference (`scripts/golden/make_wallenius_exact.py`, ≈ 170 lines)

### 3.1 Interface

```
python scripts/golden/make_wallenius_exact.py [--inputs DIR] [--out FILE] [--workers INT]
                                              [--categories CAT_a,CAT_b,...] [--log FILE]
```

- `--inputs` defaults to `tests/data/golden/r/inputs`, of which it reads only `ora_genes` and
  `ora_categories`. `--out` defaults to `tests/data/golden/wallenius/exact_tails.tsv.gz`, and
  `--workers` to 1.
- `--categories` restricts the run to a subset, for tests. The rows written for a category do not
  depend on which other categories run, or on `--workers`.
- `--log` writes the check summary of §3.4 as JSON.
- Imports: the standard library (`fractions`, `math`, `multiprocessing`, `gzip`, `csv`, `json`,
  `argparse`) and `numpy`, for long double only. **No scipy, no R, no `giap`, and no goseq or
  BiasedUrn code**, used or read (ADR-0044). The module exposes pure functions for tests:
  - `odds(...)`;
  - `r1_pmf(M, n, N, omega) -> np.ndarray[longdouble]`;
  - `r2_pmf(x, M, n, N, omega) -> Fraction`;
  - `exact_tails(k, M, n, N, omega) -> (Fraction over, Fraction under, dict method)`.

### 3.2 Parameters per category

- M = 3000 (every gene of `ora_genes`, `use_genes_without_cat = TRUE`), N = Σ detected (1,456),
  n = |C| and k = #detected in C.
- **ω.** The `p_detect` strings are parsed as exact decimals (`Fraction("0.123456")`). Then
  ω = mean_in / mean_out = (Σ_C p / n) / ((Σ_all p − Σ_C p)/(M − n)), computed exactly and then
  **correctly rounded once to float64**. The result is a dyadic rational, which R2 needs.
  - It can differ from GIAP's float ω (numpy means) by a few ulp. The effect on a tail is about
    k · 4e-16 relative, < 1e-13, against the gate's 1e-6 and the M8-13 test's 1e-12 R1 check.
    The generator logs the maximum ulp difference against the numpy formula.
  - This is question 3.

### 3.3 Algorithm

Support: lo = max(0, N − (M − n)) and hi = min(n, N). The tails are over = P(X ≥ k) and
under = P(X ≤ k), as in goseq and in `tests/helpers/golden_r_gate.py::urn_tails`.

1. **R1, for every category: the urn recursion in `np.longdouble`.**
   - It runs over the N draws; a red ball is drawn with probability ω·r/(ω·r + w). It gives the
     pmf over [lo, hi].
   - The generator asserts `np.finfo(np.longdouble).nmant >= 63`, which means 80-bit x87 on
     x86_64, and refuses to run otherwise.
   - R1 tails are **exact rational sums** of the long-double pmf terms
     (`Fraction(*v.as_integer_ratio())`), rounded once to float64.
   - It also runs in float64, for the log only.
   - Cost: N × (hi − lo + 1) vector steps. All 150 categories took 8.5 s on iron-13.
2. **Regime choice.** s = min(over_R1, under_R1).
   - If s ≥ 1e-3, **both tails are R1** (80-bit).
   - If s < 1e-3, the smaller tail is **R2** (step 3) and the larger is
     **1 − small + pmf(k)**. Here pmf(k) is the R2 term at x = k, which the small tail already
     includes, and the identity holds because Σf = 1.
   - Ties at exactly 1e-3 go to R1, as the strict `<` says. At most one tail can be < 1e-3,
     because over + under = 1 + pmf(k) ≥ 1. On the golden data 6 tails are below 1e-3:
     - CAT_003 over (1.50e-4);
     - CAT_004 over (4.17e-7);
     - CAT_010 over (6.07e-46);
     - CAT_011 over (4.44e-27);
     - CAT_012 over (2.23e-10);
     - CAT_013 under (1.06e-17).
3. **R2, the exact rational integral (ADR-0044):**
   f(x) = C(m₁,x) C(m₂,N−x) · D · Σⱼ₌₀ˣ (−1)ʲ C(x,j) B(ωj + D, N−x+1), where
   D = ω(m₁ − x) + m₂ − N + x and B(a, m+1) = m! / ∏ᵢ₌₀ᵐ (a + i).
   - **Integers only.** Write ω = W/2ᵉ. Then a + i = Aⱼ,ᵢ/2ᵉ with
     Aⱼ,ᵢ = Wj + W(m₁ − x) + 2ᵉ(m₂ − N + x + i), so each term is
     C(x,j) · m! · 2^{e(m+1)} / ∏ᵢ Aⱼ,ᵢ.
   - **Fixed point.** Tⱼ = ⌊C(x,j) · m! · 2^{e(m+1)} · 2ᴾ / ∏ᵢ Aⱼ,ᵢ⌋ and S = Σ(−1)ʲ Tⱼ. The
     absolute error is ≤ x + 1 units of 2⁻ᴾ.
   - **Precision rule (self-verifying).** Start with P = 264 + 2x. While S ≤ 0 or
     bitlen(S) < 200 + bitlen(x + 1) + 1, raise P and recompute. On exit the relative error of
     f(x) is < 2⁻¹⁹⁹. On the golden tails this needed P ≈ 3,200–3,530.
   - f(x) is the `Fraction` C(m₁,x) C(m₂,N−x) · D · S / 2ᴾ. A tail is the exact `Fraction` sum
     over its support points, rounded once to float64 (`float(Fraction)` is correctly rounded).
   - Any exact rearrangement is allowed, such as the colour swap X ↔ N − X with 1/ω, if the
     same surviving-bits rule holds.
4. **Writing.** Each value goes through `format(p, ".17g")` (round-trip exact). Rows are in
   `ora_categories` file order (CAT_001 … CAT_150), each category's `over` row before its `under`
   row. The header is `category\ttail\tp_exact`, with LF line ends, written by the same gzip
   writer as M8-12 (mtime 0, no file name, level 9). That gives 301 lines and a file of about
   5 kB.
5. **Workers.** The unit of parallel work is one (category, x) R2 point, or one category for R1.
   Results are exact integers and fractions, combined in a fixed order, so the bytes do not depend
   on `--workers` or on scheduling.

### 3.4 Checks the generator runs (logged, and summarised in PROVENANCE)

| Check | Rule | Fails the run if |
|---|---|---|
| R1 80-bit vs R1 float64, all 300 tails | max rel (ADR-0044: 3.8e-15) | > 1e-13 |
| R1 80-bit vs R2 on the 6 R2 tails | max rel (prototype: equal to all 17 digits) | > 1e-13 |
| **R2 spot checks**, ADR-0044: ≥ 10 tails ≥ 1e-3 | the 10 cheapest by Σ over the tail of (x + 1)(N − x + 1), ties by (category, tail); R2 tail vs R1 tail | rel > 1e-13 |
| **Σf = 1**, ADR-0044: ≥ 3 categories | R2 over the full support of the 3 categories with the fewest support points (ties by ID) | \|Σ − 1\| > 1e-50 |
| every exact tail | 0 < p ≤ 1 and over + under ≥ 1 | any violation; a tail < 1e-300 is refused rather than written as 0 |
| vacuity floor (ADR-0044 decision 4) | ≥ 2 tails in (0, 1e-20), ≥ 1 in [1e-20, 1e-8) | not met |

### 3.5 Runtime (measured on iron-13 with a scratch prototype of step 3)

- R1 for all categories: 8.5 s.
- R2 costs 6–12 s per support point at x ≈ 80–300 and N − x ≈ 1,150–1,450:
  - CAT_010: 1 point, 8.7 s;
  - CAT_012: 21 points, 148 s;
  - CAT_013: 41 points, 72 s;
  - CAT_004: 64 points × about 12 s, roughly 13 min;
  - CAT_011: 31 points × about 7 s, roughly 4 min;
  - CAT_003: 21 points, about 2 min.
- The 6 tails take about **25 min serial**, or **about 8 min with 4 workers**. The spot checks
  and the Σf checks are on categories with ≤ 11 support points and x ≤ 10, so they take seconds.

### 3.6 The job (`scripts/slurm/golden_wallenius.sbatch`, ≈ 50 lines)

- `#SBATCH --job-name=giap-golden-wallenius`, 4 CPUs, 4 GB, `--time=01:00:00` (ADR-0044: ≤ 4 CPUs,
  ≤ 1 h). It takes one optional argument, the checkout, and requires `SLURM_JOB_ID`.
- It runs the generator **twice** with `--workers 4` into `$TMPDIR` and `cmp`s the two files,
  about 16 min. Then it copies the file to `tests/data/golden/wallenius/` and writes
  `PROVENANCE.md` (§4.2).
- It uses `$UV_PROJECT_ENVIRONMENT/bin/python` and never `uv run`.

## 4. PROVENANCE

### 4.1 `tests/data/golden/r/PROVENANCE.md`

It is rewritten by the job (§2.1). The new fields are `Reason` and `Supersedes: Slurm job 1868639`.
Everything else keeps the M8-12 template.

### 4.2 `tests/data/golden/wallenius/PROVENANCE.md` (new, written by the job)

- **Command, job and scripts:**
  - the command and `sbatch scripts/slurm/golden_wallenius.sbatch`;
  - the Slurm job ID, host and UTC date;
  - the git commit of `scripts/golden/make_wallenius_exact.py`, with the "(scripts modified after
    this commit)" marker when the scripts are dirty.
- **Environment:** the Python and numpy versions, `np.finfo(np.longdouble)` (nmant 63) and the
  machine.
- **Inputs:** `ora_genes.tsv.gz` and `ora_categories.tsv.gz` with their sha256, which ties the
  reference to those exact ORA bytes.
- **Method:**
  - R1 is used for tails ≥ 1e-3;
  - R2 is used for tails < 1e-3, with ≥ 200 surviving bits;
  - the larger tail is 1 − small + pmf(k);
  - ω is the correctly rounded rational mean ratio;
  - the reference is ADR-0044 decision 1;
  - the generator uses no goseq, BiasedUrn or scipy code.
- **Checks:** the §3.4 numbers (the max relative differences, the 10 spot-check tails by name,
  the 3 Σf categories and \|Σ − 1\|, and the vacuity counts).
- **Reproducibility:** two runs in the job were byte-identical.
- **Regeneration:** only through this job and only with a reason. A change to the ORA inputs
  requires it.

### 4.3 CI: `scripts/ci/check_golden_provenance.py`

`GOLDEN_R` becomes a tuple, `("tests/data/golden/r/", "tests/data/golden/wallenius/")`. A changed
file under either directory needs **that** directory's `PROVENANCE.md` in the same diff.
`violations()` keeps its signature. The existing cases still hold, for example that
`tests/data/golden/rx/a.tsv` is not golden. The CI step name in `.github/workflows/ci.yml` becomes
"golden data carry a PROVENANCE.md update (M8-12, M8-12b)".

## 5. What else changes

| Item | Change | Owner |
|---|---|---|
| `tests/unit/test_golden_r_inputs.py` | **No edit needed**: `abs(median) ≤ 0.05`, the ties, zeros and planted-z checks all still pass (checked on the scratch output). The stronger rule goes in a new test file (§Test plan) | — |
| `tests/golden/test_golden_r_files.py` | no edit; must pass on the new files, including `test_committed_inputs_equal_a_fresh_generation`, `min p ≤ 1e-25`, a p in (1e-3, 1e-2], NA `log2err` ≤ 5 %, and calc_gsea_stat vs fgsea ≤ 1e-12 | — |
| `tests/unit/test_golden_r_boundary.py` | no edit; new cases for the wallenius directory in a new test | test-author |
| **A1 (#204)**, `tests/contract/test_fgsea_r_contract.py` (`test_golden_contract_strings`, `…_through_wrapper`, `test_seed_is_set_before_every_job`) | no edit; must pass on the new `fgsea.tsv.gz` against the live `giap-fgsea` env (Slurm). This is the "A1 recheck" of lead decision 1 | lead runs, records the job |
| M8-13 branch (`m8/golden-equivalence`) | after this merges, it merges main. The skipped G1 exact equality, example 9 and G3a get written. Its numbers ("9 fgsea sets with p < 1e-10 today") are re-read from the new file | M8-13 test-author |
| `tolerances.yaml` | none here (lead, M8-13 PR, `Gate-Change-Approved: ADR-0044`) | lead |
| M8-index row 12 | add "12b" (lead's edit) | lead |

## Examples (input → expected output)

| # | Input | Expected |
|---|---|---|
| 1 | toy: N = 8, two frozen zeros, free micro `[-6, -3, -2, -1, 1, 5]` | n_neg 4, n_pos 2, d = 1. −1 is flipped to +1, which collides with 1, so it steps to +2: `[-6, -3, -2, 0, 0, 1, 2, 5]`. The median is (0 + 0)/2 = 0.0 exactly |
| 2 | toy: N = 6, one zero, values `[-2, -1, 0, 1, 2, 3]` (n_neg 2, n_pos 3) | median (0 + 1)/2 = 0.5. No sign flip can fix it: N − Z = 5 is odd, so the middle pair cannot both be 0 (§1.2 general rule, Z ≥ 2 for even N). The golden data have Z = 6, so this case never arises; `_balance_signs` raises `SystemExit` if N − Z is odd |
| 3 | `make_golden_inputs.py --seed 50000` | `np.median == 0.0` (`==`); 997 negative, 6 zero and 997 positive; sorted positions 997–1002 are 0; the 8 rows of §1.3 changed; the ORA sha256 values of §2.3 |
| 4 | the same command, run twice | byte-identical four files |
| 5 | the generator with the balancing step removed (mutation) | it exits non-zero with "median −0.005705 != 0.0" or "n_neg 1005 != 997" |
| 6 | `golden_r.sbatch` with the ORA RNG stream disturbed (mutation) | the job exits 1 at the ORA guard before touching `tests/data/golden/r/` |
| 7 | `r2_pmf` with M = 4, n = 2, N = 2, ω = 2 | pmf(0) = 1/15, pmf(1) = 3/5, pmf(2) = 1/3 (urn: (4/6)(2/4) = 1/3; (2/6)(1/5) = 1/15), as exact `Fraction`s once rounded to floats. `exact_tails(k = 2)` gives over 1/3 and under 1 |
| 8 | `r2_pmf` with M = 3, n = 1, N = 1, ω = 2; and any small urn with ω = 1 | P(X = 1) = 2/(2 + 2) = 1/2. With ω = 1 it equals the central hypergeometric C(n,x)C(M−n,N−x)/C(M,N) to ≤ 2⁻¹⁹⁰ relative |
| 9 | the generator on the golden ORA inputs | `p_exact` values: CAT_010 over `6.0681368749300294e-46` and under `1`; CAT_012 over `2.2283895161517117e-10`; CAT_013 under `1.0597943893784124e-17`; CAT_004 over `4.1695578881156134e-07` (R2 = R1 to all 17 digits in the prototype; ADR-0044 lists 4.169557888115581e-07) |
| 10 | a diff touching `tests/data/golden/wallenius/exact_tails.tsv.gz` without its PROVENANCE.md | `check_golden_provenance.py` exits 1 naming the file; with the PROVENANCE.md it exits 0 |


## Degenerate and edge cases

| Case | Rule |
|---|---|
| circular or linear replicons, minus strand | not applicable: the golden inputs have no coordinates. The M8-13 fixture puts genes on linear replicons, and nothing here changes that |
| ties (2 values × 50 TUs, 6 zeros) | kept, on the same TUs. `_balance_signs` never touches frozen TUs, and no flipped value may land on a tie value or 0 (the collision step moves away from zero and skips taken values) |
| d = 0 (already balanced, another seed) | no flip; the post-conditions still run |
| \|d\| larger than the number of majority-side free TUs | impossible for N = 2000 (it needs > 1,894 on one side); `SystemExit` if it ever happens |
| flipped value collides repeatedly | step away from zero until free (TU_0807 needs one step at GOLDEN_SEED); it terminates because at most 2,000 values are taken |
| another `--seed` (tests use `GOLDEN_SEED + 1`) | the same balancing, so the median is 0.0 for every seed. The post-condition guarantees it or fails loudly |
| −0.0 | cannot occur: `micro_str` writes `0.000000`, and negation is applied to non-zero integers only |
| empty input (no categories, or a category with n = 0) | the Wallenius generator exits non-zero. It never writes an empty or partial file (it writes a temp file and renames it) |
| a category with ω undefined (mean_out = 0 or n = 0) or k outside [lo, hi] | exit non-zero, naming the category (none in the golden data) |
| k = hi or k = lo | under = 1 (k = hi) or over = 1 (k = lo) exactly; the other tail is pmf(k) |
| an exact tail < 1e-300 | refused, not written as 0 (none; the smallest is 6.07e-46) |
| `np.longdouble` not 80-bit (aarch64 quad, or macOS arm64 = float64) | refuse to run, with the reason; R1's ≥ 1e-3 values would otherwise depend on the platform |
| `golden_r.sbatch` sees changed ORA bytes or `goseq.tsv.gz` | exit 1 before writing; `GOLDEN_R_ALLOW_ORA_CHANGE=1` is the only override, with a stated reason |
| `GOLDEN_R_REASON` unset | the job exits 1 before running R |
| new fgsea outputs fail an M8-12 range test or the M8-13 A5 floor (≥ 5 sets with p < 1e-10, ≥ 1 tie-insensitive one) | **stop and report to the lead.** Do not reseed (`GOLDEN_SEED` is protected) and do not tune the planting silently |

## Acceptance criteria (each checked by a test or a recorded command)

- [ ] **B1** `make_golden_inputs.py` at GOLDEN_SEED gives `np.median(scores) == 0.0` exactly,
      with 997 / 6 / 997 signs, zeros at sorted positions 997–1002, 50 + 50 tie values on the same
      TUs as before, and exactly the 8 changed rows of §1.3 (examples 1, 3, 5).
- [ ] **B2** Two runs are byte-identical (example 4). The ORA files equal the pinned sha256
      values; `ora_*` bytes do not change for any `--seed` change confined to the GSEA stream
      (separate-stream test).
- [ ] **B3** The new files pass the unedited `tests/unit/test_golden_r_inputs.py` and
      `tests/golden/test_golden_r_files.py`, including fresh-generation equality.
- [ ] **B4** `golden_r.sbatch` ran on Slurm with the job ID in PROGRESS.md. `goseq.tsv.gz`,
      `VERSIONS.txt` (apart from the warnings caveat) and the lock are unchanged, and the ORA guard
      and `GOLDEN_R_REASON` behave as in examples 6 and §2.1.
- [ ] **B5** A1 (#204) passes on the new `fgsea.tsv.gz` on Slurm, with the job ID recorded:
      `sbatch scripts/slurm/test_branch.sbatch <worktree> tests/contract/test_fgsea_r_contract.py -m requires_r -n 2`.
- [ ] **B6** `exact_tails.tsv.gz` exists: 300 rows plus the header `category tail p_exact`,
      150 categories × {over, under} in file order, 0 < p ≤ 1, over + under ≥ 1, ≤ 5 MB, gzip. It is
      read by `golden_r_gate.read_exact_tails` without error (the M8-13 layout contract).
- [ ] **B7** The pure functions satisfy examples 7–8. They hold Σf = 1 on small urns (exact) and
      the colour symmetry with dyadic 1/ω (ω ∈ {2, 4, 0.5}). R1 and R2 agree to ≤ 1e-14 on random
      small urns (M ≤ 60), and the bytes are the same for `--workers 1` and `--workers 2` on
      `--categories CAT_010` (example 9).
- [ ] **B8** The §3.4 checks pass, and their numbers are in the wallenius PROVENANCE.md together
      with the job ID, commit, input sha256 and the longdouble description.
- [ ] **B9** `check_golden_provenance.py` covers both directories (example 10), and the existing
      boundary tests still pass.
- [ ] **B10** The generator imports only the standard library and numpy (an AST test like
      `test_numpy_only_and_no_giap_import`), and `check_imports.py` passes.
- [ ] **B11** No change to `src/`, `tolerances.yaml`, `GIAP_SPEC.md` or any existing test file.

## Test plan (test-author writes; new files only, because existing tests are not edited)

1. **`tests/unit/test_golden_inputs_median_exact.py`** (CI, ≈ 60 lines). It runs the generator
   into tmp at GOLDEN_SEED and at GOLDEN_SEED + 1 and checks B1 and B2:
   - exact `==` on the median, the sign counts and the zero positions;
   - the 8-row diff against a pinned table;
   - the ORA sha256 values;
   - `_balance_signs` on the toy arrays of example 1, imported from the script by path;
   - the mutation: monkeypatch `_balance_signs` to identity and expect `SystemExit`.
2. **`tests/unit/test_make_wallenius_exact.py`** (CI, ≈ 90 lines), for B7 and B10:
   - examples 7–8, ω = 1 against hypergeometric fractions, symmetry, R1 vs R2 on small urns;
   - the precision rule (force a tiny initial P and check that it rises);
   - the longdouble refusal (monkeypatch `finfo`);
   - the error exits;
   - the `--categories CAT_010` byte equality across worker counts (≈ 20 s; mark `slow` if
     > 60 s on CI).
3. **`tests/golden/test_golden_wallenius_files.py`** (CI, ≈ 40 lines): B6 and the PROVENANCE
   fields (regexes like `test_provenance_records_the_run`), the input sha256 in the PROVENANCE
   equal to the committed ORA files, and the vacuity floor.
4. **`tests/unit/test_golden_provenance_dirs.py`** (CI, ≈ 20 lines): B9, example 10.
5. **On Slurm (lead / implementer records the job IDs):** the two generation jobs, B5, and the
   M8-12 golden tests on the branch
   (`sbatch scripts/slurm/test_branch.sbatch <worktree> tests/unit/test_golden_r_inputs.py tests/unit/test_golden_inputs_median_exact.py tests/unit/test_make_wallenius_exact.py tests/golden -q -n 4`).
   The session allocation runs only single test files.

Feature-loop order:
1. test-author (tests 1–4, red);
2. implementer (generator, sbatch changes, the provenance check, the Wallenius generator and
   sbatch), then commit and push;
3. `golden_r.sbatch` then commit the data;
4. `golden_wallenius.sbatch` then commit the data;
5. B5 on Slurm;
6. review.

Between steps 2 and 3 the branch is red on `test_committed_inputs_equal_a_fresh_generation`; it
is a branch, so that is expected.

## Estimate

- **Code ≈ 270 lines:**
  - `make_golden_inputs.py`, about 25;
  - the `golden_r.sbatch` guard, reason and lineage, about 20;
  - `check_golden_provenance.py` and the CI step name, about 8;
  - `make_wallenius_exact.py`, about 170;
  - `golden_wallenius.sbatch`, about 50.
- **Tests:** about 200 lines.
- **Data:** two gzip inputs and two gzip outputs regenerated, plus one new file of about 5 kB
  and two PROVENANCE files.
- **Slurm:** golden_r about 1 min (2 CPUs); wallenius about 16 min (4 CPUs, two runs); one
  test_branch job.
- **PRs:** one PR on `m8/golden-median-exact`, within the ~400-line guideline excluding data
  (question 4).

## Questions for the lead (each has a reversible default)

1. **How to balance the signs.** **Default: negate the |d| = 8 free values closest to zero**
   (§1.3), which changes 8 scores by ≤ 0.0175. *Alternative:* re-centre the 1,894 free values on
   their own median before rounding, with a sign-preserving `_make_unique`. That also works
   (checked: median 0.0, 997 / 6 / 997), but it changes every free score, and it changes the
   existing `_make_unique`. Either way the change is reversible by regenerating.
2. **Plant the sets on the final scores or on the pre-balance scores?** **Default: final**, so the
   sets are planted on exactly the scores R sees, and about 110 of 120 sets change. *Alternative:*
   plant on the pre-balance `values`. `gsea_sets.tsv.gz` then stays byte-identical, and the diff
   shrinks to 8 score rows, at the cost of planting weights that are 8 values off the written
   scores.
3. **ω of the exact reference.** **Default: the exact rational mean ratio from the decimal
   strings, correctly rounded once** (the definition; independent of numpy's summation order).
   *Alternative:* the same float as GIAP's `category_odds`, recomputed with numpy means, which
   isolates the tail computation completely. The two differ by < 1e-13 relative in p; the
   generator logs the ulp gap either way.
4. **One PR or two?** **Default: one PR on this branch**, with commits in the order median →
   R rerun → Wallenius. *Alternative:* split, with the median part first, because it unblocks
   M8-13 G1/G3a while the Wallenius part unblocks A1 (Wallenius). The two parts touch disjoint
   files apart from `check_golden_provenance.py`.
5. **The ORA guard override.** **Default:** `GOLDEN_R_ALLOW_ORA_CHANGE=1`, with the reason in
   `GOLDEN_R_REASON`. *Alternative:* no override; a changed ORA input then needs an edit of the
   sbatch.
6. **The wallenius directory in the provenance check.** **Default: yes** (§4.3). Without it,
   ADR-0044's "PROVENANCE entry" is not enforced.
7. **Run the Wallenius generator twice in the job** (about 16 min of 4 CPUs) to prove byte
   determinism? **Default: yes**, as golden_r.sbatch does. *Alternative:* one run, with
   determinism covered only by the unit test on CAT_010.
8. **If the regenerated fgsea outputs miss an M8-12 range test or the M8-13 A5 floor.**
   **Default: stop and ask.** The only levers are the planting constants in the generator, which
   would be a spec change, never `GOLDEN_SEED`. Judging by the scratch planted-z spread
   (2.3–13.8, 24 sets with |z| ≥ 2), a miss is unlikely.
9. **The `VERSIONS.txt` warnings caveat.** goseq's warnings come from the ORA call, so the text
   should be identical. **Default:** any `VERSIONS.txt` difference is reported in the PR, and only
   a warnings-only difference is accepted without asking.

## Lead decisions (2026-10-23)
1. **Q1–Q9: the defaults are accepted.**
   - **Q1:** `_balance_signs` flips the 8 values.
   - **Q2:** the sets are planted on the final scores.
   - **Q3:** ω is the exactly computed mean ratio, rounded once.
   - **Q4:** one PR, with commits in the order median → R rerun → Wallenius.
   - **Q5:** the ORA guard can be overridden with `GOLDEN_R_ALLOW_ORA_CHANGE=1`.
   - **Q6:** `check_golden_provenance.py` is extended to `tests/data/golden/wallenius/`.
   - **Q7:** the Wallenius generator runs twice in the job.
   - **Q8:** a miss stops the work and is raised with the lead; never reseed.
   - **Q9:** a `VERSIONS.txt` difference is accepted only in its warnings section.
2. **The tolerance key.** The lead commits `golden_r.wallenius_goseq_abs_tol: 1.0e-5` with `Gate-Change-Approved: ADR-0044` in the M8-13 PR, not here.
3. **Order of work:** test-author, then implementer (the generators and the sbatch), then the lead runs the two Slurm jobs and commits the generated data with PROVENANCE. `test_committed_inputs_equal_a_fresh_generation` stays red between the generator commit and the data commit; that is expected.
