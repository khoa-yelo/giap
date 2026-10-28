# Variants (Module 5)

Methods notes for `gippa.variants` (GIPPA_SPEC.md §7.5). Only the class error rates section exists
so far; the pileup, call and VCF sections follow with their modules.

## Class error rates

Implements `gippa.variants.error_rates` (GIPPA_SPEC.md §7.5.4, §8.1 row "ALT call is real",
§7.5.5 `CLASS_ERR`; docs/specs/M7-10-class-error-rates.md; statistician review of PR #239).

### Estimator

- **Pairs.** A de novo row pairs with a reference row (`--reference-peptides`, via the C17 table
  `peptide_spectra.parquet`) of the same `(run, spectrum_id)` when the reference row is a target
  with `q_value ≤ 0.01` (`REF_MAX_Q`, as in M3-11). Preflight rejects a spectrum with two
  different peptides (ADR-0005), so each spectrum has at most one reference identification.
- **Observations.** Each paired de novo row whose peptide supports an ALT site is one observation
  of each class `(alt_type, mass_equivalent)` among the ALTs it supports. A peptide supporting two
  ALTs of different classes counts once in each; two sites of one class count once.
  `mass_equivalent` comes from the site's M7-05 flags; frameshifts are their own class
  `(frameshift, false)`.
- **Error.** An observation is an error when the de novo and reference peptides differ after I→L
  folding.
- **Rate.** Per class, `class_error_rate = errors / n` and `class_error_n = n`, pooled over all
  sites of the class and written on every ALT row of the class (and VCF `CLASS_ERR`). A class
  with `n = 0`, a run without `--reference-peptides`, reference rows without `spectrum_id`
  (warning `no_spectrum_pairing`) or a database-only input give `unestimated`. A class with
  observations but no error has rate 0 with its `n`.

### Estimand

The rate estimates **the share of erroneous reads among ALT-supporting de novo reads whose
spectrum a database search identifies confidently**. It is not the false-ALT rate among all
ALT-supporting reads, and its relation to that rate depends on the reference search. Three biases
act, two upward and one downward.

1. **Reference database without the variant: upward, up to degenerate (spec rule 5).** When the
   reference identifications come from a search against the same reference proteome GIPPA maps
   to, a correct ALT read can never agree with its spectrum's reference identification: the ALT
   peptide is not in that database. Such a spectrum is either not identified confidently (no
   pair) or identified as some other peptide (an "error"). Every paired observation is then an
   error, and the **estimate is 1 (up to I/L-only ALTs) whatever the true error rate**. It is
   then trivially at least the true rate, which carries no information.
   `tests/calibration/attacks/test_attack_class_error_rates.py::test_reference_db_degenerate`
   pins this. The rate is informative only when the reference search could return the ALT
   peptide: a sample-specific database (the strain's own genome or a variant-aware search), or a
   panel of known variants. Read a rate of 1 from a reference-proteome search as "not measured",
   not as "every ALT is false".
2. **Spectrum-quality selection: downward.** De novo errors are more common on poor spectra, and
   poor spectra are less often identified at q ≤ 0.01. The paired spectra are therefore cleaner
   than the unpaired ones, and the rate is biased low for all ALT-supporting reads (attack
   `test_quality_selection`: all-reads truth 0.30, paired truth and estimate 0.20).
3. **Reference false matches: upward, small.** `q_value ≤ 0.01` bounds the expected false share
   of the confident reference identifications, not zero. A correct de novo read on a spectrum
   whose confident reference ID is wrong counts as an error. A wrong de novo read is missed only
   if the wrong reference ID happens to equal it, which is rare. With de novo error rate e on the
   paired spectra of a class and false-match share φ among them, the estimate is
   e + (1 − e)·φ, so the bias is (1 − e)·φ: largest for classes with a low error rate, at most φ,
   and zero for classes whose reads are all errors. φ is a per-class share: q bounds the FDR of
   the whole reference list, and a class can sit on spectra with more false matches than average.
   Attack `test_reference_false_matches` pins the formula (e 0.2, φ 0.05, estimate 0.24).

So **the rate is not an upper bound in general**. It is an upper bound only when the reference
search lacks the variant, where it is degenerate. With a variant-aware reference it is close to
the error rate on confidently identified spectra and, if quality selection matters, low for the
rest. No correction is applied in v1; the bias is reported, not removed.

### Uncertainty and small n

Each rate carries its `n`, so a binomial interval can be computed from the row. GIPPA reports no
interval column in v1 (the §7.5.5 schema has none). The recommended reading is the Wilson 95%
interval, (p + z²/2n ± z √(p(1 − p)/n + z²/4n²)) / (1 + z²/n) with z = 1.96. At n = 5 and
p = 1 it is [0.57, 1]; at n = 20 and p = 0.2 it is [0.08, 0.42]. The interval assumes that
observations are independent spectra. One spectrum reported at two charges is one observation,
not two. Observations of one class at one site are correlated, because the same error mode
repeats. The interval is therefore optimistic when a few sites dominate a class.

No minimum n is imposed: a rate with tiny n is shown, always with its n, and `unestimated`
means exactly n = 0. A `low_n` flag (n < 20; there the Wilson half-width at p = 0.5 is about 0.2
or more) is a backlog candidate, not a v1 rule.

### Pooling

Rates are pooled per class over sites (lead decision 2). A site's own reliability is not
measured. A site-level rate needs tens of paired spectra per site, and a pilot run gives a few
paired ALT-supporting spectra per site at most. The class rate is also a mixture across sites:
a class dominated by one bad site inherits that site's rate.

**Stop-codon sites (M7-16 rule 11).** A read-through site has REF `*` at an annotated terminal
stop. It and the U-over-TGA / O-over-TAG zero-change ALTs are in class
`(substitution, mass_equivalent = false)`, and their `class_error_rate` is that class's pooled
rate. The rate is dominated by ordinary SNV sites and does not measure read-through reads: with a
reference search on Bakta proteins, which end at the stop, every read-through read is in the
degenerate case of bias 1 (the reference cannot contain the peptide). A second-pass search on a
`makedb --patch` database does not lift this, since its identification carries `X` at the stop,
which differs from every de novo residue. The rate feeds only `class_error_rate`/`CLASS_ERR`, never
`confidence` (`high` = `n_alt` ≥ k), so TE2 does not depend on it.

### Calibration

- `tests/calibration/test_variant_class_error.py`, no-variant arm: every ALT-supporting read is
  an error, so truth and estimate both sit at 0.97–1.0. This arm is a sanity check only; a
  constant estimator `rate = 1` passes it.
- The same file's planted arm uses gippasim's reference peptides, which are the sample's own
  identifications and so a variant-aware reference. Truth there is 0.23–0.73 for the
  substitution and frameshift classes. It is gated two-sided at |bias| ≤ 0.05
  (`denovo_precision_max_abs_bias`) on classes estimated in at least half the replicates; at
  least two gated classes must have truth below 0.9, so a constant estimator (1 or 0) fails.
  This is the arm with teeth; it runs nightly (fast, R = 20) and in the full suite (R = 100).
- **Measured planted-arm bias (R = 100, full run, 2026-10-23):** frameshift +0.0087 (SE 0.0031),
  substitution +0.0037 (0.0007), substitution with `mass_equivalent` +0.0021 (0.0007); block
  and swap 0 (truth 1). Every informative class is biased up, by 3 to 5 SE for the
  substitution classes and 2.8 SE for frameshift. The frameshift class is not singled out: its
  larger SE comes from n ≈ 25 paired observations per replicate against ≈ 470 and ≈ 220.
  The bias is the reference false-match term (bias 3 above), not quality selection:
  - Against the truth on paired spectra only, the bias is the same: +0.0088, +0.0042, +0.0029.
    So selection of confident spectra contributes nothing in this scenario. gippasim draws
    reference confidence independently of de novo errors.
  - In a decomposition over the first 20 planted replicates (statistician, from gippasim truth),
    the estimate minus the paired truth equals the share of paired observations where the de novo
    read is correct and the confident reference ID is wrong: 0.016 frameshift, 0.0047 and 0.0041
    substitution. No wrong de novo read was masked by an agreeing wrong reference ID.
  - The share of wrong confident reference IDs on these spectra was φ ≈ 0.011 for the
    substitution classes, near the q ≤ 0.01 threshold. It was 0.020 for frameshift spectra,
    a small class where a few false matches weigh more. (1 − e)·φ then predicts 0.015, 0.0040
    and 0.0033, matching the measured shares.

  The bias is therefore expected, has a known sign and size (≤ φ, about the reference FDR), and
  sits well inside the 0.05 bound. No correction is applied: subtracting (1 − e)·φ needs φ per
  class, which a real run does not know.
- `tests/calibration/attacks/test_attack_class_error_rates.py`:
  - `test_known_rate_mix` plants per-class rates of 0, 0.05, 0.3 and 0.8 and checks exact
    agreement with a brute-force oracle, |bias| ≤ 0.05 and Wilson coverage. Its traps are I/L,
    spectrum IDs reused across runs, q = 0.02 rows and confident decoys.
  - `test_reference_db_degenerate` and `test_quality_selection` pin the two biases above.
  - `test_duplicate_charge_rows_count_once` checks the one-spectrum-one-observation rule.
- M7-13, on M7-12's `variant_f1` scenario, adds the arm with a reference search against the
  reference genome's proteome (the degenerate case, recorded), and a known-rate mix at two de
  novo error levels.

### Sign-off

Statistician, PR #239: blocked 2026-10-22 (B1 double-counted charge states, B2 gate without
teeth, B3 rule 5 claimed an upper bound). Re-reviewed and **approved 2026-10-23** at 0585894f:

- B1 is fixed: one observation per `(run, spectrum_id, peptide)`, and preflight rejects a
  spectrum with two stripped peptides.
- B2 is fixed: the planted arm is gated two-sided with a teeth guard.
- B3 is fixed: rule 5 and acceptance 2 are reworded.
- The residual upward bias is explained above (reference false matches).

Open items (backlog, not blocking): a run flag when every class rate is 1, and an optional
`low_n` flag or interval columns.
