# M5-03: Stats engine: window scan, maximum-statistic null, FWER p-values, window merging

Milestone M5 · track S (engine) · depends on: none (reuses merged `perm_p`, M4-03) · estimate ≈ 240
code lines (+ ≈ 300 test lines) · **statistical** (label `stats`) · **can start now** · no tools, no
Slurm

## Spec sections implemented

- §8.3 "`scan_max(z, wmin, wmax, circular)` | numba; returns the window table and the maximum".
- §7.3.5(d) steps 4–6 quoted: "Scan all windows of 2–20 consecutive TUs (wrapping on circular replicons):
  S_w = Σ z_t / √n_w, for both high (S) and low (−S) islands." "For each replicate b, M_b = max over
  windows of |S_w|." "FWER. p_w = (1 + #{M_b ≥ |S_w|}) / (1 + B). Islands are windows with p_w ≤ 0.05;
  overlapping significant windows with the same direction merge into one island, represented by its
  best window."
- §8.1 "Islands | window | detrended TU scores, circular block bootstrap, maximum statistic | FWER |
  0.05".
- §7.3.6 unit "scan window enumeration with wrap-around".

## Public interface

`src/giap/stats/scan.py` (pure; numba allowed for the kernels, §6.6 deterministic order)
```python
@dataclass(frozen=True, slots=True)
class ScanResult:
    start: I64          # 0-based first position of each window
    length: I64         # n_w
    s: F64              # S_w = Σ z / √n_w (signed)
    max_abs: float      # max |S_w|; NaN when there is no window

def n_windows(n: int, wmin: int, wmax: int, circular: bool) -> int: ...
def scan_max(z: ArrayLike, wmin: int, wmax: int, circular: bool) -> ScanResult: ...
def scan_null_max(Z: F64, wmin: int, wmax: int, circular: bool) -> F64: ...   # (B, n) → (B,)
def fwer_p(s_abs: ArrayLike, null_max: ArrayLike) -> F64: ...
@dataclass(frozen=True, slots=True)
class Cluster:
    direction: Literal["high", "low"]
    best: int                    # index into the ScanResult arrays
    members: tuple[int, ...]     # window indices in the cluster (sorted)
def merge_windows(res: ScanResult, significant: ArrayLike, n: int, circular: bool) -> list[Cluster]: ...
```
- **Windows.** Lengths w = wmin … min(wmax, w_cap) with w_cap = n − 1 (circular: no window covers the
  whole circle) or n (linear). Circular: every start 0 … n − 1 for each w (wrapping mod n); linear:
  starts 0 … n − w. Window order: by length, then start.
- **S.** Σ over the window's positions of z, divided by √w; computed from prefix sums and equal to direct
  summation within 1e-9·max(1, Σ|z|). Positive S = high, negative = low; |S| is the scan statistic.
- **Does not exist:** n < wmin (no window) → empty arrays and `max_abs` NaN; NaN or ±inf in z →
  `ValueError` (the caller drops TUs without a score before scanning).
- **scan_null_max.** Row b of the bootstrap matrix → M_b = `scan_max(Z[b], …).max_abs` (vectorised; same
  windows for every row).
- **fwer_p.** p_w = `perm_p(|S_w|, M, "greater")` element-wise, i.e. (1 + #{M_b ≥ |S_w|})/(1 + B) with
  the engine's tie tolerance; NaN M_b dropped; no valid M_b → NaN.
- **merge_windows.** Among windows with `significant` True, two windows of the same direction are
  linked when they share ≥ 1 position (on circular input, cyclically); clusters are the connected
  components (transitive). Best window = largest |S|, then fewer positions, then smaller start (linear
  start; circular: start in [0, n)). Clusters are returned ordered by the best window's start. Windows of
  opposite directions never merge.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | z = [0, 0, 3, 3, 0, 0], linear, wmin 2, wmax 3 | 9 windows (5 + 4); `max_abs` 4.242641 at start 2, length 2 |
| 2 | same z, circular | 12 windows (6 + 6) |
| 3 | z = [3, 0, 0, 0, 0, 3], circular, w 2–3 | best start 5, length 2 (wraps to position 0), S 4.242641 |
| 4 | same z, linear | `max_abs` 2.121320 (windows [0, 1] and [4, 5] tie; best = start 0) |
| 5 | z = [0, −2, −2, 0], linear, w 2 | min S −2.828427 at start 1; `max_abs` 2.828427; direction low |
| 6 | n = 1, wmin 2 | no window; `max_abs` NaN; `n_windows` 0 |
| 7 | circular n = 3, wmin 2, wmax 20 | only w = 2 (cap n − 1): 3 windows |
| 8 | `fwer_p([6.5, 9, 10], M = [1 … 9])` (B = 9) | [0.4, 0.2, 0.1] (the tie at 9 counts) |
| 9 | significant windows (start, length, S) linear n = 10: (2, 2, +4.2), (3, 2, +3.0), (5, 2, +2.9), (3, 2, −3.1) | clusters: high {(2,2), (3,2)} best (2,2); high {(5,2)}; low {(3,2,−)} |
| 10 | significant high windows (0,3), (2,3), (4,3) | one cluster (chain), best = largest S |
| 11 | circular n = 10, significant high (9, 2) and (0, 2) | one cluster (share position 0) |

## Edge cases

- All z = 0 → every S = 0, `max_abs` 0.0 (exists); with any null, p_w = 1.
- wmin > wmax (bad config) → `ValueError`; wmin < 1 → `ValueError`.
- Bootstrap rows of a replicon with fewer TUs than wmin: the caller does not scan it (M5-06).
- numba first-call compile time is excluded from the benchmark; results are identical with and without
  numba (`NUMBA_DISABLE_JIT=1` test).

## Acceptance criteria

- [ ] `tests/unit/stats/test_scan.py`: examples 1–11; `n_windows` equals the window count for n ∈ [1, 50],
  w ranges, both topologies.
- [ ] Property: S on random z equals the direct sum; circular `max_abs` invariant under rotation of z;
  linear `max_abs` invariant under reversal; `scan_null_max` equals row-wise `scan_max`.
- [ ] Performance: `scan_null_max` with B = 1,000, n = 2,500, w 2–20 < 2 s on 1 CPU after compilation
  (benchmark recorded).

## Methods note (statistician sign-off)

`docs/methods/islands.md` §"Scan and maximum statistic": the max over both directions controls FWER over
all windows of one family (both signs), why S/√n is not variance-standardised for TU size and what that
does to power (TU scores of large TUs have smaller variance), merging and the best-window representative.

## Questions for the lead

1. Window cap on small circular replicons n − 1 (never the whole circle)? **Proposed default:** yes; a
   window covering the whole replicon has S = Σz/√n ≈ 0 after detrending and is not an island.
2. "Overlapping" = sharing ≥ 1 TU; adjacent (touching) windows stay separate islands? **Proposed default:**
   yes, only shared TUs merge.
