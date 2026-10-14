# M3-05: Residue masses, the generated mass-equivalence list, and the difference-operation table

Milestone M3 · depends on: none (M0 only) · creates `src/giap/data/residue_masses.tsv`,
`src/giap/data/modifications.tsv`, `src/giap/mapping/masses.py`, `src/giap/mapping/diffops.py` ·
pure Python · estimate ≈ 200 code lines (data files excluded)

## Spec sections implemented

- §4 `data/`: "bundled small data: residue masses, mass equivalents".
- Appendix B.5 residue and modification masses (monoisotopic, Da) — copied verbatim into the two
  TSVs — and "**Checks the generated list MUST pass** at the default tolerance of 0.05 Da: N = GG
  (Δ 0); Q = AG (Δ 0); GT = AS (Δ 0); AV = GL (Δ 0); AN = GQ (Δ 0); W vs AD/EG (Δ 0.0153); R vs GV
  (Δ 0.0112); K vs Q (Δ 0.0364); A + Oxidation = S; N + Deamidated = D; G + Methyl = A; M + Oxidation
  vs F (Δ 0.0329)."
- §7.2.4 the cost-1 operations: substitution; adjacent swap (XY ↔ YX); "a mass-equivalent block
  exchange between the residue compositions {N} ↔ {G,G}, {Q} ↔ {A,G}, {G,T} ↔ {A,S}, {A,V} ↔ {G,L}
  or {A,N} ↔ {G,Q} (both orders of each pair)"; "a near-equal exchange: {W} ↔ {A,D} or {E,G}, {R} ↔
  {G,V}"; "a ±1-nt frameshift split"; "**I/L** is free."
- Consumers: the DP verifier (M3-06), frameshift (M3-08), the simulator's de novo error model
  (M3-20), and Module 5's `mass_equivalent` flag (M7).

## Public interface

```python
# masses.py
RESIDUE_MASS: Mapping[str, float]          # from residue_masses.tsv (I and L both 113.08406)
MOD_DELTA: Mapping[str, float]             # from modifications.tsv
@dataclass(frozen=True, slots=True)
class MassEquivalence:
    left: str; right: str                  # residue strings, or "A+Oxidation" for modified forms
    delta_da: float                        # |m(left) − m(right)|
    kind: Literal["block", "near_block", "substitution", "modification"]
def mass_equivalents(tol_da: float = 0.05) -> list[MassEquivalence]: ...  # generated, sorted, deterministic

# diffops.py
DiffType = Literal["substitution", "swap", "block", "near_block", "frameshift"]
DIFF_TYPES: tuple[DiffType, ...]           # tie-break precedence, in this order
BLOCKS: tuple[tuple[str, str], ...]        # composition pairs of §7.2.4, sorted letters: ("N","GG"), ("Q","AG"), ("GT","AS"), ("AV","GL"), ("AN","GQ")
NEAR_BLOCKS: tuple[tuple[str, str], ...]   # ("W","AD"), ("W","EG"), ("R","GV")
def block_tables() -> BlockTables: ...      # numba-friendly uint8 arrays over M3-02 codes, I→L applied
```

`mass_equivalents` enumerates, over the 20 standard residues (I and L merged): 1↔1 residue pairs,
1↔2 and 2↔2 compositions, and residue + modification ↔ residue, keeping |Δ| ≤ `tol_da`.
`BLOCKS`/`NEAR_BLOCKS` are **fixed by the spec** (not derived from the generated list); a test
checks they are consistent with it (Δ < 0.001 for `BLOCKS`, Δ ≤ 0.05 for `NEAR_BLOCKS`).
`block_tables` expands every composition into all residue orders ("both orders of each pair"),
e.g. {A,G} → `AG`, `GA`; {G,L} with L standing for I/L.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `RESIDUE_MASS["N"] − 2·RESIDUE_MASS["G"]` | `0.00001` (B.5 "Δ 0") |
| 2 | Δ(W, AD) and Δ(W, EG) | `0.01526`, `0.01526` (reported `0.0153`) |
| 3 | Δ(R, GV) | `0.01124` |
| 4 | Δ(K, Q) | `0.03638`; present in `mass_equivalents()` as `substitution` |
| 5 | `A+Oxidation` vs `S` | Δ `0.000005` → `modification` entry |
| 6 | `M+Oxidation` vs `F` | Δ `0.0329` → entry at 0.05, absent at `tol_da=0.03` |
| 7 | `N+Deamidated` vs `D`; `G+Methyl` vs `A` | both present |
| 8 | `block_tables()` | contains `N`↔`GG`, `Q`↔`AG`, `Q`↔`GA`, `GT`↔`AS`, `TG`↔`SA`, `GT`↔`SA`, `AV`↔`GL`, `VA`↔`LG`, `AN`↔`QG`; `W`↔`DA`; not `K`↔`Q` |

## Edge cases

- I and L are one residue for every enumeration (no `I`↔`L` entry; it is free, not a difference).
- U and O are in the mass table (B.5) but excluded from enumeration (not in table-11 translations).
- The list is sorted by (`kind`, `delta_da`, `left`, `right`) and identical across runs.

## Acceptance criteria

- [ ] `tests/unit/test_masses.py`: all twelve B.5 checks, examples 1–8.
- [ ] Test: every pair in `BLOCKS` and `NEAR_BLOCKS` appears in `mass_equivalents(0.05)`.
- [ ] The two TSVs carry a header comment with the source (Appendix B.5; Unimod for Δ values).

## Questions for the lead

1. `diff_types` vocabulary for `peptide_hits.tsv` is not fixed by §9.3 (`variants.tsv` has
   `alt_type` enum{substitution, swap, block, frameshift}). **Proposed default:**
   {`substitution`, `swap`, `block`, `near_block`, `frameshift`}; Module 5 maps `near_block` →
   `block` with the `mass_equivalent` flag.
