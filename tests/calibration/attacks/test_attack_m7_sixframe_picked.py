"""Attack nulls for the M7-12 six-frame decoys under a picked class FDR (GIPPA_SPEC.md §7.6.4,
§8.1 row "Novel ORF", §8.2 "decoys mimic targets", §8.5; M7-11 rule 2; statistician, PR #225
round 3).

M7-11 scores a target ORF and its makedb decoy entry by their best PSM, keeps the winner and
estimates the class FDR as (D + 1) / T (``gippa.stats.tda.picked_fdr``). That is the knockoff+
estimator, valid when the label of each null entry's winner is a fair coin independent of its
score. These attacks check that the simulator's decoys (lead decisions 21, 24) give exactly that,
with an oracle that does not use GIPPA's novel-ORF code:

- **Fair label** (fast): per six-frame wrong peptide, the side (target or ``rev_GIPPAORF_`` decoy)
  holding its best PSM is a fair coin, pooled exact binomial p >= ``ks_min_p``. A mutant whose
  decoys do not mimic the targets (decoy scores shifted by -0.5 SD) must fail the same test, so
  the check has power.
- **Oracle picked class FDP** (full): entry-level picking within each class (wrong peptides keyed
  by the makedb entry their decoys name, planted ORFs by ``orf_id``), knockoff+ q; mean FDP
  <= q + ``se_multiplier`` SE at q in {0.01, 0.05, 0.10}, and calls exist at 0.01. The scenario
  is planted densely (200 per class) so that q = 0.01 is reachable at all; ``nterm_extension`` is
  left out because its planted peptides carry the gene's locus tag (M7-11 pairs it by locus tag).
  The same mutant must break the bound at q = 0.10 in at least one class.
- **Class q floor** (fast): with (D + 1) / T, a class whose planted
  ORFs with a true peptide number fewer than 1 / alpha (``stats.alpha.novel_orfs`` = 0.01) cannot
  call any ORF, so G3's F1 is 0 by construction. The default ``novel_orf_power`` at 1 Mb plants 40
  per class (about 30 with >= 2 peptides; oracle min class q 0.025-0.028 over 12 seeds, 0 calls).
  Statistician round-3 block, item 1. Lead decision 26: the test runs at
  ``n_per_class=NOVEL_PER_CLASS_GATE`` (300) at 1 Mb, the count M7-13 G3/G4 use. Intergenic is
  limited by gap supply (88-122 ORFs with a peptide over seeds 1-8), so the floor is asserted on
  the other three classes only (lead decision 28).
"""

from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np
import pytest
from helpers import m7_round3_oracle as Q
from helpers import m7_sim_oracle as O
from scipy import stats

import gippasim

pytestmark = pytest.mark.calsuite("novel_orfs")

ALPHA_NOVEL = 0.01  # stats.alpha.novel_orfs (src/gippa/config/defaults.yaml)
QS = (0.01, 0.05, 0.10)
CLASSES3 = ("intergenic", "cds_outframe", "antisense")
DENSE = 200
SHIFT = -0.5  # mutant: decoy scores moved down by half a null SD


def _simulate(out: Path, seed: int, size: str, per_class: int | None = None,
              monkeypatch: pytest.MonkeyPatch | None = None) -> Path:  # fmt: skip
    import gippasim.scenarios as sc

    if per_class is not None:
        assert monkeypatch is not None
        monkeypatch.setattr(
            sc, "NOVEL_PER_CLASS", MappingProxyType({**sc.NOVEL_PER_CLASS, size: per_class})
        )
    gippasim.simulate(sc.novel_orf_power(seed, size), seed, out)
    return out


def _side_wins(out: Path, shift: float = 0.0) -> tuple[int, int]:
    """(peptides whose best PSM is a target row, wrong peptides with any row)."""
    loci = Q.wrong_loci(out)
    t: dict[str, float] = defaultdict(lambda: -math.inf)
    d: dict[str, float] = defaultdict(lambda: -math.inf)
    for r in O.rows(out / "peptides.tsv"):
        if r["is_decoy"] == "true" and r["protein_ids"].startswith(Q.SIXFRAME_DECOY_PREFIX):
            p = r["peptide"][::-1]
            d[p] = max(d[p], float(r["score"]) + shift)
        elif r["is_decoy"] == "false" and r["peptide"] in loci:
            t[r["peptide"]] = max(t[r["peptide"]], float(r["score"]))
    seen = [p for p in loci if p in t or p in d]
    return sum(t[p] > d[p] for p in seen), len(seen)


def _knockoff_q(t: np.ndarray, d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Picked (D + 1) / T monotone q of each entry (NaN unless the target wins) and target wins."""
    has_t, has_d = ~np.isnan(t), ~np.isnan(d)
    dw = has_d & (~has_t | (d >= np.where(has_t, t, -np.inf)))
    tw = has_t & ~dw
    ts, ds = np.sort(t[tw])[::-1], np.sort(d[dw])
    raw = (len(ds) - np.searchsorted(ds, ts, side="left") + 1) / np.arange(1, len(ts) + 1)
    qs = np.minimum.accumulate(raw[::-1])[::-1]
    q = np.full(len(t), np.nan)
    order = np.argsort(-t[tw], kind="stable")
    q[np.flatnonzero(tw)[order]] = qs
    return q, tw


def _class_fdp(out: Path, shift: float = 0.0) -> dict[str, dict[float, tuple[int, int]]]:
    """Per class and q: (calls, false calls) of the oracle entry-level picked class FDR."""
    entries = Q.makedb_entries(out, out / "db")
    index = Q.EntryIndex(entries, Q.lengths(out))
    loci = Q.wrong_loci(out)
    nov = {r["orf_id"]: r["class"] for r in O.rows(out / "truth" / "novel_orfs.tsv")}
    bt: dict[str, float] = defaultdict(lambda: -math.inf)
    bd: dict[str, float] = defaultdict(lambda: -math.inf)
    cls: dict[str, str] = {}
    true: dict[str, bool] = {}
    for r in O.rows(out / "peptides.tsv"):
        s = float(r["score"]) if r["score"] not in ("", "NA") else math.nan
        if r["is_decoy"] == "true":
            m = Q.DECOY_ID.match(r["protein_ids"])
            if m is None:
                continue
            k = r["protein_ids"][len("rev_") :]
            cls[k], bd[k] = m.group(1), max(bd[k], s + shift)
            true.setdefault(k, False)
        elif r["peptide"] in loci:
            sp = loci[r["peptide"]]
            (e,) = index.entry_of(sp.contig, sp.s1, sp.e1, sp.strand)
            cls[e.db_id], bt[e.db_id] = e.cls, max(bt[e.db_id], s)
            true.setdefault(e.db_id, False)
        elif r["protein_ids"] in nov:
            k = r["protein_ids"]
            cls[k], bt[k], true[k] = nov[k], max(bt[k], s), True
    res: dict[str, dict[float, tuple[int, int]]] = {}
    for c in CLASSES3:
        ks = [k for k in cls if cls[k] == c]
        t = np.array([bt.get(k, math.nan) for k in ks])
        d = np.array([bd.get(k, math.nan) for k in ks])
        tr = np.array([true[k] for k in ks], dtype=bool)
        q, tw = _knockoff_q(t, d)
        res[c] = {}
        for a in QS:
            call = tw & (q <= a)
            res[c][a] = (int(call.sum()), int((call & ~tr).sum()))
    return res


def _mean_se(fdp: list[float]) -> tuple[float, float]:
    x = np.asarray(fdp)
    return float(x.mean()), float(x.std(ddof=1) / math.sqrt(len(x)))


@pytest.fixture(scope="module")
def base(tolerances: dict[str, Any]) -> int:
    return int(tolerances["seeds"]["calibration_base"])


@pytest.mark.calibration
def test_attack_picked_label_is_fair_coin(
    tmp_path_factory: pytest.TempPathFactory, base: int, tolerances: dict[str, Any]
) -> None:
    outs = [_simulate(tmp_path_factory.mktemp(f"m7pk{r}"), base + 7000 + r, "300k")
            for r in range(3)]  # fmt: skip
    p_min = float(tolerances["ks_min_p"])
    wins, n = map(sum, zip(*(_side_wins(o) for o in outs), strict=True))
    assert n >= 300, n
    assert stats.binomtest(wins, n, 0.5).pvalue >= p_min, (wins, n)
    mw, mn = map(sum, zip(*(_side_wins(o, SHIFT) for o in outs), strict=True))
    assert stats.binomtest(mw, mn, 0.5).pvalue < p_min, ("mutant not detected", mw, mn)


@pytest.mark.calibration
def test_attack_default_scenario_can_reach_class_alpha(
    tmp_path_factory: pytest.TempPathFactory, base: int
) -> None:
    import gippasim.scenarios as sc

    # Lead decision 26: G3/G4 plant NOVEL_PER_CLASS_GATE per class at 1 Mb; decision 28: the
    # gap-limited intergenic class is not asserted
    out = tmp_path_factory.mktemp("m7floor")
    cfg = sc.novel_orf_power(base + 7100, "1M", n_per_class=sc.NOVEL_PER_CLASS_GATE)
    gippasim.simulate(cfg, base + 7100, out)
    n = defaultdict(int)
    for r in O.rows(out / "truth" / "novel_orfs.tsv"):
        n[r["class"]] += int(r["n_true_peptides"]) >= 1
    floor = {c: 1 / max(n[c], 1) for c in Q.CLASSES4 if c != "intergenic"}
    assert all(f <= ALPHA_NOVEL for f in floor.values()), floor


@pytest.mark.calibration_full
def test_attack_oracle_picked_class_fdp_dense(
    tmp_path_factory: pytest.TempPathFactory, base: int, tolerances: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:  # fmt: skip
    k = float(tolerances["se_multiplier"])
    outs = [_simulate(tmp_path_factory.mktemp(f"m7dense{r}"), base + 7200 + r, "1M", DENSE,
                      monkeypatch) for r in range(15)]  # fmt: skip
    ok, mut = [_class_fdp(o) for o in outs], [_class_fdp(o, SHIFT) for o in outs]
    broken = False
    for c in CLASSES3:
        for a in QS:
            m, se = _mean_se([r[c][a][1] / max(r[c][a][0], 1) for r in ok])
            assert m <= a + k * se, (c, a, m, se)
            mm, mse = _mean_se([r[c][a][1] / max(r[c][a][0], 1) for r in mut])
            broken |= a == 0.10 and mm > a + k * mse
        assert sum(r[c][0.01][0] for r in ok) > 0, ("no calls at 0.01", c)
    assert broken, "mutant decoys (scores -0.5 SD) not detected at q = 0.10"
