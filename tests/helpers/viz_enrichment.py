"""Small synthetic ``OUT/`` trees for the V14 ``enrichment`` tests (M9-19).

Tables follow the ``enrichment_gsea`` / ``enrichment_detection`` Frictionless schemas (written
with ``helpers.viz_out_tree.write_schema_table``: schema header, absent cells ``NA``). The support
file ``plots/support/gsea_running.tsv[.gz]`` follows the M9-17 contract
(``helpers.support_analysis_data.COLUMNS["gsea_running"]``: ``collection, term, rank_position``
(1..n TUs), ``tu_id, score, in_set, running_es``), the six sets in M9-17's order (smallest q, then
largest |nes|, then collection, term). It is written from the contract, not from M9-17 code.

``running`` is fgsea's running sum (scoreType std, gseaParam 1; tests/unit/stats/
test_fgsea_running_es.py): a hit at rank i adds |score_i| / N_R, a miss subtracts 1 / (n - m).
A set's ``es`` in the table is the signed max deviation of that sum, so support and table agree.
Nothing here touches ``$GIPPA_DB_DIR`` or the network.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from helpers.support_analysis_data import COLUMNS as SUPPORT_COLUMNS
from helpers.viz_out_tree import make_out, write_events, write_resolved_config, write_schema_table

RUNNING_COLUMNS: tuple[str, ...] = SUPPORT_COLUMNS["gsea_running"]
COLLECTIONS: tuple[str, ...] = ("cog", "kegg", "go_bp", "go_mf", "go_cc", "go_slim")
N_TUS = 20
#: centred TU scores in rank order (descending), as fgsea ranked them
SCORES: tuple[float, ...] = tuple(round(2.0 - 0.2 * i, 6) for i in range(N_TUS))
ALPHA = 0.05


def running(scores: Sequence[float], members: Sequence[int]) -> list[float]:
    """fgsea's running sum after each ranked position; ``members`` are 1-based ranks."""
    n, hits = len(scores), {int(m) for m in members}
    w = np.abs(np.asarray(scores, dtype=float))
    nr = float(sum(w[i - 1] for i in hits))
    nm = n - len(hits)
    cur, out = 0.0, []
    for i in range(1, n + 1):
        if i in hits:
            cur += (w[i - 1] / nr) if nr > 0 else 1.0 / len(hits)
        else:
            cur -= 1.0 / nm
        out.append(cur)
    return out


def max_dev(r: Sequence[float]) -> float:
    lo, hi = float(min(0.0, min(r))), float(max(0.0, max(r)))
    return hi if hi > -lo else lo


def peak_rank(r: Sequence[float]) -> int:
    """1-based rank of the max deviation."""
    arr = np.abs(np.asarray(r, dtype=float))
    return int(np.argmax(arr)) + 1


def g6(x: float) -> float:
    """``x`` as a ``%.6g`` cell reads back."""
    return float(f"{x:.6g}")


def gsea_row(
    collection: str,
    term: str,
    *,
    name: str | None = None,
    es: float | None = 0.5,
    nes: float | None = 1.2,
    p: float | None = 0.01,
    q: float | None = 0.02,
    redundant: bool = False,
    flags: str | None = None,
    coverage: float = 0.8,
    size: int = 4,
) -> dict[str, object]:
    return {
        "collection": collection,
        "term": term,
        "name": name if name is not None else f"{term} name",
        "size_tus": size,
        "es": es,
        "nes": nes,
        "log2err": 0.5,
        "leading_edge": "TU0001;TU0002",
        "redundant": redundant,
        "parent_group": None,
        "coverage": coverage,
        "stat": nes,
        "effect": nes,
        "ci_low": None,
        "ci_high": None,
        "p_value": p,
        "q_value": q,
        "method": "fgsea_multilevel",
        "null_model": "tu_set_permutation",
        "n_perm": 1000,
        "flags": flags,
    }


def ora_row(
    collection: str,
    term: str,
    direction: str,
    stat: float,
    effect: float,
    q: float | None,
    *,
    name: str | None = None,
    coverage: float = 0.75,
    flags: str | None = "ci_unestimated",
) -> dict[str, object]:
    return {
        "collection": collection,
        "term": term,
        "name": name if name is not None else f"{term} name",
        "direction": direction,
        "n_universe": 200,
        "n_in_set": 20,
        "n_detected": 100,
        "k_detected_in_set": int(stat),
        "n_in_set_genes": 25,
        "k_detected_in_set_genes": int(stat),
        "odds_weight": 1.0,
        "p_classic": q,
        "p_noelim": q,
        "parent_group": None,
        "coverage": coverage,
        "stat": stat,
        "effect": effect,
        "ci_low": None,
        "ci_high": None,
        "p_value": None if q is None else min(1.0, q / 2),
        "q_value": q,
        "method": "cpb_tu_ora",
        "null_model": "conditional_poisson_binomial_tu",
        "n_perm": 0,
        "flags": flags,
    }


# ---------------------------------------------------------------- the main example

#: (collection, term, name, members (1-based ranks), q, nes factor, redundant, flags)
_SETS: tuple[tuple[str, str, str, tuple[int, ...], float, float, bool, str | None], ...] = (
    ("kegg", "ko00010", "Glycolysis", (1, 2, 4), 0.001, 2.5, False, None),
    ("go_bp", "GO:0006412", "translation", (1, 3, 5, 6), 0.002, 2.4, False, "genomic_cluster"),
    ("cog", "J", "Translation", (17, 19, 20), 0.01, 2.3, False, None),
    ("go_mf", "GO:0003735", "ribosome constituent", (2, 3, 7), 0.02, 2.2, True, None),
    ("go_cc", "GO:0005840", "ribosome", (1, 9, 12), 0.2, 2.1, False, "detection_confounded"),
    ("go_slim", "GO:0008150", "biological process", (11, 15, 18), 0.3, 2.0, False, None),
)
#: not in the top 6 (larger q) or not drawn (NA nes / q)
_OTHER: tuple[dict[str, object], ...] = (
    gsea_row("kegg", "ko00020", name="TCA cycle", es=0.3, nes=0.9, p=0.4, q=0.5),
    gsea_row("cog", "K", name="Transcription", es=0.2, nes=0.5, p=0.7, q=0.8),
    gsea_row("go_bp", "GO:0000001", name="too small", es=None, nes=None, p=None, q=None,
             flags="size_filtered", size=1),
)  # fmt: skip


def main_sets() -> list[dict[str, object]]:
    """The six top sets: members, running sums and table rows (``es`` = max deviation)."""
    out = []
    for coll, term, name, members, q, factor, redundant, flags in _SETS:
        r = running(SCORES, members)
        es = max_dev(r)
        out.append(
            {
                "collection": coll,
                "term": term,
                "name": name,
                "members": members,
                "running": r,
                "es": es,
                "row": gsea_row(coll, term, name=name, es=es, nes=round(es * factor, 4),
                                p=q / 2, q=q, redundant=redundant, flags=flags),
            }
        )  # fmt: skip
    return out


def main_gsea() -> list[dict[str, object]]:
    return [s["row"] for s in main_sets()] + [dict(r) for r in _OTHER]  # type: ignore[misc]


def main_ora() -> list[dict[str, object]]:
    """Example 1 (stat 12, effect 4 -> O/E 1.5) and example 3 (2 rows at q <= 0.05)."""
    return [
        ora_row("cog", "J", "over", 12, 4, 0.001, name="Translation"),
        ora_row("kegg", "ko00010", "under", 2, -6, 0.03, name="Glycolysis"),
        ora_row("go_bp", "GO:0006412", "over", 5, 1, 0.4, name="translation"),
    ]


def running_rows(sets: Sequence[Mapping[str, object]]) -> list[list[object]]:
    rows: list[list[object]] = []
    for s in sets:
        members = set(s["members"])  # type: ignore[call-overload]
        for i, (score, value) in enumerate(zip(SCORES, s["running"], strict=True), start=1):  # type: ignore[call-overload]
            rows.append([s["collection"], s["term"], i, f"TU{i:04d}", score, i in members, value])
    return rows


def write_running(out: Path, rows: Sequence[Sequence[object]], *, gz: bool = True) -> Path:
    from gippa.io.tsv import write_tsv

    sup = out / "plots" / "support"
    sup.mkdir(parents=True, exist_ok=True)
    path = sup / ("gsea_running.tsv.gz" if gz else "gsea_running.tsv")
    write_tsv(path, list(RUNNING_COLUMNS), rows)
    return path


def build_out(
    tmp: Path,
    gsea: Sequence[Mapping[str, object]] | None,
    ora: Sequence[Mapping[str, object]] | None,
    running_table: Sequence[Sequence[object]] | None = None,
    *,
    flags: Sequence[str] = (),
    modules: Sequence[str] | None = None,
    gz: bool = True,
) -> Path:
    """``OUT/`` with the two tables (``None``: absent; ``[]``: header-only), the support file
    (``None``: absent), run flags as ``function.run_flags`` events and a resolved config."""
    out = make_out(tmp)
    if gsea is not None:
        write_schema_table(out / "tables" / "enrichment_gsea.tsv", "enrichment_gsea", gsea)
    if ora is not None:
        write_schema_table(out / "tables" / "enrichment_detection.tsv", "enrichment_detection", ora)
    if running_table is not None:
        write_running(out, running_table, gz=gz)
    if flags:
        write_events(
            out / "logs" / "events.jsonl",
            [{"event": "function.run_flags", "data": {"flags": list(flags)}}],
        )
    if modules is not None:
        write_resolved_config(out / "resolved_config.yaml", modules)
    return out


def main_out(tmp: Path, **kw: object) -> Path:
    return build_out(tmp, main_gsea(), main_ora(), running_rows(main_sets()), **kw)  # type: ignore[arg-type]
