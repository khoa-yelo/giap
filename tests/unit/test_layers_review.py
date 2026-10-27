"""M2-20 after the code review of PR #104: lead decisions 9 and 11
(docs/specs/M2-20-context-layers.md; GIPPA_SPEC.md §6.7, §9.3, §9.4).

* Decision 9: ``layer_missing:<layer>`` is a structured event (``extra={"event": ...}``) that
  reaches ``OUT/logs/events.jsonl`` through the run's ``gippa.steps`` events handler, not only a
  log line. Accepted shapes: the event name is the flag, or the flag (or ``layer=<layer>`` with an
  event name ending in ``layer_missing``) is in the event's ``data``.
* Decision 11: tests that kill the review's surviving mutants — the ``crosses_origin`` flag, the
  circular wrap of the next same-strand gene (``terminator_conf_downstream``) and of
  ``assign_tus``, both ``tu_sources`` rules, the precomputed ``tx_rates`` exemption of
  ``layer_flags`` (pseudogene phase: ``test_layers_gff3.py``).

Inputs: the merged producers on gippasim SIM100K (``helpers.layers_chain``) and the simulator's
precomputed layers (M2-02 TSVs).
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers import layers_oracle as orc
from helpers.layers_chain import (
    Chain,
    build_chain,
    call_gene_layers,
    fresh_steps,
    read_layers_frames,
)

from gippa.config import GippaConfig

ALL_KEYS = {
    "genes", "features", "replicons", "pairs_raw", "tus", "promoter_windows", "tx_rates",
    "terminators", "terminators_crosscheck", "rbs", "codon", "regions", "region_genes", "genome",
}  # fmt: skip


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Chain:
    return build_chain(tmp_path_factory.mktemp("layers_review_chain"))


def _row(df: pd.DataFrame, tag: str) -> dict[str, Any]:
    rows = df[df["locus_tag"] == tag]
    assert len(rows) == 1, tag
    return rows.iloc[0].to_dict()


# ----------------------------------------------------------------------------- decision 9


def _matches(ev: dict[str, Any], flag: str) -> bool:
    data = ev.get("data") or {}
    layer = flag.split(":", 1)[1]
    return (
        ev.get("event") == flag
        or flag in json.dumps(data, sort_keys=True)
        or (str(ev.get("event", "")).endswith("layer_missing") and data.get("layer") == layer)
    )


def test_layer_missing_events_reach_events_jsonl(chain: Chain, tmp_path: Path) -> None:
    """``--modules core`` with Promoter Calculator off: layer_missing:genomad / antismash /
    tx_rates are WARNING events of step context_layers in logs/events.jsonl."""
    from gippa.cli import _record_step_events

    cfg = GippaConfig.model_validate(
        {"seed": 1, "modules": ["core"],
         "annotation": {"promoter_calculator": {"enabled": False}}}
    )  # fmt: skip
    inputs = {k: v for k, v in chain.inputs.items() if k != "tx_rates"}
    outs = {"annotation": tmp_path / "annotation", "layers": tmp_path / "layers.parquet"}
    logs = tmp_path / "logs"
    restore = _record_step_events(logs)
    try:
        with pytest.MonkeyPatch.context() as mp:
            fresh_steps(mp).run_step("context_layers", cfg, inputs, outs)
    finally:
        restore()
    path = logs / "events.jsonl"
    assert path.is_file(), "no structured event reached events.jsonl"
    events = [json.loads(x) for x in path.read_text(encoding="utf-8").splitlines()]
    for flag in ("layer_missing:genomad", "layer_missing:antismash", "layer_missing:tx_rates"):
        hits = [e for e in events if _matches(e, flag)]
        assert hits, (flag, events)
        assert all(e["step"] == "context_layers" and e["level"] == "WARNING" for e in hits), hits


# ----------------------------------------------------------------------------- decision 11


def test_crosses_origin_flag_follows_the_gene_table(chain: Chain) -> None:
    """A gene with ``crosses_origin`` in genes.parquet carries the ``crosses_origin`` flag;
    the others do not."""
    frames = read_layers_frames(chain)
    genes = frames["genes"].copy()
    assert genes is not None and not genes["crosses_origin"].any()
    n = int(chain.frame("replicons").set_index("contig").loc["chr1", "length"])
    k = int(genes.loc[genes["contig"] == "chr1", "start"].idxmax())
    tag = str(genes.at[k, "locus_tag"])
    genes.at[k, "end"] = n + 30  # the last chr1 gene extended across the origin
    genes.at[k, "crosses_origin"] = True
    frames["genes"] = genes
    df = call_gene_layers(frames)
    assert "crosses_origin" in [str(x) for x in _row(df, tag)["flags"]]
    others = [r for r in df.itertuples(index=False) if r.locus_tag != tag]
    assert all("crosses_origin" not in [str(x) for x in r.flags] for r in others)


def test_terminator_downstream_wraps_to_the_first_gene_on_circular(chain: Chain) -> None:
    """On a circular replicon the last ``+`` gene's next same-strand gene is the first ``+``
    gene (across the origin): a terminator just before that gene counts, although it lies
    beyond the 300-bp no-next-gene window."""
    from gippa.annotate.terminators import terminator_conf_downstream

    frames = read_layers_frames(chain)
    genes, terms = frames["genes"], frames["terminators"]
    assert genes is not None and terms is not None
    reps = orc.replicons(chain.frame("replicons"))
    rep = reps["chr1"]
    assert rep.topology == "circular"
    plus = genes[(genes["contig"] == "chr1") & (genes["strand"] == "+")].sort_values(
        ["start", "end", "locus_tag"]
    )
    last, first = plus.iloc[-1], plus.iloc[0]
    assert int(first["start"]) >= 40
    assert int(first["start"]) + rep.length - 40 >= int(last["end"]) + 300, "gap too small"
    extra = pd.DataFrame(
        [{**{c: None for c in terms.columns}, "contig": "chr1", "start": int(first["start"]) - 40,
          "end": int(first["start"]) - 10, "strand": "+", "confidence": 0.987654,
          "source": "TerminatorNet"}]
    ).astype(terms.dtypes.to_dict())  # fmt: skip
    frames["terminators"] = extra  # the only call: no other terminator decides the max
    gs = {g.locus_tag: g for g in orc.genes(genes)}
    tobjs = orc.terminators(frames["terminators"])
    want = terminator_conf_downstream(gs[last["locus_tag"]], gs[first["locus_tag"]], tobjs, rep)
    no_next = terminator_conf_downstream(gs[last["locus_tag"]], None, tobjs, rep)
    assert want is not None and want >= 0.987654 and not orc.same(want, no_next)
    got = _row(call_gene_layers(frames), str(last["locus_tag"]))["terminator_conf_downstream"]
    assert orc.same(got, want), (got, want, no_next)


def test_assign_tus_wraps_upstream_regions_across_the_origin(chain: Chain) -> None:
    """Precomputed TSSs (no tu_id) go to every same-strand TU whose upstream region holds them,
    modulo the length on circular replicons (an upstream region crossing the origin)."""
    from gippa.annotate.layers import assign_tus, tus_of
    from gippa.annotate.precomputed import load_precomputed

    reps = orc.replicons(chain.frame("replicons"))
    tus = chain.frame("tus")
    tx = load_precomputed("tx_rates", chain.sim.precomputed["tx_rates"], reps)
    wraps = [t for t in tus.itertuples(index=False)
             if reps[t.contig].topology == "circular"
             and (t.upstream_start < 0 or t.upstream_end > reps[t.contig].length)]  # fmt: skip
    assert wraps, "no TU upstream region crosses the origin on sim100k"
    extra = []
    for t in wraps:
        n = reps[t.contig].length
        for pos in (int(t.upstream_start), int(t.upstream_end) - 1):
            extra.append({"contig": t.contig, "strand": t.strand, "tss": pos % n, "tx_rate": 2.5})
    tx = pd.concat([tx, pd.DataFrame(extra)], ignore_index=True).astype(tx.dtypes.to_dict())
    got = assign_tus(tx, tus_of(tus), reps)
    want = set()
    for t in tus.itertuples(index=False):
        rep = reps[t.contig]
        ks = (-1, 0, 1) if rep.topology == "circular" else (0,)
        sel = tx[(tx["contig"] == t.contig) & (tx["strand"] == t.strand)]
        for s, r in zip(sel["tss"], sel["tx_rate"], strict=True):
            if any(t.upstream_start <= int(s) + k * rep.length < t.upstream_end for k in ks):
                want.add((t.contig, t.strand, int(s), float(r), t.tu_id))
    cols = got[["contig", "strand", "tss", "tx_rate", "tu_id"]].itertuples(index=False)
    have = {(str(c), str(st), int(s), float(r), str(u)) for c, st, s, r, u in cols}
    assert have == want
    assert {t.tu_id for t in wraps} <= set(got["tu_id"])


def _pairs(chain: Chain, source: Any) -> list[Any]:
    from gippa.annotate.uniop_adapter import PairProb

    p = chain.frame("pairs_raw")
    return [PairProb(str(a), str(b), float(x), source(str(a), str(b)))
            for a, b, x in p[["gene_a", "gene_b", "uniop_p"]].itertuples(index=False)]  # fmt: skip


def _tu_pairs(tu: Any) -> list[tuple[str, str]]:
    """The TU's adjacent pairs as pairs_raw keys (genome order: gene_a before gene_b)."""
    order = tu.genes if tu.strand == "+" else tu.genes[::-1]
    return list(itertools.pairwise(order))


@pytest.mark.parametrize("strand", ["+", "-"])
def test_tu_sources_any_distance_model_link_marks_the_tu(chain: Chain, strand: str) -> None:
    """A multi-gene TU with one distance-model link among UniOP links is ``distance_model``;
    pairs_raw keys are in genome order, a ``-`` TU lists its genes in transcription order."""
    from gippa.annotate.layers import tu_sources, tus_of

    tus = tus_of(chain.frame("tus"))
    big = next(t for t in tus if len(t.genes) >= 3 and t.strand == strand)
    keys = {(str(a), str(b)) for a, b in chain.frame("pairs_raw")[["gene_a", "gene_b"]].values}
    assert set(_tu_pairs(big)) <= keys
    dm = _tu_pairs(big)[-1]
    pairs = _pairs(chain, lambda a, b: "distance_model" if (a, b) == dm else "uniop")
    src = tu_sources(tus, pairs)
    assert src[big.tu_id] == "distance_model"
    assert all(v == "UniOP" for k, v in src.items() if k != big.tu_id)


def test_tu_sources_singletons_follow_the_run(chain: Chain) -> None:
    """A singleton TU has no link: ``distance_model`` when every pair of the run comes from the
    distance model (UniOP absent), ``UniOP`` otherwise."""
    from gippa.annotate.layers import tu_sources, tus_of

    tus = tus_of(chain.frame("tus"))
    singles = [t.tu_id for t in tus if len(t.genes) == 1]
    assert singles
    all_dm = tu_sources(tus, _pairs(chain, lambda a, b: "distance_model"))
    assert all(all_dm[t] == "distance_model" for t in singles)
    p0 = chain.frame("pairs_raw").iloc[0]
    one = (str(p0["gene_a"]), str(p0["gene_b"]))
    src = tu_sources(
        tus, _pairs(chain, lambda a, b: "distance_model" if (a, b) == one else "uniop")
    )
    assert all(src[t] == "UniOP" for t in singles)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_tu_sources_all_distance_model_run(chain: Chain, strand: str) -> None:
    """Every pair from the distance model (UniOP absent): every TU of either strand is
    ``distance_model``."""
    from gippa.annotate.layers import tu_sources, tus_of

    tus = tus_of(chain.frame("tus"))
    multi = [t.tu_id for t in tus if len(t.genes) > 1 and t.strand == strand]
    assert multi
    src = tu_sources(tus, _pairs(chain, lambda a, b: "distance_model"))
    assert {src[t] for t in multi} == {"distance_model"}


def test_layer_flags_precomputed_tx_rates_with_promoter_calculator_off(chain: Chain) -> None:
    """``promoter_calculator.enabled=false`` with ``annotation.precomputed.tx_rates`` set: the
    layer is present (precomputed), so no ``layer_missing:tx_rates``."""
    from gippa.annotate.layers import layer_flags

    pre = {"tx_rates": str(chain.sim.precomputed["tx_rates"])}
    cfg = GippaConfig.model_validate(
        {"seed": 1, "annotation": {"promoter_calculator": {"enabled": False}, "precomputed": pre}}
    )
    assert "layer_missing:tx_rates" not in layer_flags(cfg, set(ALL_KEYS))
    off = GippaConfig.model_validate(
        {"seed": 1, "annotation": {"promoter_calculator": {"enabled": False}}}
    )
    assert "layer_missing:tx_rates" in layer_flags(off, set(ALL_KEYS))
