"""Codon usage, lead decisions after code review (PR #86), items 4-8 of
docs/specs/M2-17-codon-usage.md.

4. Run flags are stored in the ``codon.parquet`` file metadata, key ``gippa_flags`` (JSON list).
5. Translation tables 4, 11 and 25 are accepted (M11-16 item 6.1, superseding "only 11"); any
   other table (1, 2, 12) raises ``ValueError``.
6. ``gc3`` returns ``None`` when no third position is A/C/G/T.
7. codon-bias RuntimeWarnings (genes with nothing to weight) are suppressed and logged once as
   a count.
8. Step-level tests that kill the mutants that survived the reviewer's mutation run
   (local/scratch/review-m2-17/mut/mutants.py): ``step_cairef_unscored``,
   ``step_gene_partial_false``, ``step_table_hard1``, ``any_table``, ``shares_no_upper``.
   (``shares_keep_stops`` is an equivalent mutant: stop counts never enter a sense family.)

Resolutions (test author):
* ``gippa_flags`` is always present; it is ``[]`` when the run has no flag, so readers need no
  special case. Reruns stay byte-identical.
* Item 7's single log record is emitted on a ``gippa.*`` logger, mentions "warning" and the count.
* Lower-case nucleotide sequences give exactly the same ``codon.parquet`` as upper-case ones.
"""

from __future__ import annotations

import dataclasses
import importlib
import json
import logging
import random
import re
import warnings
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from helpers.codon_excerpt import write_inputs

from gippa.config import GippaConfig
from gippa.model.coords import Interval
from gippa.model.genes import Gene

PARSE_OUT = {
    "genes": "genes.parquet",
    "features": "features.parquet",
    "replicons": "replicons.parquet",
    "rrna16s": "rrna16s.fa",
    "trna_counts": "trna_counts.tsv",
    "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json",
}
RPLB = "KOIAIB_03524"
RPSA = "KOIAIB_00982"
FLAG = "cai_reference_missing"
TRNA = {"CAT": 3, "CAG": 4, "GCC": 2, "TTC": 2, "GTC": 3, "TGC": 2, "GAA": 1, "CGT": 2}  # no CCA


def cu() -> ModuleType:
    return importlib.import_module("gippa.annotate.codon_usage")


_POS = [0]


def make_gene(locus: str, nt: str, **kw: Any) -> Gene:
    start = _POS[0]
    _POS[0] += len(nt) + 100
    base: dict[str, Any] = {"gene": None, "product": "hypothetical protein",
                            "feature_type": "cds", "partial": False}  # fmt: skip
    base.update(kw)
    return Gene(
        locus_tag=locus,
        interval=Interval("c1", start, start + len(nt), "+"),
        aa="M",
        nt=nt,
        start_codon=nt[:3],
        **base,
    )


def random_cds(rnd: random.Random, n: int) -> str:
    sense = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
    sense = [c for c in sense if c not in ("TAA", "TAG", "TGA")]
    return "ATG" + "".join(rnd.choices(sense, k=n - 2)) + "TAA"


def gene_set(seed: int = 3) -> list[Gene]:
    rnd = random.Random(seed)
    genes = [
        make_gene(f"R_{k}", random_cds(rnd, 120), gene=f"rpl{chr(65 + k)}",
                  product=f"50S ribosomal protein L{k + 1}")
        for k in range(6)
    ]  # fmt: skip
    genes += [make_gene(f"O_{k}", random_cds(rnd, rnd.randint(40, 300))) for k in range(20)]
    return genes


# ------------------------------------------------------------------ step harness


def _fresh_steps(mp: pytest.MonkeyPatch) -> ModuleType:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    names = load_step_modules()
    mp.setattr(reg, "_REGISTRY", {})
    for n in names:
        if n != "registry":
            importlib.reload(importlib.import_module(f"gippa.steps.{n}"))
    return reg


@pytest.fixture
def steps() -> Iterator[ModuleType]:
    with pytest.MonkeyPatch.context() as mp:
        yield _fresh_steps(mp)


@pytest.fixture
def work(steps: ModuleType, tmp_path: Path) -> dict[str, Path]:
    js, fa = write_inputs(tmp_path / "in")
    outs = {k: tmp_path / "work" / v for k, v in PARSE_OUT.items()}
    steps.run_step("parse_bakta", GippaConfig(seed=1), {"bakta_json": js, "genome": fa}, outs)
    return outs


def _cfg(table: int = 11) -> GippaConfig:
    return GippaConfig(seed=1, annotation={"translation_table": table})


def _run(
    reg: ModuleType, genes: Path, trna: Path, out: Path, cfg: GippaConfig | None = None
) -> Path:
    reg.run_step(
        "codon_usage", cfg or _cfg(), {"genes": genes, "trna_counts": trna}, {"codon": out}
    )
    return out


def _edit_genes(src: Path, dst: Path, edit: Any) -> Path:
    """Copy genes.parquet with ``edit(row)`` applied to every row (same schema)."""
    table = pq.read_table(src)
    rows = [edit(dict(r)) for r in table.to_pylist()]
    pq.write_table(pa.Table.from_pylist(rows, schema=table.schema), dst)
    return dst


def _rows(path: Path) -> dict[str, dict[str, Any]]:
    return {r["locus_tag"]: r for r in pq.read_table(path).to_pylist()}


def _flags(path: Path) -> list[str]:
    meta = pq.read_metadata(path).metadata or {}
    assert b"gippa_flags" in meta, f"codon.parquet metadata keys: {sorted(meta)}"
    value = json.loads(meta[b"gippa_flags"])
    assert isinstance(value, list)
    return value


# ------------------------------------------------------------------ 4. run flags in metadata


def test_flags_metadata_lists_cai_reference_missing(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path
) -> None:
    from gippa.annotate.bakta_layers import is_ribosomal_protein

    def strip(r: dict[str, Any]) -> dict[str, Any]:
        g = make_gene("x", "ATG", gene=r["gene"], product=r["product"] or "")
        if is_ribosomal_protein(g):
            r["gene"], r["product"] = None, "hypothetical protein"
        return r

    genes = _edit_genes(work["genes"], tmp_path / "g.parquet", strip)
    out = _run(steps, genes, work["trna_counts"], tmp_path / "c.parquet")
    assert _flags(out) == [FLAG]
    assert all(r["cai"] is None for r in _rows(out).values())


def test_flags_metadata_is_empty_list_without_flags(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path
) -> None:
    a = _run(steps, work["genes"], work["trna_counts"], tmp_path / "a.parquet")
    b = _run(steps, work["genes"], work["trna_counts"], tmp_path / "b.parquet")
    assert _flags(a) == []
    assert a.read_bytes() == b.read_bytes()


# ------------------------------------------------------------------ 5. tables 4, 11 and 25 only


@pytest.mark.parametrize("table", [1, 2, 12])
def test_codon_usage_rejects_unsupported_tables(table: int) -> None:
    with pytest.raises(ValueError, match="11"):
        cu().codon_usage(gene_set(), TRNA, table)
    with pytest.raises(ValueError):
        cu().synonymous_shares(["ATGGCCTAA"], table)


def test_table_11_is_accepted() -> None:
    feats, _ = cu().codon_usage(gene_set(), TRNA, 11)
    assert feats


@pytest.mark.parametrize("table", [1, 12])
def test_step_reads_the_configured_table(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path, table: int
) -> None:
    """Kills ``step_table_hard1``: the step must pass ``annotation.translation_table``."""
    from gippa.steps.registry import StepError

    with pytest.raises(StepError) as ei:
        _run(steps, work["genes"], work["trna_counts"], tmp_path / "c.parquet", _cfg(table))
    assert isinstance(ei.value.__cause__, ValueError)
    assert not (tmp_path / "c.parquet").exists()


# ------------------------------------------------------------------ 6. gc3 without ACGT thirds


@pytest.mark.parametrize("seq", ["", "NNN", "ACNGTNAAR", "acn"])
def test_gc3_none_without_acgt_third_positions(seq: str) -> None:
    assert cu().gc3(seq) is None


def test_gc3_lower_case_and_partial_ambiguity() -> None:
    assert cu().gc3("gccgcggca") == pytest.approx(2 / 3)
    assert cu().gc3("GCNGCG") == 1.0


def test_codon_usage_gene_without_acgt_thirds_gets_null_gc3() -> None:
    genes = [*gene_set(), make_gene("N_1", "ATG" + "AAN" * 12 + "TAA")]
    feats, _ = cu().codon_usage(genes, TRNA)
    assert feats["N_1"].gc3 is None
    assert feats["O_0"].gc3 is not None


# ------------------------------------------------------------------ 7. RuntimeWarnings


def _nothing_to_weight() -> list[Gene]:
    # TRNA has no CCA anticodon: an all-Trp gene has no tAI weight (codon-bias divides 0 by 0)
    return [*gene_set(), make_gene("W_1", "ATG" + "TGG" * 12 + "TAA"),
            make_gene("W_2", "ATG" + "TGG" * 15 + "TAA")]  # fmt: skip


def test_codon_usage_emits_no_runtime_warning() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        feats, _ = cu().codon_usage(_nothing_to_weight(), TRNA)
    assert feats["W_1"].tai is None and feats["W_2"].tai is None
    assert feats["O_0"].tai is not None


def test_step_logs_suppressed_warnings_once_as_a_count(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    rows = pq.read_table(work["genes"]).to_pylist()
    victims = [r["locus_tag"] for r in rows if r["feature_type"] == "cds"][-3:]

    def trp(r: dict[str, Any]) -> dict[str, Any]:
        if r["locus_tag"] in victims:
            r["nt"] = "ATG" + "TGG" * 12 + "TAA"
        return r

    genes = _edit_genes(work["genes"], tmp_path / "g.parquet", trp)
    trna = tmp_path / "trna_counts.tsv"
    trna.write_text(
        "anticodon\tamino_acid\tn_copies\n"
        + "".join(f"{a}\tX\t{n}\n" for a, n in sorted(TRNA.items()))
    )
    caplog.set_level(logging.DEBUG, logger="gippa")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        out = _run(steps, genes, trna, tmp_path / "c.parquet")
    got = _rows(out)
    assert all(got[v]["tai"] is None for v in victims)
    hits = [
        r for r in caplog.records
        if r.name.startswith("gippa") and re.search("warning", r.getMessage(), re.I)
    ]  # fmt: skip
    assert len(hits) == 1, [r.getMessage() for r in hits]
    assert re.search(r"\b3\b", hits[0].getMessage()), hits[0].getMessage()


# ------------------------------------------------------------------ 8. surviving mutants


@pytest.mark.parametrize(
    "edit",
    [{"partial": True}, {"feature_type": "pseudogene"}, {"nt": "ATG" + "GCC" * 7 + "TAA"}],
    ids=["partial", "pseudogene", "too_short"],
)
def test_step_unscored_ribosomal_gene_is_not_cai_reference(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path, edit: dict[str, Any]
) -> None:
    """Kills ``step_cairef_unscored`` and (partial) ``step_gene_partial_false``."""

    def apply(r: dict[str, Any]) -> dict[str, Any]:
        if r["locus_tag"] == RPLB:
            r.update(edit)
        return r

    genes = _edit_genes(work["genes"], tmp_path / "g.parquet", apply)
    got = _rows(_run(steps, genes, work["trna_counts"], tmp_path / "c.parquet"))
    assert got[RPLB]["cai_reference"] is False
    assert got[RPLB]["cai"] is None and got[RPLB]["enc"] is None
    assert got[RPSA]["cai_reference"] is True and got[RPSA]["cai"] is not None


def test_step_reads_the_real_partial_value(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path
) -> None:
    """Kills ``step_gene_partial_false``: a partial non-ribosomal gene gets no scores."""
    rows = pq.read_table(work["genes"]).to_pylist()
    lt = next(r["locus_tag"] for r in rows
              if r["feature_type"] == "cds" and r["gene"] == "lacZ")  # fmt: skip
    base = _rows(_run(steps, work["genes"], work["trna_counts"], tmp_path / "a.parquet"))
    assert base[lt]["cai"] is not None

    def part(r: dict[str, Any]) -> dict[str, Any]:
        if r["locus_tag"] == lt:
            r["partial"] = True
        return r

    genes = _edit_genes(work["genes"], tmp_path / "g.parquet", part)
    got = _rows(_run(steps, genes, work["trna_counts"], tmp_path / "b.parquet"))
    assert all(got[lt][k] is None for k in ("cai", "tai", "enc", "fop", "gc3", "cpb"))
    assert got[lt]["at_2_10"] == base[lt]["at_2_10"]


def test_synonymous_shares_ignore_case() -> None:
    """Kills ``shares_no_upper``."""
    seqs = ["CTGCTGCTACTCGCCGCA", "TTAAAAAAG"]
    up = cu().synonymous_shares(seqs)
    assert up and cu().synonymous_shares([s.lower() for s in seqs]) == up
    assert cu().synonymous_shares(["CtGcTa"]) == cu().synonymous_shares(["CTGCTA"])


def test_codon_usage_ignores_case() -> None:
    genes = gene_set(seed=5)
    lower = [dataclasses.replace(g, nt=(g.nt or "").lower()) for g in genes]
    a, fa = cu().codon_usage(genes, TRNA)
    b, fb = cu().codon_usage(lower, TRNA)
    assert a == b and fa == fb


def test_step_lower_case_sequences_give_the_same_codon_parquet(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path
) -> None:
    """Kills ``shares_no_upper`` at step level (rare codons come from the genome-wide shares)."""

    def low(r: dict[str, Any]) -> dict[str, Any]:
        r["nt"] = r["nt"].lower() if r["nt"] else r["nt"]
        return r

    genes = _edit_genes(work["genes"], tmp_path / "g.parquet", low)
    a = _run(steps, work["genes"], work["trna_counts"], tmp_path / "a.parquet")
    b = _run(steps, genes, work["trna_counts"], tmp_path / "b.parquet")
    ra, rb = _rows(a), _rows(b)
    assert sum(r["rare_codons_2_10"] or 0 for r in ra.values()) > 0
    assert ra == rb
