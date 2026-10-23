"""Invariants of ``scripts/dev/mvp04_summary.py`` (MVP-04 part 2; interface pinned in
tests/unit/test_mvp04_summary_gates.py).

* status: ``fail`` iff the gate is not met, whatever the expected range says;
* gate 9: a directory compared with its copy is identical; any one-byte change in a compared
  table is reported, and only that table;
* gate 10: the three shares lie in [0, 1] and sum to 1; ``exact`` counts each peptide once;
* gate 5: decoy wins <= paired eligible genes <= eligible genes, and dropping a decoy row never
  adds a win;
* gate 11: precision and recall lie in [0, 1]; adding a false ``high`` row never raises precision.
"""

from __future__ import annotations

import importlib.util
import shutil
import sys
import tempfile
from pathlib import Path
from types import ModuleType
from typing import Any

from helpers import mvp04_summary_data as S
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "dev" / "mvp04_summary.py"
SETTINGS = settings(max_examples=40, deadline=None,
                    suppress_health_check=[HealthCheck.function_scoped_fixture])  # fmt: skip


def mod() -> ModuleType:
    spec = importlib.util.spec_from_file_location("mvp04_summary_props", SCRIPT)
    assert spec is not None and spec.loader is not None, SCRIPT
    m = importlib.util.module_from_spec(spec)
    sys.modules["mvp04_summary_props"] = m
    spec.loader.exec_module(m)
    return m


@given(st.booleans(), st.booleans())
def test_status_fail_iff_gate_not_met(ok: bool, inside: bool) -> None:
    assert (mod().status_of(ok, inside) == "fail") is (not ok)


# ---------------------------------------------------------------- gate 9


def _inputs(root: Path, rows: list[str]) -> Path:
    root.mkdir(parents=True)
    (root / "peptides.tsv").write_text(
        "peptide\tq_value\n" + "".join(f"{r}\t0.001\n" for r in rows)
    )
    (root / "truth").mkdir()
    for name in S.TRUTH_TABLES:
        (root / "truth" / f"{name}.tsv").write_text(
            f"col\n{name}\n" + "".join(r + "\n" for r in rows)
        )
    return root


@SETTINGS
@given(rows=st.lists(st.text("ACDEFGHIKLMNPQRSTVWY", min_size=6, max_size=12), min_size=1,
                     max_size=8),
       which=st.sampled_from(["peptides.tsv", *(f"truth/{n}.tsv" for n in S.TRUTH_TABLES)]),
       data=st.data())  # fmt: skip
def test_gate9_reports_exactly_the_changed_table(rows: list[str], which: str, data: Any) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        a = _inputs(Path(tmp) / "a", rows)
        b = Path(tmp) / "b"
        shutil.copytree(a, b)
        assert mod().gate9(a, b)["value"] == []
        p = b / which
        raw = bytearray(p.read_bytes())
        i = data.draw(st.integers(0, len(raw) - 1))
        raw[i] = (raw[i] + 1) % 256
        p.write_bytes(bytes(raw))
        res = mod().gate9(a, b)
        assert res["value"] == [which] and res["status"] == "fail"


# ---------------------------------------------------------------- gate 10

ROUTES = st.lists(st.sampled_from(["exact", "tolerant"]), min_size=0, max_size=3)


@SETTINGS
@given(hits=st.lists(ROUTES, min_size=1, max_size=15), n_unmapped=st.integers(0, 5))
def test_gate10_shares_partition_the_peptides(hits: list[list[str]], n_unmapped: int) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        out, inp = Path(tmp) / "run", Path(tmp) / "inputs"
        rows, unmapped = [], []
        for k, routes in enumerate(hits):
            pep = S.peptide(k)
            if not routes:
                unmapped.append([pep, pep, "peptide", "unmapped", "no_hit", "NA"])
            rows += [[pep, r, 0, "cds_inframe", "unique", "G"] for r in routes]
        unmapped += [[S.peptide(100 + k), S.peptide(100 + k), "peptide", "unmapped", "no_hit",
                      "NA"] for k in range(n_unmapped)]  # fmt: skip
        S.tsv(
            out / "mapping" / "peptide_hits.tsv",
            ["peptide", "route", "n_diff", "locus_class", "uniqueness", "gene"],
            rows,
        )
        S.tsv(
            out / "mapping" / "unmapped.tsv",
            ["sequence_id", "sequence", "kind", "label", "reason", "contaminant_match"],
            unmapped,
        )
        S.tsv(inp / "peptides.tsv", ["peptide", "protein_ids", "q_value", "is_decoy"], [])
        res = mod().gate10(out, inp)
        sh = res["shares"]
        n = len(hits) + n_unmapped
        assert res["n_peptides"] == n
        assert all(0.0 <= sh[k] <= 1.0 for k in ("exact", "tolerant", "unmapped"))
        assert abs(sh["exact"] + sh["tolerant"] + sh["unmapped"] - 1.0) < 1e-12
        assert sh["exact"] * n == sum(1 for r in hits if "exact" in r)
        assert (res["gate_ok"] is True) == any(r and "exact" not in r for r in hits)


# ---------------------------------------------------------------- gate 5

GENE = st.tuples(st.sampled_from(["partial_excluded", "detected", "not_detected"]),
                 st.booleans(),  # gene q NA
                 st.sampled_from([None, 0.001, 0.05]))  # decoy q (None = no decoy)  # fmt: skip


def _gate5_dirs(root: Path, genes: list[tuple[str, bool, float | None]]) -> tuple[Path, Path]:
    out, inp = root / "run", root / "inputs"
    S.tsv(
        out / "tables" / "genes.tsv",
        ["locus_tag", "detection_status", "gene_q_value"],
        [[S.locus(i), s, "NA" if na else 0.001] for i, (s, na, _) in enumerate(genes)],
    )
    S.tsv(inp / "peptides.tsv", ["peptide", "protein_ids", "q_value", "is_decoy"],
          [[S.peptide(i), "rev_" + S.locus(i), q, "true"]
           for i, (_, _, q) in enumerate(genes) if q is not None])  # fmt: skip
    return out, inp


@SETTINGS
@given(genes=st.lists(GENE, min_size=1, max_size=20), drop=st.integers(0, 19))
def test_gate5_counts_are_nested(genes: list[tuple[str, bool, float | None]], drop: int) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        out, inp = _gate5_dirs(Path(tmp) / "a", genes)
        res = mod().gate5(out, inp)
        eligible = sum(1 for s, _, _ in genes if s != "partial_excluded")
        assert 0 <= res["value"] <= res["n_paired_genes"] <= eligible
        want = sum(1 for s, na, q in genes if s != "partial_excluded" and na and q == 0.001)
        assert res["value"] == want
        fewer = [g if i != drop else (g[0], g[1], None) for i, g in enumerate(genes)]
        out2, inp2 = _gate5_dirs(Path(tmp) / "b", fewer)
        assert mod().gate5(out2, inp2)["value"] <= res["value"]


# ---------------------------------------------------------------- gate 11

SITE = st.tuples(st.integers(0, 30), st.sampled_from("VLT"), st.integers(0, 3))


@SETTINGS
@given(truth=st.lists(SITE, min_size=1, max_size=10, unique_by=lambda s: s[0]),
       calls=st.lists(st.tuples(st.integers(0, 30), st.sampled_from("VLT")), min_size=1,
                      max_size=10, unique=True))  # fmt: skip
def test_gate11_rates_are_bounded_and_a_false_call_never_helps(
    truth: list[tuple[int, str, int]], calls: list[tuple[int, str]]
) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        out, inp = Path(tmp) / "run", Path(tmp) / "inputs"
        S.tsv(
            inp / "truth" / "alt_sites.tsv",
            [
                "contig",
                "codon_start",
                "strand",
                "gene",
                "ref_aa",
                "alt_aa",
                "n_true_peptides",
                "frameshift_upstream",
            ],
            [[S.CONTIG, 1 + 3 * c, "+", "G", "A", a, n, "false"] for c, a, n in truth],
        )

        def write(rows: list[tuple[int, str]]) -> dict[str, Any]:
            S.tsv(out / "tables" / "variants.tsv",
                  ["variant_id", "contig", "codon_start", "strand", "gene", "ref_aa", "alt_aa",
                   "call", "confidence", "flags"],
                  [[f"V{k}", S.CONTIG, 1 + 3 * c, "+", "G", "A", a, "ALT", "high",
                    "snv;strain_difference;database_derived"]
                   for k, (c, a) in enumerate(rows)])  # fmt: skip
            (out / "tables" / "variants.vcf").write_text(S.VCF_HEADER)
            res: dict[str, Any] = mod().gate11(out, inp)
            return res

        res = write(calls)
        for k in ("precision", "recall"):
            assert res[k] is None or 0.0 <= res[k] <= 1.0
        more = write([*calls, (31, "W")])  # codon 94: never a truth site
        if res["precision"] is not None:
            assert more["precision"] <= res["precision"]
