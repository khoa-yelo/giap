"""Builders for the M10-07 free-entrapment tests (docs/specs/M10-07-free-entrapment.md).

Everything is synthetic and hand-made from the spec's examples: K-12-like proteins, a
contaminant protein, a Nissle-like genome back-translated from the spec's translation texts with
fixed table-11 codons and stop-codon spacers, a ``strain_db`` protein list and A1-like discovery
rows. No network, nothing from ``$GIPPA_ROOT`` or ``$GIPPA_DB_DIR``, nothing derived from a DIA-NN
report.

The genome of the shared fixture (examples 1-2) holds ``MSAGLLKVEPTLDEKAAGGYWR`` and
``AGGGVLTEPWR`` in frame +1 of one contig, each flanked by ``TAA TAG TGA``; the other frames of
the spacers read stops, ``N``/``I`` repeats and similar short texts that match no fixture
peptide (checked against ``match_exact`` / ``match_tolerant`` when the fixture was written).
"""

from __future__ import annotations

import functools
import gzip
import hashlib
import importlib
import json
import os
import subprocess
import sys
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
ENTRAPMENT_PY = REPO / "validation" / "pilot" / "entrapment.py"
METHODS_NOTE = REPO / "docs" / "methods" / "free-entrapment.md"

# ------------------------------------------------------------------ the shared fixture (ex. 1-2)

P1 = "MSAGLLKVEPTIDEKAAGGWWR"
P2 = "AGNVLTEPWR"
P3 = "QQWWLLDR"
C1 = "QQWWLLDRPEPTIDES"
K12 = (P1, P2, P3)
CONTAMINANTS = (C1,)
#: Nissle's translation texts of the shared fixture
NISSLE = ("MSAGLLKVEPTLDEKAAGGYWR", "AGGGVLTEPWR")

SHARED_1 = {"MSAGLLK", "VEPTLDEK", "MSAGLLKVEPTLDEK", "AGNVLTEPWR"}
K12_ONLY_1 = {"AAGGWWR", "VEPTLDEKAAGGWWR"}
#: example 1's precursors per universe element
N_PREC_1 = {"MSAGLLK": 4, "VEPTLDEK": 3, "MSAGLLKVEPTLDEK": 8, "AGNVLTEPWR": 3,
            "AAGGWWR": 2, "VEPTLDEKAAGGWWR": 4}  # fmt: skip

#: one fixed codon per residue (table 11); I is never written (the genome reads L)
CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
         "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG",
         "R": "CGT", "S": "AGC", "T": "ACC", "V": "GTG", "W": "TGG", "Y": "TAT"}  # fmt: skip
SPACER = "TAATAGTGA"
_COMP = str.maketrans("ACGT", "TGCA")

CAVEATS = (
    "near-homolog entrapment: the estimate errs high (§15.5)",
    "an estimated upper bound in expectation (Wen 2025, Suppl. Note 2.2), not a bound on this "
    "sample",
    "assumes the Nissle genome is complete, is the sample's strain, and that no K-12 material is "
    "in the sample",
    "assessment is made at precursor level, the unit DIA-NN's global q controls",
)
PEPTIDE_NOT_ASSESSED = (
    "peptide-level FDR is not controlled by the search (q_value is DIA-NN's precursor-level "
    "Global.Q.Value)"
)
#: the statistician's erratum (DIA-NN 2.3.1 also filters the main report at run-level Q.Value)
PER_RUN_NOT_ASSESSED = (
    "runs were selected by the global and the run-level q-value together; no FDR claim attaches "
    "to a run's intersection"
)
CI_METHOD = "binomial_inversion_equal_chance"
CI_ASSUMES = (
    "equal chance and independent false matches; not valid under clustered (protein-level) errors"
)


def back_translate(protein: str) -> str:
    return "".join(CODON[c] for c in protein)


def coding_text(*proteins: str) -> str:
    """The proteins in one reading frame, each flanked by stop spacers."""
    return SPACER + SPACER.join(back_translate(p) for p in proteins) + SPACER


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def replicons(genome: Mapping[str, str], circular: Iterable[str] = ()) -> dict[str, Any]:
    from gippa.model.coords import Replicon

    circ = set(circular)
    return {
        c: Replicon(c, len(s), "circular" if c in circ else "linear") for c, s in genome.items()
    }


def origin_split(protein: str, cut: int) -> str:
    """A contig whose end holds ``protein[:cut]`` and whose start holds ``protein[cut:]``."""
    head, tail = back_translate(protein[:cut]), back_translate(protein[cut:])
    return tail + SPACER + "GCTGCT" + SPACER * 3 + head


# ------------------------------------------------------------------ discoveries (example 2)

ROW_COLS = ("run", "peptide", "modified_peptide", "charge", "q_value", "is_decoy", "protein_ids")


def row(run: str, pep: str, charge: int, q: float | None, decoy: bool = False, *,
        mod: str | None = None, ids: str = "sp|P0K12|X_ECOLI") -> dict[str, Any]:  # fmt: skip
    return {"run": run, "peptide": pep, "modified_peptide": mod or pep, "charge": charge,
            "q_value": q, "is_decoy": decoy, "protein_ids": ids}  # fmt: skip


EXAMPLE2_ROWS: tuple[dict[str, Any], ...] = (
    row("WT_1", "MSAGLLK", 2, 0.001),
    row("WT_2", "MSAGLLK", 2, 0.003, mod="M(UniMod:35)SAGLLK"),
    row("WT_1", "VEPTIDEK", 2, 0.004),
    row("WT_3", "AAGGWWR", 2, 0.008),
    row("WT_1", "VEPTIDEKAAGGWWR", 3, 0.02),
    row("WT_2", "QQWWLLDR", 2, 0.001, ids="Cont_P00761"),
    row("WT_1", "AGNVLTEPWR", 2, 0.01),
    row("WT_2", "KEDITPEPR", 2, 0.001, True, ids="rev_sp|P0K12|X_ECOLI"),
)
#: example 8's extra accepted target row (a Met-excised N-terminal product)
EXAMPLE8_ROW = row("WT_1", "SAGLLKVEPTIDEK", 2, 0.002)


def frame(rows: Sequence[Mapping[str, Any]]) -> Any:
    import pandas as pd

    return pd.DataFrame(list(rows), columns=list(ROW_COLS))


def accepted_strings(rows: Sequence[Mapping[str, Any]]) -> set[str]:
    """Every peptide and modified-peptide string of the accepted target rows, literal and I->L."""
    out: set[str] = set()
    for r in rows:
        q = r["q_value"]
        if r["is_decoy"] or q is None or q > 0.01:
            continue
        for s in (str(r["peptide"]), str(r["modified_peptide"])):
            out |= {s, s.replace("I", "L")}
    return out


# ------------------------------------------------------------------ module and CLI


@functools.cache
def _load() -> Any:
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    return importlib.import_module("validation.pilot.entrapment")


def load_entrapment() -> Any:
    """``validation.pilot.entrapment``; fails the calling test (not an import error at
    collection) while the module does not exist."""
    import pytest

    if not ENTRAPMENT_PY.is_file():
        pytest.fail(f"{ENTRAPMENT_PY.relative_to(REPO)} does not exist (M10-07 not implemented)")
    return _load()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_entrapment(
    *args: str, cwd: Path | None = None, env_extra: Mapping[str, str] | None = None
) -> subprocess.CompletedProcess[str]:
    """``python validation/pilot/entrapment.py ARGS`` with ``PYTHONPATH=$REPO/src:$REPO``
    (``env_extra`` is added to the environment)."""
    import pytest

    if not ENTRAPMENT_PY.is_file():
        pytest.fail(f"{ENTRAPMENT_PY.relative_to(REPO)} does not exist (M10-07 not implemented)")
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(REPO / "src"), str(REPO), str(REPO / "tests" / "sim")])
    env |= dict(env_extra or {})
    return subprocess.run(
        [sys.executable, str(ENTRAPMENT_PY), *args],
        capture_output=True,
        text=True,
        env=env,
        cwd=str(cwd) if cwd else None,
        check=False,
    )


def write_faa(path: Path, proteins: Sequence[str], prefix: str = "sp|P") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(
        f">{prefix}{i:05d}|X_ECOLI test protein {i}\n{p}\n" for i, p in enumerate(proteins)
    )
    path.write_text(text, encoding="utf-8")
    return path


def write_genome(path: Path, contigs: Mapping[str, str], circular: Iterable[str]) -> Path:
    circ = set(circular)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for c, s in contigs.items():
        tag = " [topology=circular]" if c in circ else ""
        lines += [f">{c} synthetic Nissle-like contig{tag}", s]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_a1(path: Path, rows: Sequence[Mapping[str, Any]]) -> Path:
    """``a1.tsv.gz`` in the generic peptide-table layout (M10-03), gzip with a fixed mtime."""
    from helpers.pilot_metrics_m10 import pep_row, write_peptides

    out = []
    for r in rows:
        d = pep_row(str(r["peptide"]), str(r["protein_ids"]), q=r["q_value"],
                    decoy=bool(r["is_decoy"]), run=str(r["run"]), search_db="other")  # fmt: skip
        d["modified_peptide"] = r["modified_peptide"]
        d["charge"] = r["charge"]
        out.append(d)
    return write_peptides(path, out)


def write_runs_manifest(path: Path, runs: Sequence[str] = ("WT_1", "WT_2", "WT_3")) -> Path:
    cols = ["run", "file_name", "url", "url_record", "bytes", "checksum_record", "sha256",
            "sample", "selected_because"]  # fmt: skip
    lines = ["\t".join(cols)]
    for r in runs:
        lines.append(
            "\t".join([r, f"{r}.d.zip", "NA", "NA", "1", "NA", "0" * 64, "NA", "token:WT"])
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def diann_step1_tokens(cfg: Path, k12: Path, cont: Path, **override: str) -> list[str]:
    """``command_step1.txt`` tokens built from ``diann_A1.cfg`` ``[step1]`` (one per line),
    with the two ``--fasta`` paths expanded and ``override`` replacing a flag's value."""
    tokens = ["apptainer", "exec", "diann_2.3.1.sif", "/diann-2.3.1/diann-linux"]
    section = None
    fastas = iter([str(k12), str(cont)])
    for line in cfg.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if s.startswith("["):
            section = s
            continue
        if section != "[step1]":
            continue
        flag, _, value = s.partition(" ")
        if flag == "--fasta":
            value = next(fastas)
        if flag in override:
            value = override[flag]
        tokens += [flag, value] if value else [flag]
    return [*tokens, "--threads", "16", "--out-lib", "/scratch/A1.parquet"]


def write_out_a1(root: Path, exact: Sequence[str], gi_count: int, gi_n: int) -> Path:
    """An A1 run directory with ``tables/peptide_hits.tsv`` (one ``route = exact`` target row per
    peptide of ``exact``), an empty ``tables/unmapped.tsv`` and ``tables/detection_summary.tsv``
    with ``genome_inconsistent_rate`` = ``gi_count / gi_n``."""
    from types import SimpleNamespace

    from gippa.model.coords import Interval
    from helpers.pilot_metrics_m10 import hit_row, write_table

    g = SimpleNamespace(interval=Interval("nissle", 9, 300, "+"), locus_tag="NISSLE_00001")
    write_table(root / "tables" / "peptide_hits.tsv", "peptide_hits",
                [hit_row(p, g) for p in exact])  # fmt: skip
    write_table(root / "tables" / "unmapped.tsv", "unmapped", [])
    ds = {"section": "qc", "group_type": "all", "group": "all",
          "metric": "genome_inconsistent_rate", "value": gi_count / gi_n, "n": gi_n}  # fmt: skip
    write_table(root / "tables" / "detection_summary.tsv", "detection_summary", [ds])
    return root


def read_json(path: Path) -> dict[str, Any]:
    doc: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return doc


def gz_text(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return fh.read()


def codes(warnings: Sequence[Any]) -> list[str]:
    """Warning codes (a warning is a string or a dict with ``code``; M10-05's convention)."""
    return [w if isinstance(w, str) else str(w.get("code")) for w in warnings]


def reasons(block: Mapping[str, Any]) -> set[str]:
    """Every string value under a key containing ``reason`` in ``block`` or one level below."""
    out: set[str] = set()
    for k, v in block.items():
        if "reason" in k and isinstance(v, str):
            out.add(v)
        elif isinstance(v, Mapping):
            out |= {w for kk, w in v.items() if "reason" in kk and isinstance(w, str)}
    return out


def flagged(block: Mapping[str, Any], flag: str, warnings: Sequence[Any] = ()) -> bool:
    """``flag`` is set on ``block`` (``block[flag] is True`` or listed in ``block["flags"]``)
    or is a warning code."""
    flags = block.get("flags") or []
    return block.get(flag) is True or flag in flags or flag in codes(warnings)
