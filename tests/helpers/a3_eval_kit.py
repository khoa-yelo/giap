"""Synthetic ``.ssl`` files, synthetic genomes and an in-process driver for the M10-08b tests.

Spec: docs/specs/M10-08-a3-feasibility.md (part 08b: ``a3_feasibility.py evaluate`` and
``early-fail``, the verdict JSON, the failure rule L6, examples 1-10 and 12, the
``denovo_pilot.sbatch`` Cascadia branch). Builds on ``helpers.a3_kit`` (08a).

Nothing here reaches the network, a GPU, Cascadia, the real Nissle genome or ``$GIPPA_ROOT``. Every
``.ssl`` is written at test time (M10-09 guard: no ``.ssl`` is ever committed, not even a
synthetic one). Genomes are a few kb of random DNA with planted peptide codings.

Conventions settled by the test author (the spec names the subcommand but not its options; see
the report to the lead):

* ``a3_feasibility.py evaluate SSL --run RUN --out J`` is the lead-pinned minimal form
  (2026-10-24; M10-10's ``cascadia_convert.sbatch`` calls it so). Optional: ``--genome
  FNA[.gz] --bakta-json JSON --contaminants FASTA --cutoff X``. ``--genome`` and
  ``--bakta-json`` default to the Nissle files (spec step 8), ``--contaminants`` to
  ``$GIPPA_DB_DIR/contaminants/contaminants.fasta``. The tests pass all three except in the
  two minimal-form tests (which point ``$GIPPA_DB_DIR`` at ``tmp_path``). Other provenance
  options (route, weights, mzML check, job IDs) are the implementer's; the tests only require
  their keys in ``J``.
* ``early-fail --reason R --evidence TEXT --out J`` (exactly the spec's interface).
* ``main(argv) -> int`` in-process (the ``helpers.a3_kit`` convention); ``J`` is written for a
  verdict (exit 0 or 3) and not required on exit 1 or 2.
* Replicons come from the Bakta JSON's ``sequences`` records (``id``, ``type``, ``topology``,
  ``length``) with the nucleotides from the genome FASTA, as ``validation/proxies/te_proxies.py``
  reads them (``parse_bakta_json(json, genome=seqs).model.replicons``).
"""

from __future__ import annotations

import contextlib
import gzip
import importlib.util
import io
import json
import random
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType
from typing import Any

from helpers.a3_kit import A3_FEASIBILITY, require

SSL_COLUMNS = ("file", "scan", "charge", "sequence", "score-type", "score", "retention-time",
               "start-time", "end-time")  # fmt: skip

# Every field of validation/pilot/metrics/a3_feasibility.json (spec "aggregates only").
JSON_FIELDS = frozenset({
    "schema_version", "run", "cascadia_version", "weights_sha256",
    "route", "route_version", "mzml_sha256", "mzml_check",
    "genome_sha256", "cutoff", "cutoff_source",
    "n_rows", "n_below_cutoff", "n_len_7_30", "n_pass", "n_top", "n_match", "n_nonstandard",
    "share", "share_ci95", "decoy_share", "contaminant_share",
    "matches_by_replicon", "length_hist", "rule", "fails", "verdict", "reason",
    "evidence", "job_ids", "date",
})  # fmt: skip
COUNT_FIELDS = ("n_rows", "n_below_cutoff", "n_len_7_30", "n_pass", "n_top", "n_match",
                "n_nonstandard")  # fmt: skip
RULE = {"min_peptides": 500, "min_share": 0.20, "top_n": 1000}
REASONS = ("route_unavailable", "structurally_incompatible", "tool_incompatible")
VERDICT_REASONS = (None, "too_few_peptides", "low_genome_match", "both", *REASONS)

RESIDUES = "ACDEFGHKLMNPQRSTVWY"  # I folded into L; B J O U X Z never match
SEQ_LIKE = re.compile(r"^[A-Z]{7,}$")

# One table-11 codon per residue (none is a stop).
CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
    "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG",
    "R": "CGT", "S": "AGC", "T": "ACC", "V": "GTG", "W": "TGG", "Y": "TAT",
}  # fmt: skip
_COMP = str.maketrans("ACGTN", "TGCAN")


def coding(peptide: str) -> str:
    return "".join(CODON[a] for a in peptide)


def revcomp(nt: str) -> str:
    return nt.translate(_COMP)[::-1]


def fold(seq: str) -> str:
    """The spec's stripping rule (step 2) and I -> L, as an independent oracle."""
    s = re.sub(r"\[[^\]]*\]|\([^)]*\)", "", seq)
    s = re.sub(r"[+-]\d+(\.\d+)?", "", s)
    s = s.replace("-", "")
    s = re.sub(r"[^A-Z]", "", s)
    return s.upper().replace("I", "L")


def random_peptides(
    n: int, rng: random.Random, *, lengths: tuple[int, int] = (8, 14), avoid: Iterable[str] = ()
) -> list[str]:
    """``n`` distinct random peptides over ``RESIDUES`` (no I, so folding is the identity)."""
    seen = set(avoid)
    out: list[str] = []
    while len(out) < n:
        p = "".join(rng.choice(RESIDUES) for _ in range(rng.randint(*lengths)))
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def random_dna(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


# ------------------------------------------------------------------ .ssl files


@dataclass(frozen=True)
class Row:
    sequence: str
    score: float | str


def write_ssl(
    path: Path,
    rows: Sequence[Row],
    *,
    columns: Sequence[str] = SSL_COLUMNS,
    run: str = "WT_1",
) -> Path:
    """A Cascadia 0.0.7 ``.ssl`` (tab-separated, ``score-type`` ``UNKNOWN``) with ``rows``."""
    lines = ["\t".join(columns)]
    for i, r in enumerate(rows):
        score = r.score if isinstance(r.score, str) else f"{r.score:.6f}"
        vals = {"file": f"{run}.mzML", "scan": str(i), "charge": str(2 + i % 3),
                "sequence": r.sequence, "score-type": "UNKNOWN", "score": score,
                "retention-time": f"{1.0 + i * 0.01:.4f}", "start-time": "0.0",
                "end-time": "0.0"}  # fmt: skip
        lines.append("\t".join(vals.get(c, "") for c in columns))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def ladder(n: int, top: float = 0.99, step: float = 2e-5) -> list[float]:
    """``n`` distinct descending scores, all strictly above 0.8 for ``n`` up to 9,000."""
    out = [round(top - i * step, 6) for i in range(n)]
    assert out[-1] > 0.8, "ladder too long for the cut-off"
    return out


# ------------------------------------------------------------------ synthetic genomes


@dataclass
class Contig:
    id: str
    seq: str
    topology: str = "circular"
    kind: str = "chromosome"


@dataclass
class Genome:
    contigs: list[Contig]

    def write(self, d: Path) -> tuple[Path, Path]:
        """``genome.fna.gz`` (mtime 0) and a minimal Bakta JSON with the replicon records."""
        d.mkdir(parents=True, exist_ok=True)
        fasta = "".join(f">{c.id} synthetic\n{c.seq}\n" for c in self.contigs)
        fna = d / "genome.fna.gz"
        fna.write_bytes(gzip.compress(fasta.encode(), mtime=0))
        bakta = {
            "version": {"bakta": "1.11.4", "db": {"version": "6.0", "type": "full"}},
            "sequences": [
                {"id": c.id, "type": c.kind, "topology": c.topology, "length": len(c.seq)}
                for c in self.contigs
            ],
            "features": [],
        }
        js = d / "bakta.json"
        js.write_text(json.dumps(bakta))
        return fna, js

    def as_inputs(self) -> tuple[dict[str, str], dict[str, Any]]:
        """``genome`` and ``replicons`` mappings for a direct ``build_sixframe`` call."""
        from gippa.model.coords import Replicon

        genome = {c.id: c.seq for c in self.contigs}
        reps = {c.id: Replicon(c.id, len(c.seq), c.topology, c.kind) for c in self.contigs}  # type: ignore[arg-type]
        return genome, reps


def planted_contig(
    cid: str,
    peptides: Iterable[str],
    rng: random.Random,
    *,
    topology: str = "circular",
    minus: Iterable[str] = (),
    pad: int = 60,
) -> Contig:
    """Random DNA with each of ``peptides`` coded on the plus strand (and ``minus`` on the minus
    strand), separated by random spacers."""
    parts = [random_dna(pad, rng)]
    for p in peptides:
        parts += [coding(p), random_dna(rng.randint(5, 13), rng)]
    for p in minus:
        parts += [revcomp(coding(p)), random_dna(rng.randint(5, 13), rng)]
    parts.append(random_dna(pad, rng))
    return Contig(cid, "".join(parts), topology)


def write_contaminants(path: Path, proteins: Mapping[str, str] | None = None) -> Path:
    """A small protein FASTA (keratin-like stand-ins by default)."""
    proteins = proteins or {
        "CONT_KRT1": "MSRQFSSRSGYRSGGGFSSGSAGLLNYQRRTTSSSTRRSGGGGGRFS",
        "CONT_TRYP": "IVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKS",
    }
    path.write_text("".join(f">{k}\n{v}\n" for k, v in proteins.items()))
    return path


# ------------------------------------------------------------------ driver

_MODULE: ModuleType | None = None


def module() -> ModuleType:
    """``validation/pilot/a3_feasibility.py`` imported once (it is not a package module)."""
    global _MODULE
    if _MODULE is None:
        require(A3_FEASIBILITY, "deliverable 9")
        spec = importlib.util.spec_from_file_location("a3_feasibility_08b", A3_FEASIBILITY)
        assert spec is not None and spec.loader is not None
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _MODULE = mod
    return _MODULE


def has_subcommand(name: str) -> None:
    """Fail clearly (not with an argparse usage error) while a subcommand is missing."""
    text = require(A3_FEASIBILITY, "deliverable 9").read_text()
    assert re.search(rf"""["']{re.escape(name)}["']""", text), (
        f"validation/pilot/a3_feasibility.py has no {name!r} subcommand yet (M10-08b)"
    )


@dataclass
class Result:
    rc: int
    out: Path
    stdout: str
    stderr: str
    j: dict[str, Any] | None = None
    raw: str | None = None


def call(argv: Sequence[str], out: Path) -> Result:
    buf_out, buf_err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(buf_out), contextlib.redirect_stderr(buf_err):
        try:
            rc = module().main(list(argv))
        except SystemExit as exc:  # argparse inside main
            rc = int(exc.code or 0) if not isinstance(exc.code, str) else 2
    raw = out.read_text() if out.is_file() else None
    return Result(rc=int(rc), out=out, stdout=buf_out.getvalue(), stderr=buf_err.getvalue(),
                  j=json.loads(raw) if raw else None, raw=raw)  # fmt: skip


@dataclass
class Case:
    """One ``evaluate`` input set in a directory."""

    d: Path
    ssl: Path
    genome: Path
    bakta: Path
    contaminants: Path
    peptides: list[str] = field(default_factory=list)

    def argv(self, out: Path, *extra: str, run: str = "WT_1") -> list[str]:
        return ["evaluate", str(self.ssl), "--run", run, "--genome", str(self.genome),
                "--bakta-json", str(self.bakta), "--contaminants", str(self.contaminants),
                "--out", str(out), *extra]  # fmt: skip

    def evaluate(self, *extra: str, out: Path | None = None) -> Result:
        has_subcommand("evaluate")
        o = out or (self.d / "a3_feasibility.json")
        if o.exists():
            o.unlink()
        return call(self.argv(o, *extra), o)


def make_case(
    d: Path,
    rows: Sequence[Row],
    genome: Genome,
    *,
    contaminants: Mapping[str, str] | None = None,
    columns: Sequence[str] = SSL_COLUMNS,
) -> Case:
    d.mkdir(parents=True, exist_ok=True)
    ssl = write_ssl(d / "WT_1.ssl", rows, columns=columns)
    fna, js = genome.write(d / "genome")
    cont = write_contaminants(d / "contaminants.fasta", contaminants)
    return Case(d=d, ssl=ssl, genome=fna, bakta=js, contaminants=cont)


def ranked_case(
    d: Path,
    n_pass: int,
    match_ranks: Iterable[int],
    *,
    seed: int,
    planted_below_top: int = 0,
    duplicates: int = 0,
) -> Case:
    """``n_pass`` distinct passing peptides with descending scores; the 0-based ranks in
    ``match_ranks`` are planted in a circular contig, plus ``planted_below_top`` peptides ranked
    after 1,000 (they must not count). ``duplicates`` extra rows repeat random peptides with a
    lower score (one peptide on several charges)."""
    rng = random.Random(seed)
    peps = random_peptides(n_pass, rng)
    scores = ladder(n_pass)
    ranks = set(match_ranks)
    below = [r for r in range(1000, n_pass)][:planted_below_top]
    plant = [peps[r] for r in sorted(ranks | set(below))]
    genome = Genome([planted_contig("chr1", plant, rng)])
    rows = [Row(p, s) for p, s in zip(peps, scores, strict=True)]
    for _ in range(duplicates):
        i = rng.randrange(n_pass)
        rows.append(Row(peps[i], round(max(0.8001, scores[i] - 0.01), 6)))
    rng.shuffle(rows)
    case = make_case(d, rows, genome)
    case.peptides = peps
    return case


# ------------------------------------------------------------------ JSON checks


def string_values(obj: Any) -> list[str]:
    if isinstance(obj, str):
        return [obj]
    if isinstance(obj, Mapping):
        return [s for v in obj.values() for s in string_values(v)]
    if isinstance(obj, list | tuple):
        return [s for v in obj for s in string_values(v)]
    return []


def sequence_like_strings(obj: Any) -> list[str]:
    return [s for s in string_values(obj) if SEQ_LIKE.match(s)]
