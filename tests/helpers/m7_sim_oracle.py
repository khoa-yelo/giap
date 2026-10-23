"""Independent oracle for the M7-12 giapsim scenarios (docs/specs/M7-12-*.md).

Brute-force recounts of the truth denominators of rule 2 and coordinate checks of rule 1, written
from the spec and the M3-18/M3-19/M3-20 file formats only: nothing here imports giapsim, so a bug
in the simulator cannot hide in the oracle. Shared by tests/unit/test_sim_m7_scenarios.py,
tests/property/test_sim_m7_props.py and tests/calibration/test_sim_m7_calibration.py.

Conventions (pinned by the tests; see the unit-test docstring):

* every truth table is tab-separated with a header and 1-based closed coordinates;
* ``truth/novel_orfs.tsv`` keeps its seven M3-18 columns first (``orf_id, class, contig, start,
  end, strand, aa``; ``start..end`` includes the stop codon and, for ``nterm_extension``, runs
  from the true start) and adds ``stop`` (0-based lowest coordinate of the stop codon, as M7-06's
  ``Orf.stop``), ``true_start`` (1-based first base of the true start codon in transcription
  direction: ``start`` on ``+``, ``end`` on ``-``), ``n_true_peptides`` and ``rbs_planted``;
* ``n_true_peptides`` counts distinct I/L-folded peptides of ``truth/peptides.tsv`` with
  ``is_wrong = false`` that cover the object: for a variant site, peptides of the site's gene
  whose sample interval contains the SNV; for a novel ORF of class intergenic, cds_outframe or
  antisense, peptides placed in that ORF (``locus_tag = orf_id``); for ``nterm_extension``,
  peptides of the extended gene with at least one residue upstream of the annotated start.
"""

from __future__ import annotations

import csv
import json
import re
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

_BASES = "TCAG"
_AA11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE11: dict[str, str] = {
    a + b + c: _AA11[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
STOPS = frozenset(c for c, a in CODE11.items() if a == "*")
STARTS = frozenset({"ATG", "GTG", "TTG"})
NOVEL_COLUMNS_M318 = ["orf_id", "class", "contig", "start", "end", "strand", "aa"]
NOVEL_NEW_COLUMNS = ["stop", "true_start", "n_true_peptides", "rbs_planted"]
NOVEL_CLASSES_M7 = ("intergenic", "cds_outframe", "antisense", "nterm_extension")
ALT_COLUMNS = ["contig", "codon_start", "strand", "gene", "ref_aa", "alt_aa", "n_true_peptides"]
INDEL_COLUMNS = ["contig", "pos_reference", "ref", "alt", "gene", "aa_change"]
ISOBARIC_TYPES = frozenset({"block", "near_block", "modification", "swap"})
RBS = "AGGAGG"
_COMP = str.maketrans("ACGT", "TGCA")
_SNV_CHANGE = re.compile(r"^([A-Z*])(\d+)([A-Z*])$")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def il(s: str) -> str:
    return s.replace("I", "L")


def truthy(v: str | None) -> bool:
    return str(v).strip().lower() in ("true", "1")


def rows(path: Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def header(path: Path) -> list[str]:
    return Path(path).read_text().splitlines()[0].split("\t")


def read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    cur = ""
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif line:
            out[cur].append(line.strip().upper())
    return {k: "".join(v) for k, v in out.items()}


def bakta_features(path: Path) -> dict[str, dict[str, Any]]:
    """CDS/sORF features of a giapsim Bakta-like JSON by locus (``contig`` or ``sequence``)."""
    data = json.loads(Path(path).read_text())
    out = {}
    for f in data["features"]:
        if f["type"] in ("cds", "sorf"):
            f = dict(f)
            f["contig"] = f.get("contig", f.get("sequence"))
            out[str(f["locus"])] = f
    return out


def strand_slice(seq: str, s0: int, e0: int, strand: str) -> str:
    """``seq[s0:e0]`` (wrapping past the end of a circular contig) in ``strand`` orientation."""
    n = len(seq)
    nt = seq[s0:e0] if e0 <= n else seq[s0:] + seq[: e0 - n]
    return revcomp(nt) if strand == "-" else nt


def translate_orf(nt: str) -> str:
    """Start ... stop codons → protein with the initiator read as M and the stop dropped."""
    cod = [nt[i : i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]
    return "M" + "".join(CODE11[c] for c in cod[1:-1])


def covers(s1: int, e1: int, p0: int, n: int) -> bool:
    """1-based closed ``[s1, e1]`` (``e1`` may exceed ``n``) contains 0-based ``p0``."""
    return any(s1 - 1 <= q < e1 for q in (p0, p0 + n))


def correct_peptide_rows(outdir: Path) -> list[dict[str, str]]:
    return [
        r
        for r in rows(outdir / "truth" / "peptides.tsv")
        if r["origin"] != "decoy" and not truthy(r["is_wrong"])
    ]


# ---------------------------------------------------------------- variants


def nonsyn_snvs(outdir: Path) -> list[dict[str, Any]]:
    """Non-synonymous SNV rows of ``truth/strain_differences.tsv`` with the parsed change."""
    out = []
    for r in rows(outdir / "truth" / "strain_differences.tsv"):
        if r["type"] != "snv" or not r["gene"] or not r["aa_change"]:
            continue
        m = _SNV_CHANGE.match(r["aa_change"])
        assert m, r
        out.append({**r, "ref_aa": m[1], "codon_no": int(m[2]), "alt_aa": m[3]})
    return out


def inframe_indels(outdir: Path) -> list[tuple[str, str, str, str, str, str]]:
    """In-frame indels inside a gene: ``(contig, pos_reference, ref, alt, gene, aa_change)``."""
    out = []
    for r in rows(outdir / "truth" / "strain_differences.tsv"):
        if r["type"] in ("ins", "del") and r["gene"] and (len(r["ref"]) - len(r["alt"])) % 3 == 0:
            out.append((r["contig"], r["pos_reference"], r["ref"], r["alt"], r["gene"],
                        r["aa_change"]))  # fmt: skip
    return out


def ref_codon_start(gene: Mapping[str, Any], pos_reference: int, length: int) -> int:
    """1-based lowest coordinate of the reference codon of ``gene`` holding ``pos_reference``."""
    s0, e0 = int(gene["start"]) - 1, int(gene["stop"])
    if e0 <= s0:  # Bakta-style start > stop across the origin
        e0 += length
    p0 = pos_reference - 1
    if p0 < s0:
        p0 += length
    if gene["strand"] == "+":
        c0 = s0 + 3 * ((p0 - s0) // 3)
    else:
        k = (e0 - 1 - p0) // 3
        c0 = e0 - 3 * (k + 1)
    return c0 % length + 1


def recount_site_peptides(
    outdir: Path, gene: str, contig: str, pos_sample: int, length: int,
    correct: Iterable[Mapping[str, str]] | None = None,
) -> int:  # fmt: skip
    """Distinct I/L-folded correct peptides of ``gene`` whose sample interval holds the SNV."""
    pool = correct if correct is not None else correct_peptide_rows(outdir)
    peps = {
        il(r["peptide"])
        for r in pool
        if r["locus_tag"] == gene
        and r["contig"] == contig
        and covers(int(r["start"]), int(r["end"]), pos_sample - 1, length)
    }
    return len(peps)


# ---------------------------------------------------------------- novel ORFs


def receiving_genome(outdir: Path) -> dict[str, str]:
    """The genome GIAP is given (rule 1): ``reference.fasta`` in strain scenarios."""
    ref = outdir / "reference.fasta"
    return read_fasta(ref if ref.exists() else outdir / "genome.fasta")


def nterm_feature(outdir: Path, row: Mapping[str, str]) -> dict[str, str]:
    """The ``truth/features.tsv`` CDS whose true start is the row's ``true_start``."""
    hits = [
        f
        for f in rows(outdir / "truth" / "features.tsv")
        if f["type"] == "cds"
        and f["true_start"]
        and f["contig"] == row["contig"]
        and f["strand"] == row["strand"]
        and int(f["true_start"]) == int(row["true_start"])
    ]
    assert len(hits) == 1, (row, hits)
    return hits[0]


def recount_orf_peptides(
    outdir: Path, row: Mapping[str, str], correct: Iterable[Mapping[str, str]] | None = None
) -> int:
    pool = list(correct if correct is not None else correct_peptide_rows(outdir))
    if row["class"] != "nterm_extension":
        return len({il(r["peptide"]) for r in pool if r["locus_tag"] == row["orf_id"]})
    f = nterm_feature(outdir, row)
    a1, b1 = int(f["start"]), int(f["end"])  # annotated CDS, 1-based closed
    out = set()
    for r in pool:
        if r["locus_tag"] != f["feature_id"]:
            continue
        s1, e1 = int(r["start"]), int(r["end"])
        upstream = s1 < a1 if row["strand"] == "+" else e1 > b1
        if upstream:
            out.add(il(r["peptide"]))
    return len(out)


def has_rbs_motif(seq: str, strand: str, first0: int) -> bool:
    """``AGGAGG`` ending 5-10 nt before the start codon whose first transcribed base is the
    0-based ``first0`` (spacer = nt between the motif's last G and the start codon)."""
    n = len(seq)
    for sp in range(5, 11):
        if strand == "+":
            a = first0 - sp - len(RBS)
            got = "".join(seq[(a + i) % n] for i in range(len(RBS)))
        else:
            a = first0 + sp + 1
            got = revcomp("".join(seq[(a + i) % n] for i in range(len(RBS))))
        if got == RBS:
            return True
    return False


def tryptic(protein: str, lo: int = 7, hi: int = 30, missed: int = 1) -> set[str]:
    """Trypsin/P products (cut after every K/R) of length ``lo..hi`` with ≤ ``missed`` missed
    cleavages."""
    cuts = [0, *(i + 1 for i, c in enumerate(protein[:-1]) if c in "KR"), len(protein)]
    out = set()
    for i in range(len(cuts) - 1):
        for j in range(i + 1, min(i + 2 + missed, len(cuts))):
            if lo <= cuts[j] - cuts[i] <= hi:
                out.add(protein[cuts[i] : cuts[j]])
    return out
