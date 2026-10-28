"""Hand-made tables for the fixed MVP-04 gate 11 scorer (tests/unit/test_mvp04_gate11_scorer.py
pins the contract). Nothing here imports gippa or gippasim."""

from __future__ import annotations

import importlib.util
import sys
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "dev" / "mvp04_summary.py"
CONTIG = "U00096.3"
N_SITES = 40
GATE_TEXT = "n_high >= 30; precision >= 0.95; recall >= 0.90; flags; VCF"
FLAGS = "snv;strain_difference;database_derived"
VAR_COLS = ["variant_id", "contig", "codon_start", "strand", "gene", "aa_pos", "ref_codon",
            "ref_aa", "alt_aa", "alt_type", "call", "confidence", "n_ref", "n_alt", "nt_changes",
            "alt_codon_vcf", "flags"]  # fmt: skip
ALT_COLS = ["contig", "codon_start", "strand", "gene", "ref_aa", "alt_aa", "n_true_peptides",
            "frameshift_upstream"]  # fmt: skip
NOVEL_ALT_COLS = ["contig", "codon_start", "strand", "orf_id", "ref_aa", "alt_aa",
                  "n_true_peptides"]  # fmt: skip
NOVEL_COLS = ["orf_id", "class", "contig", "start", "end", "strand", "aa", "stop", "true_start",
              "n_true_peptides", "rbs_planted"]  # fmt: skip
PEP_COLS = ["peptide", "is_wrong", "origin", "locus_tag", "contig", "start", "end", "strand"]
VCF = (
    "##fileformat=VCFv4.3\n"
    f"##contig=<ID={CONTIG},length=4641652>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    f"{CONTIG}\t1001\tVAR_1\tC\tT\t.\tPASS\t.\n"
)
_AA = "ACDEFGHMNPQSTVWY"  # no I/L (folding stays explicit), no K/R


def mod() -> ModuleType:
    spec = importlib.util.spec_from_file_location("mvp04_summary_gate11", SCRIPT)
    assert spec is not None and spec.loader is not None, SCRIPT
    m = importlib.util.module_from_spec(spec)
    sys.modules["mvp04_summary_gate11"] = m
    spec.loader.exec_module(m)
    return m


def pep(i: int, n: int = 12) -> str:
    """A distinct peptide of length ``n`` (``i`` written in base 16 over ``_AA``), ending in K."""
    body = []
    for _ in range(n - 1):
        body.append(_AA[i % 16])
        i //= 16
    return "".join(body) + "K"


def codon(k: int) -> int:
    return 1000 + 600 * k


def gene(k: int) -> str:
    return f"K12_{k:05d}"


def tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)] + ["\t".join(str(x) for x in r) for r in rows]
    path.write_text("\n".join(lines) + "\n")


@dataclass
class Case:
    """Editable tables; ``write`` puts them on disk and ``gate`` scores them."""

    tmp: Path
    variants: list[list[Any]]
    alt: list[list[Any]]
    novel_alt: list[list[Any]] | None
    novel: list[list[Any]]
    peptides: list[list[Any]]
    min_len: int = 10
    _n: int = 0

    @property
    def out(self) -> Path:
        return self.tmp / "run"

    @property
    def inp(self) -> Path:
        return self.tmp / "inputs"

    def next_pep(self, n: int = 12) -> str:
        self._n += 1
        return pep(10_000 + self._n, n)

    def true_peps(self, locus: str, cs: int, strand: str = "+", n: int = 2, length: int = 12,
                  **kw: Any) -> None:  # fmt: skip
        """``n`` true peptides of ``length`` on ``locus`` whose interval holds codon ``cs``."""
        for _ in range(n):
            s = cs - 3 * (length // 2)
            row = {"peptide": self.next_pep(length), "is_wrong": "false", "origin": "gene",
                   "locus_tag": locus, "contig": CONTIG, "start": s,
                   "end": s + 3 * length - 1, "strand": strand} | kw  # fmt: skip
            self.peptides.append([row[c] for c in PEP_COLS])

    def var(self, cs: int, alt: str, *, conf: str = "high", gene: str = "NA", strand: str = "+",
            ref_codon: str = "GCC", ref: str = "A", flags: str = FLAGS) -> None:  # fmt: skip
        i = len(self.variants) + 1
        self.variants.append([f"VAR_{i}", CONTIG, cs, strand, gene, "NA" if gene == "NA" else 5,
                              ref_codon, ref, alt, "substitution", "ALT", conf, 0, 2, 1, "NNN",
                              flags])  # fmt: skip

    def write(self) -> None:
        tsv(self.out / "tables" / "variants.tsv", VAR_COLS, self.variants)
        (self.out / "tables" / "variants.vcf").write_text(VCF)
        (self.out / "resolved_config.yaml").write_text(
            "mapping:\n  exact:\n    min_len: 7\n    max_len: 40\n"
            f"  tolerant:\n    enabled: true\n    min_len: {self.min_len}\n"
        )
        truth = self.inp / "truth"
        tsv(truth / "alt_sites.tsv", ALT_COLS, self.alt)
        nas = truth / "novel_alt_sites.tsv"
        if self.novel_alt is None:
            nas.unlink(missing_ok=True)
        else:
            tsv(nas, NOVEL_ALT_COLS, self.novel_alt)
        tsv(truth / "novel_orfs.tsv", NOVEL_COLS, self.novel)
        tsv(truth / "peptides.tsv", PEP_COLS, self.peptides)

    def gate(self) -> dict[str, Any]:
        self.write()
        res: dict[str, Any] = mod().gate11(self.out, self.inp)
        return res


def novel_orf(orf_id: str, start: int, end: int, strand: str, cls: str = "intergenic") -> list[Any]:
    """A truth novel ORF row (1-based closed ``start..end``; ``true_start`` = first base of the
    start codon in transcription direction, ``stop`` = lowest base of the stop codon)."""
    ts, stop = (start, end - 2) if strand == "+" else (end, start)
    return [orf_id, cls, CONTIG, start, end, strand, "M" + "A" * ((end - start + 1) // 3 - 2),
            stop, ts, 2, "false"]  # fmt: skip


def base_case(tmp: Path) -> Case:
    """40 gene sites ``A -> V`` on ``+``, each with 2 true 12-mers covering it and one matching
    ``high`` row, plus two truth novel ORFs (one per strand) far from them: precision 1, recall
    1, no initiator rows."""
    c = Case(tmp, [], [], [], [], [])
    for k in range(N_SITES):
        c.alt.append([CONTIG, codon(k), "+", gene(k), "A", "V", 2, "false"])
        c.true_peps(gene(k), codon(k))
        c.var(codon(k), "V", gene=gene(k))
    c.novel.append(novel_orf("NOVEL_0001", 100_001, 100_300, "+"))
    c.novel.append(novel_orf("NOVEL_0002", 200_001, 200_300, "-"))
    return c
