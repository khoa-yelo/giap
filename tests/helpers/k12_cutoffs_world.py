"""Toy K-12-shaped input worlds for the M6-12 driver ``validation/calibrate_cutoffs.py``.

Every file the driver reads is written here under ``tmp_path`` (never RegulonDB bytes, never
``$GIAP_ROOT/recordings`` or ``$GIAP_DB_DIR``):

- ``regulondb/``: a release directory as ``validation/fetch_regulondb.py`` writes it (three
  RegulonDB-*shaped* tables + ``provenance.json``) and ``truth/`` as M6-11 writes it
  (``pair_truth.tsv``, ``terminator_truth.tsv``, ``truth_summary.json``);
- ``bakta.json`` (no ``nt``) and ``genome.fa`` (written with ``giap.io.fasta.write_fasta``, so its
  bytes are already the normalised genome);
- ``uniop/``: ``uniop_cds.gff`` (``giap.annotate.uniop_adapter.write_prodigal_gff``),
  ``uniop.pred`` and a deliberately **shifted** ``uniop_map.tsv`` that the driver must ignore;
- ``terminatornet/terms.tsv.gz`` (TerminatorNet ``-s`` rows, 1-based closed) and its
  ``PROVENANCE.md`` carrying ``tool commit `<sha>```.

Coordinates in the specs below are 0-based half-open, as everywhere inside GIAP.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import random
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from giap.annotate.operons import adjacent_pairs
from giap.annotate.uniop_adapter import uniop_gene_order, write_prodigal_gff
from giap.io.fasta import write_fasta
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

CONTIG = "U00096.3"
TN_COMMIT = "72180309545988dba63b8bd5db5d0b09bf4dba70"  # ADR-0008 pin
UNIOP_PIN = "3ed14dfc2a5d69b993d99a26e50c1f8e9723c087"  # ADR-0008 pin
DATASETS = ("TUSet", "TerminatorSet", "GeneProductSet")


@dataclass(frozen=True)
class G:
    """A toy gene: locus tag, 0-based half-open span, strand, contig."""

    tag: str
    start: int
    end: int
    strand: str
    contig: str = CONTIG


@dataclass(frozen=True)
class Call:
    """A TerminatorNet call (0-based half-open) with its ``pred`` confidence."""

    start: int
    end: int
    strand: str
    conf: float
    contig: str = CONTIG


@dataclass(frozen=True)
class Term:
    """A RegulonDB terminator (0-based half-open) with its id and ``type``."""

    tid: str
    start: int
    end: int
    strand: str
    kind: str = "rho-independent"
    contig: str = CONTIG


@dataclass
class World:
    root: Path
    regulondb_dir: Path
    bakta_json: Path
    genome: Path
    uniop_dir: Path
    terminatornet: Path
    genes: list[Gene]
    replicons: dict[str, Replicon]
    pairs: list[tuple[str, str, str]]  # pair_truth.tsv rows (gene_a, gene_b, label), file order
    probs: dict[tuple[str, str], float]
    calls: list[Call]
    truth: list[Term]
    pins: dict[str, str] = field(default_factory=dict)

    def run_args(self, out_dir: Path) -> list[Any]:
        return [
            self.regulondb_dir,
            self.bakta_json,
            self.genome,
            self.uniop_dir,
            self.terminatornet,
            out_dir,
        ]

    def cli_args(self, out_dir: Path) -> list[str]:
        return [
            "--regulondb-dir",
            str(self.regulondb_dir),
            "--bakta-json",
            str(self.bakta_json),
            "--genome",
            str(self.genome),
            "--uniop-dir",
            str(self.uniop_dir),
            "--terminatornet",
            str(self.terminatornet),
            "--out-dir",
            str(out_dir),
        ]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gene_obj(g: G) -> Gene:
    return Gene(
        locus_tag=g.tag,
        feature_type="cds",
        interval=Interval(g.contig, g.start, g.end, g.strand),  # type: ignore[arg-type]
        aa="M" * max(1, (g.end - g.start) // 3 - 1),
        nt=None,
        start_codon="ATG",
    )


def _bakta_json(path: Path, genes: Sequence[G], replicons: Mapping[str, Replicon]) -> None:
    features = [
        {
            "type": "cds",
            "locus": g.tag,
            "sequence": g.contig,
            "start": g.start + 1,
            "stop": g.end,
            "strand": g.strand,
            "aa": "M" * max(1, (g.end - g.start) // 3 - 1),
        }
        for g in genes
    ]
    seqs = [
        {"id": r.id, "length": r.length, "type": "chromosome", "topology": r.topology}
        for r in replicons.values()
    ]
    doc = {
        "version": {"bakta": "1.12.2", "db": {"version": "6.0", "type": "full"}},
        "sequences": seqs,
        "features": features,
    }
    path.write_text(json.dumps(doc), encoding="utf-8")


def _genome(path: Path, replicons: Mapping[str, Replicon], seed: int = 11) -> None:
    rnd = random.Random(seed)
    write_fasta(
        path,
        ((r.id, "".join(rnd.choice("ACGT") for _ in range(r.length))) for r in replicons.values()),
    )


def _release(regdir: Path, truth: Sequence[Term], version: str = "14.5.0") -> dict[str, str]:
    regdir.mkdir(parents=True, exist_ok=True)
    rows = [
        f"{t.tid}\t{t.start + 1}\t{t.end}\t{'forward' if t.strand == '+' else 'reverse'}\tacgt\t\t"
        f"{t.kind}"
        for t in truth
    ]
    tables = {
        "TUSet": "1)id\t2)name\t3)operonId\t4)operonName\t5)tuGenes\n",
        "TerminatorSet": "\n".join(
            ["1)id\t2)leftEndPos\t3)rightEndPos\t4)strand\t5)sequence\t6)relatedTus\t7)type", *rows]
        )
        + "\n",
        "GeneProductSet": (
            "1)geneId\t2)geneName\t3)bnumber\t4)leftEndPos\t5)rightEndPos\t6)strand\n"
        ),
    }
    files, shas = [], {}
    for ds in DATASETS:
        data = tables[ds].encode("utf-8")
        (regdir / f"{ds}.txt").write_bytes(data)
        shas[ds] = hashlib.sha256(data).hexdigest()
        files.append({"dataset": ds, "file": f"{ds}.txt", "bytes": len(data), "sha256": shas[ds]})
    meta = {
        "source": "RegulonDB",
        "regulondb_version": version,
        "release_date": "2026-01-28",
        "files": files,
    }
    (regdir / "provenance.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return shas


def _truth_dir(
    tdir: Path,
    pairs: Sequence[tuple[str, str, str]],
    truth: Sequence[Term],
    shas: Mapping[str, str],
    contig: str,
) -> None:
    tdir.mkdir(parents=True, exist_ok=True)
    lines = ["gene_a\tgene_b\tlabel\ttu_ids"]
    lines += [f"{a}\t{b}\t{lab}\t{'TU' + a if lab == 'positive' else ''}" for a, b, lab in pairs]
    (tdir / "pair_truth.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rows = sorted(truth, key=lambda t: (t.contig, t.start, t.end, t.strand, t.tid))
    tl = ["contig\tstart\tend\tstrand\tsource_id"]
    tl += [f"{t.contig}\t{t.start}\t{t.end}\t{t.strand}\t{t.tid}" for t in rows]
    (tdir / "terminator_truth.tsv").write_text("\n".join(tl) + "\n", encoding="utf-8")
    n = {lab: sum(1 for p in pairs if p[2] == lab) for lab in ("positive", "negative", "ambiguous")}
    summary = {
        "pairs_total": len(pairs),
        "pairs_positive": n["positive"],
        "pairs_negative": n["negative"],
        "pairs_ambiguous": n["ambiguous"],
        "pairs_labelled": n["positive"] + n["negative"],
        "terminators": len(truth),
        "contig": contig,
        "regulondb_version": "14.5.0",
        "release_date": "2026-01-28",
        "sha256": dict(shas),
    }
    (tdir / "truth_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def write_uniop(
    udir: Path,
    genes: Sequence[Gene],
    replicons: Mapping[str, Replicon],
    probs: Mapping[tuple[str, str], float | str],
) -> None:
    """GFF via the real adapter, a pred holding the pairs in ``probs`` and a shifted map."""
    udir.mkdir(parents=True, exist_ok=True)
    gff = udir / "uniop_cds.gff"
    id_map = write_prodigal_gff(gff, genes, list(replicons.values()))
    tag_of = dict(id_map)
    order = uniop_gene_order(gff)
    lines = ["operon prediction scores", "Gene A\tGene B\tPrediction"]
    for k, (pid, contig, _strand) in enumerate(order, start=1):
        nxt = next((o for o in order[k:] if o[1] == contig), None)
        if nxt is None:
            continue
        key = (tag_of[pid], tag_of[nxt[0]])
        if key in probs:
            lines.append(f"{k}\t{k + 1}\t{probs[key]} ")
    (udir / "uniop.pred").write_text("\n".join(lines) + "\n", encoding="utf-8")
    # the recorded map's defect: every ID labelled with its neighbour's tag (must be ignored)
    tags = [t for _, t in id_map]
    shifted = ["prodigal_id\tlocus_tag"]
    shifted += [f"{pid}\t{tags[(i + 1) % len(tags)]}" for i, (pid, _) in enumerate(id_map)]
    (udir / "uniop_map.tsv").write_text("\n".join(shifted) + "\n", encoding="utf-8")


def tn_line(c: Call, k: int) -> str:
    cols = [
        c.contig,
        f"TERM_{k}",
        "GGAAACAC",
        "AGAAAAAAGCCCGCAC",
        "TCGACCAA",
        "((((....))))",
        "-12.5",
        repr(c.conf),
        str(c.start + 1),
        str(c.end),
        c.strand,
    ]
    return "\t".join(cols)


def write_terminatornet(tdir: Path, calls: Sequence[Call], commit: str = TN_COMMIT) -> Path:
    tdir.mkdir(parents=True, exist_ok=True)
    lines = [tn_line(c, k) for k, c in enumerate(calls, start=1)]
    lines.append("plasmid_absent\tNO TERMINATORS FOUND")  # skipped by the parser
    path = tdir / "terms.tsv.gz"
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    (tdir / "PROVENANCE.md").write_text(
        "# terminatornet on toy\n\n"
        "- Env lock: `workflow/envs/locks/terminatornet.explicit.txt` sha256 `00`; "
        f"tool commit `{commit}`\n",
        encoding="utf-8",
    )
    return path


def build_world(
    root: Path,
    genes: Sequence[G],
    replicons: Iterable[Replicon],
    labels: Mapping[tuple[str, str], str],
    probs: Mapping[tuple[str, str], float],
    calls: Sequence[Call] = (),
    truth: Sequence[Term] = (),
    *,
    tn_commit: str = TN_COMMIT,
    default_label: str = "ambiguous",
) -> World:
    """Write every input; ``pair_truth.tsv`` holds every adjacent same-strand pair of the model on
    ``CONTIG`` (wrap pair included), labelled from ``labels`` (else ``default_label``)."""
    reps = {r.id: r for r in replicons}
    gobjs = [gene_obj(g) for g in genes]
    on_contig = [g for g in gobjs if g.interval.contig == CONTIG]
    pairs = [
        (a.locus_tag, b.locus_tag, labels.get((a.locus_tag, b.locus_tag), default_label))
        for a, b in adjacent_pairs(on_contig, {CONTIG: reps[CONTIG]})
    ]
    regdir = root / "regulondb"
    shas = _release(regdir, truth)
    _truth_dir(regdir / "truth", pairs, [t for t in truth if t.contig == CONTIG], shas, CONTIG)
    bakta = root / "bakta.json"
    _bakta_json(bakta, genes, reps)
    genome = root / "genome.fa"
    _genome(genome, reps)
    udir = root / "uniop"
    write_uniop(udir, gobjs, reps, probs)
    terms = write_terminatornet(root / "terminatornet", calls, tn_commit)
    pins = {
        "terminatornet_commit": TN_COMMIT,
        "genome_sha256": sha256(genome),
        "uniop_pred_sha256": sha256(udir / "uniop.pred"),
        "uniop_gff_sha256": sha256(udir / "uniop_cds.gff"),
        "terms_sha256": sha256(terms),
    }
    return World(
        root,
        regdir,
        bakta,
        genome,
        udir,
        terms,
        gobjs,
        reps,
        pairs,
        dict(probs),
        list(calls),
        list(truth),
        pins,
    )


def refresh_pins(w: World) -> None:
    """Recompute the UniOP / genome / terms pins after a test rewrote a recording file."""
    w.pins["genome_sha256"] = sha256(w.genome)
    w.pins["uniop_pred_sha256"] = sha256(w.uniop_dir / "uniop.pred")
    w.pins["uniop_gff_sha256"] = sha256(w.uniop_dir / "uniop_cds.gff")
    w.pins["terms_sha256"] = sha256(w.terminatornet)


# ----------------------------------------------------------------- spec examples 6 and 7

LIN10K = Replicon(CONTIG, 10_000, "linear")
CIRC10K = Replicon(CONTIG, 10_000, "circular")

EX6_GENES = [
    G("g1", 100, 1000, "+"),
    G("g2", 1500, 2500, "+"),
    G("g3", 3000, 4000, "-"),
    G("g4", 5000, 6000, "-"),
]
EX6_CALLS = [
    Call(1100, 1140, "+", 0.95),  # c1
    Call(2460, 2490, "+", 0.90),  # c2
    Call(2600, 2640, "-", 0.85),  # c3
    Call(3500, 3540, "+", 0.99),  # c4: in no + window
    Call(4100, 4130, "-", 0.60),  # c5
    Call(1460, 1520, "+", 0.99),  # c6: crosses g2's start
]
EX6_TRUTH = [
    Term("T1", 1120, 1150, "+", "rho-independent"),
    Term("T2", 2560, 2580, "+", "rho-independent"),
    Term("T3", 2680, 2700, "-", "rho-dependent"),
    Term("T4", 4100, 4130, "+", "rho-independent"),
]
EX6_LABELS = {("g1", "g2"): "positive", ("g3", "g4"): "negative"}
EX6_PROBS = {("g1", "g2"): 0.9, ("g3", "g4"): 0.2}


def example6(root: Path) -> World:
    return build_world(root, EX6_GENES, [LIN10K], EX6_LABELS, EX6_PROBS, EX6_CALLS, EX6_TRUTH)


# Example 7 on the circular 10 kb replicon, with four extra genes between gA and gB so the UniOP
# half has both classes; gB's window still wraps to gA. The wrap pair gB->gA is ambiguous (as on
# K-12), so every labelled pair has a UniOP p.
EX7_GENES = [
    G("gA", 200, 800, "+"),
    G("gC", 2000, 2600, "+"),
    G("gD", 2700, 3300, "+"),
    G("gE", 5000, 5600, "-"),
    G("gF", 5650, 6200, "-"),
    G("gB", 9500, 9900, "+"),
]
EX7_LABELS = {("gA", "gC"): "negative", ("gC", "gD"): "positive", ("gE", "gF"): "positive"}
EX7_PROBS = {("gA", "gC"): 0.3, ("gC", "gD"): 0.8, ("gE", "gF"): 0.7}
EX7_CALLS = [Call(20, 60, "+", 0.8)]
EX7_TRUTH = [Term("TW", 80, 100, "+")]


def example7(root: Path) -> World:
    return build_world(root, EX7_GENES, [CIRC10K], EX7_LABELS, EX7_PROBS, EX7_CALLS, EX7_TRUTH)


# ----------------------------------------------------------------- a K-12-sized-down world


@dataclass
class BigSpec:
    genes: list[G]
    labels: dict[tuple[str, str], str]
    probs: dict[tuple[str, str], float]
    calls: list[Call]
    truth: list[Term]
    directon_of: dict[str, int]
    length: int


def big_spec(seed: int = 7, n_directons: int = 81) -> BigSpec:
    """About 640 same-strand pairs in 81 directons (first and last on ``+``, so the wrap joins
    them) on a circular replicon. Labels cluster by TU inside directons; UniOP p is tied
    (2 decimals) and informative; distances overlap between classes so the refit converges.
    Calls sit in downstream windows (some matched within 50 bp, some not) plus some inside genes.
    """
    rnd = random.Random(seed)
    genes: list[G] = []
    labels: dict[tuple[str, str], str] = {}
    probs: dict[tuple[str, str], float] = {}
    directon_of: dict[str, int] = {}
    pos = 300
    k = 0
    for di in range(n_directons):
        strand = "+" if di % 2 == 0 else "-"
        size = rnd.randint(4, 14)
        run: list[G] = []
        for gi in range(size):
            if gi:
                prev = run[-1]
                tag_a = prev.tag
                same_tu = rnd.random() > 0.3
                d = rnd.randint(-8, 90) if same_tu else rnd.randint(20, 380)
                start = prev.end + d
            else:
                start = pos
            length = rnd.randint(200, 400) * 3
            k += 1
            g = G(f"KT_{k:05d}", start, start + length, strand)
            run.append(g)
            directon_of[g.tag] = di
            if gi:
                # in genome order the pair is (prev, g) whatever the strand
                key = (tag_a, g.tag)
                u = rnd.random()
                label = "ambiguous" if u < 0.07 else ("positive" if same_tu else "negative")
                labels[key] = label
                if label == "positive":
                    p = rnd.betavariate(4.0, 1.4)
                elif label == "negative":
                    p = rnd.betavariate(1.4, 3.0)
                else:
                    p = rnd.random()
                probs[key] = round(min(max(p, 0.01), 0.99), 2)
        genes.extend(run)
        pos = run[-1].end + rnd.randint(150, 420)
    length = pos + 500
    # wrap join: the last directon is on +, like the first
    first_d = directon_of[genes[0].tag]
    last_d = directon_of[genes[-1].tag]
    for tag, di in list(directon_of.items()):
        if di == last_d:
            directon_of[tag] = first_d
    labels[(genes[-1].tag, genes[0].tag)] = "ambiguous"

    calls: list[Call] = []
    truth: list[Term] = []
    nt = 0
    for i, g in enumerate(genes):
        if g.strand == "+":
            nxt = genes[i + 1] if i + 1 < len(genes) else None
            room = (nxt.start if nxt else length + genes[0].start) - g.end
            if room < 70 or rnd.random() < 0.35:
                continue
            c = Call(g.end + 8, g.end + 40, "+", round(rnd.random(), 4))
        else:
            prv = genes[i - 1] if i else None
            room = g.start - (prv.end if prv else 0)
            if room < 70 or rnd.random() < 0.35:
                continue
            c = Call(g.start - 40, g.start - 8, "-", round(rnd.random(), 4))
        calls.append(c)
        u = rnd.random()
        if u < 0.45 + 0.3 * c.conf - 0.15:  # matched truth, more often at high confidence
            nt += 1
            kind = "rho-independent" if rnd.random() < 0.8 else "rho-dependent"
            if c.strand == "+":
                truth.append(Term(f"TM{nt:04d}", c.end + 3, c.end + 25, "+", kind))
            else:
                truth.append(Term(f"TM{nt:04d}", c.start - 25, c.start - 3, "-", kind))
        elif u > 0.93:  # an opposite-strand terminator on top of the call: never a match
            nt += 1
            other = "-" if c.strand == "+" else "+"
            truth.append(Term(f"TM{nt:04d}", c.start, c.end, other, "rho-independent"))
    # calls inside genes: never eligible
    for g in genes[::9]:
        mid = (g.start + g.end) // 2
        calls.append(Call(mid, mid + 30, g.strand, round(rnd.random(), 4)))
    calls.sort(key=lambda c: (c.start, c.end, c.strand))
    return BigSpec(genes, labels, probs, calls, truth, directon_of, length)


def big_world(root: Path, seed: int = 7) -> tuple[World, BigSpec]:
    spec = big_spec(seed)
    rep = Replicon(CONTIG, spec.length, "circular")
    w = build_world(root, spec.genes, [rep], spec.labels, spec.probs, spec.calls, spec.truth)
    return w, spec
