"""miniprot route: lead decisions after code review of PR #47 (items 1, 4, 5; issue M3-12).

docs/specs/M3-12-miniprot-route.md, "Lead decisions after code review":
1. Alternative starts: the ``cs`` tag of the ``##PAF`` line preceding an alignment is read; a lone
   ``*[acgt]tgM`` at query position 0 counts as a match for identity and sets
   ``ProteinAlignment.alt_start``.
4. ``assign_alignments`` raises ``ValueError`` on non-canonical intervals
   (``check_interval``: 0 <= start < L, unwrapped end <= start + L; linear contigs never wrap).
5. The parser raises ``ValueError`` when a ``Target`` ID contains ``;`` or ``=``.

The fixture ``tests/data/miniprot/alt_start_sim100k.gff`` is an excerpt of real miniprot 0.18
output (``-S -N 5 --outs=0.95``, image pinned in workflow/conf/containers.config) on the 100-kb
simulated genome: GSIM_00005 (ATG, exact), GSIM_00025 (GTG, minus), GSIM_00085 (TTG, minus),
GSIM_00090 (GTG, plus). miniprot reports the start codon as a substitution (``cs:Z:*gtgM:258``,
``Identity=0.9961`` = 258/259).

Choices made here (spec silent):
- The corrected identity is the miniprot identity with the start column counted as a match,
  i.e. (matches + 1) / aligned columns; asserted to 5e-4 (``Identity`` is printed to 4 decimals).
- ``alt_start`` is False when there is no ``##PAF`` line, when the first ``cs`` op is not a
  ``*..tgM`` substitution, when the substitution is not at query position 0 (PAF query start > 0
  or preceded by matches) or when the query residue is not ``M``.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from gippa.config.model import ProteinsConfig
from gippa.model import Interval, Replicon
from gippa.model.genes import Gene, GeneModel

DATA = Path(__file__).resolve().parents[1] / "data" / "miniprot" / "alt_start_sim100k.gff"
L = 100_000
REPS = {"chr1": Replicon("chr1", L, "circular", "chromosome")}
LIN = {"ctg1": Replicon("ctg1", 20_000, "linear", "unknown")}
CFG = ProteinsConfig()
PLEN = {"GSIM_00005": 647, "GSIM_00025": 259, "GSIM_00085": 188, "GSIM_00090": 154}


def _mp():  # type: ignore[no-untyped-def]
    from gippa.mapping import miniprot

    return miniprot


def _pa():  # type: ignore[no-untyped-def]
    from gippa.mapping import protein_assign

    return protein_assign


def _parse(path: Path, plen: dict[str, int] = PLEN, reps: dict[str, Replicon] = REPS):  # type: ignore[no-untyped-def]
    return _mp().parse_miniprot_gff(path, reps, plen)


def _by_id(alns) -> dict:  # type: ignore[no-untyped-def,type-arg]
    return {a.protein_id: a for a in alns}


def _block(pid: str) -> list[str]:
    """The ``##PAF`` line and the three GFF lines of ``pid`` in the fixture."""
    lines = DATA.read_text().splitlines()
    i = next(k for k, line in enumerate(lines) if line.startswith(f"##PAF\t{pid}\t"))
    return lines[i : i + 4]


def _write(tmp_path: Path, lines: list[str], name: str = "mp.gff") -> Path:
    p = tmp_path / name
    p.write_text("##gff-version 3\n" + "".join(line + "\n" for line in lines))
    return p


def _with_cs(lines: list[str], cs: str, identity: str | None = None) -> list[str]:
    paf = lines[0].split("\t")
    paf[-1] = f"cs:Z:{cs}"
    out = ["\t".join(paf), *lines[1:]]
    if identity is not None:
        out = [
            line.replace("Identity=0.9961", f"Identity={identity}") if i > 0 else line
            for i, line in enumerate(out)
        ]
    return out


def _cds(tag: str, aln_iv: Interval) -> Gene:
    n_aa = aln_iv.length_bp // 3
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=aln_iv,
        aa="M" + "A" * max(0, n_aa - 2),
        nt=None,
        start_codon="GTG",
        partial=False,
    )


# ---------------------------------------------------------------- item 1: alternative starts


def test_fixture_is_real_miniprot_output() -> None:
    text = DATA.read_text()
    assert text.count("##PAF\t") == 4
    assert "cs:Z:*gtgM:258" in text and "cs:Z:*ttgM:187" in text


def test_alt_start_field_exists_and_defaults_false() -> None:
    fields = {f.name: f for f in dataclasses.fields(_mp().ProteinAlignment)}
    assert "alt_start" in fields
    assert fields["alt_start"].default is False


def test_atg_start_exact_alignment_is_not_alt_start() -> None:
    a = _by_id(_parse(DATA))["GSIM_00005"]
    assert a.alt_start is False
    assert a.identity == pytest.approx(1.0)


@pytest.mark.parametrize(
    ("pid", "strand"),
    [("GSIM_00025", "-"), ("GSIM_00090", "+")],
    ids=["gtg-minus", "gtg-plus"],
)
def test_gtg_start_counts_as_match(pid: str, strand: str) -> None:
    a = _by_id(_parse(DATA))[pid]
    assert a.interval.strand == strand
    assert a.alt_start is True
    assert a.identity == pytest.approx(1.0, abs=5e-4)
    assert (a.q_start, a.q_end) == (0, PLEN[pid])


def test_ttg_start_counts_as_match() -> None:
    a = _by_id(_parse(DATA))["GSIM_00085"]
    assert a.alt_start is True
    assert a.identity == pytest.approx(1.0, abs=5e-4)


def test_all_four_fixture_alignments_parsed_in_file_order() -> None:
    alns = _parse(DATA)
    assert [a.protein_id for a in alns] == ["GSIM_00005", "GSIM_00025", "GSIM_00085", "GSIM_00090"]
    assert [a.alt_start for a in alns] == [False, True, True, True]


@pytest.mark.parametrize("pid", ["GSIM_00025", "GSIM_00085", "GSIM_00090"])
def test_alt_start_otherwise_exact_is_identical(pid: str) -> None:
    a = _by_id(_parse(DATA))[pid]
    gm = GeneModel(REPS, (_cds("LT_1", a.interval),), ())
    df, leftovers = _pa().assign_alignments([a], gm, CFG)
    assert leftovers == set()
    assert len(df) == 1
    r = df.iloc[0]
    assert r["gene"] == "LT_1"
    assert r["class"] == "identical"
    assert r["identity"] == pytest.approx(1.0, abs=5e-4)


def test_alt_start_with_mismatch_elsewhere_start_still_counts(tmp_path: Path) -> None:
    """``*gtgM:100*gcaP:157`` (259 aa): start counted, one real mismatch -> 258/259, near."""
    lines = _with_cs(_block("GSIM_00025"), "*gtgM:100*gcaP:157", identity="0.9923")
    a = _by_id(_parse(_write(tmp_path, lines)))["GSIM_00025"]
    assert a.alt_start is True
    assert a.identity == pytest.approx(258 / 259, abs=5e-4)
    gm = GeneModel(REPS, (_cds("LT_1", a.interval),), ())
    df, _ = _pa().assign_alignments([a], gm, CFG)
    assert df.iloc[0]["class"] == "near"


def test_ttg_start_with_mismatch_elsewhere(tmp_path: Path) -> None:
    lines = _with_cs(_block("GSIM_00085"), "*ttgM:50*gcaP:136")
    lines = [line.replace("Identity=0.9947", "Identity=0.9894") for line in lines]
    a = _by_id(_parse(_write(tmp_path, lines)))["GSIM_00085"]
    assert a.alt_start is True
    assert a.identity == pytest.approx(187 / 188, abs=5e-4)


@pytest.mark.parametrize(
    "cs",
    [":10*gtgM:248", "*gtgV:258", "*gcgA:258"],
    ids=["not-at-query-0", "query-not-M", "not-a-start-codon"],
)
def test_substitution_not_an_alt_start(tmp_path: Path, cs: str) -> None:
    a = _by_id(_parse(_write(tmp_path, _with_cs(_block("GSIM_00025"), cs))))["GSIM_00025"]
    assert a.alt_start is False
    assert a.identity == pytest.approx(0.9961, abs=1e-6)


def test_paf_query_start_not_zero_is_not_alt_start(tmp_path: Path) -> None:
    """Partial alignment from residue 2: the leading ``*gtgM`` is at query position 1."""
    lines = _block("GSIM_00025")
    paf = lines[0].split("\t")
    paf[3] = "1"
    paf[-1] = "cs:Z:*gtgM:257"
    lines = [
        "\t".join(paf),
        *(line.replace("Target=GSIM_00025 1 259", "Target=GSIM_00025 2 259") for line in lines[1:]),
    ]
    a = _by_id(_parse(_write(tmp_path, lines)))["GSIM_00025"]
    assert a.q_start == 1
    assert a.alt_start is False
    assert a.identity == pytest.approx(0.9961, abs=1e-6)


def test_no_paf_line_means_no_alt_start(tmp_path: Path) -> None:
    a = _by_id(_parse(_write(tmp_path, _block("GSIM_00025")[1:])))["GSIM_00025"]
    assert a.alt_start is False
    assert a.identity == pytest.approx(0.9961, abs=1e-6)


def test_paf_cs_belongs_to_the_following_alignment_only(tmp_path: Path) -> None:
    """A PAF line applies to the next mRNA only: an ATG alignment after a GTG one stays exact."""
    lines = _block("GSIM_00090") + _block("GSIM_00005")
    alns = _by_id(_parse(_write(tmp_path, lines)))
    assert alns["GSIM_00090"].alt_start is True
    assert alns["GSIM_00005"].alt_start is False


# ---------------------------------------------------------------- item 4: canonical intervals


def _aln(iv: Interval, pid: str = "P1", plen: int = 99):  # type: ignore[no-untyped-def]
    return _mp().ProteinAlignment(
        protein_id=pid,
        protein_len_aa=plen,
        route="miniprot",
        interval=iv,
        rank=1,
        identity=1.0,
        positive=1.0,
        q_start=0,
        q_end=plen,
        frameshift=0,
        stop_codon=0,
    )


def _negative_start(contig: str = "chr1") -> Interval:
    iv = Interval(contig, 10, 310, "+")
    object.__setattr__(iv, "start", -10)  # bypass the constructor check on purpose
    return iv


@pytest.mark.parametrize(
    "iv",
    [
        Interval("chr1", L, L + 300, "+"),
        Interval("chr1", L + 50, L + 350, "+"),
        Interval("chr1", 99_000, 99_000 + L + 3, "+"),
        _negative_start(),
    ],
    ids=["start-eq-L", "start-gt-L", "longer-than-replicon", "negative-start"],
)
def test_assign_rejects_non_canonical_interval(iv: Interval) -> None:
    gm = GeneModel(REPS, (_cds("LT_1", Interval("chr1", 100, 400, "+")),), ())
    with pytest.raises(ValueError):
        _pa().assign_alignments([_aln(iv)], gm, CFG)


def test_assign_rejects_wrapped_interval_on_linear_contig() -> None:
    gm = GeneModel(LIN, (_cds("LT_1", Interval("ctg1", 100, 400, "+")),), ())
    with pytest.raises(ValueError):
        _pa().assign_alignments([_aln(Interval("ctg1", 19_900, 20_200, "+"))], gm, CFG)


def test_assign_rejects_even_when_alignment_fails_thresholds() -> None:
    """Validation is on input, not only on rows that pass identity and coverage."""
    gm = GeneModel(REPS, (), ())
    bad = dataclasses.replace(_aln(Interval("chr1", L + 5, L + 305, "+")), identity=0.1)
    with pytest.raises(ValueError):
        _pa().assign_alignments([bad], gm, CFG)


def test_assign_rejects_one_bad_among_good() -> None:
    gm = GeneModel(REPS, (_cds("LT_1", Interval("chr1", 100, 400, "+")),), ())
    good = _aln(Interval("chr1", 100, 400, "+"), pid="P1")
    bad = _aln(Interval("chr1", L, L + 300, "+"), pid="P2")
    with pytest.raises(ValueError):
        _pa().assign_alignments([good, bad], gm, CFG)


@pytest.mark.parametrize(
    "iv",
    [
        Interval("chr1", 99_900, 100_197, "+"),  # canonical origin-crossing (unwrapped)
        Interval("chr1", 0, L, "+"),  # whole replicon
        Interval("chr1", L - 1, 2 * L - 1, "+"),  # end == start + L
    ],
    ids=["origin-crossing", "whole-replicon", "end-eq-start-plus-L"],
)
def test_assign_accepts_canonical_intervals(iv: Interval) -> None:
    gm = GeneModel(REPS, (), ())
    plen = iv.length_bp // 3
    df, leftovers = _pa().assign_alignments([_aln(iv, plen=plen)], gm, CFG)
    assert len(df) == 1
    assert leftovers == {"P1"}


# ---------------------------------------------------------------- item 5: unsafe Target IDs


def _mrna_line(pid: str) -> str:
    attrs = f"ID=MP000001;Rank=1;Identity=1.0000;Positive=1.0000;Target={pid} 1 99"
    return "\t".join(["chr1", "miniprot", "mRNA", "101", "400", "500", "+", ".", attrs])


@pytest.mark.parametrize("pid", ["P=1", "a=b=c", "P;1", "P;Rank=2"])
def test_parser_rejects_target_id_with_separator(tmp_path: Path, pid: str) -> None:
    paf_cols = ["##PAF", pid, "99", "0", "99", "+", "chr1", str(L), "100", "397", "297", "297"]
    paf = "\t".join([*paf_cols, "0", "AS:i:500", "cg:Z:99M", "cs:Z::99"])
    path = _write(tmp_path, [paf, _mrna_line(pid)])
    with pytest.raises(ValueError):
        _parse(path, plen={pid: 99, pid.split(";")[0]: 99})


def test_parser_accepts_safe_ids_with_other_punctuation(tmp_path: Path) -> None:
    pid = "sp|P0A7G6|RECA_ECOLI.1-x_y"
    path = _write(tmp_path, [_mrna_line(pid)])
    (a,) = _parse(path, plen={pid: 99})
    assert a.protein_id == pid


# ---------------------------------------------------------------- item 3 (Python side)


def test_main_nf_has_no_hard_coded_miniprot_args_default() -> None:
    """``mapping.proteins.miniprot_args`` reaches the process via ``ext.args`` only."""
    main_nf = Path(__file__).resolve().parents[2] / "workflow/modules/local/miniprot/main.nf"
    code = "\n".join(
        line for line in main_nf.read_text().splitlines() if not line.lstrip().startswith("//")
    )
    assert "task.ext.args" in code
    assert "--outs=0.95" not in code, "default miniprot_args hard-coded in main.nf"
