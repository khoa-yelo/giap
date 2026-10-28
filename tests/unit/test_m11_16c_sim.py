"""M11-16c: gippasim at genetic codes 4 and 25 (docs/specs/M11-16-translation-tables-4-25.md
item 8, AC 14).

Interface assumed (spec item 8):

- ``gippasim.codons.code(table) -> dict[str, str]`` (64 codons, ``*`` for stops) and
  ``gippasim.codons.stops(table) -> frozenset[str]``; ``CODE_11`` / ``STOPS`` stay the table-11
  values; any table outside {4, 11, 25} raises ``ValueError``;
- ``gippasim.codons.sense_usage(table, tga_sense_share=None) -> dict[str, float]``: today's
  ``SENSE_USAGE`` at 11; at 4 TGA joins Trp with ``tga_sense_share`` of Trp's total (default
  0.6), at 25 TGA joins Gly with ``tga_sense_share`` of Gly's total (default 0.05); the total
  weight of each amino acid is unchanged (the other codons of the family share the rest in
  today's proportions);
- ``SimConfig.genetic_code: int = 11`` in {4, 11, 25}, else ``ValueError``;
  ``SimConfig.tga_sense_share`` (``None`` = the table's default) must lie in (0, 1) at 4/25;
- default stop weights at 4/25: TAA/TAG renormalised (0.61/0.09 -> 0.871/0.129); a user
  ``stop_codons`` containing TGA at 4/25 raises ``ValueError``;
- ``n_selenoproteins > 0`` with ``genetic_code != 11`` raises ``ValueError``; ``n_readthrough``
  works (terminal stops TAA/TAG, read as Q);
- the simulated Bakta JSON has ``genome.translation_table`` only when ``genetic_code != 11``.

The codon-usage defaults at 4/25 are placeholders (spec Q-L1), so the statistical tests check
the structure (TGA never a stop; TGA's share of its family equals ``tga_sense_share``; the amino
acid composition is that of table 11; stops come only from the table's stop set) and use the
configured share, whatever its default. Proportions are pooled over 20 seeds
(``calibration_base + i`` of tests/calibration/tolerances.yaml) and must lie within 3 binomial
SEs (AC 14). Codons within 7 of either CDS end are left out of the share counts: the intra-TU
overlap of up to 20 bp rewrites them outside the sense-usage draw.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from collections.abc import Iterator
from functools import cache
from pathlib import Path
from typing import Any

import pytest
import yaml

TOL = Path(__file__).resolve().parents[1] / "calibration" / "tolerances.yaml"
BASE = int(yaml.safe_load(TOL.read_text())["seeds"]["calibration_base"])
N_SEEDS = 20
EDGE = 7  # codons at each CDS end excluded from share counts (overlaps up to 20 bp)
NCBI = {
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    25: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}
FAMILY = {4: "W", 25: "G"}
DEFAULT_SHARE = {4: 0.6, 25: 0.05}
UNSUPPORTED = (0, 1, 2, 12, 26, -1)
_COMP = str.maketrans("ACGT", "TGCA")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _ncbi(table: int) -> dict[str, str]:
    b = "TCAG"
    aa = NCBI[table]
    return {
        x + y + z: aa[16 * i + 4 * j + k]
        for i, x in enumerate(b)
        for j, y in enumerate(b)
        for k, z in enumerate(b)
    }


def _codons(nt: str) -> list[str]:
    return [nt[i : i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]


def _sim(tmp: Path, seed: int, **kw: Any) -> Any:
    from gippasim.simulate import SimConfig, simulate

    return simulate(SimConfig(**kw), seed, tmp)


def _json(path: Path) -> dict[str, Any]:
    data: dict[str, Any] = json.loads(path.read_text())
    return data


def _cds(data: dict[str, Any]) -> Iterator[dict[str, Any]]:
    """Intact annotated CDSs (no pseudogenes) of a simulated Bakta JSON."""
    for f in data["features"]:
        if f["type"] == "cds" and not f.get("pseudo"):
            yield f


def _tree_digest(root: Path) -> str:
    h = hashlib.sha256()
    for p in sorted(q for q in root.rglob("*") if q.is_file()):
        h.update(p.relative_to(root).as_posix().encode())
        h.update(p.read_bytes().replace(str(root).encode(), b"<OUT>"))
    return h.hexdigest()


def _within_3se(k: int, n: int, p: float) -> None:
    se = math.sqrt(p * (1 - p) / n)
    assert abs(k / n - p) <= 3 * se, f"{k}/{n} = {k / n:.4f}, expected {p:.4f} +- 3 SE ({se:.4f})"


# --------------------------------------------------------------------- codons module


@pytest.mark.parametrize("table", [4, 11, 25])
def test_code_matches_ncbi(table: int) -> None:
    from gippasim.codons import code

    assert code(table) == _ncbi(table)


def test_code_11_and_stops_unchanged() -> None:
    from gippasim.codons import CODE_11, STOPS, code, stops

    assert code(11) == CODE_11
    assert stops(11) == STOPS == frozenset({"TAA", "TAG", "TGA"})


@pytest.mark.parametrize("table", [4, 25])
def test_tga_is_sense_at_4_and_25(table: int) -> None:
    from gippasim.codons import CODE_11, code, stops

    t = code(table)
    assert t["TGA"] == FAMILY[table]
    assert stops(table) == frozenset({"TAA", "TAG"})
    assert isinstance(stops(table), frozenset)
    assert {c for c in t if t[c] != CODE_11[c]} == {"TGA"}  # the only differing codon


@pytest.mark.parametrize("table", UNSUPPORTED)
def test_codon_functions_refuse_other_tables(table: int) -> None:
    from gippasim import codons

    for fn in (codons.code, codons.stops, codons.sense_usage):
        with pytest.raises(ValueError):
            fn(table)


def test_sense_usage_11_is_today() -> None:
    from gippasim.codons import SENSE_USAGE, sense_usage

    assert sense_usage(11) == SENSE_USAGE
    assert sense_usage(11, tga_sense_share=0.3) == SENSE_USAGE  # ignored at 11


def _family_totals(usage: dict[str, float], table: int) -> dict[str, float]:
    t = _ncbi(table)
    out: Counter[str] = Counter()
    for c, w in usage.items():
        out[t[c]] += w
    return dict(out)


@pytest.mark.parametrize("table", [4, 25])
@pytest.mark.parametrize("share", [None, 0.3])
def test_sense_usage_tga_joins_family(table: int, share: float | None) -> None:
    from gippasim.codons import SENSE_USAGE, sense_usage

    u = sense_usage(table) if share is None else sense_usage(table, tga_sense_share=share)
    s = DEFAULT_SHARE[table] if share is None else share
    aa = FAMILY[table]
    t11, tt = _ncbi(11), _ncbi(table)
    assert set(u) == set(SENSE_USAGE) | {"TGA"}
    assert not set(u) & {"TAA", "TAG"}
    fam11 = [c for c in SENSE_USAGE if t11[c] == aa]
    total = sum(SENSE_USAGE[c] for c in fam11)
    assert u["TGA"] == pytest.approx(s * total)
    for c in fam11:  # the rest of the family keeps today's proportions
        assert u[c] == pytest.approx((1 - s) * SENSE_USAGE[c])
    for c in SENSE_USAGE:
        if tt[c] != aa:
            assert u[c] == SENSE_USAGE[c]
    # protein composition is the same in all three tables
    want = _family_totals(SENSE_USAGE, 11)
    got = _family_totals(u, table)
    assert got.keys() == want.keys()
    for a in want:
        assert got[a] == pytest.approx(want[a])


# --------------------------------------------------------------------- SimConfig validation


@pytest.mark.parametrize("table", [4, 11, 25])
def test_simconfig_accepts_supported_tables(table: int) -> None:
    from gippasim.simulate import SimConfig

    assert SimConfig(genetic_code=table).genetic_code == table


@pytest.mark.parametrize("table", UNSUPPORTED)
def test_simconfig_refuses_other_tables(table: int) -> None:
    from gippasim.simulate import SimConfig

    with pytest.raises(ValueError):
        SimConfig(genetic_code=table)


@pytest.mark.parametrize("table", [4, 25])
@pytest.mark.parametrize("share", [0.0, 1.0, -0.1, 1.5])
def test_tga_sense_share_must_be_in_open_unit_interval(table: int, share: float) -> None:
    from gippasim.simulate import SimConfig

    with pytest.raises(ValueError):
        SimConfig(genetic_code=table, tga_sense_share=share)


@pytest.mark.parametrize("table", [4, 25])
def test_user_stop_codons_with_tga_refused(table: int) -> None:
    from gippasim.simulate import SimConfig

    with pytest.raises(ValueError):
        SimConfig(genetic_code=table, stop_codons=(("TAA", 0.5), ("TGA", 0.5)))
    cfg = SimConfig(genetic_code=table, stop_codons=(("TAA", 0.5), ("TAG", 0.5)))
    assert dict(cfg.stop_codons) == {"TAA": 0.5, "TAG": 0.5}


@pytest.mark.parametrize("table", [4, 25])
def test_selenoproteins_refused_without_code_11(table: int, tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        _sim(tmp_path, 1, genetic_code=table, n_selenoproteins=1)


# --------------------------------------------------------------------- outputs


@pytest.mark.parametrize("table", [4, 25])
def test_bakta_json_records_table(table: int, tmp_path: Path) -> None:
    res = _sim(tmp_path / "m0", 1, genetic_code=table)
    assert _json(res.bakta_json)["genome"]["translation_table"] == table
    v2 = _sim(tmp_path / "v2", 1, genetic_code=table, genome_length=60_000, write_precomputed=True)
    assert _json(v2.bakta_json)["genome"]["translation_table"] == table


def test_bakta_json_has_no_genome_block_at_11(tmp_path: Path) -> None:
    res = _sim(tmp_path, 1, genetic_code=11)
    assert "genome" not in _json(res.bakta_json)


@pytest.mark.parametrize("table", [4, 25])
def test_deterministic_through_the_seed(table: int, tmp_path: Path) -> None:
    from gippasim.proteomics import ProteomicsConfig

    kw = {"genetic_code": table, "n_readthrough": 2, "proteomics": ProteomicsConfig()}
    _sim(tmp_path / "a", 5, **kw)
    _sim(tmp_path / "b", 5, **kw)
    _sim(tmp_path / "c", 6, **kw)
    assert _tree_digest(tmp_path / "a") == _tree_digest(tmp_path / "b")
    assert _tree_digest(tmp_path / "a") != _tree_digest(tmp_path / "c")


def _check_cds(f: dict[str, Any], table: int) -> None:
    from gippasim.codons import stops

    cod = _codons(f["nt"])
    st = stops(table)
    assert len(f["nt"]) % 3 == 0
    assert cod[-1] in st, f"{f['locus']} ends in {cod[-1]}"
    assert not st.intersection(cod[1:-1]), f"{f['locus']} has an internal stop"
    t = _ncbi(table)
    assert f["aa"] == "M" + "".join(t[c] for c in cod[1:-1]), f"{f['locus']} aa != nt"


@pytest.mark.parametrize("table", [4, 25])
def test_cds_structure_and_translation(table: int, tmp_path: Path) -> None:
    """No CDS has an internal TAA/TAG or ends in TGA; ``aa`` reads TGA as W (4) or G (25)."""
    res = _sim(
        tmp_path,
        3,
        genetic_code=table,
        genome_length=80_000,
        n_plasmids=1,
        plasmid_length=(6_000, 6_000),
        linear_contig_length=5_000,
        origin_crossing_gene=True,
        n_sorfs=2,
        n_nterm_extensions=2,
    )
    feats = list(_cds(_json(res.bakta_json)))
    assert len(feats) > 50
    for f in feats:
        _check_cds(f, table)
    assert any("TGA" in _codons(f["nt"])[1:-1] for f in feats)  # TGA is used as sense


@pytest.mark.parametrize("table", [4, 25])
def test_strain_sample_and_reference_follow_table(table: int, tmp_path: Path) -> None:
    from gippasim.strain import StrainConfig

    res = _sim(tmp_path, 2, genetic_code=table, strain=StrainConfig(snv_rate_per_kb=3.0))
    assert res.reference_json is not None
    for path in (res.bakta_json, res.reference_json):
        feats = list(_cds(_json(path)))
        assert feats
        for f in feats:
            _check_cds(f, table)


@pytest.mark.parametrize("table", [4, 25])
def test_readthrough_at_4_and_25(table: int, tmp_path: Path) -> None:
    """Terminal stops are TAA/TAG read as Q; the new stop is the next in-frame TAA/TAG (TGA is
    read through as sense)."""
    from gippasim.codons import stops

    res = _sim(tmp_path, 4, genetic_code=table, n_readthrough=5)
    seqs: dict[str, str] = {}
    cid = None
    for line in res.genome_fasta.read_text().splitlines():
        if line.startswith(">"):
            cid = line[1:].split()[0]
            seqs[cid] = ""
        else:
            assert cid is not None
            seqs[cid] += line
    with (res.outdir / "truth" / "readthrough.tsv").open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 5
    st = stops(table)
    for r in rows:
        assert r["codon"] in st and r["residue"] == "Q"
        seq, a = seqs[r["contig"]], int(r["codon_start"])
        n, step = len(seq), 3 if r["strand"] == "+" else -3

        def at(p: int, seq: str = seq, n: int = n, r: dict[str, str] = r) -> str:
            c = "".join(seq[(p + i) % n] for i in range(3))
            return c if r["strand"] == "+" else _revcomp(c)

        assert at(a) == r["codon"]
        p = a + step
        for _ in range(61):
            if at(p) in st:
                break
            p += step
        assert p % n == int(r["ext_stop"]), r


# --------------------------------------------------------------------- calibration (AC 14)


@cache
def _pool(table: int, share: float | None) -> tuple[Counter[str], int]:
    """Interior codon counts (EDGE codons trimmed) and terminal-TAA count over N_SEEDS genomes;
    also checks the per-CDS structure of every genome."""
    import tempfile

    from gippasim.simulate import SimConfig, simulate

    counts: Counter[str] = Counter()
    n_taa = 0
    kw: dict[str, Any] = {"genetic_code": table}
    if share is not None:
        kw["tga_sense_share"] = share
    with tempfile.TemporaryDirectory() as d:
        for i in range(N_SEEDS):
            out = Path(d) / str(i)
            res = simulate(SimConfig(**kw), BASE + i, out)
            data = _json(res.bakta_json)
            free = {f["locus"] for f in _free_stops(data)}
            for f in _cds(data):
                _check_cds(f, table)
                cod = _codons(f["nt"])
                counts.update(cod[1 + EDGE : -1 - EDGE])
                if f["locus"] in free:
                    counts["__stop__"] += 1
                    n_taa += cod[-1] == "TAA"
    return counts, n_taa


def _free_stops(data: dict[str, Any]) -> Iterator[dict[str, Any]]:
    """Coding features whose stop codon overlaps no other coding feature: an intra-TU overlap
    rewrites the upstream stop from the downstream start, outside the stop-weight draw."""
    feats = [f for f in data["features"] if f["type"] in ("cds", "sorf")]
    iv = [(f["contig"], f["start"] - 1, f["stop"]) for f in feats]
    for j, f in enumerate(feats):
        c, a, b = iv[j]
        s0, s1 = (b - 3, b) if f["strand"] == "+" else (a, a + 3)
        if not any(k != j and x == c and y0 < s1 and s0 < y1 for k, (x, y0, y1) in enumerate(iv)):
            yield f


SHARE_CASES = [(4, None), (4, 0.3), (25, None), (25, 0.4)]


@pytest.mark.calibration
@pytest.mark.parametrize(("table", "share"), SHARE_CASES)
def test_tga_share_of_family(table: int, share: float | None) -> None:
    """TGA's share of its family (W at 4, G at 25) in CDSs is ``tga_sense_share``."""
    counts, _ = _pool(table, share)
    t = _ncbi(table)
    fam = [c for c in t if t[c] == FAMILY[table]]
    n = sum(counts[c] for c in fam)
    s = DEFAULT_SHARE[table] if share is None else share
    _within_3se(counts["TGA"], n, s)


@pytest.mark.calibration
@pytest.mark.parametrize(("table", "share"), SHARE_CASES)
def test_family_frequency_is_table_11_composition(table: int, share: float | None) -> None:
    """The amino-acid weight of the family is unchanged: W (or G) per interior codon is the
    table-11 per-mille weight over the total sense weight."""
    from gippasim.codons import SENSE_USAGE

    counts, _ = _pool(table, share)
    t, t11 = _ncbi(table), _ncbi(11)
    aa = FAMILY[table]
    n = sum(v for k, v in counts.items() if k != "__stop__")
    k = sum(counts[c] for c in t if t[c] == aa)
    p = sum(w for c, w in SENSE_USAGE.items() if t11[c] == aa) / sum(SENSE_USAGE.values())
    _within_3se(k, n, p)


@pytest.mark.calibration
@pytest.mark.parametrize("table", [4, 25])
def test_default_stop_weights_renormalised(table: int) -> None:
    """Default stops at 4/25: TAA 0.61 / (0.61 + 0.09) = 0.871, TAG the rest, never TGA
    (``_check_cds``). Counted over stops that overlap no other gene (``_free_stops``; at table 11
    today this gives TAA 0.607 over 822 stops against 0.61)."""
    counts, n_taa = _pool(table, None)
    _within_3se(n_taa, counts["__stop__"], 0.61 / 0.70)
