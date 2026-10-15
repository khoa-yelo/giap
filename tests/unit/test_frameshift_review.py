"""M3-08 review pins (docs/specs/M3-08-frameshift-splits.md, "Lead decisions after code review").

1. Tests that kill the mutants surviving the PR #60 review: the hashed seed prefilter (forced on a
   small genome through the keyword-only knob ``match_frameshift(..., prefilter_min_hits=N)``; the
   prefilter runs when a query has more than ``N`` q-gram hits, default 256; a ~1 Mb equivalence
   test lives in tests/unit/test_frameshift_review_slow.py), the > 16-placement serial path,
   reverse-decoy ``frame_id``/``aa_start``/``aa_end``, the precedence boundaries (``smax``,
   ``emin``, ``lo`` slack) and the 2↔1 operation of the non-frameshift overlap check, and the
   circular shifts at the origin.
2. ``mapping.tolerant.frameshift.min_side_aa <= 0`` is rejected when the config is parsed.
3. Public helper names (test-author choice, the implementer adopts them):
   ``giap.mapping.tolerant.Search``, ``giap.mapping.tolerant.budget``,
   ``giap.mapping.tolerant.parse_max_diff``, ``giap.mapping.tolerant_dp.rank_weight``; the two
   M3-08 re-implementations of M3-06 operations are importable as
   ``giap.mapping.frameshift.band_row`` and ``giap.mapping.frameshift.nonfs_overlap`` (current
   signatures) and are cross-checked against ``tolerant_dp.verify`` on random cases.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "property"))

from fs_genomes import background, codons_for, join_split, rc
from oracle_frameshift import frameshift_oracle
from oracle_tolerant import align_window

from giap.config.model import FrameshiftConfig, TolerantConfig
from giap.mapping.decoys import DecoyType, build_decoys, mirror_aa
from giap.mapping.diffops import block_tables
from giap.mapping.frameshift import match_frameshift
from giap.mapping.sixframe import (
    Frame,
    SixFrameIndex,
    build_sixframe,
    decode,
    encode_query,
    encode_target,
)
from giap.mapping.tolerant_dp import decode_key, exchange_lookup, verify
from giap.model import Interval, Replicon
from giap.rng import for_step

CFG = TolerantConfig()
FS_SRC = Path(__file__).resolve().parents[2] / "src" / "giap" / "mapping" / "frameshift.py"
KEY = ["contig", "start", "end", "strand", "frame", "n_diff", "diff_types", "blocks"]
BIG = 1 << 62
UNIT = 1 << 45  # key of one difference (B**3, B = 2**15; M3-06)
NEVER = 10**12  # prefilter_min_hits that never triggers the prefilter


def _index(genome: dict[str, str], topology: str = "linear") -> SixFrameIndex:
    reps = {c: Replicon(c, len(s), topology) for c, s in genome.items()}  # type: ignore[arg-type]
    return build_sixframe(genome, reps)


def _decoys(idx: SixFrameIndex, *kinds: DecoyType) -> dict[DecoyType, tuple[Frame, ...]]:
    return {k: build_decoys(idx, k, 50, for_step("sixframe_index", 0)) for k in kinds}


def _run(idx: SixFrameIndex, queries: list[str], **kw: int) -> pd.DataFrame:
    dec = _decoys(idx, "reverse")
    return match_frameshift(
        queries, idx, dec, dict.fromkeys(queries, "reverse"), set(queries), CFG, **kw
    )


def _target(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["search"] == "target"].reset_index(drop=True)


def _blocks(s: object) -> list[tuple[int, int]]:
    return [(int(a), int(b)) for a, b in (x.split("-") for x in str(s).split(";"))]


def _ex() -> np.ndarray:
    return exchange_lookup(block_tables())


# =============================================================================================
# item 2: min_side_aa <= 0 is a config error
# =============================================================================================


@pytest.mark.parametrize("value", [0, -1, -4])
def test_min_side_aa_not_positive_is_rejected_by_the_config_model(value: int) -> None:
    with pytest.raises(ValueError, match="min_side_aa"):
        FrameshiftConfig(min_side_aa=value)
    with pytest.raises(ValueError, match="min_side_aa"):
        TolerantConfig.model_validate({"frameshift": {"enabled": True, "min_side_aa": value}})


@pytest.mark.parametrize("value", [0, -2])
def test_min_side_aa_not_positive_is_rejected_by_resolve_config(value: int) -> None:
    from giap.config import resolve_config

    with pytest.raises(ValueError, match="min_side_aa"):
        resolve_config(overrides=[f"mapping.tolerant.frameshift.min_side_aa={value}"])


def test_min_side_aa_one_is_still_accepted() -> None:
    assert FrameshiftConfig(min_side_aa=1).min_side_aa == 1
    assert TolerantConfig().frameshift.min_side_aa == 4


# =============================================================================================
# item 3: public helper names
# =============================================================================================


def test_tolerant_helpers_are_public_with_docstrings() -> None:
    from giap.mapping import tolerant, tolerant_dp

    for mod, names in (
        (tolerant, ("Search", "budget", "parse_max_diff")),
        (tolerant_dp, ("rank_weight",)),
    ):
        for name in names:
            assert name in mod.__all__, (mod.__name__, name)
            obj = getattr(mod, name)
            doc = getattr(obj, "__doc__", None) or getattr(
                getattr(obj, "py_func", None), "__doc__", None
            )
            assert doc and doc.strip(), (mod.__name__, name)
            assert not hasattr(mod, "_" + name), f"private alias {mod.__name__}._{name} left behind"


def test_frameshift_imports_no_private_names_from_tolerant_modules() -> None:
    tree = ast.parse(FS_SRC.read_text(encoding="utf-8"))
    bad = [
        (node.module, a.name)
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
        and node.module in ("giap.mapping.tolerant", "giap.mapping.tolerant_dp")
        for a in node.names
        if a.name.startswith("_")
    ]
    assert bad == []


def test_frameshift_dp_helpers_are_importable_by_public_name() -> None:
    from giap.mapping import frameshift

    for name in ("band_row", "nonfs_overlap"):
        assert callable(getattr(frameshift, name)), name


def test_parse_max_diff_and_budget() -> None:
    from giap.mapping.tolerant import budget, parse_max_diff

    ranges = parse_max_diff({"12-40": 2, "10-11": 1})
    assert ranges == [(10, 11, 1), (12, 40, 2)]
    assert [budget(L, CFG, ranges) for L in (4, 9, 10, 11, 12, 40, 41)] == [
        None,
        None,
        1,
        1,
        2,
        2,
        None,
    ]
    for bad in ({"10": 1}, {"10-600": 1}, {"10-40": 3}, {"10-40": -1}):
        with pytest.raises(ValueError):
            parse_max_diff(bad)


def test_rank_weight_is_the_m306_key_increment() -> None:
    from giap.mapping.tolerant_dp import rank_weight

    b = 1 << 15
    assert [rank_weight(r, b) for r in range(4)] == [
        b**3,
        b**3 + b**2,
        b**3 + b**2 + b,
        b**3 + b**2 + b + 1,
    ]
    names = ["substitution", "swap", "block", "near_block"]
    for r in range(4):
        assert decode_key(rank_weight(r, b)) == (1, (names[r],))


def test_search_run_returns_frame_local_windows() -> None:
    from giap.mapping.tolerant import Search

    idx = _index({"fs1": "GGCCGGAACCGACCATTGATAGAAAAACTGATGAACGCGTAA"})
    text = decode(idx.frames[2].aa)[:12]  # PEPTLDRKTDER (+3)
    out = Search(idx.frames).run([text], {text: 1}, _ex())
    assert [(key, fid, a, b) for key, fid, a, b in out[text]] == [(0, 2, 0, 12)]


# =============================================================================================
# item 3: band_row and nonfs_overlap against tolerant_dp.verify (random cases)
# =============================================================================================

ALPHA = "GNAQSTVLDWERKP"  # residues of the mass-equivalent and near-equal blocks + fillers
TEXT_ALPHA = ALPHA + "*X"


def _band_cells(q: np.ndarray, t: np.ndarray, d0: int, k: int) -> dict[int, int]:
    """Row-L keys (end column -> key) of ``band_row`` run with a free start inside the
    ``±2`` band around ``d0`` (the M3-06 band: every cell has ``|j - i - d0| <= 2``)."""
    from giap.mapping.frameshift import band_row

    L, n = len(q), len(t)
    nd = 5
    base = d0 - 2
    key = np.full((L + 1) * nd, BIG, dtype=np.int64)
    jn = np.full((L + 1) * nd, -1, dtype=np.int64)
    for dd in range(nd):
        if 0 <= base + dd <= n:
            key[dd] = 0
            jn[dd] = base + dd
    limit = (k + 1) * UNIT
    for i in range(1, L + 1):
        band_row(q, i, base, t, 0, n, _ex(), key, jn, nd, limit, 0)
    out = {}
    for dd in range(nd):
        j = base + L + dd
        if 0 <= j <= n and key[L * nd + dd] < limit:
            out[j] = int(key[L * nd + dd])
    return out


@settings(max_examples=400, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    q=st.text(alphabet=ALPHA, min_size=5, max_size=14),
    t=st.text(alphabet=TEXT_ALPHA, min_size=1, max_size=24),
    d0=st.integers(-3, 22),
    k=st.integers(0, 2),
)
def test_band_row_matches_verify(q: str, t: str, d0: int, k: int) -> None:
    qc, tc = encode_query(q), encode_target(t)
    cells = _band_cells(qc, tc, d0, k)
    got = verify(qc, tc, d0, k, block_tables())
    if not got:
        assert cells == {}
        return
    best = min(cells.values())
    assert decode_key(best) == (got[0].cost, got[0].diff_types)
    assert cells.get(got[0].t_end) == best


def test_band_row_matches_verify_on_planted_blocks() -> None:
    """Deterministic: each M3-06 operation (incl. 2↔1 and 1↔2 blocks) inside the band."""
    rng = np.random.default_rng(8)
    cases = [
        ("KDNWTEHFPY", "KDGGWTEHFPY"),
        ("KDGGWTEHFPY", "KDNWTEHFPY"),
        ("KDQWTEHFPY", "KDAGWTEHFPY"),
        ("KDAGWTEHFPY", "KDQWTEHFPY"),
        ("KDGTWTEHFPY", "KDASWTEHFPY"),
        ("KDWTEHFPYM", "KDWETHFPYM"),
        ("KDWTEHFPYM", "KDADTEHFPYM"),
    ]
    for q, w in cases:
        for pad in range(3):
            t = "".join(rng.choice(list(ALPHA), size=pad)) + w + "SS"
            for d0 in range(-1, pad + 3):
                qc, tc = encode_query(q), encode_target(t)
                got = verify(qc, tc, d0, 2, block_tables())
                cells = _band_cells(qc, tc, d0, 2)
                if got:
                    assert decode_key(min(cells.values())) == (got[0].cost, got[0].diff_types), (
                        q,
                        t,
                        d0,
                    )
                else:
                    assert cells == {}, (q, t, d0)


def _nonfs(q: str, t: str, o: int, a: int, b: int, c: int) -> bool:
    from giap.mapping.frameshift import nonfs_overlap

    tc = encode_target(t)
    return bool(nonfs_overlap(encode_query(q), tc, 0, len(tc), o, a, b, c, _ex()))


def _verify_windows(q: str, t: str, c: int) -> set[tuple[int, int, int]]:
    """(cost, start, end) of every verify() window over all diagonals."""
    qc, tc = encode_query(q), encode_target(t)
    out = set()
    for d0 in range(-2 - c, len(t) + 3):
        for al in verify(qc, tc, d0, c, block_tables()):
            out.add((al.cost, al.t_start, al.t_end))
    return out


@settings(max_examples=300, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    q=st.text(alphabet=ALPHA, min_size=5, max_size=12),
    t=st.text(alphabet=TEXT_ALPHA, min_size=3, max_size=30),
    c=st.integers(0, 2),
    o=st.integers(0, 2),
    data=st.data(),
)
def test_nonfs_overlap_agrees_with_verify(
    q: str, t: str, c: int, o: int, data: st.DataObject
) -> None:
    n = len(t)
    windows = _verify_windows(q, t, c)
    # whole frame: any window of cost <= c
    assert _nonfs(q, t, o, o, o + 3 * n, c) == bool(windows)
    a = data.draw(st.integers(o, o + 3 * n - 1))
    b = data.draw(st.integers(a + 1, o + 3 * n))
    got = _nonfs(q, t, o, a, b, c)
    # soundness from verify: a verify window overlapping [a, b) forces True
    if any(o + 3 * s < b and o + 3 * e > a for _, s, e in windows):
        assert got
    # completeness towards verify: True needs some window of cost <= c in the frame
    if got:
        assert windows


def _nonfs_brute(q: str, t: str, o: int, a: int, b: int, c: int) -> bool:
    """Exact: some window ``t[s:e]`` overlapping strand nt ``[a, b)`` aligns at cost <= c."""
    n = len(t)
    for s in range(n):
        for e in range(s + 1, n + 1):
            if not (o + 3 * s < b and o + 3 * e > a and abs((e - s) - len(q)) <= c):
                continue
            if align_window(q.replace("I", "L"), t.replace("I", "L"), s, e, None, c) is not None:
                return True
    return False


def test_nonfs_overlap_matches_exact_windows_on_random_intervals() -> None:
    rng = np.random.default_rng(306)
    for _ in range(250):
        L = int(rng.integers(5, 10))
        n = int(rng.integers(L - 2, L + 8))
        w = "".join(rng.choice(list(ALPHA), size=L))
        if rng.random() < 0.7:  # plant q (with blocks) in t so that hits exist
            t = "".join(rng.choice(list(TEXT_ALPHA), size=max(n - L, 0) // 2)) + w
            t = t + "".join(rng.choice(list(TEXT_ALPHA), size=max(n - len(t), 0)))
            q = w
            for src, dst in (("GG", "N"), ("N", "GG"), ("AG", "Q"), ("Q", "AG"), ("AS", "GT")):
                if src in q and rng.random() < 0.5:
                    q = q.replace(src, dst, 1)
        else:
            q = w
            t = "".join(rng.choice(list(TEXT_ALPHA), size=n))
        if len(q) < 5:
            continue
        c = int(rng.integers(0, 3))
        o = int(rng.integers(0, 3))
        a = int(rng.integers(o, o + 3 * len(t)))
        b = int(rng.integers(a + 1, o + 3 * len(t) + 1))
        assert _nonfs(q, t, o, a, b, c) == _nonfs_brute(q, t, o, a, b, c), (q, t, o, a, b, c)


# --- nonfs_overlap boundaries (mutants smax, emin, lo, ok21) ----------------------------------

Q = "KDNWTEHFPYMK"  # 12 aa, Q[0] == Q[-1]


@pytest.mark.parametrize("o", [0, 1, 2])
def test_nonfs_overlap_end_boundary_emin(o: int) -> None:
    """A window ending exactly at ``a`` does not overlap ``[a, b)``; one nt more does."""
    t = "SSSSS" + "KDNWTEHFPYMR" + "SSSSSSSSSSSSSSSS"  # Q with one substitution at aa [5, 17)
    end = o + 3 * 17
    assert not _nonfs(Q, t, o, end, end + 30, 1)
    assert _nonfs(Q, t, o, end - 1, end + 30, 1)
    assert _nonfs(Q, t, o, end - 3, end + 30, 1)


@pytest.mark.parametrize("o", [0, 1, 2])
def test_nonfs_overlap_start_boundary_smax(o: int) -> None:
    """A window starting exactly at ``b`` does not overlap ``[a, b)``; one nt earlier does."""
    t = "SSSSSSSSSSSSSSSS" + "KDNWTEHFPYMR" + "SSSSS"  # aa [16, 28)
    start = o + 3 * 16
    assert not _nonfs(Q, t, o, start - 30, start, 1)
    assert _nonfs(Q, t, o, start - 30, start + 1, 1)
    assert _nonfs(Q, t, o, start - 30, start + 3, 1)


@pytest.mark.parametrize("c", [1, 2])
def test_nonfs_overlap_long_window_reaching_one_codon_lo(c: int) -> None:
    """A 13-aa window (N -> GG block) whose last residue is the first codon of ``[a, b)``: its
    start lies L + 1 residues before ``a`` (the ``lo`` slack must cover target insertions)."""
    w = "KDGGWTEHFPYMK"  # 13 aa at [20, 33)
    t = "S" * 20 + w + "S" * 20
    a = 3 * 32  # codon of the window's last residue
    assert _nonfs(Q, t, 0, a, a + 36, c)
    assert not _nonfs(Q, t, 0, a + 3, a + 36, c)
    w2 = "KDGGWTEHFPYAGK"  # two insertions (N -> GG, Q -> AG): 14 aa, cost 2
    q2 = "KDNWTEHFPYQK"
    t2 = "S" * 20 + w2 + "S" * 20
    a2 = 3 * 33
    assert _nonfs(q2, t2, 0, a2, a2 + 36, 2)
    assert not _nonfs(q2, t2, 0, a2, a2 + 36, 1)


def test_nonfs_overlap_two_to_one_block() -> None:
    """Query GG against target N (2↔1, cost 1) is a non-frameshift alignment."""
    q = "KDGGWTEHFPYMK"
    t = "S" * 10 + "KDNWTEHFPYMK" + "S" * 10
    assert _nonfs(q, t, 0, 0, 3 * len(t), 1)
    assert _nonfs(q, t, 0, 3 * 21, 3 * 22, 1)  # last codon only
    assert not _nonfs(q, t, 0, 0, 3 * len(t), 0)
    q2 = "KDAGWTEHFPYMK"  # AG -> Q
    t2 = "S" * 10 + "KDQWTEHFPYMK" + "S" * 10
    assert _nonfs(q2, t2, 0, 3 * 21, 3 * 22, 1)


# =============================================================================================
# item 1: precedence at genome level, and circular shifts at the origin
# =============================================================================================


def _split_nt(q: str, p: int, delta: int) -> str:
    nt = join_split(codons_for(q[:p]), codons_for(q[p:]), delta)
    assert nt is not None
    return str(nt)


# configurations on the + strand: (name, nt, split locus start offset, suppressed?)
def _configs() -> list[tuple[str, str, int, bool]]:
    s1 = _split_nt(Q, 6, 1)  # 37 nt
    out = []
    # (emin) Q with one substitution ending exactly where the split starts: kept
    w = "".join(codons_for("KDNWTEHFPYMR"))
    out.append(("adjacent_before", w + s1, len(w), False))
    # (smax) the same window starting exactly where the split ends: kept
    out.append(("adjacent_after", s1 + w, 0, False))
    # (lo) N -> GG window (13 aa) sharing its last codon with the split's first codon: suppressed
    w13 = "".join(codons_for("KDGGWTEHFPYMK"))
    out.append(("block_before", w13[:-3] + s1, len(w13) - 3, True))
    # window sharing its first codon with the split's last codon (suffix frame): suppressed
    sm = _split_nt(Q, 5, -1)  # 35 nt: suffix frame = prefix frame - 1 nt; last codon at 32
    out.append(("block_after", sm + w13[3:], 0, True))
    # (ok21) query GG against target N, window sharing the split's first codon: suppressed
    q2 = "KDGGWTEHFPYMK"
    s2 = _split_nt(q2, 6, 1)
    w12 = "".join(codons_for("KDNWTEHFPYMK"))
    out.append(("two_to_one_before", w12[:-3] + s2, len(w12) - 3, True))
    return out


def _query_of(name: str) -> str:
    return "KDGGWTEHFPYMK" if name.startswith("two_to_one") else Q


CONFIGS = _configs()


@pytest.mark.parametrize(("name", "nt", "off", "suppressed"), CONFIGS, ids=[c[0] for c in CONFIGS])
def test_precedence_boundaries_on_a_linear_genome(
    name: str, nt: str, off: int, suppressed: bool
) -> None:
    pos = 300
    g = background(900, 11)
    g = g[:pos] + nt + g[pos + len(nt) :]
    idx = _index({"c": g})
    q = _query_of(name)
    tgt = _target(_run(idx, [q]))
    oracle = frameshift_oracle(q, idx, idx.frames, None, 2, 4)
    want = sorted((p.start, p.end, p.cost) for p in oracle.placements)
    assert sorted(zip(tgt["start"], tgt["end"], tgt["n_diff"], strict=True)) == want
    if suppressed:
        assert tgt.empty, tgt[KEY]
    else:
        assert (pos + off) in set(tgt["start"]), tgt[KEY]
        assert set(tgt["n_diff"]) == {1}


def _rotated_rows(g: str, cut: int, q: str) -> list[tuple[int, int, int, str]]:
    """Target rows of ``q`` on the circular genome rotated by ``cut``, in the coordinates of
    the unrotated genome (start mod length, end and blocks shifted alike)."""
    n = len(g)
    idx = _index({"c": g[cut:] + g[:cut]}, "circular")
    tgt = _target(_run(idx, [q]))
    out = []
    for r in tgt.itertuples():
        s = (int(r.start) + cut) % n
        sh = s - int(r.start)
        bl = ";".join(f"{a + sh}-{b + sh}" for a, b in _blocks(r.blocks))
        out.append((s, int(r.end) + sh, int(r.n_diff), bl))
    return sorted(out)


@pytest.mark.parametrize(("name", "nt", "off", "suppressed"), CONFIGS, ids=[c[0] for c in CONFIGS])
def test_precedence_is_invariant_under_circular_rotation(
    name: str, nt: str, off: int, suppressed: bool
) -> None:
    """The configuration is moved across the origin (and next to it, where only one copy of a
    window lies in the 3'-extension): the rows are those of the unrotated genome."""
    n = 420
    pos = 150
    g = background(n, 12)
    g = g[:pos] + nt + g[pos + len(nt) :]
    q = _query_of(name)
    ref = _rotated_rows(g, 0, q)
    if suppressed:
        assert ref == []
    else:
        assert [r[0] for r in ref].count(pos + off) == 1
    cuts = sorted(
        {
            *range(pos - 140, pos - 40, 3),
            *range(pos - 40, pos + len(nt) + 3),
            *range(pos + len(nt), pos + len(nt) + 140, 3),
        }
    )
    bad = [cut for cut in cuts if _rotated_rows(g, cut % n, q) != ref]
    assert bad == [], f"rows change under rotation by {bad[:10]}"


# =============================================================================================
# item 1: > 16 placements per query (serial recomputation, no truncation)
# =============================================================================================


def test_more_than_sixteen_split_loci_are_all_reported() -> None:
    q = "PEPTWDEKLMNAHR"  # 14 aa, k = 2
    nt = _split_nt(q, 7, 1)
    g = background(24_000, 21)
    pos = [500 + 1000 * j + (j % 3) for j in range(22)]
    for x in pos:
        g = g[:x] + nt + g[x + len(nt) :]
    idx = _index({"c": g})
    tgt = _target(_run(idx, [q]))
    starts = sorted(set(tgt.loc[tgt["n_diff"] == 1, "start"]))
    assert starts == pos


def test_more_than_sixteen_split_loci_minus_strand_and_decoys() -> None:
    q = "PEPTWDEKLMNAHR"
    nt = rc(_split_nt(q, 7, -1) or "")
    g = background(24_000, 22)
    pos = [700 + 1000 * j for j in range(20)]
    for x in pos:
        g = g[:x] + nt + g[x + len(nt) :]
    idx = _index({"c": g})
    rev = q[::-1]
    dec = _decoys(idx, "reverse")
    df = match_frameshift([q, rev], idx, dec, {q: "reverse", rev: "reverse"}, {q}, CFG)
    tgt = df[(df["search"] == "target") & (df["peptide"] == q)]
    assert sorted(tgt["start"]) == pos
    assert (tgt["strand"] == "-").all()
    d = df[(df["search"] == "decoy") & (df["peptide"] == rev)]
    assert sorted(d["start"]) == pos  # reverse-decoy splits of reversed(Q) = target splits of Q


# =============================================================================================
# item 1: reverse-decoy frame_id and aa_start/aa_end
# =============================================================================================


def _assert_prefix_window(
    row: pd.Series, frames: tuple[Frame, ...], idx: SixFrameIndex, kind: str | None
) -> None:
    """``frames[frame_id].aa[aa_start:aa_end]`` reads the query's prefix (exact cost-1 split) and
    maps (mirrored for decoys) onto one of the two genome blocks, labelled ``frame``."""
    fid, a, b = int(row["frame_id"]), int(row["aa_start"]), int(row["aa_end"])
    fr = frames[fid]
    seg = decode(fr.aa[a:b])
    L = len(str(row["query"]))
    assert 4 <= len(seg) <= L - 4, (seg, row["query"])
    assert str(row["query"]).startswith(seg), (seg, row["query"])
    ta, tb = (a, b) if kind is None else mirror_aa(kind, len(fr.aa), a, b)  # type: ignore[arg-type]
    iv = idx.to_interval(fid, ta, tb)
    assert (iv.contig, iv.strand) == (row["contig"], row["strand"])
    assert (iv.start, iv.end) in [tuple(x) for x in _blocks(row["blocks"])], (iv, row["blocks"])
    assert idx.frame_label(Interval(iv.contig, iv.start, iv.end, iv.strand)) == row["frame"]


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("delta", [1, -1])
@pytest.mark.parametrize("p", [4, 6, 9])
def test_reverse_decoy_rows_point_at_the_prefix_window(strand: str, delta: int, p: int) -> None:
    pep = "PEPTSDSKLSNAHR"  # S (any third base) ends every planted prefix: delta = -1 fits
    rev = pep[::-1]
    nt = _split_nt(rev, len(rev) - p, delta)  # reversed(pep) planted: pep's prefix has p aa
    nt = nt if strand == "+" else rc(nt)
    g = background(3000, 31)
    g = g[:1201] + nt + g[1201 + len(nt) :]
    idx = _index({"c": g})
    dec = _decoys(idx, "reverse")
    df = match_frameshift([pep, rev], idx, dec, {pep: "reverse", rev: "reverse"}, {rev}, CFG)
    drows = df[(df["search"] == "decoy") & (df["peptide"] == pep)]
    trows = df[(df["search"] == "target") & (df["peptide"] == rev)]
    assert len(drows) >= 1 and len(trows) >= 1
    assert sorted(drows["start"]) == sorted(trows["start"])
    for _, row in drows.iterrows():
        assert row["n_diff"] == 1
        _assert_prefix_window(row, dec["reverse"], idx, "reverse")
    for _, row in trows.iterrows():
        _assert_prefix_window(row, idx.frames, idx, None)


# =============================================================================================
# item 1: hashed seed prefilter (forced on a small genome)
# =============================================================================================


def _planted_genome(n: int, seed: int, n_plant: int) -> tuple[str, list[str]]:
    """Background with ``n_plant`` split peptides of 26..40 aa carrying one extra difference
    (k = 2: frameshift + substitution or block) at varied phases, both strands."""
    rng = np.random.default_rng(seed)
    g = background(n, seed)
    residues = "ACDEFGHKMNPQRSTVWY"
    qs: list[str] = []
    span = (n - 400) // n_plant
    for j in range(n_plant):
        L = int(rng.integers(26, 41))
        p = int(rng.integers(4, L - 3))
        s = "".join(rng.choice(list(residues), size=L))
        gp, gs = s[:p], s[p:]
        q = s
        i = int(rng.integers(L))
        if j % 3 == 1:
            q = q[:i] + ("W" if q[i] != "W" else "Y") + q[i + 1 :]
        elif j % 3 == 2 and "N" not in q:
            q = q[:i] + "N" + q[i + 1 :]
            if i < p:
                gp = gp[:i] + "GG" + gp[i + 1 :]
            else:
                gs = gs[: i - p] + "GG" + gs[i - p + 1 :]
        nt = join_split(codons_for(gp, rng), codons_for(gs, rng), (1, -1)[j % 2])
        if nt is None:
            continue
        nt = nt if (j // 2) % 2 == 0 else rc(nt)
        x = 200 + span * j + j % 8
        g = g[:x] + nt + g[x + len(nt) :]
        qs.append(q)
    return g, qs


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_forced_prefilter_equals_unfiltered_seeds(seed: int) -> None:
    g, qs = _planted_genome(30_000, seed, 24)
    idx = _index({"c": g})
    on = _run(idx, qs, prefilter_min_hits=0)
    off = _run(idx, qs, prefilter_min_hits=NEVER)
    pd.testing.assert_frame_equal(on, off)
    found = set(_target(on)["peptide"])
    assert found == set(qs)


def test_default_prefilter_threshold_is_256_and_keyword_only() -> None:
    import inspect

    sig = inspect.signature(match_frameshift)
    par = sig.parameters["prefilter_min_hits"]
    assert par.kind is inspect.Parameter.KEYWORD_ONLY
    assert par.default == 256


def test_forced_prefilter_matches_oracle_on_planted_splits() -> None:
    g, qs = _planted_genome(6_000, 5, 6)
    idx = _index({"c": g})
    df = _target(_run(idx, qs, prefilter_min_hits=0))
    for q in qs:
        res = frameshift_oracle(q, idx, idx.frames, None, 2, 4)
        want = sorted({(p.start, p.end, p.strand) for p in res.placements})
        sub = df[df["peptide"] == q]
        assert sorted(zip(sub["start"], sub["end"], sub["strand"], strict=True)) == want, q
