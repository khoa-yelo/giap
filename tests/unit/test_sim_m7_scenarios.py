"""giapsim scenarios for variants and novel ORFs (issue M7-12; GIAP_SPEC.md §11.3, §7.5.6, §7.6.6,
§8.4 suite "Novel-ORF class FDR", §8.5 attack 6, §11.6 TE2).

Examples 1-5, rules 1-3 and the acceptance items of docs/specs/M7-12-sim-variant-novel-scenarios.md
with lead decisions 1-4 (2026-10-22). The spec names ``tests/sim/test_m7_scenarios.py``; like the
M4-07/M5-04 simulator tests this file lives in tests/unit (tests/sim holds only the giapsim package
and tests/unit/sim/ is reserved by M6-03 on another branch). The 1-Mb examples are marked ``slow``;
the rule-4 runtime budget belongs in ``benchmarks/`` (see the report to the lead). Brute-force
recounts come from tests/helpers/m7_sim_oracle.py, which imports nothing from giapsim.

Interface choices pinned here (beyond the spec's four signatures):

* ``giapsim.scenarios.variant_f1(seed, size)``, ``novel_orf_power(seed, size)``,
  ``same_strain_insilico(seed)`` and ``isobaric_attack(seed)`` return a ``SimConfig`` that is
  simulated with ``giapsim.simulate(cfg, seed, outdir)``; ``size`` "100k" / "1M" sets
  ``genome_length`` 100_000 / 1_000_000; any other size raises ``ValueError``. The same seed gives
  an equal config.
* Every M7 scenario writes ``truth/alt_sites.tsv`` (header only when no SNV is planted) and lists
  it as ``SimResult.truth["alt_sites"]``. Columns, in order: ``contig, codon_start, strand, gene,
  ref_aa, alt_aa, n_true_peptides``. ``codon_start`` is the lowest 1-based coordinate of the codon
  on ``reference.fasta`` (as ``variants.tsv``, M7-03 rule 10); ``strand`` and ``gene`` are the
  reference gene's; ``ref_aa`` is the reference residue, ``alt_aa`` the sample residue (the two
  letters of ``strain_differences.tsv``'s ``aa_change``). Rows in (contig, codon_start, strand)
  order. Every non-synonymous SNV whose two residues are not ``*`` has a row; nonsense SNVs may
  be left out; nothing else appears.
* In-frame indels inside genes (lead decision 2) go to ``truth/indel_sites.tsv``
  (``SimResult.truth["indel_sites"]``) with at least the columns ``contig, pos_reference, ref,
  alt, gene, aa_change`` copied from ``strain_differences.tsv``; written whenever a strain is set.
  ``variant_f1`` plants indels (``indel_rate_per_kb > 0``).
* ``truth/novel_orfs.tsv`` in an M7 scenario keeps the seven M3-18 columns first and adds
  ``stop, true_start, n_true_peptides, rbs_planted`` (definitions in the oracle docstring).
  ``novel_orf_power`` lists ``nterm_extension`` rows there too (class ``nterm_extension``,
  ``start..end`` = true start to stop, ``aa`` = the true protein), so G3's truth set is one table.
  ``rbs_planted`` is true for every intergenic and nterm_extension row of ``novel_orf_power``
  and false for cds_outframe/antisense rows (found, not planted).
* ``novel_orf_power``, ``same_strain_insilico`` and ``isobaric_attack`` set no strain (truth and
  peptides then share ``genome.fasta``); ``novel_orf_power`` writes the precomputed layers
  (``rbs``, ``rbs_novel``) the RBS calibration reads.
* Truth readers (the spec's "≈ 40 lines of loaders the calibration tests share"):
  ``giapsim.truth.read_alt_sites(outdir)``, ``read_indel_sites(outdir)`` and
  ``read_novel_orfs(outdir)`` return pandas DataFrames with integer coordinates and counts and a
  boolean ``rbs_planted``.
* Example 5 is pinned on this branch's base (f6a8bf8, which includes M7-08's ``rbs_novel`` layer
  merged after c4d8841); ``config.yaml`` is left out of the hashes because it holds absolute paths.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import pytest
from helpers import m7_sim_oracle as O

import giapsim

SEED = 20261022
SEED_B = 4242

#: Pre-M7-12 outputs (f6a8bf8) of configs that set no M7-12 field (example 5).
SIM100K_SHA256 = {
    "annotation.gff3": "fa95ab8fe14014e3dc65af7003beed7647b561ebb6e81f5b40189097b4224e2e",
    "bakta_like.json": "ccf7c77f2a5b6b985ce28229c06f81eac8dda17897112a27cf7f379d893358b2",
    "genome.fasta": "dd10216ee21cb0ac13940b2f7d5a26bf0fdc8c2af54af8baf26b380fc65a8d69",
    "precomputed/operon_pairs.tsv": "84d2f05a42eefab68bc21569ecc2e6270edb797af86a3786fc7244710354098f",  # noqa: E501
    "precomputed/promoters.tsv": "14e304b6ff7536ec0a842c56d8de3f56062143f953149a369d744b1a2d9222b6",
    "precomputed/rbs.tsv": "3dbe7d15a7421ccc02d466847460b6aaf4cc749ed4696065e7ef389a9c1d0fa6",
    "precomputed/rbs_novel.tsv": "24f10c068fdfe924c19a845a970d54a5ca3615c22c4336a0336c13b25eb21025",
    "precomputed/regions.tsv": "7ae7a89e972819ee81f5390566dd5b8c97990780822536e2774d1a4bef3fd334",
    "precomputed/terminators.tsv": "ad8c479d19c0dae7494bcf93bbcac816e7740a182a3cc7d3de8fa666c47e3bdf",  # noqa: E501
    "precomputed/tx_rates.tsv": "4afd83ff497e0deac1f342b8545101e8a6b2b7c4e8fa12940e8251444b0b60d2",
    "truth/features.tsv": "7ce5191ca875a1ea02d0b58f2108cd7277a725f47b19477e92dca07cf1692f56",
    "truth/genes.tsv": "d2ec79f1ccae930b52dc3fa86430f71b35f9a9cfa1a934cc4cf70b4f2cff404a",
    "truth/novel_orfs.tsv": "57b8bf862b97cc98740f605389ee2c51d5b7da03f2e03ce6bb42827c7a7da9f2",
    "truth/promoters.tsv": "b85564f428548e8f2c1425dbd8c23bb9bfd59f10014388bd26582cfb191475fe",
    "truth/rbs.tsv": "c82822bd4836abb92e0a5823177af0e790cdb9a02b448f8b04f0f4b5b98b3a84",
    "truth/regions.tsv": "59b12822a325309aa6b08cc9a2d7771e8535a5ec2a0abe93d5522256e9491471",
    "truth/terminators.tsv": "084d949118b861b7c65f7f611474f61101b704e2fd28896491ff5ccc66faca83",
    "truth/tus.tsv": "b73ae5ce50c7a114c806e2fd189997e26d44a39fc2f5de5ab19914e719fa6721",
}  # fmt: skip
MIXED_SEED = 7
MIXED_SHA256 = {
    "abundance.tsv": "6e8e869668805ea15637e953b7850ef2f17ba01878b64800581bcd871a0955d8",
    "annotation.gff3": "bf1da82987863569db251f6386d17a3ac6454be6e2403fa0ce7f017b0cc0bd1f",
    "bakta_like.json": "ffd43720276e5e68deeac0d40ac4862889f033fe51195e2e7107dc6232010e48",
    "genome.fasta": "4ba749a484ce201a1ff6af659ead1ea7c2306fe63c56c9c063eaf31a2a845c62",
    "peptides.tsv": "d33e6722017a81789c3dd5e2e7db77c60e439e3e84041980a856773247c44369",
    "peptides_denovo.tsv": "dda6fa4327610cdbb38081d64c58576b06d7622e999c86c2370431324ae75bd0",
    "precomputed/operon_pairs.tsv": "666495c5b5d0b97a9c9f80007a7f659beb6cb9748d18d3bbd4831b85866c6446",  # noqa: E501
    "precomputed/promoters.tsv": "dc18ac15361ae0d058d0898eda8bc96a1512bea9f4a44826c4e37ac803fb658e",
    "precomputed/rbs.tsv": "62f1708ac7437eee26fdab6099de2da08910838056546f6483629210477ef6ce",
    "precomputed/rbs_novel.tsv": "587975169cb9019926ebba5f5d359fbec437d49576bbbb18b9f7ee198f5ffcd3",
    "precomputed/regions.tsv": "f5718b0001a9bc220c9151b1ec3b193091530fcd472731d9659925e937733cb1",
    "precomputed/terminators.tsv": "764c5c32adbaa3a952998c2cb960cd359dff5fe317bbb44480abf12264973912",  # noqa: E501
    "precomputed/tx_rates.tsv": "130025031d32c61b968ef63159fa4e3162d78767a9be98894f0ed9356421abde",
    "proteins.faa": "d8a02fe249ec1c59e7744418896f75fec842cc59a8437286d4a5f989949897d0",
    "reference.fasta": "fee1cb3a6a6d9d5ac10f3699be8e6f8676b9e316d5c475f44d7c2774b232c32b",
    "reference_bakta_like.json": "87756f79b84967a30301d0f400f78cbd8460d683eaa3cbb6faca4c046b845b09",
    "reference_peptides.tsv": "2612ab7d72ee5aa204a3b406ab64a89cca24ceef5af67cfc41ba7b526b1884d5",
    "truth/denovo.tsv": "6d31377e13b368e184ad207789c7a0b0b1230fa0289988da7707244092592117",
    "truth/denovo_precision.tsv": "2a7b79288af6ff53f5d66de7e2433ec12e35b9fa42b02b54e991177a3beb0879",  # noqa: E501
    "truth/features.tsv": "7ee8d28f13b9d46a3da754ea50806b5bce2b9b33fc4821f4f4dc9e1a97aa2454",
    "truth/genes.tsv": "10d8c1d25c25aca9354785c699ed00df172a06e3514ae3470f524ecfd18d93b7",
    "truth/novel_orfs.tsv": "9ea3bbe89cab7b7177417a774992349b2ec81f84d10b03c248d8a7c7d1488c42",
    "truth/peptides.tsv": "3724da73aa8a7ab74bffd45e4941b6cc26b6a143d02678f6cef76c084fe7f512",
    "truth/promoters.tsv": "eb0d45316c2bcf064b47826d3b15ab826134264ca2628460d1349bc5b8e22576",
    "truth/proteomics.json": "e0a9f04827f08d0b4ff573513c37e56bb58df3fbd9c769babdb16e0627e0d5ba",
    "truth/rbs.tsv": "a42d5c5255322ee4595d28bf7831ef7ab982c6c18f3115103036373d96c73632",
    "truth/regions.tsv": "8c31320cff1a63dbf46cc5bed7dff2d354035ba8f869175f3142036e51039fc1",
    "truth/strain_differences.tsv": "2580abb16875e259de9c79e96578b26337449171e299283f1e5d3b2dc84aea74",  # noqa: E501
    "truth/terminators.tsv": "b704a3813ea33d487d406892718151331b207674b9ed9a0ebb6a02eed11eecbf",
    "truth/tus.tsv": "16d116933f8ec6bc05da894d6fcfcbaace264103a6cfdcc382761a06d83cc1b2",
}  # fmt: skip
M7_SCENARIOS = ("variant_f1", "novel_orf_power", "same_strain_insilico", "isobaric_attack")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sc() -> Any:
    import giapsim.scenarios

    return giapsim.scenarios


def config(name: str, seed: int, size: str | None = None) -> Any:
    fn = getattr(sc(), name)
    return fn(seed) if size is None else fn(seed, size)


def run(name: str, seed: int, out: Path, size: str | None = None) -> Any:
    return giapsim.simulate(config(name, seed, size), seed, out)


def files(out: Path) -> dict[str, str]:
    return {
        p.relative_to(out).as_posix(): sha(p)
        for p in sorted(out.rglob("*"))
        if p.is_file() and p.name != "config.yaml"
    }


def mixed_config() -> Any:
    from giapsim.denovo import DenovoConfig
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.strain import StrainConfig

    return giapsim.SimConfig(
        genome_length=60_000,
        n_novel_orfs=(("intergenic", 3), ("cds_outframe", 2), ("antisense", 2)),
        n_nterm_extensions=3,
        strain=StrainConfig(snv_rate_per_kb=2.0, indel_rate_per_kb=0.5),
        proteomics=ProteomicsConfig(pi_wrong=0.05),
        denovo=DenovoConfig(),
        rbs_strength=0.8,
        write_precomputed=True,
    )


@pytest.fixture(scope="module")
def outs(tmp_path_factory: pytest.TempPathFactory) -> Any:
    cache: dict[tuple[str, int, str | None], tuple[Path, Any]] = {}

    def make(name: str, seed: int = SEED, size: str | None = None) -> tuple[Path, Any]:
        key = (name, seed, size)
        if key not in cache:
            out = tmp_path_factory.mktemp(f"{name}_{seed}_{size}")
            cache[key] = (out, run(name, seed, out, size))
        return cache[key]

    return make


@pytest.fixture(scope="module")
def vf(outs: Any) -> Any:
    return lambda: outs("variant_f1", SEED, "100k")


@pytest.fixture(scope="module")
def nop(outs: Any) -> Any:
    return lambda: outs("novel_orf_power", SEED, "100k")


@pytest.fixture(scope="module")
def ssi(outs: Any) -> Any:
    return lambda: outs("same_strain_insilico", SEED)


@pytest.fixture(scope="module")
def iso(outs: Any) -> Any:
    return lambda: outs("isobaric_attack", SEED)


def scenario_out(request: pytest.FixtureRequest, name: str) -> tuple[Path, Any]:
    fixture = {"variant_f1": "vf", "novel_orf_power": "nop", "same_strain_insilico": "ssi",
               "isobaric_attack": "iso"}[name]  # fmt: skip
    out: tuple[Path, Any] = request.getfixturevalue(fixture)()
    return out


# ---------------------------------------------------------------- the interface


@pytest.mark.parametrize("name", M7_SCENARIOS)
def test_scenarios_return_simconfig_deterministically(name: str) -> None:
    a, b = config(name, SEED), config(name, SEED)
    assert isinstance(a, giapsim.SimConfig)
    assert a == b


@pytest.mark.parametrize("name", ["variant_f1", "novel_orf_power"])
@pytest.mark.parametrize(("size", "length"), [("100k", 100_000), ("1M", 1_000_000)])
def test_size_sets_genome_length(name: str, size: str, length: int) -> None:
    assert config(name, SEED, size).genome_length == length
    assert config(name, SEED).genome_length == 100_000  # default "100k"


@pytest.mark.parametrize("name", ["variant_f1", "novel_orf_power"])
def test_unknown_size_raises(name: str) -> None:
    with pytest.raises(ValueError):
        config(name, SEED, "10M")


@pytest.mark.parametrize("name", M7_SCENARIOS)
def test_every_scenario_writes_alt_sites(name: str, request: pytest.FixtureRequest) -> None:
    out, res = scenario_out(request, name)
    path = out / "truth" / "alt_sites.tsv"
    assert path.exists()
    assert Path(res.truth["alt_sites"]).resolve() == path.resolve()
    head = O.header(path)  # decision 11 appends frameshift_upstream (lead edit)
    assert head[: len(O.ALT_COLUMNS)] == O.ALT_COLUMNS and "frameshift_upstream" in head


# ---------------------------------------------------------------- example 1: variant_f1


def test_example1_variant_f1_files(vf: Any) -> None:
    out, res = vf()
    cfg = config("variant_f1", SEED, "100k")
    assert cfg.strain is not None and cfg.proteomics is not None and cfg.denovo is not None
    assert cfg.strain.snv_rate_per_kb > 0 and cfg.strain.indel_rate_per_kb > 0
    for f in ("reference.fasta", "reference_bakta_like.json", "peptides_denovo.tsv",
              "reference_peptides.tsv", "truth/strain_differences.tsv", "truth/alt_sites.tsv",
              "truth/indel_sites.tsv", "truth/denovo.tsv"):  # fmt: skip
        assert (out / f).exists(), f
    assert res.reference_fasta is not None


def test_example1_reference_peptides_pair_denovo_spectra(vf: Any) -> None:
    out, _ = vf()
    dn = {(r["run"], r["spectrum_id"]) for r in O.rows(out / "peptides_denovo.tsv")}
    ref = [(r["run"], r["spectrum_id"]) for r in O.rows(out / "reference_peptides.tsv")]
    assert ref and dn
    assert len(ref) == len(set(ref)), "one reference identification per spectrum"
    assert set(ref) <= dn
    assert {r["source"] for r in O.rows(out / "peptides_denovo.tsv")} == {"denovo"}
    assert {r["source"] for r in O.rows(out / "reference_peptides.tsv")} == {"database"}


def test_example1_peptides_are_sample_derived(vf: Any) -> None:
    """Rule 1: peptides come from the sample; a correct peptide over a non-synonymous SNV carries
    the sample residue, so it occurs in the sample protein and not at that reference codon."""
    out, _ = vf()
    sample = O.bakta_features(out / "bakta_like.json")
    correct = O.correct_peptide_rows(out)
    assert correct
    for r in correct[:500]:
        if r["locus_tag"] in sample:
            assert O.il(r["peptide"]) in O.il(sample[r["locus_tag"]]["aa"]), r


def test_alt_sites_match_strain_differences(vf: Any) -> None:
    out, _ = vf()
    alt = O.rows(out / "truth" / "alt_sites.tsv")
    # decision 11 (lead edit): genes with a frameshifting indel are checked by the round-2 lift
    fs = {r["gene"] for r in O.rows(out / "truth" / "indel_sites.tsv") if O.truthy(r["frameshift"])}
    alt = [r for r in alt if r["gene"] not in fs]
    snvs = [s for s in O.nonsyn_snvs(out) if s["gene"] not in fs]
    ref_genes = O.bakta_features(out / "reference_bakta_like.json")
    ref_seq = O.read_fasta(out / "reference.fasta")
    want = set()
    allowed = set()
    for s in snvs:
        g = ref_genes[s["gene"]]
        cs = O.ref_codon_start(g, int(s["pos_reference"]), len(ref_seq[s["contig"]]))
        key = (s["contig"], cs, g["strand"], s["gene"], s["ref_aa"], s["alt_aa"])
        allowed.add(key)
        if "*" not in (s["ref_aa"], s["alt_aa"]):
            want.add(key)
    got = {(r["contig"], int(r["codon_start"]), r["strand"], r["gene"], r["ref_aa"], r["alt_aa"])
           for r in alt}  # fmt: skip
    assert len(got) == len(alt), "one row per codon"
    assert want <= got, sorted(want - got)[:5]
    assert got <= allowed, sorted(got - allowed)[:5]
    assert want, "the 100-kb variant_f1 plants non-synonymous SNVs"


def test_alt_sites_codons_read_on_reference(vf: Any) -> None:
    """Rule 1: ``codon_start``/``strand`` address the reference codon encoding ``ref_aa``."""
    out, _ = vf()
    ref_seq = O.read_fasta(out / "reference.fasta")
    ref_genes = O.bakta_features(out / "reference_bakta_like.json")
    for r in O.rows(out / "truth" / "alt_sites.tsv"):
        seq = ref_seq[r["contig"]]
        c0 = int(r["codon_start"]) - 1
        codon = O.strand_slice(seq, c0, c0 + 3, r["strand"])
        g = ref_genes[r["gene"]]
        first = (int(g["start"]) - 1 if r["strand"] == "+" else int(g["stop"]) - 3) % len(seq)
        aa = "M" if c0 == first and codon in O.STARTS else O.CODE11[codon]
        assert aa == r["ref_aa"], (r, codon)
        assert r["alt_aa"] != r["ref_aa"], r


def test_alt_sites_sorted(vf: Any) -> None:
    out, _ = vf()
    order = {c: k for k, c in enumerate(O.read_fasta(out / "reference.fasta"))}
    keys = [(order[r["contig"]], int(r["codon_start"]), r["strand"])
            for r in O.rows(out / "truth" / "alt_sites.tsv")]  # fmt: skip
    assert keys == sorted(keys)


def test_rule2_alt_site_counts_match_bruteforce(vf: Any) -> None:
    out, _ = vf()
    sample = O.read_fasta(out / "genome.fasta")
    ref_genes = O.bakta_features(out / "reference_bakta_like.json")
    ref_seq = O.read_fasta(out / "reference.fasta")
    by_codon = {}
    for s in O.nonsyn_snvs(out):
        g = ref_genes[s["gene"]]
        cs = O.ref_codon_start(g, int(s["pos_reference"]), len(ref_seq[s["contig"]]))
        by_codon[(s["contig"], cs, s["gene"])] = s
    correct = O.correct_peptide_rows(out)
    for r in O.rows(out / "truth" / "alt_sites.tsv"):
        s = by_codon[(r["contig"], int(r["codon_start"]), r["gene"])]
        n = O.recount_site_peptides(out, s["gene"], s["contig"], int(s["pos_sample"]),
                                    len(sample[s["contig"]]), correct)  # fmt: skip
        assert int(r["n_true_peptides"]) == n, (r, n)


def test_decision2_indels_excluded_and_reported_separately(vf: Any) -> None:
    out, res = vf()
    path = out / "truth" / "indel_sites.tsv"
    assert Path(res.truth["indel_sites"]).resolve() == path.resolve()
    assert set(O.INDEL_COLUMNS) <= set(O.header(path))
    # decision 11 (lead edit): frameshifts are listed too, flagged; the in-frame ones are as before
    rows = [r for r in O.rows(path) if not O.truthy(r["frameshift"])]
    got = {tuple(r[c] for c in O.INDEL_COLUMNS) for r in rows}
    assert got == set(O.inframe_indels(out))
    # no alt_sites row sits on an indel
    alt = O.rows(out / "truth" / "alt_sites.tsv")
    for contig, pos, ref, _alt, _gene, _ch in got:
        p = int(pos)
        for r in alt:
            if r["contig"] == contig:
                c = int(r["codon_start"])
                assert not c <= p < c + max(3, len(ref)), (r, pos)


@pytest.mark.slow
def test_example1_1mb_has_20_sites_with_two_peptides(outs: Any) -> None:
    out, _ = outs("variant_f1", SEED, "1M")
    alt = O.rows(out / "truth" / "alt_sites.tsv")
    assert sum(int(r["n_true_peptides"]) >= 2 for r in alt) >= 20
    assert O.rows(out / "truth" / "indel_sites.tsv"), "1 Mb plants in-frame indels in genes"


# ---------------------------------------------------------------- example 2: novel_orf_power


def novel_rows(out: Path) -> list[dict[str, str]]:
    return O.rows(out / "truth" / "novel_orfs.tsv")


def test_novel_truth_columns(nop: Any) -> None:
    out, _ = nop()
    h = O.header(out / "truth" / "novel_orfs.tsv")
    assert h[:7] == O.NOVEL_COLUMNS_M318
    assert set(O.NOVEL_NEW_COLUMNS) <= set(h[7:])


def test_novel_orf_power_classes_and_layers(nop: Any) -> None:
    out, _ = nop()
    cfg = config("novel_orf_power", SEED, "100k")
    assert cfg.strain is None and cfg.proteomics is not None and cfg.write_precomputed
    assert (out / "precomputed" / "rbs.tsv").exists()
    assert (out / "precomputed" / "rbs_novel.tsv").exists()
    classes = {r["class"] for r in novel_rows(out)}
    assert classes == set(O.NOVEL_CLASSES_M7)


def test_novel_ids_unique(nop: Any) -> None:
    out, _ = nop()
    ids = [r["orf_id"] for r in novel_rows(out)]
    assert len(ids) == len(set(ids))


def test_novel_coordinates_stop_true_start_and_protein(nop: Any) -> None:
    """Rule 1: ``start..end`` (1-based closed) spans true start to stop on the receiving genome;
    ``stop`` is the 1-based lowest base of the stop codon (decision 16); ``true_start`` is the
    transcription-first base of the start codon; ``aa`` is the translation (initiator as M,
    stop dropped)."""
    out, _ = nop()
    genome = O.receiving_genome(out)
    for r in novel_rows(out):
        seq = genome[r["contig"]]
        s1, e1, st = int(r["start"]), int(r["end"]), r["strand"]
        nt = O.strand_slice(seq, s1 - 1, e1, st)
        assert len(nt) % 3 == 0, r
        assert nt[:3] in O.STARTS and nt[-3:] in O.STOPS, (r["orf_id"], nt[:3], nt[-3:])
        assert O.translate_orf(nt) == r["aa"], r["orf_id"]
        assert int(r["stop"]) == ((e1 - 2) if st == "+" else s1), r  # 1-based (decision 16)
        assert int(r["true_start"]) == (s1 if st == "+" else e1), r


def test_novel_nterm_rows_match_features(nop: Any) -> None:
    out, _ = nop()
    feats = {f["feature_id"]: f for f in O.rows(out / "truth" / "features.tsv")}
    nterm = [r for r in novel_rows(out) if r["class"] == "nterm_extension"]
    assert nterm
    for r in nterm:
        f = O.nterm_feature(out, r)
        assert f["feature_id"] in feats
        # the annotated CDS shares the stop and sits strictly inside the true span
        if r["strand"] == "+":
            assert int(f["end"]) == int(r["end"]) and int(f["start"]) > int(r["start"])
        else:
            assert int(f["start"]) == int(r["start"]) and int(f["end"]) < int(r["end"])


def test_novel_intergenic_length_range(nop: Any) -> None:
    """Scenarios use 30-80 aa planted ORFs (≥ 2 tryptic peptides typical)."""
    out, _ = nop()
    inter = [r for r in novel_rows(out) if r["class"] == "intergenic"]
    assert inter
    for r in inter:
        assert 30 <= len(r["aa"]) <= 80, (r["orf_id"], len(r["aa"]))


def test_rule2_novel_counts_match_bruteforce(nop: Any) -> None:
    out, _ = nop()
    correct = O.correct_peptide_rows(out)
    for r in novel_rows(out):
        assert int(r["n_true_peptides"]) == O.recount_orf_peptides(out, r, correct), r["orf_id"]


def test_rbs_planted_flags_and_motifs(nop: Any) -> None:
    """RBS (AGGAGG, spacer 5-10) before every planted intergenic start and N-terminal true
    start; found ORFs (cds_outframe, antisense) are not planted."""
    out, _ = nop()
    genome = O.receiving_genome(out)
    for r in novel_rows(out):
        planted = O.truthy(r["rbs_planted"])
        if r["class"] in ("intergenic", "nterm_extension"):
            assert planted, r["orf_id"]
            first0 = int(r["true_start"]) - 1
            assert O.has_rbs_motif(genome[r["contig"]], r["strand"], first0), r["orf_id"]
        else:
            assert not planted, r["orf_id"]


def test_rbs_novel_layer_lists_planted_starts(nop: Any) -> None:
    out, _ = nop()
    keys = {(r["contig"], r["strand"], int(r["start"]))
            for r in O.rows(out / "precomputed" / "rbs_novel.tsv")}  # fmt: skip
    for r in novel_rows(out):
        if O.truthy(r["rbs_planted"]):
            assert (r["contig"], r["strand"], int(r["true_start"])) in keys, r["orf_id"]


@pytest.mark.slow
def test_example2_1mb_ten_per_class_most_with_two_peptides(outs: Any) -> None:
    out, _ = outs("novel_orf_power", SEED, "1M")
    rws = novel_rows(out)
    for cls in O.NOVEL_CLASSES_M7:
        assert sum(r["class"] == cls for r in rws) >= 10, cls
    two = sum(int(r["n_true_peptides"]) >= 2 for r in rws)
    assert two > len(rws) / 2, (two, len(rws))
    genome = O.receiving_genome(out)
    for r in rws:
        if r["class"] == "intergenic":
            assert O.has_rbs_motif(genome[r["contig"]], r["strand"], int(r["true_start"]) - 1)


# ---------------------------------------------------------------- example 3: same_strain_insilico


def test_example3_same_strain_no_errors(ssi: Any) -> None:
    out, _ = ssi()
    cfg = config("same_strain_insilico", SEED)
    assert cfg.strain is None and cfg.proteomics is not None
    assert not (out / "reference.fasta").exists()
    db = O.rows(out / "peptides.tsv")
    assert db and {r["search_db"] for r in db} == {"same"}
    for t in O.rows(out / "truth" / "peptides.tsv"):
        if t["origin"] != "decoy":
            assert not O.truthy(t["is_wrong"]), t
    if (out / "truth" / "denovo.tsv").exists():
        assert {t["error_type"] for t in O.rows(out / "truth" / "denovo.tsv")} == {"none"}
    assert O.rows(out / "truth" / "alt_sites.tsv") == []


def test_example3_every_tryptic_peptide_of_every_annotated_protein(ssi: Any) -> None:
    out, _ = ssi()
    feats = O.bakta_features(out / "bakta_like.json")
    want: set[str] = set()
    for f in feats.values():
        if not f.get("pseudo"):
            want |= O.tryptic(f["aa"])
    db = O.rows(out / "peptides.tsv")
    got = {r["peptide"] for r in db if not O.truthy(r["is_decoy"])}
    assert got == want, (len(got - want), len(want - got))


def test_example3_non_atg_starts_give_m_nterm_peptides(ssi: Any) -> None:
    out, _ = ssi()
    feats = O.bakta_features(out / "bakta_like.json")
    db = {r["peptide"] for r in O.rows(out / "peptides.tsv") if not O.truthy(r["is_decoy"])}
    alt_truth = {r["peptide"] for r in O.rows(out / "truth" / "peptides.tsv")
                 if O.truthy(r["alt_start_nterm"])}  # fmt: skip
    seen = 0
    for f in feats.values():
        if f.get("pseudo") or f["nt"][:3] == "ATG":
            continue
        nterm = sorted(p for p in O.tryptic(f["aa"]) if f["aa"].startswith(p))
        for p in nterm:
            assert p.startswith("M") and p in db, (f["locus"], p)
            assert p in alt_truth, (f["locus"], p)
            seen += 1
    assert seen > 0, "at least one non-ATG gene with an N-terminal peptide"


# ---------------------------------------------------------------- example 4: isobaric_attack


def test_example4_isobaric_attack(iso: Any) -> None:
    out, _ = iso()
    cfg = config("isobaric_attack", SEED)
    assert cfg.strain is None and cfg.denovo is not None
    assert O.rows(out / "truth" / "alt_sites.tsv") == []
    assert not (out / "truth" / "strain_differences.tsv").exists()
    w = {k for k, v in dict(cfg.denovo.error_weights).items() if v > 0}
    assert w == O.ISOBARIC_TYPES
    kinds = {t["error_type"] for t in O.rows(out / "truth" / "denovo.tsv")} - {"none"}
    assert kinds == O.ISOBARIC_TYPES, kinds


# ---------------------------------------------------------------- contract round trips


@pytest.mark.parametrize("name", M7_SCENARIOS)
def test_genomes_and_annotations_load(name: str, request: pytest.FixtureRequest) -> None:
    from giap.annotate.bakta import parse_bakta_json
    from giap.io.fasta import read_genome_fasta

    out, _ = scenario_out(request, name)
    pairs = [("genome.fasta", "bakta_like.json")]
    if (out / "reference.fasta").exists():
        pairs.append(("reference.fasta", "reference_bakta_like.json"))
    for fa, js in pairs:
        recs = read_genome_fasta(out / fa)
        assert recs
        ann = parse_bakta_json(out / js)
        assert ann is not None


@pytest.mark.parametrize("name", M7_SCENARIOS)
def test_peptide_tables_pass_m106(name: str, request: pytest.FixtureRequest) -> None:
    from giap.config.model import InputConfig
    from giap.io.peptides import read_peptide_table, validate_peptides

    out, _ = scenario_out(request, name)
    tables = [p for p in ("peptides.tsv", "peptides_denovo.tsv", "reference_peptides.tsv")
              if (out / p).exists()]  # fmt: skip
    assert tables
    for t in tables:
        raw, issues = read_peptide_table(out / t)
        df, issues2 = validate_peptides(raw, InputConfig())
        errors = [i for i in [*issues, *issues2] if i.severity == "error"]
        assert not errors, (t, errors)
        assert len(df) == len(O.rows(out / t))


def test_truth_readers(vf: Any, nop: Any) -> None:
    from giapsim import truth

    out, _ = vf()
    alt = truth.read_alt_sites(out)
    assert list(alt.columns)[: len(O.ALT_COLUMNS)] == O.ALT_COLUMNS
    assert len(alt) == len(O.rows(out / "truth" / "alt_sites.tsv"))
    for c in ("codon_start", "n_true_peptides"):
        assert str(alt[c].dtype).lower().startswith("int"), c
    ind = truth.read_indel_sites(out)
    assert len(ind) == len(O.rows(out / "truth" / "indel_sites.tsv"))
    out2, _ = nop()
    nov = truth.read_novel_orfs(out2)
    assert len(nov) == len(novel_rows(out2))
    for c in ("start", "end", "stop", "true_start", "n_true_peptides"):
        assert str(nov[c].dtype).lower().startswith("int"), c
    assert str(nov["rbs_planted"].dtype).lower() in ("bool", "boolean")


def test_truth_reader_empty_alt_sites(iso: Any) -> None:
    from giapsim import truth

    out, _ = iso()
    alt = truth.read_alt_sites(out)
    assert len(alt) == 0
    assert list(alt.columns)[: len(O.ALT_COLUMNS)] == O.ALT_COLUMNS


# ---------------------------------------------------------------- rule 3 / example 5: determinism


@pytest.mark.parametrize("name", M7_SCENARIOS)
def test_same_seed_byte_identical(
    name: str, request: pytest.FixtureRequest, tmp_path: Path
) -> None:
    out, _ = scenario_out(request, name)
    size = "100k" if name in ("variant_f1", "novel_orf_power") else None
    run(name, SEED, tmp_path, size)
    assert files(tmp_path) == files(out)


@pytest.mark.parametrize(
    ("name", "truth_file"),
    [("variant_f1", "alt_sites.tsv"), ("novel_orf_power", "novel_orfs.tsv")],
)
def test_other_seed_changes_truth(name: str, truth_file: str, outs: Any) -> None:
    a, _ = outs(name, SEED, "100k")
    b, _ = outs(name, SEED_B, "100k")
    assert sha(a / "truth" / truth_file) != sha(b / "truth" / truth_file)


def test_example5_sim100k_byte_identical(tmp_path: Path) -> None:
    giapsim.simulate(giapsim.SIM100K, giapsim.SIM100K_SEED, tmp_path)
    assert files(tmp_path) == SIM100K_SHA256


def test_example5_mixed_config_byte_identical(tmp_path: Path) -> None:
    """Strain with indels, novel ORFs, N-terminal extensions, proteomics, de novo and layers, with
    no M7-12 field set: same files, same bytes, no alt_sites/indel_sites written."""
    giapsim.simulate(mixed_config(), MIXED_SEED, tmp_path)
    assert files(tmp_path) == MIXED_SHA256


def test_default_novel_length_range_stays_10_30(tmp_path: Path) -> None:
    """The configurable range defaults to M3-18's 10-30 aa for planted intergenic ORFs."""
    cfg = giapsim.SimConfig(genome_length=80_000, n_novel_orfs=(("intergenic", 12),))
    giapsim.simulate(cfg, 3, tmp_path)
    inter = [r for r in O.rows(tmp_path / "truth" / "novel_orfs.tsv") if r["class"] == "intergenic"]
    assert len(inter) == 12
    assert all(10 <= len(r["aa"]) <= 30 for r in inter)
    assert O.header(tmp_path / "truth" / "novel_orfs.tsv") == O.NOVEL_COLUMNS_M318
