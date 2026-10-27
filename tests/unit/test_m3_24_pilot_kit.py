"""The pilot kit of GIPPA_SPEC.md §15: cfg files, README, ``db_a2`` and the three sbatch scripts
(M3-24 and its lead decisions 1-6 in docs/specs/M3-24-pilot-kit.md).

Static checks only: files are read, the scripts are parsed with ``bash -n``. Nothing is submitted,
no container runs and nothing reaches the network (the stubbed runs live in
``test_m3_24_diann_pilot_run.py``). Decision 6: M3 prepares the kit, it does not run it.

Conventions settled by the test author (the cfg grammar is documented in
``tests/helpers/pilot_kit.py``):

* Each cfg holds the flags of both DIA-NN 2.x steps under ``[step1]`` and ``[step2]`` headers;
  the two arms differ only in ``--fasta`` (spec acceptance criterion 2).
* Digest flags are compared against ``gippa.mapping.digest`` itself, not against a second copy of
  the numbers: ``--cut`` is derived from ``cut_sites``, the lengths and missed cleavages from the
  ``digest`` signature, and ``--met-excision`` must be absent because GIPPA digests the protein as
  written (digest.py lead decision 1).
* Every flag of the §15.2 template that the cfg files do not use must be named in the README
  (decision 5: "any flag ... that 2.3.1 does not accept is listed in the README").
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest
from helpers.pilot_kit import (
    BINDS,
    CFG,
    DIANN_BIN,
    DIANN_SIF,
    DIANN_SIF_SHA256,
    DIANN_VERSION,
    ID_MAP,
    ID_MAP_COLUMNS,
    MANIFEST,
    PILOT,
    README,
    RECORDED_NISSLE,
    REPO,
    SBATCH_SUPPLIED,
    SCRIPTS,
    SPEC_VERSION,
    STRAIN_DB,
    gippa_digest_flags,
    gippa_excises_initiator_met,
    names,
    read_cfg,
    spec_template_flags,
    values,
)

ARMS = ("A1", "A2")
MAX_SMALL_BYTES = 1_000_000  # decision 3: nothing big is committed under validation/
R11_MAX_BYTES = 5 * 1024 * 1024  # R11: no single committed file exceeds 5 MB


def _text(path: Path) -> str:
    assert path.is_file(), f"{path.relative_to(REPO)} does not exist (M3-24)"
    return path.read_text()


def _code(path: Path) -> str:
    """A shell script without comment lines (``#SBATCH`` directives kept)."""
    return "\n".join(ln for ln in _text(path).splitlines()
                     if not ln.lstrip().startswith("#") or ln.startswith("#SBATCH"))  # fmt: skip


def _directive(path: Path, key: str) -> str | None:
    m = re.search(rf"^#SBATCH\s+--{key}(?:[= ](\S+))?", _text(path), re.M)
    return (m.group(1) or "") if m else None


def _mem_gb(value: str) -> float:
    m = re.fullmatch(r"(\d+(?:\.\d+)?)\s*([KMGT]?)B?", value.strip(), re.I)
    assert m, f"unparsable --mem {value!r}"
    return float(m.group(1)) * {"": 1 / 1024, "K": 1 / 1024**2, "M": 1 / 1024, "G": 1,
                                "T": 1024}[m.group(2).upper()]  # fmt: skip


# ---------------------------------------------------------------- deliverables exist


@pytest.mark.parametrize(
    "path",
    [README, CFG["A1"], CFG["A2"], STRAIN_DB, MANIFEST, *SCRIPTS.values()],
    ids=lambda p: str(p.relative_to(REPO)),
)
def test_deliverable_exists(path: Path) -> None:
    assert path.is_file(), f"{path.relative_to(REPO)} does not exist (M3-24 deliverables)"


def test_db_a2_has_an_id_map() -> None:
    """Decision 10: the name ``gippa makedb`` writes, ``db_id_map.tsv`` (the spec's is a typo)."""
    assert ID_MAP.is_file(), f"{ID_MAP.relative_to(REPO)} does not exist"
    header = ID_MAP.read_text().splitlines()[0].split("\t")
    assert header == list(ID_MAP_COLUMNS), header


# ---------------------------------------------------------------- the two cfg files


@pytest.mark.parametrize("arm", ARMS)
def test_cfg_has_both_steps(arm: str) -> None:
    """Decision 4: the predict step and the search step both live in the cfg."""
    cfg = read_cfg(CFG[arm])
    assert set(cfg) == {"step1", "step2"}, sorted(cfg)
    assert cfg["step1"] and cfg["step2"], {k: len(v) for k, v in cfg.items()}


@pytest.mark.parametrize("arm", ARMS)
def test_step1_predicts_a_library_and_step2_searches(arm: str) -> None:
    """Decision 4: ``--gen-spec-lib --predictor`` in step 1, no raw files anywhere in the cfg."""
    cfg = read_cfg(CFG[arm])
    step1, step2 = names(cfg["step1"]), names(cfg["step2"])
    assert {"--gen-spec-lib", "--predictor", "--fasta-search"} <= step1, sorted(step1)
    assert "--f" not in step1, "step 1 takes no raw files (decision 4)"
    assert "--f" not in step2, "the raw files come from $RAW_DIR, not from the cfg"
    assert "--predictor" not in step2, (
        "step 2 searches the predicted library; --predictor belongs to step 1 only"
    )
    assert "--lib" not in step1, "step 1 builds the library it must not consume"


@pytest.mark.parametrize("arm", ARMS)
def test_digest_flags_match_gippa_digest(arm: str) -> None:
    """Decision 5: the DIA-NN digest is GIPPA's digest (M3-22, src/gippa/mapping/digest.py)."""
    cfg = read_cfg(CFG[arm])
    expected = gippa_digest_flags()
    got = {flag: vals for flag in expected if (vals := values(cfg["step1"], flag))}
    missing = sorted(set(expected) - set(got))
    assert not missing, f"{CFG[arm].name} [step1]: digest flags missing: {missing}"
    for flag, want in expected.items():
        assert got[flag] == [want], f"{CFG[arm].name} [step1] {flag}: {got[flag]} != [{want!r}]"
    for flag, want in expected.items():  # if repeated in step 2, with the same value
        assert values(cfg["step2"], flag) in ([], [want]), (
            f"{CFG[arm].name} [step2] {flag} disagrees with [step1]"
        )


@pytest.mark.parametrize("arm", ARMS)
def test_no_met_excision(arm: str) -> None:
    """GIPPA digests the protein as written (digest.py decision 1), so DIA-NN must not excise Met."""
    cfg = read_cfg(CFG[arm])
    used = any("--met-excision" in names(section) for section in cfg.values())
    assert used == gippa_excises_initiator_met(), (
        f"{CFG[arm].name}: --met-excision is "
        f"{'absent' if not used else 'used'} while gippa.mapping.digest "
        f"{'excises' if gippa_excises_initiator_met() else 'keeps'} the initiator Met"
    )


@pytest.mark.parametrize("arm", ARMS)
def test_search_step_reports_decoys_and_fixes_qvalue(arm: str) -> None:
    """§15.2 step 5: decoy reporting, so GIPPA can use the picked tier."""
    step2 = read_cfg(CFG[arm])["step2"]
    assert "--report-decoys" in names(step2), sorted(names(step2))
    assert values(step2, "--qvalue") == ["0.01"], values(step2, "--qvalue")


def test_a1_searches_the_reference_proteome() -> None:
    cfg = read_cfg(CFG["A1"])
    fastas = values(cfg["step1"], "--fasta")
    assert any("UP000000625" in f for f in fastas), fastas
    assert any("contaminant" in f.lower() for f in fastas), fastas
    assert not any("db_a2" in f for f in fastas), fastas


def test_a2_searches_the_makedb_strain_database() -> None:
    cfg = read_cfg(CFG["A2"])
    fastas = values(cfg["step1"], "--fasta")
    assert any(re.search(r"db_a2/strain_db\.faa", f) for f in fastas), fastas
    assert any("contaminant" in f.lower() for f in fastas), fastas
    assert not any("UP000000625" in f for f in fastas), fastas


def test_the_two_arms_differ_only_in_fasta() -> None:
    """Spec acceptance criterion: same flag set, only ``--fasta`` differs."""
    a1, a2 = (read_cfg(CFG[arm]) for arm in ARMS)
    assert set(a1) == set(a2)
    for step in sorted(a1):
        left = [f for f in a1[step] if f[0] != "--fasta"]
        right = [f for f in a2[step] if f[0] != "--fasta"]
        assert left == right, f"[{step}] differs beyond --fasta:\n{left}\n{right}"


def test_template_flags_are_used_or_explained_in_the_readme() -> None:
    """Decision 5: a §15.2 flag that 2.3.1 does not take is named in the README."""
    used = {flag for arm in ARMS for section in read_cfg(CFG[arm]).values()
            for flag in names(section)}  # fmt: skip
    text = _text(README)
    unexplained = sorted(
        flag for flag in spec_template_flags() - SBATCH_SUPPLIED - used if flag not in text
    )
    assert not unexplained, (
        f"§15.2 template flags neither used in the cfg files nor named in the README: {unexplained}"
    )


# ---------------------------------------------------------------- the README


def test_readme_records_the_version_the_image_and_its_sha256() -> None:
    """Decision 2: the exact version string, the image path and its sha256."""
    text = _text(README)
    assert DIANN_VERSION in text, f"README does not record {DIANN_VERSION!r}"
    assert DIANN_SIF in text, f"README does not record the image path {DIANN_SIF}"
    assert DIANN_BIN in text, f"README does not record the binary inside the image {DIANN_BIN}"
    assert re.search(r"\b[0-9a-f]{64}\b", text), "README records no sha256 of the image"


def test_readme_has_a_version_actually_used_section() -> None:
    """Decision 2: note 2.3.1 against the spec's 2.7.0; do not rewrite the spec."""
    text = _text(README)
    m = re.search(r"^#{1,6}\s*Version actually used\s*$", text, re.M | re.I)
    assert m, "README has no 'Version actually used' section (decision 2)"
    section = re.split(r"^#{1,6}\s", text[m.end() :], maxsplit=1, flags=re.M)[0]
    assert "2.3.1" in section and SPEC_VERSION in section, section


def test_readme_gives_the_pilot_commands() -> None:
    """Deliverable 6: the commands for §15.2 steps 4-9 and who does what."""
    text = _text(README)
    for script in SCRIPTS:
        assert f"scripts/slurm/{script}.sbatch" in text, f"README does not run {script}.sbatch"
    assert re.search(r"\bsbatch\s", text), (
        "README must submit the work with sbatch, never a login node"
    )
    for outdir in ("out_a1", "out_a2", "out_a3"):
        assert re.search(rf"gippa run[^\n]*--outdir\s+{outdir}\b", text), (
            f"README misses the §15.2 `gippa run … --outdir {outdir}` line"
        )
    assert "--reference-peptides" in text, "the A3 run passes --reference-peptides a2.tsv.gz"
    assert "Khoa" in text, "deliverable 6: who does what"


def test_readme_warns_that_the_image_is_never_copied() -> None:
    """Decision 3: ``DIANN_SIF`` is documented as overridable; the image stays where it is."""
    text = _text(README)
    assert "DIANN_SIF" in text, "README does not document the DIANN_SIF override (decision 3)"


def test_readme_records_the_pinned_sha256() -> None:
    """Decision 9: the README is the record, and the checksum is pinned, not re-derived."""
    assert DIANN_SIF_SHA256 in _text(README), (
        f"README does not record sha256 {DIANN_SIF_SHA256} of {DIANN_SIF}"
    )


@pytest.mark.slow
def test_the_pinned_sha256_is_still_the_image_on_disk() -> None:
    """The pin of decision 9 against the real file: a replaced image fails loudly."""
    image = Path(DIANN_SIF)
    if not image.is_file():
        pytest.skip(f"{DIANN_SIF} is not readable here")
    digest = hashlib.sha256()
    with image.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            digest.update(chunk)
    assert digest.hexdigest() == DIANN_SIF_SHA256, (
        f"{DIANN_SIF} now hashes to {digest.hexdigest()}, not the pinned {DIANN_SIF_SHA256}; "
        "re-pin it in tests/helpers/pilot_kit.py and in validation/pilot/README.md"
    )


# ---------------------------------------------------------------- db_a2 (gippa makedb)


def _recorded_counts() -> dict[str, int]:
    with gzip.open(RECORDED_NISSLE, "rt") as fh:
        data = json.load(fh)
    out: dict[str, int] = {"cds": 0, "sorf": 0}
    for feature in data["features"]:
        if feature["type"] in out:
            out[feature["type"]] += 1
    return out


def test_db_a2_manifest_records_the_recorded_nissle_json() -> None:
    """Acceptance criterion: the manifest lists the annotation's SHA-256 (§15.2 step 2)."""
    manifest = json.loads(_text(MANIFEST))
    annotation = manifest["inputs"]["annotation"]
    assert "nissle" in annotation["path"].lower(), annotation["path"]
    raw = RECORDED_NISSLE.read_bytes()
    with gzip.open(RECORDED_NISSLE, "rb") as fh:
        plain = fh.read()
    allowed = {hashlib.sha256(raw).hexdigest(), hashlib.sha256(plain).hexdigest()}
    assert annotation["sha256"] in allowed, (
        f"{MANIFEST.name} annotation sha256 {annotation['sha256']} is not that of "
        f"{RECORDED_NISSLE.relative_to(REPO)}"
    )


def test_db_a2_holds_the_bakta_proteins_including_sorfs() -> None:
    """§15.1 A2: Bakta proteins with sORFs, no six-frame ORFs, contaminants added by DIA-NN."""
    manifest = json.loads(_text(MANIFEST))
    counts, options = manifest["counts"], manifest["options"]
    expected = _recorded_counts()
    # ``gippa makedb`` drops the 3 CDS Bakta marks as pseudogenes, so the count is just under the
    # annotated one; a database built from a different genome would be far off.
    assert 0.99 * expected["cds"] <= counts["cds"] <= expected["cds"], (counts, expected)
    assert counts["sorf"] == expected["sorf"], (counts, expected)
    assert options["six_frame_orfs"] is False, "A2 is the annotated database (§15.1), no ORFs"
    assert counts.get("orf", 0) == 0, counts
    assert options["contaminants"] is False, (
        "M3-23 Q2 / deliverable 2: db_a2 is built with --no-contaminants; DIA-NN adds "
        "contaminants.fasta as a second --fasta"
    )
    assert counts.get("contaminant", 0) == 0, counts


def test_strain_db_is_a_gzipped_protein_fasta_matching_the_manifest() -> None:
    with gzip.open(STRAIN_DB, "rt") as fh:
        headers = [ln for ln in fh if ln.startswith(">")]
    counts = json.loads(_text(MANIFEST))["counts"]
    assert len(headers) == sum(counts.values()), (len(headers), counts)
    assert not any("*" in h for h in headers), "headers must not carry stop characters"


# ---------------------------------------------------------------- R11 and decision 3


def _tracked(prefix: str) -> list[Path]:
    out = subprocess.run(["git", "-C", str(REPO), "ls-files", prefix],
                         capture_output=True, text=True, check=True).stdout.split()  # fmt: skip
    return [REPO / p for p in out]


def test_no_container_image_is_committed() -> None:
    """Decision 3 and R11: the 532 MB image is never copied into the repository."""
    images = [p for p in _tracked(".") if p.suffix in {".sif", ".simg", ".img"}]
    assert not images, images
    big = [(str(p.relative_to(REPO)), p.stat().st_size) for p in _tracked(".")
           if p.is_file() and p.stat().st_size > R11_MAX_BYTES]  # fmt: skip
    assert not big, f"R11: committed files over 5 MB: {big}"


def test_validation_holds_only_small_files() -> None:
    oversized = {
        str(p.relative_to(REPO)): p.stat().st_size
        for p in _tracked("validation")
        if p.is_file() and p.stat().st_size > MAX_SMALL_BYTES and p != STRAIN_DB
    }
    assert not oversized, (
        f"only db_a2/strain_db.faa.gz may exceed 1 MB under validation/: {oversized}"
    )
    if STRAIN_DB.is_file():
        assert STRAIN_DB.stat().st_size <= R11_MAX_BYTES, STRAIN_DB.stat().st_size


def test_no_script_copies_the_image_into_the_repository() -> None:
    for name, path in SCRIPTS.items():
        code = _code(path)
        for m in re.finditer(r"^\s*(cp|rsync|ln|mv|curl|wget)\b[^\n]*", code, re.M):
            line = m.group(0)
            assert ".sif" not in line, f"{name}: the image is never copied ({line.strip()})"
        assert not re.search(r"\$\{?GIPPA_REPO\}?[^\n]*>\s", code), (
            f"{name}: writes into the repository checkout"
        )


# ---------------------------------------------------------------- the three sbatch scripts


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_script_parses(name: str) -> None:
    """The repo has no shellcheck gate; ``bash -n`` is the acceptance criterion."""
    path = SCRIPTS[name]
    _text(path)
    bash = shutil.which("bash")
    assert bash is not None
    res = subprocess.run([bash, "-n", str(path)], capture_output=True, text=True, timeout=30)
    assert res.returncode == 0, res.stderr
    assert _text(path).splitlines()[0].strip() in ("#!/bin/bash", "#!/usr/bin/env bash")


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_script_passes_shellcheck_when_available(name: str) -> None:
    sc = shutil.which("shellcheck")
    if sc is None:
        pytest.skip("shellcheck is not installed")
    res = subprocess.run([sc, "-S", "error", str(SCRIPTS[name])], capture_output=True, text=True,
                         timeout=120)  # fmt: skip
    assert res.returncode == 0, res.stdout + res.stderr


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_sbatch_house_rules(name: str) -> None:
    """R14 and §3.8: ``gippa-*`` job name, resources, the shared environment, no login-node work."""
    path = SCRIPTS[name]
    job = _directive(path, "job-name")
    assert job and job.startswith("gippa-"), f"{name}: job name {job!r} (R14: gippa-*)"
    assert _directive(path, "partition"), f"{name}: no --partition"
    assert _directive(path, "time"), f"{name}: no --time"
    cpus = _directive(path, "cpus-per-task")
    assert cpus and 1 <= int(cpus) <= 128, f"{name}: --cpus-per-task {cpus}"
    mem = _directive(path, "mem")
    assert mem, f"{name}: no --mem"
    assert _mem_gb(mem) <= 256, f"{name}: --mem {mem}"
    code = _code(path)
    assert "scripts/farmshare_env.sh" in code, f"{name}: does not source scripts/farmshare_env.sh"
    assert re.search(r"^set -euo pipefail$", code, re.M), f"{name}: no `set -euo pipefail`"
    out = _directive(path, "output")
    assert out and ("$GIPPA_ROOT" in out or "/local/logs/" in out), f"{name}: --output {out}"


@pytest.mark.parametrize("name", sorted(SCRIPTS))
def test_outputs_live_outside_the_repository(name: str) -> None:
    code = _code(SCRIPTS[name])
    assert re.search(r"\$\{?GIPPA_(ROOT|DATA_DIR|LOG_DIR)\b", code), (
        f"{name}: writes nothing under $GIPPA_ROOT (R14: large files only there)"
    )
    out_in_repo = re.compile(r"(OUT|OUTPUT|WORK|DEST)[A-Z_]*=\S*(validation/|\$\{?GIPPA_REPO)")
    bad = [ln for ln in code.splitlines() if out_in_repo.search(ln)]
    assert not bad, f"{name}: output directory inside the repository: {bad}"


def test_the_kit_fits_the_r14_cpu_budget() -> None:
    total = sum(int(_directive(p, "cpus-per-task") or 1) for p in SCRIPTS.values())
    total += int(_directive(SCRIPTS["diann_pilot"], "cpus-per-task") or 1)  # A1 and A2 in parallel
    assert total <= 128, f"the pilot kit asks for {total} CPUs at once (R14 cap 128)"
    gres = _directive(SCRIPTS["denovo_pilot"], "gres") or ""
    gpus = int(m.group(1)) if (m := re.search(r"gpu:(?:\w+:)?(\d+)", gres)) else 0
    assert gpus <= 4, f"R14 GPU cap: {gres}"


# ---------------------------------------------------------------- diann_pilot.sbatch (static)


def test_diann_pilot_uses_the_image_from_decision_1() -> None:
    code = _code(SCRIPTS["diann_pilot"])
    m = re.search(r"DIANN_SIF:-(\S+?)\}", code)
    assert m, "diann_pilot.sbatch: DIANN_SIF has no ${DIANN_SIF:-…} default (decision 3)"
    # M10-02 lead decision 1: the pilot runs the verified copy; Khoa's image is the source of
    # diann_image.sbatch
    assert m.group(1).strip("\"'") == "$GIPPA_ROOT/images/diann_2.3.1.sif", m.group(1)
    image_job = _code(REPO / "scripts" / "slurm" / "diann_image.sbatch")
    shared = re.search(r"DIANN_SIF_SHARED:-(\S+?)\}", image_job)
    assert shared and shared.group(1).strip("\"'") == DIANN_SIF, image_job[:200]
    assert DIANN_BIN in code, f"the binary inside the image is {DIANN_BIN}"
    for bind in BINDS:
        assert re.search(rf"(-B|--bind[= ])\s*\S*{re.escape(bind)}\b", code), (
            f"diann_pilot.sbatch: missing bind {bind} (decision 1)"
        )
    assert re.search(r"\bmodule\s+load\s+apptainer", code), "module load apptainer/1.4.4"


def test_diann_pilot_reads_the_cfg_files() -> None:
    code = _code(SCRIPTS["diann_pilot"])
    assert "diann_A" in code or "diann_${" in code or "diann_$" in code, (
        "diann_pilot.sbatch does not read validation/pilot/diann_<arm>.cfg"
    )
    assert "validation/pilot" in code


def test_diann_pilot_names_both_arms() -> None:
    code = _code(SCRIPTS["diann_pilot"])
    assert "A1" in code and "A2" in code


def test_denovo_pilot_is_a_gpu_job_with_the_documented_fallback() -> None:
    """Deliverable 4: Cascadia on the gpu partition, Casanovo/UTI89 fallback documented."""
    path = SCRIPTS["denovo_pilot"]
    assert _directive(path, "partition") == "gpu", _directive(path, "partition")
    assert re.fullmatch(r"gpu(:\w+)?:1", _directive(path, "gres") or ""), _directive(path, "gres")
    assert (_directive(path, "job-name") or "").startswith("gippa-denovo")
    text = _text(path)
    assert "Cascadia" in text and "0.0.7" in text, "the pinned Cascadia version (§15.1)"
    for token in ("Casanovo", "5.2.1", "PXD053727", "UTI89"):
        assert token in text, f"denovo_pilot.sbatch does not document the fallback: {token}"
    assert "diann" not in _code(path).lower(), "the de novo arm never runs DIA-NN"


def test_fetch_script_downloads_outside_the_repository() -> None:
    """Deliverable 5: ``$GIPPA_DATA_DIR/pxd078459/`` plus ``SHA256SUMS``, never under the repo."""
    path = SCRIPTS["fetch_pxd078459"]
    code = _code(path)
    assert re.search(r"\$\{?GIPPA_DATA_DIR\}?/pxd078459", code, re.I), code[:400]
    assert "SHA256SUMS" in code
    assert re.search(r"PXD078459", code, re.I), "the accession of §15.1 is not named"
    assert re.search(r"\b(curl|wget)\b", code), (
        "the fetch script must download with curl or wget (the only tools the tests stub)"
    )
    assert not re.search(r"--no-check-certificate|--insecure|(?<![-\w])-k(?![-\w])", code), (
        "R6: never disable TLS verification"
    )
    assert not re.search(r"http://(?!127\.0\.0\.1)", code), "R6: downloads use https (or ftps)"


def test_only_the_fetch_script_reaches_the_network() -> None:
    """Decision 6: the other two scripts never download anything."""
    for name in ("diann_pilot", "denovo_pilot"):
        code = _code(SCRIPTS[name])
        assert not re.search(r"\b(curl|wget|aria2c|git clone)\b", code), (
            f"{name}.sbatch reaches the network"
        )


def test_pilot_directory_holds_no_mass_spectrometry_data() -> None:
    """R11 and decision 6: no raw data, no DIA-NN report is committed in M3."""
    suspicious = [
        str(p.relative_to(REPO))
        for p in _tracked("validation")
        if p.suffix.lower() in {".raw", ".d", ".mzml", ".parquet", ".speclib", ".tdf"}
    ]
    assert not suspicious, suspicious
    assert PILOT.is_dir()
