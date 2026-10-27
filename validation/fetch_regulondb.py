"""Download the RegulonDB transcription units and terminators (GIPPA_SPEC.md §15.3, §15.4).

Usage: ``python validation/fetch_regulondb.py`` (no arguments; network, run by hand).

Writes the official RegulonDB datasets ``TUSet``, ``TerminatorSet`` and ``GeneProductSet`` into
``$GIPPA_DATA_DIR/regulondb/<release>/`` together with ``provenance.json`` (URL, release, dataset
version, byte count, SHA-256 and every retrieval date). Nothing is committed (R11) and no mirror
is used (R6): the bytes come from RegulonDB's own GraphQL endpoint.

TLS. Since 2026-10-14 ``regulondb.ccg.unam.mx`` serves a certificate issued by "GlobalSign GCC R46
OV TLS CA 2025" but **omits that intermediate** from the handshake, so a normal client fails with
``unable to get local issuer certificate``. This script completes the chain itself: it downloads
the intermediate from GlobalSign's own CA endpoint, refuses to continue unless its SHA-256 equals
``INTERMEDIATE_SHA256``, checks with ``openssl verify`` against Python's own default trust store
(``ssl.get_default_verify_paths()``, the store the handshake uses) that it chains to "GlobalSign
Root R46", and only then adds it to the CA bundle used for RegulonDB, with
``VERIFY_X509_PARTIAL_CHAIN`` cleared so the handshake must still reach that root. Verification is
never disabled, and the extra certificate is scoped to this script's requests. When the
intermediate expires (2029-06-23) or RegulonDB publishes a complete chain, the failure names the
pin so the cause is obvious; see docs/backlog.md.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import shutil
import ssl
import subprocess
import sys
import tempfile
import time
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

GRAPHQL_URL = "https://regulondb.ccg.unam.mx/graphql"
HOST = "regulondb.ccg.unam.mx"

# The intermediate RegulonDB omits, and the root it must chain to (R6: GlobalSign's own endpoint).
INTERMEDIATE_URL = "http://secure.globalsign.com/cacert/gsgccr46ovtlsca2025.crt"
INTERMEDIATE_SHA256 = "d160e2de4e56cb10b66c1cb0adcb79cf93c78dcd9db30c2018220262c04063f9"
INTERMEDIATE_CN = "GlobalSign GCC R46 OV TLS CA 2025"
INTERMEDIATE_NOT_AFTER = "2029-06-23"
ROOT_CN = "GlobalSign Root R46"

# §15.4 needs the transcription units and the terminators; the gene set carries the coordinates
# that map RegulonDB gene names onto the recorded K-12 Bakta locus tags (regulondb_truth.py).
DATASETS = ("TUSet", "TerminatorSet", "GeneProductSet")

ATTEMPTS = 3
TIMEOUT_S = 300
OPENSSL_TIMEOUT_S = 60
USER_AGENT = "gippa-validation (FarmShare)"

_FILE_QUERY = """query GetDataOfFile($fileName: String) {
  getDataOfFile(fileName: $fileName) {
    fileName
    fileFormat
    version
    rdbVersion
    creationDate
    license
    citation
    columnsDetails
    content
  }
}"""

_INFO_QUERY = """query GetDatabaseInfo {
  getDatabaseInfo {
    regulonDBVersion
    ecocycVersion
    releaseDate
    genomeVersion
  }
}"""


class FetchError(RuntimeError):
    """The download could not be completed; the message names the cause."""


@dataclass(frozen=True, slots=True)
class Dataset:
    """One RegulonDB download: its bytes and the metadata the endpoint reports."""

    name: str
    content: bytes
    file_format: str | None
    version: str | None
    rdb_version: str | None
    creation_date: str | None
    license: str | None
    citation: str | None
    columns_details: str | None

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.content).hexdigest()


# --------------------------------------------------------------------- TLS


def sha256_hex(data: bytes) -> str:
    """SHA-256 of ``data`` as lower-case hex."""
    return hashlib.sha256(data).hexdigest()


def _openssl() -> str:
    """The ``openssl`` binary that performs the chain check; absent means refuse, not downgrade."""
    found = shutil.which("openssl")
    if found is None:
        raise FetchError(
            "fetch_regulondb: openssl is not on PATH, so the pinned intermediate cannot be "
            f"checked against {ROOT_CN!r}; refusing to continue"
        )
    return found


def _run_openssl(*args: str) -> subprocess.CompletedProcess[str]:
    """Run ``openssl``; a hang becomes a :class:`FetchError`, never a traceback."""
    argv = [_openssl(), *args]
    try:
        return subprocess.run(  # fixed argv, on a file we just wrote
            argv, capture_output=True, text=True, check=False, timeout=OPENSSL_TIMEOUT_S
        )
    except subprocess.TimeoutExpired:
        raise FetchError(
            f"fetch_regulondb: `openssl {args[0]}` timed out after {OPENSSL_TIMEOUT_S} s while "
            f"checking the pinned intermediate {INTERMEDIATE_SHA256}; refusing to continue"
        ) from None


def _store_args() -> list[str]:
    """``-CAfile`` / ``-CApath`` for Python's default trust store, so both checks agree."""
    paths = ssl.get_default_verify_paths()
    args: list[str] = []
    if paths.cafile:
        args += ["-CAfile", paths.cafile]
    if paths.capath:
        args += ["-CApath", paths.capath]
    if not args:
        raise FetchError(
            "fetch_regulondb: Python has no default trust store (ssl.get_default_verify_paths()), "
            f"so the pinned intermediate cannot be checked against {ROOT_CN!r}; refusing"
        )
    return args


def verify_intermediate(der: bytes) -> str:
    """Pin ``der`` by SHA-256, check it chains to :data:`ROOT_CN`, and return it as PEM.

    The chain check is ``openssl verify`` against Python's default trust store, named explicitly
    with ``-CAfile`` / ``-CApath`` from ``ssl.get_default_verify_paths()`` (a bundle file, a hashed
    directory or both), so it trusts exactly what the TLS handshake will trust; the issuer is then
    read back and must be :data:`ROOT_CN`, so the anchor is the root we expect.

    Raises :class:`FetchError` naming :data:`INTERMEDIATE_SHA256` when the pin does not match or
    the chain does not verify, so an expiry or a re-issue is obvious rather than mysterious.
    """
    got = sha256_hex(der)
    if got != INTERMEDIATE_SHA256:
        raise FetchError(
            f"fetch_regulondb: the GlobalSign intermediate changed: sha256 {got} != the pin "
            f"{INTERMEDIATE_SHA256} ({INTERMEDIATE_CN}, pinned valid to {INTERMEDIATE_NOT_AFTER}). "
            "Re-pin it only after checking it against GlobalSign; never disable verification."
        )
    pem = ssl.DER_cert_to_PEM_cert(der)
    with tempfile.TemporaryDirectory() as tmp:
        inter = Path(tmp) / "intermediate.pem"
        inter.write_text(pem, encoding="ascii")
        verified = _run_openssl("verify", *_store_args(), str(inter))
        issuer = _run_openssl("x509", "-in", str(inter), "-noout", "-issuer")
    if verified.returncode != 0:
        raise FetchError(
            f"fetch_regulondb: the pinned intermediate {INTERMEDIATE_SHA256} does not chain to a "
            f"trusted root ({ROOT_CN!r} should be in Python's default trust store): "
            f"{verified.stdout.strip()} {verified.stderr.strip()}"
        )
    if ROOT_CN not in issuer.stdout:
        raise FetchError(
            f"fetch_regulondb: the pinned intermediate is anchored to {issuer.stdout.strip()!r}, "
            f"not to {ROOT_CN!r}; refusing to continue"
        )
    return pem


def _get(url: str, timeout: int = TIMEOUT_S) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # fixed http(s) URL
        data: bytes = resp.read()
    return data


def regulondb_context(fetch: Callable[[str], bytes] = _get) -> ssl.SSLContext:
    """Build the TLS context for RegulonDB: Python's default roots plus the pinned intermediate.

    ``fetch`` is the one network call (injected so tests stay offline). Verification stays on:
    the intermediate is only added after the SHA-256 pin and the ``openssl verify`` check pass.
    """
    der = fetch(INTERMEDIATE_URL)
    if der[:1] == b"-":  # the pin covers the DER GlobalSign serves, not a PEM re-encoding
        raise FetchError(
            f"fetch_regulondb: {INTERMEDIATE_URL} returned PEM, not the DER the pin covers"
        )
    pem = verify_intermediate(der)
    ctx = ssl.create_default_context()
    # Python >= 3.13 sets VERIFY_X509_PARTIAL_CHAIN by default, which would let the added
    # intermediate act as a trust anchor on its own; clear it so the chain must reach the root.
    ctx.verify_flags &= ~ssl.VERIFY_X509_PARTIAL_CHAIN
    ctx.load_verify_locations(cadata=pem)
    return ctx


# ----------------------------------------------------------------- GraphQL


def graphql(query: str, variables: dict[str, str], context: ssl.SSLContext) -> dict[str, Any]:
    """POST ``query`` to RegulonDB and return the ``data`` object; raise on GraphQL errors."""
    body = json.dumps({"query": query, "variables": variables}).encode()
    req = urllib.request.Request(
        GRAPHQL_URL,
        data=body,
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
        method="POST",
    )
    with urllib.request.urlopen(  # GRAPHQL_URL is a constant https URL
        req, timeout=TIMEOUT_S, context=context
    ) as resp:
        payload = json.loads(resp.read().decode("utf-8"))
    if payload.get("errors"):
        raise FetchError(f"fetch_regulondb: GraphQL error from {HOST}: {payload['errors']!r:.400}")
    data = payload.get("data")
    if not isinstance(data, dict):
        raise FetchError(f"fetch_regulondb: GraphQL response without data from {HOST}")
    return data


def _one(value: Any, field: str) -> dict[str, Any]:
    """RegulonDB returns either an object or a one-element list; take the first object."""
    if isinstance(value, list):
        value = value[0] if value else None
    if not isinstance(value, dict):
        raise FetchError(f"fetch_regulondb: {field} is missing from the GraphQL response")
    return value


def latest_release(info: Any) -> dict[str, Any]:
    """The newest release row of ``getDatabaseInfo`` (the list is newest first, but sort anyway)."""
    rows = info if isinstance(info, list) else [info]
    clean = [r for r in rows if isinstance(r, dict) and r.get("regulonDBVersion")]
    if not clean:
        raise FetchError("fetch_regulondb: getDatabaseInfo returned no release")
    return max(clean, key=lambda r: str(r.get("releaseDate") or ""))


def fetch_dataset(name: str, context: ssl.SSLContext) -> Dataset:
    """Download one RegulonDB dataset by name (``TUSet``, ``TerminatorSet``, ...)."""
    data = graphql(_FILE_QUERY, {"fileName": name}, context)
    rec = _one(data.get("getDataOfFile"), f"getDataOfFile({name})")
    content = rec.get("content")
    if not isinstance(content, str) or not content.strip():
        raise FetchError(f"fetch_regulondb: {name} came back empty")
    return Dataset(
        name=name,
        content=content.encode("utf-8"),
        file_format=rec.get("fileFormat"),
        version=rec.get("version"),
        rdb_version=rec.get("rdbVersion"),
        creation_date=rec.get("creationDate"),
        license=rec.get("license"),
        citation=rec.get("citation"),
        columns_details=rec.get("columnsDetails"),
    )


# -------------------------------------------------------------- the writer


def release_dir(root: Path, release: str) -> Path:
    """``<root>/regulondb/<release>``; the release string is sanitised for a path component."""
    safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in release) or "unknown"
    return root / "regulondb" / safe


def write_atomically(path: Path, data: bytes) -> None:
    """Write ``data`` to ``path`` through a ``.part`` sibling, so readers never see a half file."""
    part = path.with_suffix(path.suffix + ".part")
    part.write_bytes(data)
    part.replace(path)


def provenance(
    datasets: list[Dataset],
    release: dict[str, Any],
    now: dt.datetime,
    previous: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build the provenance record; ``previous`` keeps earlier retrieval dates (example 6)."""
    today = now.date().isoformat()
    retrieved = list((previous or {}).get("retrieved_dates") or [])
    if today not in retrieved:
        retrieved.append(today)
    return {
        "source": "RegulonDB",
        "url": GRAPHQL_URL,
        "host": HOST,
        "regulondb_version": release.get("regulonDBVersion"),
        "release_date": release.get("releaseDate"),
        "ecocyc_version": release.get("ecocycVersion"),
        "genome_version": release.get("genomeVersion"),
        "retrieved": today,
        "retrieved_at": now.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "retrieved_dates": retrieved,
        "tls": {
            "note": f"{HOST} omits its TLS intermediate; GIPPA completes the chain",
            "intermediate_url": INTERMEDIATE_URL,
            "intermediate_sha256": INTERMEDIATE_SHA256,
            "intermediate_cn": INTERMEDIATE_CN,
            "intermediate_not_after": INTERMEDIATE_NOT_AFTER,
            "anchored_to": ROOT_CN,
        },
        "files": [
            {
                "dataset": d.name,
                "file": f"{d.name}.txt",
                "bytes": len(d.content),
                "sha256": d.sha256,
                "file_format": d.file_format,
                "dataset_version": d.version,
                "rdb_version": d.rdb_version,
                "creation_date": d.creation_date,
                "license": d.license,
                "citation": d.citation,
                "columns": d.columns_details,
            }
            for d in datasets
        ],
    }


def _read_previous(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        loaded = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return loaded if isinstance(loaded, dict) else None


def main() -> int:
    data_dir = os.environ.get("GIPPA_DATA_DIR")
    if not data_dir:
        print("fetch_regulondb: GIPPA_DATA_DIR is not set", file=sys.stderr)
        return 1
    try:
        context = regulondb_context()
    except (FetchError, OSError) as exc:
        print(f"fetch_regulondb: cannot reach {HOST} over a verified chain: {exc}", file=sys.stderr)
        return 1

    release: dict[str, Any] = {}
    datasets: list[Dataset] = []
    for attempt in range(1, ATTEMPTS + 1):
        try:
            release = latest_release(graphql(_INFO_QUERY, {}, context).get("getDatabaseInfo"))
            datasets = [fetch_dataset(name, context) for name in DATASETS]
            break
        except (FetchError, OSError, json.JSONDecodeError) as exc:
            print(f"fetch_regulondb: attempt {attempt}: {exc}", file=sys.stderr)
            if attempt == ATTEMPTS:
                return 1
            time.sleep(30 * attempt)

    out = release_dir(Path(data_dir), str(release.get("regulonDBVersion") or "unknown"))
    out.mkdir(parents=True, exist_ok=True)
    for d in datasets:
        write_atomically(out / f"{d.name}.txt", d.content)
    sidecar = out / "provenance.json"
    meta = provenance(datasets, release, dt.datetime.now(dt.UTC), _read_previous(sidecar))
    write_atomically(sidecar, (json.dumps(meta, indent=2) + "\n").encode("utf-8"))
    sizes = ", ".join(f"{d.name} {len(d.content)} B" for d in datasets)
    print(f"fetch_regulondb: RegulonDB {meta['regulondb_version']} -> {out} ({sizes})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
