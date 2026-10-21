"""Download the RegulonDB transcription units and terminators (GIAP_SPEC.md §15.3, §15.4).

Usage: ``python validation/fetch_regulondb.py`` (no arguments; network, run by hand).

Writes the official RegulonDB datasets ``TUSet``, ``TerminatorSet`` and ``GeneProductSet`` into
``$GIAP_DATA_DIR/regulondb/<release>/`` together with ``provenance.json`` (URL, release, dataset
version, byte count, SHA-256 and every retrieval date). Nothing is committed (R11) and no mirror
is used (R6): the bytes come from RegulonDB's own GraphQL endpoint.

TLS. Since 2026-10-14 ``regulondb.ccg.unam.mx`` serves a certificate issued by "GlobalSign GCC R46
OV TLS CA 2025" but **omits that intermediate** from the handshake, so a normal client fails with
``unable to get local issuer certificate``. This script completes the chain itself: it downloads
the intermediate from GlobalSign's own CA endpoint, refuses to continue unless its SHA-256 equals
``INTERMEDIATE_SHA256``, checks with ``openssl verify`` that it is signed by "GlobalSign Root R46"
(already in the system trust store), and only then adds it to the CA bundle used for RegulonDB.
Verification is never disabled, and the extra anchor is scoped to this script's requests. When the
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
from collections.abc import Callable, Mapping
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
USER_AGENT = "giap-validation (FarmShare)"

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


def _common_names(cert: Mapping[str, Any]) -> set[str]:
    """The subject common names of a certificate as ``ssl`` decodes it."""
    out: set[str] = set()
    subject = cert.get("subject") or ()
    for rdn in subject:
        for attribute in rdn:
            if isinstance(attribute, tuple) and len(attribute) == 2:
                key, value = attribute
                if key == "commonName":
                    out.add(str(value))
    return out


def trusted_root_pem(common_name: str) -> str:
    """Return the PEM of the trusted root whose subject CN is ``common_name``.

    Reads the system trust store only (``ssl`` default certificates); raises :class:`FetchError`
    when the root is absent, because completing the chain would then rest on nothing.
    """
    ctx = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    ctx.load_default_certs(ssl.Purpose.SERVER_AUTH)
    decoded = ctx.get_ca_certs()
    binary = ctx.get_ca_certs(binary_form=True)
    for cert, der in zip(decoded, binary, strict=True):
        if common_name in _common_names(cert):
            return ssl.DER_cert_to_PEM_cert(der)
    raise FetchError(
        f"fetch_regulondb: {common_name!r} is not in the system trust store, so the "
        f"{INTERMEDIATE_CN!r} intermediate cannot be anchored; refusing to continue"
    )


def verify_intermediate(der: bytes, root_pem: str) -> str:
    """Pin ``der`` by SHA-256, check it is signed by ``root_pem``, and return it as PEM.

    Raises :class:`FetchError` naming :data:`INTERMEDIATE_SHA256` when the pin does not match or
    the signature does not verify, so an expiry or a re-issue is obvious rather than mysterious.
    """
    got = sha256_hex(der)
    if got != INTERMEDIATE_SHA256:
        raise FetchError(
            f"fetch_regulondb: the GlobalSign intermediate changed: sha256 {got} != the pin "
            f"{INTERMEDIATE_SHA256} ({INTERMEDIATE_CN}, pinned valid to {INTERMEDIATE_NOT_AFTER}). "
            "Re-pin it only after checking it against GlobalSign; never disable verification."
        )
    openssl = shutil.which("openssl")
    if openssl is None:
        raise FetchError(
            "fetch_regulondb: openssl is not on PATH, so the pinned intermediate cannot be "
            f"checked against {ROOT_CN!r}; refusing to continue"
        )
    pem = ssl.DER_cert_to_PEM_cert(der)
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "root.pem"
        inter = Path(tmp) / "intermediate.pem"
        root.write_text(root_pem, encoding="ascii")
        inter.write_text(pem, encoding="ascii")
        proc = subprocess.run(  # fixed argv, on files we just wrote
            [openssl, "verify", "-CAfile", str(root), str(inter)],
            capture_output=True,
            text=True,
            check=False,
            timeout=60,
        )
    if proc.returncode != 0:
        raise FetchError(
            f"fetch_regulondb: the pinned intermediate {INTERMEDIATE_SHA256} does not verify "
            f"against {ROOT_CN!r}: {proc.stdout.strip()} {proc.stderr.strip()}"
        )
    return pem


def _get(url: str, timeout: int = TIMEOUT_S) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # fixed http(s) URL
        data: bytes = resp.read()
    return data


def regulondb_context(fetch: Callable[[str], bytes] = _get) -> ssl.SSLContext:
    """Build the TLS context for RegulonDB: the system roots plus the pinned intermediate.

    ``fetch`` is the one network call (injected so tests stay offline). Verification stays on:
    the intermediate is only added after the SHA-256 pin and the ``openssl verify`` check pass.
    """
    der = fetch(INTERMEDIATE_URL)
    if der[:1] == b"-":  # the pin covers the DER GlobalSign serves, not a PEM re-encoding
        raise FetchError(
            f"fetch_regulondb: {INTERMEDIATE_URL} returned PEM, not the DER the pin covers"
        )
    pem = verify_intermediate(der, trusted_root_pem(ROOT_CN))
    ctx = ssl.create_default_context()
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
            "note": f"{HOST} omits its TLS intermediate; GIAP completes the chain",
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
    data_dir = os.environ.get("GIAP_DATA_DIR")
    if not data_dir:
        print("fetch_regulondb: GIAP_DATA_DIR is not set", file=sys.stderr)
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
