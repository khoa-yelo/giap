"""Download the UniProt E. coli K-12 proteome UP000000625 (GIAP_SPEC.md §11.4, §11.6, §15.3).

Usage: ``python validation/fetch_uniprot_k12.py`` (no arguments; network, weekly job only).

Writes ``$GIAP_DATA_DIR/uniprot/UP000000625.fasta`` from the official UniProt REST stream (R6)
and ``UP000000625.provenance.json`` next to it with the URL, UniProt release, entry count,
SHA-256 and retrieval date (YYYY-MM-DD). Never committed (R11). Exit 1 when the download is
incomplete (entry count differs from ``X-Total-Results``) or the source is unreachable.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import sys
import time
import urllib.request
from pathlib import Path

URL = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000000625%29"
NAME = "UP000000625"
ATTEMPTS = 3
TIMEOUT_S = 600


def _download(url: str, dst: Path) -> dict[str, str]:
    """Stream ``url`` into ``dst``; returns the response headers (lower-case keys)."""
    req = urllib.request.Request(url, headers={"User-Agent": "giap-validation (FarmShare)"})
    with urllib.request.urlopen(req, timeout=TIMEOUT_S) as resp, dst.open("wb") as fh:
        while chunk := resp.read(1 << 20):
            fh.write(chunk)
        return {k.lower(): v for k, v in resp.headers.items()}


def main() -> int:
    data = os.environ.get("GIAP_DATA_DIR")
    if not data:
        print("fetch_uniprot_k12: GIAP_DATA_DIR is not set", file=sys.stderr)
        return 1
    out = Path(data) / "uniprot"
    out.mkdir(parents=True, exist_ok=True)
    fasta, part = out / f"{NAME}.fasta", out / f"{NAME}.part"
    headers: dict[str, str] = {}
    for attempt in range(1, ATTEMPTS + 1):
        try:
            headers = _download(URL, part)
            break
        except OSError as exc:
            print(f"fetch_uniprot_k12: attempt {attempt}: {exc}", file=sys.stderr)
            if attempt == ATTEMPTS:
                part.unlink(missing_ok=True)
                return 1
            time.sleep(30 * attempt)
    raw = part.read_bytes()
    n = raw.count(b"\n>") + int(raw.startswith(b">"))
    total = headers.get("x-total-results")
    if n == 0 or (total is not None and int(total) != n):
        print(f"fetch_uniprot_k12: incomplete download ({n} of {total} entries)", file=sys.stderr)
        part.unlink(missing_ok=True)
        return 1
    part.replace(fasta)
    now = dt.datetime.now(dt.UTC)
    meta = {
        "file": fasta.name,
        "url": URL,
        "uniprot_release": headers.get("x-uniprot-release"),
        "entries": n,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "retrieved": now.date().isoformat(),
        "retrieved_at": now.strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    sidecar = out / f"{NAME}.provenance.json"
    sidecar.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"fetch_uniprot_k12: {n} entries, release {meta['uniprot_release']} -> {fasta}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
