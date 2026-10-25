#!/usr/bin/env python3
"""Resolve an image tag to its manifest digest (registry v2 API; GIAP_SPEC.md §3.1, §12).

Usage: python scripts/resolve_digest.py <image:tag> [...]
Prints one `<image>@sha256:<64 hex>` per argument (bare, no scheme; ADR-0059); a `docker://`
prefix on the input is accepted. Only the registries listed in §3.2
(Docker Hub, ghcr.io, quay.io) are used; tokens are anonymous pull tokens.
"""

from __future__ import annotations

import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

REGISTRIES = {"docker.io": "registry-1.docker.io", "ghcr.io": "ghcr.io", "quay.io": "quay.io"}
ACCEPT = ", ".join(
    [
        "application/vnd.oci.image.index.v1+json",
        "application/vnd.docker.distribution.manifest.list.v2+json",
        "application/vnd.oci.image.manifest.v1+json",
        "application/vnd.docker.distribution.manifest.v2+json",
    ]
)
DIGEST = re.compile(r"^sha256:[0-9a-f]{64}$")


def split_ref(ref: str) -> tuple[str, str, str]:
    """`quay.io/biocontainers/x:1` → (registry, repository, tag); Docker Hub is the default."""
    ref = ref.removeprefix("docker://")
    name, _, tag = ref.rpartition(":")
    if not name or "/" in tag:
        name, tag = ref, "latest"
    first, _, rest = name.partition("/")
    if rest and ("." in first or ":" in first):
        registry, repo = first, rest
    else:
        registry, repo = "docker.io", name
    if registry == "docker.io" and "/" not in repo:
        repo = f"library/{repo}"
    if registry not in REGISTRIES:
        raise ValueError(f"registry {registry!r} is not an allowed source (§3.2)")
    return registry, repo, tag


def _token(challenge: str) -> str:
    params = dict(re.findall(r'(\w+)="([^"]*)"', challenge))
    url = params.pop("realm") + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as resp:
        body = json.load(resp)
    return str(body.get("token") or body.get("access_token"))


def resolve(ref: str) -> str:
    registry, repo, tag = split_ref(ref)
    url = f"https://{REGISTRIES[registry]}/v2/{repo}/manifests/{tag}"
    headers = {"Accept": ACCEPT}
    for _attempt in range(2):
        req = urllib.request.Request(url, headers=headers, method="HEAD")
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                digest = resp.headers.get("Docker-Content-Digest", "")
                break
        except urllib.error.HTTPError as err:
            challenge = err.headers.get("WWW-Authenticate", "")
            if err.code != 401 or not challenge.startswith("Bearer") or "Authorization" in headers:
                raise
            headers["Authorization"] = f"Bearer {_token(challenge)}"
    else:  # pragma: no cover - loop always breaks or raises
        raise RuntimeError(url)
    if not DIGEST.match(digest):
        raise RuntimeError(f"{ref}: no digest in registry response")
    image = repo.removeprefix("library/") if registry == "docker.io" else f"{registry}/{repo}"
    return f"{image}@{digest}"


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__, file=sys.stderr)
        return 2
    for ref in argv:
        print(resolve(ref))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
