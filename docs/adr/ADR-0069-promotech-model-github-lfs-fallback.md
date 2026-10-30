# ADR-0069: A second official source for the Promotech RF-HOT model (GitHub LFS fallback)
- Status: accepted-by-khoa (2026-10-30, in session: "Accept as fallback")
- Date: 2026-10-30
- Spec sections: §3.2 (source table), §0.2 R6; ADR-0014; issue #508
- Needs Khoa: yes (adds a download source; R6 allows only the sources in §3.2)

## Context
`gippa setup` downloads `RF-HOT.zip` only from `www.cs.mun.ca/~lourdes/public/PromoTech_models/`
(§3.2). On 2026-10-30 that host timed out from Stanford Sherlock (compute and login nodes) and from
FarmShare, while every other source worked (#508). Promotech is required by default, so a fresh
setup cannot complete while the host is down.

The Promotech authors' repository `github.com/BioinformaticsLabAtMUN/Promotech` added the model
as a Git LFS object in commit `bdf2f3ce06ffd0e342607eed4bfa32e160816143` (2020-08-06, "Add RF-HOT
model") and removed it from the tip in `f6c4f922` (2021-04-28). The LFS pointer at `bdf2f3ce`
records `oid sha256:fcd96ffe5a43b579af965839ee2bb896de94de3601978d88bcabd32c4fb5093f`, size
1006105109: the same SHA-256 as the ADR-0014 pin. On 2026-10-30
`https://media.githubusercontent.com/media/BioinformaticsLabAtMUN/PromoTech/bdf2f3ce06ffd0e342607eed4bfa32e160816143/models/RF-HOT.zip`
answered directly (no redirect), and its first and last KiB equal the pinned copy. No other
official copy exists (Zenodo 4737459 lacks the model; the paper's supplements and bioconda have none).

## Decision
1. That GitHub LFS URL is a second official source of `promotech_models`, tried only when
   `www.cs.mun.ca` raises a source error (unreachable, HTTP error). cs.mun.ca stays first: the
   authors advertise it, and LFS bandwidth is billed to the authors' organisation.
2. Both sources are checked against the same ADR-0014 SHA-256. A checksum mismatch from either is
   final and never falls through to the other source. No other item gains a fallback.
3. The record and the "downloaded from" message name the URL that served the file.
4. When every source fails, the error says how to install the model by hand: put `RF-HOT.zip` and
   its unzipped `RF-HOT.model` in `<db-dir>/promotech/`, then rerun `gippa setup`.

## Evidence (before → after: metrics, tests, simulations)
Before: one source, setup fails while cs.mun.ca is down. After: unit tests pin the fallback order,
no fallback on a checksum mismatch, and the recorded URL; the LFS head/tail bytes match the pin.

## Consequences
- §3.2's source table gains the row `media.githubusercontent.com` (Promotech RF-HOT, fallback).
- If the authors purge the LFS object or its quota runs out, the fallback fails like cs.mun.ca and
  the manual route remains.
