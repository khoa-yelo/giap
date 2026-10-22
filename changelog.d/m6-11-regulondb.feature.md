`validation/fetch_regulondb.py` downloads the RegulonDB transcription units, terminators and gene
coordinates into `$GIAP_DATA_DIR` with SHA-256 provenance, completing the TLS chain the server omits
with a pinned GlobalSign intermediate (verification is never disabled), and
`validation/regulondb_truth.py` builds the K-12 pair and terminator truth sets of GIAP_SPEC.md §15.4.
