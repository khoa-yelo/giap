<!-- m06b -->
Operon and terminator calls can be checked against RegulonDB: `validation/fetch_regulondb.py` (in the source repository)
downloads the RegulonDB transcription units, terminators and gene coordinates with SHA-256
provenance (TLS verification stays on), and `validation/regulondb_truth.py` builds the E. coli
K-12 pair and terminator truth sets used for validation.
