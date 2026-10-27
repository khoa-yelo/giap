# ADR-0038: GSEA runs R fgsea as an external tool
- Status: accepted-by-khoa (2026-10-22)
- Date: 2026-10-22
- Spec sections: §7.7.4, §7.7.5, §11.5, G15, R5
- Needs Khoa: no (Khoa decided the dependency, §2 decision G15 and item 5 on 2026-10-22)

## Context
G15 and §7.7.4 make the abundance GSEA a Python and numba port of fgsea 1.36.2. R fgsea was only a fallback,
used when the port misses the §11.5 golden tolerances. Port part I (`fgsea_simple`, M8-02) is merged. Part II,
the multilevel stage (M8-11, PR #193), was written without the fgsea source: the R6 question about
bioconductor.org is still open. So its sign correction and p-value indexing are unverified against fgsea.
On 2026-10-22 Khoa decided that the pipeline calls R fgsea instead of a port.

## Decision
1. **Engine.** `enrichment_gsea.tsv` comes from `fgsea::fgseaMultilevel` 1.36.2 run in R 4.5.
   - It runs as a separate process in its own environment, `workflow/envs/fgsea.yml`: conda-forge and
     bioconda only, `r-base`, `bioconductor-fgsea=1.36.2` and `r-data.table`, with an explicit lock in
     `workflow/envs/locks/`.
   - It needs no Bioconductor data packages, so nothing is downloaded from bioconductor.org.
   - Apptainer profiles use the BioContainers image of the same version, pinned by digest.
2. **What stays in Python.** `gippa.function.gsea` keeps its existing jobs and hands each collection to R as two
   TSVs, the TU scores and the sets in long form:
   - the TU collapse, the set construction and the size limits (`function.gsea_size`);
   - the redundancy re-test sets;
   - BH within collections, the leading edge mapped back to genes, flags and the table writer.
3. **The R script.** It lives in `workflow/bin/gippa_fgsea.R` and never in `src/` (R5).
   - It calls `fgseaMultilevel(sampleSize, eps, nPermSimple, minSize, maxSize, scoreType = "std",
     nproc = 1)` with the `function.gsea` defaults, which are unchanged: 101, 0, 1000.
   - It writes ES, NES, pval, log2err, size and leadingEdge, plus `R.version.string` and
     `packageVersion("fgsea")` for provenance.
   - The seed is one integer from `gippa.rng.for_step("function.gsea", seed)` passed to `set.seed`. With
     `nproc = 1` (SerialParam), a rerun is byte-identical.
4. **Interfaces.** `METHOD = "fgsea_multilevel"`. The table columns, flags and schema are unchanged. The
   engine is still an injectable callable, so unit tests use fakes or recorded outputs (§11.4) and never R.
   A slow integration test and an nf-test run the real environment, on Slurm only.
5. **ORA.** The detection ORA stays on
   `scipy.stats.nchypergeom_wallenius`. It is analytical and exact, not a port of an approximation, and it is
   still checked against goseq by the §11.5 golden data (#192). Calling goseq in the pipeline would pull in
   `bioconductor-genelendatabase` and `bioconductor-go.db`, whose install scripts download from
   bioconductor.org. The elim GO redundancy also needs many dependent re-tests, so each one would cost an R
   round trip.
6. **The port.** PR #193 (multilevel port) is not merged and stays open as a draft for reference. The merged
   `fgsea_simple` stays until the R engine lands; the hardening pass then removes `gippa.stats.fgsea` and its
   kernels, unless Khoa wants a no-R fallback.
7. **Golden data.** The §11.5 fgsea golden outputs (#192) become a contract test of the wrapper. On the
   frozen inputs with `GOLDEN_SEED`, `gippa_fgsea.R` must reproduce `make_golden_r.R`'s fgsea table exactly.
   The goseq golden data keep validating the scipy ORA.

## Evidence (before → after)
- Before: the multilevel p-values come from an unverified reimplementation (#193 notes the missing source
  check).
- After: by construction they are fgsea's own numbers; the remaining risk is the wrapper, which the contract
  test in item 7 covers.
- Licence: fgsea is MIT and R is GPL-2. Both run only as separate processes, as R5 requires;
  `check_imports.py` and `check_no_rscript.py` still hold.

## Consequences
- A new tool environment (about 1 GB) and an R start-up of a few seconds per `function` step. Runtime is
  otherwise the same order as fgsea itself.
- `function` gains an external-tool dependency. `gippa check`/preflight must report a missing fgsea
  environment the same way as other tools do.
- G15 and §7.7.4 are superseded on the engine. The spec text stays read-only; this ADR overrides it.
- The M8-11 multilevel calibration tests (#193) do not land. GSEA calibration (null uniformity, the
  operon-packing attack null) runs against the R engine on Slurm instead.
