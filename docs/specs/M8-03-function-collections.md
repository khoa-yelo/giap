# M8-03: Function resources → gene-set collections (GO, GO slim, KEGG, COG) and coverage

Milestone M8 · track R · **pass 1** · depends on: nothing new (merged M2-09 setup items, M4-15
gene-table function columns) · estimate ≈ 300 code lines (+ ≈ 380 test lines) · not statistical

## Spec sections implemented

- §7.7.1 "Source. Bakta full-database `db_xrefs` are the only annotation source (decision 10):
  `COG:` IDs (and category letters if present), `GO:`, `EC:` and `KEGG:` KOs. Genes without a close
  match stay unannotated. With the light DB, Module 7 runs on whatever cross-references exist,
  reports low coverage, and the manifest notes `bakta_db = light`."
- §7.7.1 "Resources … COG 2020 `cog-20.def.tab` (COG → category) and `fun-20.tab` (category names);
  `go-basic.obo` and `goslim_generic.obo`; KEGG KO → pathway links and pathway names."
- §7.7.1 "GO propagation. Terms are propagated to their ancestors with goatools."
- §7.7.1 "Coverage. The share of eligible genes with ≥ 1 term, per ontology, reported on every
  enrichment row and in `category_summary.tsv`."
- §7.7.4 collections "(COG, KEGG, GO BP, GO MF, GO CC, GO slim)"; §5.2 `function.go_slim: generic`;
  §7.9.1 / ADR-0017 (KEGG unversioned, recorded).

## Public interface

`src/giap/function/collections.py` (pure; readers in `src/giap/io/function_resources.py`)
```python
Collection = Literal["cog", "kegg", "go_bp", "go_mf", "go_cc", "go_slim"]

@dataclass(frozen=True, slots=True)
class GeneSets:
    collection: Collection
    sets: Mapping[str, frozenset[str]]       # term -> locus tags (eligible genes only)
    names: Mapping[str, str]                 # term -> human name
    parents: Mapping[str, frozenset[str]]    # GO only: term -> direct parents (is_a, part_of)
    coverage: float                          # share of eligible genes with >= 1 term here
    n_unknown_ids: int                       # annotation IDs not found in the resource

@dataclass(frozen=True, slots=True)
class Resources:                             # every field None when not installed
    go: GoDag | None; slim: GoDag | None
    kegg_links: pd.DataFrame | None; kegg_names: pd.DataFrame | None
    cog_names: pd.DataFrame | None

def read_resources(paths: Mapping[str, Path | None]) -> Resources: ...      # io/
def build_collections(genes: pd.DataFrame, eligible: B8, res: Resources
                      ) -> tuple[list[GeneSets], list[str]]: ...         # + run flags
```

1. **Gene annotations come from the gene table**, not from `xrefs.parquet`: `cog_category`
   (letters from COG 2020 via M4-15), `go_terms`, `kegg_ko`. One source of truth; the §6.2 row's
   "xrefs" is satisfied through the gene table (M8-06 Q2).
2. **GO.** `go-basic.obo` read with goatools (`GODag(optional_attrs={"relationship"})`). Each gene's
   terms are mapped to their alt-ID primary, obsolete terms dropped (counted in `n_unknown_ids`),
   then propagated to all ancestors over **`is_a` and `part_of`** (decision 3), and split by
   namespace into `go_bp`, `go_mf`, `go_cc`. The three roots (GO:0008150, GO:0003674, GO:0005575)
   are removed from the sets (they would contain every annotated gene of their namespace).
   `parents` holds direct `is_a`/`part_of` parents within the namespace, for elim and redundancy.
3. **GO slim** (`function.go_slim = generic` → `goslim_generic.obo`): a gene belongs to slim term T
   when any of its propagated terms is T (membership through all ancestors, goatools `mapslim`'s
   `all_ancestors`). One collection `go_slim` across namespaces; no `parents` (classic tests only,
   decision 5). The slim roots are removed as in step 2.
4. **KEGG.** `kegg_ko_pathway.tsv` rows `ko:Kxxxxx\tpath:mapNNNNN` only (the `path:koNNNNN` rows
   are exact duplicates; decision 6); `KEGG:` KOs from the gene table joined to pathways; names
   from `kegg_pathway_names.tsv`. Global and overview maps (`map011xx`, `map012xx`, e.g.
   map01100 "Metabolic pathways") are **excluded** and counted (decision 6, ADR-C).
5. **COG.** One set per category letter (a gene with `cog_category = "KL"` is in K and L);
   letters from `fun-20.tab` with their names; `-`, blank and NA are not a category. The
   category letters "R" (general function) and "S" (unknown) are kept (they are real categories).
6. **Eligibility.** Sets contain eligible genes only (§7.3: non-partial rows with finite
   `p_detect`; the same rule as `regions_table._eligible`, which this issue moves to a shared
   `giap.genetable.eligible()` helper). Coverage = (# eligible genes with ≥ 1 term in the
   collection) / (# eligible genes).
7. **Run flags.** A resource that is absent gives `layer_missing:go`, `layer_missing:kegg` or
   `layer_missing:cog` and no collection; `annotation_source = bakta_light` gives the §9.4 run flag
   `bakta_db_light`. Nothing here applies size filters (M8-04 and M8-05 do, in genes or TUs).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| GO ID absent from go-basic (wrong release, typo, simulator ID) | dropped, counted in `n_unknown_ids`, never raised |
| alt ID (`alt_id:` in the OBO) | mapped to its primary term |
| obsolete term | dropped, counted |
| gene with GO terms in only one namespace | contributes to that namespace's coverage only |
| KO with no pathway | the gene is in no KEGG set and is **not** KEGG-covered (coverage counts genes with ≥ 1 term of the collection, not with ≥ 1 KO) |
| a KEGG pathway with no name row | name `NA`, kept |
| `cog_category` of a COG ID missing from cog-20.def | already NA in the gene table (M4-15); not a category |
| no eligible gene | every collection empty, coverage NA |
| resources present but the gene table has no function columns (user annotation) | every collection empty, coverage 0, no error |

## Examples

Toy OBO (written by the test): `GO:0008150` (BP root) ← `GO:B1` ← `GO:B2` (is_a); `GO:B3`
part_of `GO:B1`; `GO:0003674` (MF root) ← `GO:M1`.

| # | Input | Expected |
|---|---|---|
| 1 | gene g1 `go_terms = [GO:B2]` | g1 in `go_bp` sets `GO:B2` and `GO:B1`; not in the root; `parents[GO:B2] = {GO:B1}` |
| 2 | g2 `[GO:B3]` | g2 in `GO:B3` and `GO:B1` (part_of propagates) |
| 3 | g3 `[GO:M1, GO:B2]` | g3 in `go_mf` `GO:M1` and in `go_bp` `GO:B2`, `GO:B1` |
| 4 | 4 eligible genes, g1–g3 as above, g4 none | `go_bp` coverage **0.75**, `go_mf` coverage **0.25** |
| 5 | g5 `[GO:9999999]` (absent) | dropped; `n_unknown_ids = 1`; g5 unannotated |
| 6 | KEGG link rows `ko:K00001 path:map00010`, `ko:K00001 path:ko00010`, `ko:K00001 path:map01100`; g6 `kegg_ko = [K00001]` | g6 in `map00010` only (ko row deduplicated, map01100 excluded) |
| 7 | `cog_category = "KL"` | gene in COG sets K and L, names from fun-20 |
| 8 | partial gene with `go_terms` | in no set; not in the coverage denominator |
| 9 | `go_basic` not installed | no `go_*` collections, run flag `layer_missing:go`; `go_slim` also absent (it needs the DAG) |
| 10 | `annotation_source = bakta_light` | run flag `bakta_db_light`; collections built from whatever exists |

## Acceptance criteria

- [ ] `tests/unit/function/test_collections.py` covers examples 1–10 and every degenerate row on a
  toy OBO, a toy slim OBO and a **synthetic** two-pathway KEGG table (never real KEGG content: the
  KEGG terms of use, index §KEGG).
- [ ] Contract test on the installed resources (`needs_go_db` / setup-manifest capability, skipped
  in GitHub CI): go-basic 2026-05-19 parses; the K-12 full-DB Bakta fixture gives `go_bp` coverage
  within [0.80, 0.95] (90 % of CDS carry GO xrefs, counted 2026-10-21) and KEGG coverage ≤ 0.40.
- [ ] Property: every gene of a child set is in every ancestor set (propagation is closed); no
  set contains a non-eligible gene.
- [ ] Performance: the K-12 fixture (4,295 CDS) builds all six collections in < 10 s on 1 CPU.
- [ ] `docs/methods/function.md` §"Collections": sources, versions, relations, the excluded maps,
  the coverage definition, KEGG's academic-use restriction.

## Hardening owed (pass 2)

- The backlog GO-header item (M2-09: `verify_go` compares against a pinned `data_version`) does
  not block M8 — the installed release is 2026-05-19 — but a GO release bump must re-run the
  contract test.
- `n_unknown_ids` > 5 % of annotated IDs on a real run should become a report warning (M9).

## Questions for the lead

1. GO relations for propagation: `is_a` only (goatools' default) or `is_a` + `part_of` (GO.db,
   topGO)? **Proposed default:** `is_a` + `part_of`.
2. KEGG: drop the duplicate `path:ko` rows and the overview maps `map011xx`/`map012xx`?
   **Proposed default:** yes to both; the overview maps overlap nearly every metabolic set and,
   where they survive the 500 cap, dominate BH.
3. GO slim: one mixed-namespace collection (as §7.7.4 lists "GO slim") with classic tests only?
   **Proposed default:** yes.

## Lead decisions

1. **Questions 1–3:** the proposed defaults are accepted.
   - Propagate over `is_a` + `part_of`.
   - Drop the `path:ko` duplicate rows and the overview maps (map011xx/map012xx).
   - GO slim is one mixed-namespace collection, tested with classic tests only.
2. **The test author's choices** (commit f850ba6) are accepted as the contract.
   - **Output order:** `cog, kegg, go_bp, go_mf, go_cc, go_slim`. Set keys are sorted.
   - **Empty sets** are left out, so every set is a non-empty frozenset.
   - **No eligible gene:** all six collections are returned empty, with coverage NaN. When the function columns are absent, coverage is 0.0.
   - **`n_unknown_ids`:** the largest value across go_bp/mf/cc equals the number of dropped GO IDs. The implementer may put the count in each namespace.
   - **`parents`:** a term's direct non-root `is_a`/`part_of` parents. Only GO terms have them, so `parents` is empty for go_slim, kegg and cog.
   - **Slim membership** uses the full propagation.
   - **KEGG keys** are bare `mapNNNNN` IDs. A pathway with no name row is kept, and its name is NA.
   - **Run flags:** an exact word, or `word: msg`.
     - Missing GO gives `layer_missing:go` and drops all four GO collections.
     - A missing slim alone drops go_slim only.
     - `bakta_db_light` is raised from `annotation_source == bakta_light`.
   - **Inputs:** numpy-array list cells behave like lists. Inputs are never mutated.
3. **Owed to pass 2:** the contract test on the installed resources (it needs a `needs_go_db` marker) and the < 10 s performance test.
4. **After the code review (#159):**
   - The overview-map exclusion also covers map01310 and map01320, KEGG's full "Global and overview maps" group (ADR-0023).
   - "Excluded and counted" becomes "excluded and documented": there is no count field in pass 1.
   - The methods section is docs/methods/function-collections.md, a separate file so it doesn't collide with the M8-01 `function.md`.
   - The implementer's choices 1–7 are accepted.

