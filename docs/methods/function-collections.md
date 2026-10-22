# Function collections (M8-03)

Gene sets for ORA (§7.7.3) and GSEA (§7.7.4) are built by `giap.function.collections.build_collections`. The installed resources are read by `giap.io.function_resources.read_resources`.

## Members
- **Who can be a member:** only genes eligible for detection-aware testing (`giap.genetable.eligible`: not partial, and a finite `p_detect`). Empty sets are left out.
- **Coverage:** the share of eligible genes in at least one set. It is NaN when no gene is eligible.

## GO (go_bp, go_mf, go_cc)
- **Annotations:** gene-table `go_terms` IDs are mapped to their primary term through go-basic `alt_id`. Obsolete or absent IDs are dropped and counted in `n_unknown_ids`, once per (eligible gene, ID).
- **Propagation:** annotations are propagated to every ancestor over `is_a` and `part_of`, as GO.db and topGO do; goatools' default would be `is_a` only. The three namespace roots are then removed.
- **Namespaces:** go-basic has no cross-namespace edges, which was checked on release 2026-05-19.
- **Parents:** `parents` holds each term's direct non-root parents in the same namespace. GO elim uses it.

## GO slim (go_slim)
- **Membership:** a gene belongs to a slim term when the term is in the gene's full propagated ancestor set.
- **Testing:** the collection mixes namespaces and is tested with classic tests only.

## KEGG
- **Sets:** one per `mapNNNNN`, built from the `path:map` rows of `link/pathway/ko`.
- **Excluded:** the global and overview maps (ADR-0023).
- **Unknown KOs:** a KO absent from the link table counts as unknown. That includes valid KOs that have no pathway.

## COG
One set per COG 2020 functional-category letter from the gene table's `cog_category`, with names from `fun-20.tab`. `-` is not a category.

## Run flags
- `layer_missing:go` drops all four GO collections. A missing slim alone drops only go_slim.
- `layer_missing:kegg` and `layer_missing:cog` are raised when those resources are absent.
- `bakta_db_light` is raised when the annotation came from the Bakta light database.
