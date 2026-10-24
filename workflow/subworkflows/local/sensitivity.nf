// SENSITIVITY_RERUNS (GIAP_SPEC.md §6.2 "sensitivity" row, §8.6; M4-19, M6-10, M8-16): the gene
// table plus C1/C5/C7/C8 -> tables/sensitivity.tsv. Out of CONTEXT since M8-16 (lead decision 4,
// gap 12): main.nf calls it after DETERMINANTS (which needs CONTEXT's gene context) and beside
// FUNCTION, so the determinants family gets DETERMINANTS' features table and the enriched
// categories family the function resources and GENETABLE's TU folds. The operon_flags family
// (M6-10) takes ANNOTATE's layers_for_pairs and, only beside pairs_raw, the C2 genes / C3
// features. Every M8-16 input is optional (the step skips a family without its inputs). Without
// C7/C8 the step does not run. The py step keeps the GIAP_STEP alias SENSITIVITY.
include { GIAP_STEP as SENSITIVITY } from '../../modules/local/giap_step/main'
include { stepIn; kv; pick } from './annotate'

workflow SENSITIVITY_RERUNS {
    take:
    gene_table            // genetable/genes.parquet
    annotation            // [genes.parquet, features.parquet, replicons.parquet]
    contract              // [xrefs.parquet, annotation_meta.json, tus.parquet, layers.parquet]
    regions               // [regions.parquet, region_genes.parquet] (region_layer)
    layers_for_pairs      // [key, path]: pairs_raw and the evidence layers (M6-10), or channel.empty()
    tu_folds              // genetable/tu_folds.parquet, or channel.empty()
    determinants_features // tables/determinants_features.tsv (DETERMINANTS), or channel.empty()
    function_resources    // [setup item name, path] (main.nf functionResources()), or channel.empty()

    main:
    def se = stepIn('sensitivity', [kv('gene_table', gene_table), kv('replicons', annotation.map { a -> a[2] }),
                                    kv('regions', regions.map { r -> r[0] }), kv('region_genes', regions.map { r -> r[1] }),
                                    kv('tus', contract.map { c -> c[2] }), layers_for_pairs,
                                    layers_for_pairs.filter { k, f -> k == 'pairs_raw' }.combine(annotation)
                                        .flatMap { e -> [['genes', e[2]], ['features', e[3]]] },
                                    kv('tu_folds', tu_folds), kv('determinants_features', determinants_features),
                                    function_resources],
                    ['gene_table', 'replicons', 'regions', 'region_genes', 'tus'])
    SENSITIVITY(se.files, se.keys, ['sensitivity=tables/sensitivity.tsv'])

    emit:
    sensitivity = pick(SENSITIVITY.out.outputs, 'sensitivity.tsv')
    versions    = SENSITIVITY.out.versions
    events      = SENSITIVITY.out.events
}
