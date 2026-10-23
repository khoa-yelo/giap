// CONTEXT subworkflow (GIAP_SPEC.md §6.2, §9.4; M4-18, MVP-01, M5-09): the final gene table plus
// the M4-06 contract tables -> tables/regions_stats.tsv, tables/detection_summary.tsv,
// tables/autocorrelation.tsv, tables/islands.tsv and tables/genes.tsv (context_stats, through a
// GIAP_STEP alias). The §8.6 reruns (sensitivity) moved to SENSITIVITY_RERUNS (this file,
// below, M8-16), which main.nf runs after DETERMINANTS.
// context_stats also writes two intermediate parquet files under context/: the genes.tsv
// context group it left-joins into genes.tsv, and the ADR-0021 Amendment 1 island diagnostics.
// The optional C15 features table (the oriC of the ori-ter rows) is the annotation channel's
// second file, so CONTEXT needs no further take: channel (M4-20 note).
// C7/C8 (regions, region_genes) are required; ANNOTATE's region_layer writes them, also under a
// precomputed regions layer (MVP-02). If the `regions` channel is empty context_stats does not
// run, the way MAPPING skips its protein route. The M3 tables come out of MAPPING's loci directory.
// `layers_for_pairs` is no longer used here (SENSITIVITY_RERUNS takes it); the take: is kept.
include { GIAP_STEP as CONTEXT_STATS; GIAP_STEP as SENSITIVITY } from '../../modules/local/giap_step/main'
include { stepIn; kv; pick; inDir } from './annotate'

// One named file of the loci output: a file of the mapping/ directory in a real run, the
// matching staged file in a stub run (inDir handles both).
def ofMapping(ch, String name) { inDir(ch, 'mapping').filter { f -> f.name == name } }

workflow CONTEXT {
    take:
    gene_table            // genetable/genes.parquet
    detectability_summary // detectability diagnostics parquet, or channel.empty()
    annotation            // [genes.parquet, features.parquet, replicons.parquet]
    contract              // [xrefs.parquet, annotation_meta.json, tus.parquet, layers.parquet]
    regions               // [regions.parquet, region_genes.parquet] (region_layer)
    mapping               // LOCI outputs (the mapping/ directory)
    stratum0_classes      // ADR-0019 class table, or channel.empty()
    peptides              // peptides.norm.parquet, or channel.empty()
    abundance             // abundance.norm.parquet, or channel.empty()
    replicon_depth        // §5.3 depth TSV, or channel.empty()
    tier_summary          // gene_table's tier counts / q-threshold curve (M4-20 decision 3), or channel.empty()
    layers_for_pairs      // [key, path]: pairs_raw and promoter_windows / terminators / tx_rates (M6-10), or channel.empty()

    main:
    def gt = kv('gene_table', gene_table)
    def reps = kv('replicons', annotation.map { a -> a[2] })
    def c7 = kv('regions', regions.map { r -> r[0] })
    def c8 = kv('region_genes', regions.map { r -> r[1] })
    def c5 = kv('tus', contract.map { c -> c[2] })
    def tables = [kv('genes', annotation.map { a -> a[0] }), reps,
                  kv('xrefs', contract.map { c -> c[0] }), kv('annotation_meta', contract.map { c -> c[1] }),
                  c5, kv('layers', contract.map { c -> c[3] }), c7, c8]
    def m3 = ['peptide_hits', 'protein_hits', 'unmapped'].collect { n -> kv(n, ofMapping(mapping, "${n}.tsv")) }
    def optional = [kv('peptides', peptides), kv('abundance', abundance),
                    kv('detectability_summary', detectability_summary),
                    kv('features', annotation.map { a -> a[1] }),
                    kv('stratum0_classes', stratum0_classes), kv('replicon_depth', replicon_depth),
                    kv('tier_summary', tier_summary)]
    def needed = ['gene_table', 'genes', 'replicons', 'xrefs', 'annotation_meta', 'tus', 'layers',
                  'regions', 'region_genes', 'peptide_hits', 'protein_hits', 'unmapped']

    def cs = stepIn('context_stats', [gt] + tables + m3 + optional, needed)
    CONTEXT_STATS(cs.files, cs.keys, ['regions_stats=tables/regions_stats.tsv',
        'detection_summary=tables/detection_summary.tsv',
        'autocorrelation=tables/autocorrelation.tsv', 'islands=tables/islands.tsv',
        'genes_tsv=tables/genes.tsv', 'gene_context=context/gene_context.parquet',
        'island_diagnostics=context/island_diagnostics.parquet'])

    emit:
    regions_stats    = pick(CONTEXT_STATS.out.outputs, 'regions_stats.tsv')
    detection_summary = pick(CONTEXT_STATS.out.outputs, 'detection_summary.tsv')
    autocorrelation  = pick(CONTEXT_STATS.out.outputs, 'autocorrelation.tsv')
    islands          = pick(CONTEXT_STATS.out.outputs, 'islands.tsv')
    genes_tsv        = pick(CONTEXT_STATS.out.outputs, 'genes.tsv')
    gene_context     = pick(CONTEXT_STATS.out.outputs, 'gene_context.parquet')
    island_diagnostics = pick(CONTEXT_STATS.out.outputs, 'island_diagnostics.parquet')
    versions         = CONTEXT_STATS.out.versions
    events           = CONTEXT_STATS.out.events
}

// SENSITIVITY_RERUNS (GIAP_SPEC.md §6.2 "sensitivity" row, §8.6; M4-19, M6-10, M8-16): the gene
// table plus C1/C5/C7/C8 -> tables/sensitivity.tsv. Out of CONTEXT since M8-16 (lead decision 4,
// gap 12): main.nf calls it after DETERMINANTS (which needs CONTEXT's gene context) and beside
// FUNCTION, so the determinants family gets DETERMINANTS' features table and the enriched
// categories family the function resources and GENETABLE's TU folds. The operon_flags family
// (M6-10) takes ANNOTATE's layers_for_pairs and, only beside pairs_raw, the C2 genes / C3
// features. Every M8-16 input is optional (the step skips a family without its inputs). Without
// C7/C8 the step does not run.
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
