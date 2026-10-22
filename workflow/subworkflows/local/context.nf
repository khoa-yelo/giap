// CONTEXT subworkflow (GIAP_SPEC.md §6.2, §8.6, §9.4; M4-18, M4-19, MVP-01, M5-09): the final
// gene table plus the M4-06 contract tables -> tables/regions_stats.tsv,
// tables/detection_summary.tsv, tables/autocorrelation.tsv, tables/islands.tsv and
// tables/genes.tsv (context_stats), and the §8.6 reruns -> tables/sensitivity.tsv
// (sensitivity, registered by M4-19 decision 16). Both py steps run through GIAP_STEP aliases.
// context_stats also writes two intermediate parquet files under context/: the genes.tsv
// context group it left-joins into genes.tsv, and the ADR-0021 Amendment 1 island diagnostics.
// The optional C15 features table (the oriC of the ori-ter rows) is the annotation channel's
// second file, so CONTEXT needs no further take: channel (M4-20 note).
// C7/C8 (regions, region_genes) are required by both steps; ANNOTATE's region_layer writes them,
// also under a precomputed regions layer (MVP-02). If the `regions` channel is empty neither step
// runs, the way MAPPING skips its protein route. The M3 tables come out of MAPPING's loci directory.
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

    def se = stepIn('sensitivity', [gt, reps, c7, c8, c5],
                    ['gene_table', 'replicons', 'regions', 'region_genes', 'tus'])
    SENSITIVITY(se.files, se.keys, ['sensitivity=tables/sensitivity.tsv'])

    emit:
    regions_stats    = pick(CONTEXT_STATS.out.outputs, 'regions_stats.tsv')
    detection_summary = pick(CONTEXT_STATS.out.outputs, 'detection_summary.tsv')
    autocorrelation  = pick(CONTEXT_STATS.out.outputs, 'autocorrelation.tsv')
    islands          = pick(CONTEXT_STATS.out.outputs, 'islands.tsv')
    genes_tsv        = pick(CONTEXT_STATS.out.outputs, 'genes.tsv')
    gene_context     = pick(CONTEXT_STATS.out.outputs, 'gene_context.parquet')
    island_diagnostics = pick(CONTEXT_STATS.out.outputs, 'island_diagnostics.parquet')
    sensitivity      = pick(SENSITIVITY.out.outputs, 'sensitivity.tsv')
    versions         = CONTEXT_STATS.out.versions.mix(SENSITIVITY.out.versions)
}
