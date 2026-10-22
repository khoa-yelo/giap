// GENETABLE subworkflow (GIAP_SPEC.md §6.2, §7.3, §9.3; M4-15, MVP-01): the M4-06 contract
// tables plus MAPPING's loci directory -> genetable/genes_pre.parquet and the M4-12 tier rows of
// detection_summary (gene_table), then the cross-fitted p_detect -> genetable/genes.parquet, the
// detectability diagnostics and the TU folds (detectability, genetable/tu_folds.parquet: M8-06). Both py steps run through GIAP_STEP aliases, as in annotate.nf / mapping.nf.
// Every required input must arrive or neither step runs (stepIn's `need`), the way MAPPING skips
// its protein route; the optional ones (peptides, abundance, preflight.json, the COG 2020 table)
// are staged only when `giap run` passed them.
include { GIAP_STEP as GENE_TABLE; GIAP_STEP as DETECTABILITY } from '../../modules/local/giap_step/main'
include { stepIn; kv; pick } from './annotate'

// The loci output as one path for the `mapping` key: the mapping/ directory itself in a real
// run (loci declares `mapping=mapping`), and in a stub run — where loci declares the seven
// §7.2.9 files one by one — the directory those files were written into.
def mappingDir(ch) {
    ch.map { s, fs ->
        def files = fs instanceof List ? fs : [fs]
        files.find { f -> f.name == 'mapping' && f.isDirectory() } ?: files[0].parent
    }
}

workflow GENETABLE {
    take:
    annotation  // [genes.parquet, features.parquet, replicons.parquet]
    contract    // [xrefs.parquet, annotation_meta.json, tus.parquet, layers.parquet]
    mapping     // LOCI outputs (the mapping/ directory)
    peptides    // peptides.norm.parquet, or channel.empty()
    abundance   // abundance.norm.parquet, or channel.empty()
    preflight   // preflight.json, or channel.empty()
    cog2020_def // cog-20.def.tab, or channel.empty()

    main:
    def required = ['genes', 'replicons', 'xrefs', 'annotation_meta', 'tus', 'layers', 'mapping']
    def tables = [kv('genes', annotation.map { a -> a[0] }), kv('replicons', annotation.map { a -> a[2] }),
                  kv('xrefs', contract.map { c -> c[0] }), kv('annotation_meta', contract.map { c -> c[1] }),
                  kv('tus', contract.map { c -> c[2] }), kv('layers', contract.map { c -> c[3] })]
    def optional = [kv('peptides', peptides), kv('abundance', abundance),
                    kv('preflight', preflight), kv('cog2020_def', cog2020_def)]

    def gt = stepIn('gene_table', tables + [kv('mapping', mappingDir(mapping))] + optional, required)
    GENE_TABLE(gt.files, gt.keys, ['genes_pre=genetable/genes_pre.parquet',
        'tier_summary=genetable/detection_summary_parts/tier.parquet'])

    // `need` even with one input: stepIn emits an empty task when every channel is empty
    def dt = stepIn('detectability', [kv('genes_pre', pick(GENE_TABLE.out.outputs, 'genes_pre.parquet'))], ['genes_pre'])
    DETECTABILITY(dt.files, dt.keys, ['gene_table=genetable/genes.parquet',
        'detectability_summary=genetable/detection_summary_parts/detectability.parquet',
        'tu_folds=genetable/tu_folds.parquet'])
    def D = DETECTABILITY.out.outputs

    emit:
    gene_table            = pick(D, 'genes.parquet')
    detectability_summary = pick(D, 'detectability.parquet')
    tu_folds              = pick(D, 'tu_folds.parquet')  // M8-06: the cross-fit folds (FUNCTION)
    tier_summary          = pick(GENE_TABLE.out.outputs, 'tier.parquet')
    versions              = GENE_TABLE.out.versions.mix(DETECTABILITY.out.versions)
}
