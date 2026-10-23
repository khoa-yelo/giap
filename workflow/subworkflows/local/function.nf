// FUNCTION subworkflow (GIAP_SPEC.md §6.2 step row "function", §6.3, §7.7, §9.1; M8-06b): the
// final gene table, the C1 replicons, the TU folds and the installed function resources ->
// tables/category_summary.tsv, tables/enrichment_detection.tsv and tables/enrichment_gsea.tsv
// (the `function` py step, through a GIAP_STEP alias). main.nf calls it only when 'function' is
// in params.modules. gene_table, replicons and tu_folds are required (stepIn's `need`); each
// resource ([setup item name, path]) is staged only when installed, else the step records
// layer_missing:<layer>. Only the three tables/*.tsv are published (M6-08 decision 7).
include { GIAP_STEP as FUNCTION_STEP } from '../../modules/local/giap_step/main'
include { stepIn; kv; pick } from './annotate'

workflow FUNCTION {
    take:
    gene_table  // genetable/genes.parquet
    annotation  // [genes.parquet, features.parquet, replicons.parquet]
    tu_folds    // genetable/tu_folds.parquet
    resources   // [key, path] per installed function resource, or channel.empty()

    main:
    def chans = [kv('gene_table', gene_table), kv('replicons', annotation.map { a -> a[2] }),
                 kv('tu_folds', tu_folds), resources]
    def v = stepIn('function', chans, ['gene_table', 'replicons', 'tu_folds'])
    FUNCTION_STEP(v.files, v.keys, ['category_summary=tables/category_summary.tsv',
        'enrichment_detection=tables/enrichment_detection.tsv',
        'enrichment_gsea=tables/enrichment_gsea.tsv'])
    def out = FUNCTION_STEP.out.outputs

    emit:
    category_summary     = pick(out, 'category_summary.tsv')
    enrichment_detection = pick(out, 'enrichment_detection.tsv')
    enrichment_gsea      = pick(out, 'enrichment_gsea.tsv')
    versions             = FUNCTION_STEP.out.versions
    events               = FUNCTION_STEP.out.events
}
