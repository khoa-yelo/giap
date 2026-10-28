// OPERONS subworkflow (GIPPA_SPEC.md §6.2 step row "operons", §7.4, §9.1, §9.3; M6-08): the final
// gene table, the C2/C15/C1 gene model, the C5 TUs and ANNOTATE's `layers_for_pairs` ->
// tables/pairs.tsv, tables/operons.tsv and tables/operon_flags.tsv (the `operons` py step,
// through a GIPPA_STEP alias). It runs beside CONTEXT and needs no C7/C8.
// Publishing rule (M6-08 item 7, for later milestones to copy): only the user-facing
// tables/*.tsv are named in the out spec, so GIPPA_STEP's publishDir copies nothing else; work
// Parquet stays in work/. The required keys go to stepIn's `need` explicitly, so the step does
// not run when one of them is absent (e.g. no gene table without peptides); the evidence layers
// promoter_windows / terminators / tx_rates are optional (NA evidence + layer_missing:<layer>).
include { GIPPA_STEP as OPERONS_STEP } from '../../modules/local/gippa_step/main'
include { stepIn; kv; pick } from './annotate'

workflow OPERONS {
    take:
    gene_table        // genetable/genes.parquet
    annotation        // [genes.parquet, features.parquet, replicons.parquet]
    contract          // [xrefs.parquet, annotation_meta.json, tus.parquet, layers.parquet]
    layers_for_pairs  // [key, path]: pairs_raw, and promoter_windows / terminators / tx_rates /
                      // peptide_hits (C11, main.nf; M6-05H) when present

    main:
    def chans = [kv('gene_table', gene_table), kv('genes', annotation.map { a -> a[0] }),
                 kv('features', annotation.map { a -> a[1] }), kv('replicons', annotation.map { a -> a[2] }),
                 kv('tus', contract.map { c -> c[2] }), layers_for_pairs]
    def need = ['gene_table', 'genes', 'replicons', 'tus', 'pairs_raw', 'features']
    def v = stepIn('operons', chans, need)
    OPERONS_STEP(v.files, v.keys, ['pairs=tables/pairs.tsv', 'operons=tables/operons.tsv',
        'operon_flags=tables/operon_flags.tsv'])
    def out = OPERONS_STEP.out.outputs

    emit:
    pairs        = pick(out, 'pairs.tsv')
    operons      = pick(out, 'operons.tsv')
    operon_flags = pick(out, 'operon_flags.tsv')
    versions     = OPERONS_STEP.out.versions
    events       = OPERONS_STEP.out.events
}
