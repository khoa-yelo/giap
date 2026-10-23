// DETERMINANTS subworkflow (GIAP_SPEC.md §6.2, §7.8.3, §9.1; M8-10): the final gene table plus
// C1/C2/C5/C6, and optionally C7/C8 and context_stats' gene context, ->
// tables/determinants_{features,model,attribution}.tsv. The py step runs through the GIAP_STEP
// alias DETERMINANTS_STEP (an alias cannot be the subworkflow's own name); base.config gives it
// process_high resources so the permutation null runs in parallel (lead decision 2).
// Without a gene table the step does not run; without C7/C8 or the gene context (CONTEXT
// skipped) it runs on the required inputs.
include { GIAP_STEP as DETERMINANTS_STEP } from '../../modules/local/giap_step/main'
include { stepIn; kv; pick } from './annotate'

workflow DETERMINANTS {
    take:
    gene_table   // genetable/genes.parquet (GENETABLE's final gene table)
    annotation   // [genes.parquet, features.parquet, replicons.parquet]
    contract     // [xrefs.parquet, annotation_meta.json, tus.parquet, layers.parquet]
    regions      // [regions.parquet, region_genes.parquet], or channel.empty()
    gene_context // context/gene_context.parquet (CONTEXT), or channel.empty()

    main:
    def chans = [kv('gene_table', gene_table), kv('genes', annotation.map { a -> a[0] }),
                 kv('replicons', annotation.map { a -> a[2] }),
                 kv('tus', contract.map { c -> c[2] }), kv('layers', contract.map { c -> c[3] }),
                 kv('regions', regions.map { r -> r[0] }), kv('region_genes', regions.map { r -> r[1] }),
                 kv('gene_context', gene_context)]
    def d = stepIn('determinants', chans, ['gene_table', 'genes', 'replicons', 'tus', 'layers'])
    DETERMINANTS_STEP(d.files, d.keys, ['determinants_features=tables/determinants_features.tsv',
        'determinants_model=tables/determinants_model.tsv',
        'determinants_attribution=tables/determinants_attribution.tsv'])

    emit:
    determinants_features    = pick(DETERMINANTS_STEP.out.outputs, 'determinants_features.tsv')
    determinants_model       = pick(DETERMINANTS_STEP.out.outputs, 'determinants_model.tsv')
    determinants_attribution = pick(DETERMINANTS_STEP.out.outputs, 'determinants_attribution.tsv')
    versions                 = DETERMINANTS_STEP.out.versions
}
