// GIAP entry workflow (GIAP_SPEC.md §6.3). Parameters come only from OUT/work/params.json
// written by `giap run`. With a genome and replicons.tsv it runs ANNOTATE (M2-07), then MAPPING
// (M3-17), VARIANTS (M7-04), GENETABLE (M4-15), CONTEXT (M4-18, M4-19), beside CONTEXT, OPERONS
// (M6-08), NOVEL (M7-08) and DETERMINANTS (M8-10, after CONTEXT: it takes the gene context) when
// params.inputs.peptides (peptides.norm.parquet) is given, and FUNCTION (M8-06) when 'function' is
// in params.modules; without inputs it runs one GIAP_STEP smoke task (skeleton check, M2-06).
// SENSITIVITY_RERUNS (M8-16) runs after DETERMINANTS and FUNCTION.
// GENETABLE, CONTEXT, OPERONS and DETERMINANTS skip themselves when an input they require is absent
// (MVP-01). Under a precomputed regions layer region_layer still writes C7/C8 (MVP-02). REPORT runs
// last on every published [key, path] pair (M9-02).
include { GIAP_STEP } from './modules/local/giap_step/main'
include { ANNOTATE; kv; pick } from './subworkflows/local/annotate'
include { MAPPING } from './subworkflows/local/mapping'
include { GENETABLE } from './subworkflows/local/genetable'
include { CONTEXT } from './subworkflows/local/context'
include { SENSITIVITY_RERUNS } from './subworkflows/local/sensitivity'
include { VARIANTS } from './subworkflows/local/variants'
include { NOVEL } from './subworkflows/local/novel'
include { DETERMINANTS } from './subworkflows/local/determinants'
include { OPERONS } from './subworkflows/local/operons'
include { REPORT } from './subworkflows/local/report'
include { FUNCTION } from './subworkflows/local/function'

// A path input as a value channel, or channel.empty() when absent.
def opt(path) { path ? channel.value(file(path, checkIfExists: true)) : channel.empty() }

// <db_dir>/<relative path> when the db_dir setup manifest lists `name` and the file is there,
// else null (giap.setupdb.installed / giap.runner.local, which resolve the same two items).
// Manifest shapes as giap.preflight.environment.read_manifest: the `giap setup` form
// {"schema_version", "items": [{"name", ..., "status"}]} (only `ok` items count) or the
// provisional {"<name>": ...} map (M3-17 lead decision 18).
def installed(String name, String relative) {
    def db = params.db_dir ?: System.getenv('GIAP_DB_DIR')
    def manifest = db ? file("${db}/db_manifest.json") : null
    if (manifest == null || !manifest.exists()) return null
    def path = file("${db}/${relative}")
    def data = new groovy.json.JsonSlurper().parse(manifest.toFile())
    if (!(data instanceof Map)) return null
    def listed = data.items instanceof List
        ? data.items.any { it instanceof Map && it.name == name && it.status == 'ok' }
        : data.containsKey(name)
    return listed && path.exists() ? path : null
}

// params.inputs.contaminants, else the installed M3-15 library (giap.runner.local.contaminant_library).
def contaminantLibrary(inputs) {
    return inputs.contaminants ?: installed('contaminants', 'contaminants/contaminants.fasta')
}

// The raw NCBI COG 2020 table for gene_table's cog_category (M4-15 decision 5; the path is
// db_sources.yaml cog2020_def.path, as giap.setupdb.installed.cog_table resolves it).
def cogTable() { installed('cog2020_def', 'cog/cog-20.def.tab') }

// [key, path] pairs of a step's outputs channel ([step, path | list]): every path, under `key`.
def outs(String key, ch) { ch.flatMap { s, fs -> (fs instanceof List ? fs : [fs]).collect { f -> [key, f] } } }
// The installed function resources (M8-06 item 7) as a [setup item name, path] channel; the names
// are giap.setupdb.installed.FUNCTION_RESOURCES and the paths db_sources.yaml's.
def functionResources() {
    def found = [go_basic: installed('go_basic', 'go/go-basic.obo'),
                 goslim_generic: installed('goslim_generic', 'go/goslim_generic.obo'),
                 kegg_ko_pathway: installed('kegg_ko_pathway', 'kegg/ko_pathway.tsv'),
                 kegg_pathway_names: installed('kegg_pathway_names', 'kegg/pathway_names.tsv'),
                 cog2020_fun: installed('cog2020_fun', 'cog/fun-20.tab')]
    return channel.fromList(found.findAll { k, p -> p != null }.collect { k, p -> [k, p] })
}

workflow {
    def inputs = (params.inputs ?: [:]).findAll { k, v -> v }
    if (inputs.genome && inputs.replicons) {
        def genome = channel.value(file(inputs.genome, checkIfExists: true))
        ANNOTATE(genome, channel.value(file(inputs.replicons, checkIfExists: true)))
        def gm = ANNOTATE.out.gene_model
        def pub = [kv('preflight', opt(inputs.preflight)), kv('genome', genome),
                   kv('annotation', pick(ANNOTATE.out.context, 'annotation')),
                   kv('genes', gm.map { a -> a[0] }), kv('features', gm.map { a -> a[1] }),
                   kv('replicons', gm.map { a -> a[2] })]
        if (inputs.peptides) {
            def gene_model = ANNOTATE.out.gene_model
            MAPPING(genome, opt(inputs.peptides), opt(inputs.proteins_long), gene_model,
                    opt(inputs.reference_peptides), opt(contaminantLibrary(inputs)))
            // Module 5 (M7-04): independent of GENETABLE; CONTEXT takes its variants_qc (M7-14)
            VARIANTS(genome, gene_model, MAPPING.out.mapping, MAPPING.out.hits_tolerant,
                     opt(inputs.peptides), opt(inputs.reference_peptides), opt(inputs.peptide_spectra))
            GENETABLE(gene_model, ANNOTATE.out.contract, MAPPING.out.mapping, opt(inputs.peptides),
                      opt(inputs.abundance), opt(inputs.preflight), opt(cogTable()))
            CONTEXT(GENETABLE.out.gene_table, GENETABLE.out.detectability_summary, gene_model,
                    ANNOTATE.out.contract, ANNOTATE.out.regions, MAPPING.out.mapping,
                    MAPPING.out.stratum0, opt(inputs.peptides), opt(inputs.abundance),
                    opt(inputs.replicon_depth), GENETABLE.out.tier_summary,
                    ANNOTATE.out.layers_for_pairs, VARIANTS.out.variants_qc)
            NOVEL(genome, gene_model, ANNOTATE.out.contract, MAPPING.out.mapping, opt(inputs.peptides),
                  opt(contaminantLibrary(inputs)), GENETABLE.out.gene_table, ANNOTATE.out.layer_pairs,
                  MAPPING.out.hits_exact.combine(MAPPING.out.hits_tolerant).combine(MAPPING.out.placement_costs).first())
            OPERONS(GENETABLE.out.gene_table, gene_model, ANNOTATE.out.contract, ANNOTATE.out.layers_for_pairs)
            DETERMINANTS(GENETABLE.out.gene_table, gene_model, ANNOTATE.out.contract,
                         ANNOTATE.out.regions, CONTEXT.out.gene_context)
            // M8-16: after DETERMINANTS and beside FUNCTION (its enriched categories family)
            SENSITIVITY_RERUNS(GENETABLE.out.gene_table, gene_model, ANNOTATE.out.contract, ANNOTATE.out.regions,
                               ANNOTATE.out.layers_for_pairs, GENETABLE.out.tu_folds,
                               DETERMINANTS.out.determinants_features,
                               'function' in (params.modules ?: []) ? functionResources() : channel.empty())
            pub += [outs('mapping', MAPPING.out.mapping),
                    kv('pileup', VARIANTS.out.pileup), kv('variants', VARIANTS.out.variants),
                    kv('variants_vcf', VARIANTS.out.variants_vcf), kv('variants_qc', VARIANTS.out.variants_qc),
                    kv('gene_table', GENETABLE.out.gene_table), kv('detectability_summary', GENETABLE.out.detectability_summary),
                    kv('tier_summary', GENETABLE.out.tier_summary), kv('genes_pre', GENETABLE.out.genes_pre),
                    kv('tu_folds', GENETABLE.out.tu_folds),
                    kv('gene_context', CONTEXT.out.gene_context),
                    kv('island_diagnostics', CONTEXT.out.island_diagnostics),
                    kv('regions_stats', CONTEXT.out.regions_stats), kv('detection_summary', CONTEXT.out.detection_summary),
                    kv('autocorrelation', CONTEXT.out.autocorrelation), kv('islands', CONTEXT.out.islands),
                    kv('genes_tsv', CONTEXT.out.genes_tsv), kv('sensitivity', SENSITIVITY_RERUNS.out.sensitivity),
                    kv('pairs', OPERONS.out.pairs), kv('operons', OPERONS.out.operons),
                    kv('operon_flags', OPERONS.out.operon_flags),
                    kv('novel_orfs', NOVEL.out.novel_orfs), kv('novel_orfs_gff3', NOVEL.out.novel_orfs_gff3),
                    kv('novel_orfs_faa', NOVEL.out.novel_orfs_faa), kv('annotation_patch', NOVEL.out.annotation_patch),
                    kv('novel_class_fdr', NOVEL.out.novel_class_fdr),
                    kv('determinants_features', DETERMINANTS.out.determinants_features),
                    kv('determinants_model', DETERMINANTS.out.determinants_model),
                    kv('determinants_attribution', DETERMINANTS.out.determinants_attribution),
                    // M9-17: the analysis steps' optional plot support, staged for the plots step
                    kv('plots_support', CONTEXT.out.context_stats_support),
                    kv('plots_support', DETERMINANTS.out.determinants_support)]
            if ('function' in (params.modules ?: [])) {  // §6.3: the optional module switch
                FUNCTION(GENETABLE.out.gene_table, gene_model, GENETABLE.out.tu_folds, functionResources())
                pub += [kv('category_summary', FUNCTION.out.category_summary),
                        kv('enrichment_detection', FUNCTION.out.enrichment_detection),
                        kv('enrichment_gsea', FUNCTION.out.enrichment_gsea),
                        kv('plots_support', FUNCTION.out.function_support)]
            }
        }
        // the step event logs, merged in a content order, under the key 'events' as the local
        // executor's live logs/events.jsonl: the plots read their run flags from it (M9-13)
        def ev = inputs.peptides ? ANNOTATE.out.events.mix(MAPPING.out.events, VARIANTS.out.events,
            GENETABLE.out.events, CONTEXT.out.events, SENSITIVITY_RERUNS.out.events, NOVEL.out.events,
            OPERONS.out.events) : ANNOTATE.out.events
        if (inputs.peptides && 'function' in (params.modules ?: [])) ev = ev.mix(FUNCTION.out.events)
        pub += [kv('events', ev.collectFile(name: 'events.jsonl', sort: { f -> f.text }))]
        REPORT(pub.inject(channel.empty()) { acc, c -> acc.mix(c) })
    } else {
        def keys = inputs.keySet().sort()
        def files = keys.collect { k -> file(inputs[k], checkIfExists: true) }
        GIAP_STEP(channel.of(['smoke', files]), keys, ['marker=smoke.done'])
    }
}
