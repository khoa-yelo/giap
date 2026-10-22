// GIAP entry workflow (GIAP_SPEC.md §6.3). Parameters come only from OUT/work/params.json
// written by `giap run`. With a genome and replicons.tsv it runs ANNOTATE (M2-07), then MAPPING
// (M3-17), GENETABLE (M4-15) and CONTEXT (M4-18, M4-19) when params.inputs.peptides
// (peptides.norm.parquet) is given; without inputs it runs one GIAP_STEP smoke task (skeleton
// check, M2-06). GENETABLE and CONTEXT skip themselves when an input they require is absent
// (MVP-01). Under a precomputed regions layer region_layer still writes C7/C8 (MVP-02).
include { GIAP_STEP } from './modules/local/giap_step/main'
include { ANNOTATE } from './subworkflows/local/annotate'
include { MAPPING } from './subworkflows/local/mapping'
include { GENETABLE } from './subworkflows/local/genetable'
include { CONTEXT } from './subworkflows/local/context'

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

workflow {
    def inputs = (params.inputs ?: [:]).findAll { k, v -> v }
    if (inputs.genome && inputs.replicons) {
        def genome = channel.value(file(inputs.genome, checkIfExists: true))
        ANNOTATE(genome, channel.value(file(inputs.replicons, checkIfExists: true)))
        if (inputs.peptides) {
            def gene_model = ANNOTATE.out.gene_model
            MAPPING(genome, opt(inputs.peptides), opt(inputs.proteins_long), gene_model,
                    opt(inputs.reference_peptides), opt(contaminantLibrary(inputs)))
            GENETABLE(gene_model, ANNOTATE.out.contract, MAPPING.out.mapping, opt(inputs.peptides),
                      opt(inputs.abundance), opt(inputs.preflight), opt(cogTable()))
            CONTEXT(GENETABLE.out.gene_table, GENETABLE.out.detectability_summary, gene_model,
                    ANNOTATE.out.contract, ANNOTATE.out.regions, MAPPING.out.mapping,
                    MAPPING.out.stratum0, opt(inputs.peptides), opt(inputs.abundance),
                    opt(inputs.replicon_depth), GENETABLE.out.tier_summary)
        }
    } else {
        def keys = inputs.keySet().sort()
        def files = keys.collect { k -> file(inputs[k], checkIfExists: true) }
        GIAP_STEP(channel.of(['smoke', files]), keys, ['marker=smoke.done'])
    }
}
