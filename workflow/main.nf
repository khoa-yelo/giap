// GIAP entry workflow (GIAP_SPEC.md §6.3). Parameters come only from OUT/work/params.json
// written by `giap run`. With a genome and replicons.tsv it runs ANNOTATE (M2-07), then MAPPING
// (M3-17) when params.inputs.peptides (peptides.norm.parquet) is given; without inputs it runs one
// GIAP_STEP smoke task (skeleton check, M2-06).
include { GIAP_STEP } from './modules/local/giap_step/main'
include { ANNOTATE } from './subworkflows/local/annotate'
include { MAPPING } from './subworkflows/local/mapping'

// A path input as a value channel, or channel.empty() when absent.
def opt(path) { path ? channel.value(file(path, checkIfExists: true)) : channel.empty() }

// params.inputs.contaminants, else <db_dir>/contaminants/contaminants.fasta when the db_dir
// manifest lists `contaminants` (M3-15; giap.runner.local.contaminant_library), else null.
// Manifest shapes as giap.preflight.environment.read_manifest: the `giap setup` form
// {"schema_version", "items": [{"name", ..., "status"}]} (only `ok` items count) or the
// provisional {"<name>": ...} map (M3-17 lead decision 18).
def contaminantLibrary(inputs) {
    if (inputs.contaminants) return inputs.contaminants
    def db = params.db_dir ?: System.getenv('GIAP_DB_DIR')
    def manifest = db ? file("${db}/db_manifest.json") : null
    if (manifest == null || !manifest.exists()) return null
    def lib = file("${db}/contaminants/contaminants.fasta")
    def data = new groovy.json.JsonSlurper().parse(manifest.toFile())
    if (!(data instanceof Map)) return null
    def listed = data.items instanceof List
        ? data.items.any { it instanceof Map && it.name == 'contaminants' && it.status == 'ok' }
        : data.containsKey('contaminants')
    return listed && lib.exists() ? lib : null
}

workflow {
    def inputs = (params.inputs ?: [:]).findAll { k, v -> v }
    if (inputs.genome && inputs.replicons) {
        def genome = channel.value(file(inputs.genome, checkIfExists: true))
        ANNOTATE(genome, channel.value(file(inputs.replicons, checkIfExists: true)))
        if (inputs.peptides) {
            MAPPING(genome, opt(inputs.peptides), opt(inputs.proteins_long), ANNOTATE.out.gene_model,
                    opt(inputs.reference_peptides), opt(contaminantLibrary(inputs)))
        }
    } else {
        def keys = inputs.keySet().sort()
        def files = keys.collect { k -> file(inputs[k], checkIfExists: true) }
        GIAP_STEP(channel.of(['smoke', files]), keys, ['marker=smoke.done'])
    }
}
