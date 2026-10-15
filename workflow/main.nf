// GIAP entry workflow (GIAP_SPEC.md §6.3). Parameters come only from OUT/work/params.json
// written by `giap run`. With a genome and replicons.tsv it runs ANNOTATE (M2-07); without
// inputs it runs one GIAP_STEP smoke task (skeleton check, M2-06).
include { GIAP_STEP } from './modules/local/giap_step/main'
include { ANNOTATE } from './subworkflows/local/annotate'

workflow {
    def inputs = (params.inputs ?: [:]).findAll { k, v -> v }
    if (inputs.genome && inputs.replicons) {
        ANNOTATE(channel.value(file(inputs.genome, checkIfExists: true)),
                 channel.value(file(inputs.replicons, checkIfExists: true)))
    } else {
        def keys = inputs.keySet().sort()
        def files = keys.collect { k -> file(inputs[k], checkIfExists: true) }
        GIAP_STEP(channel.of(['smoke', files]), keys, ['marker=smoke.done'])
    }
}
