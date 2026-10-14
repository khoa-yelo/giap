// GIAP entry workflow (GIAP_SPEC.md §6.3). Parameters come only from OUT/work/params.json
// written by `giap run`. Until ANNOTATE (M2-07) lands, the workflow runs one GIAP_STEP smoke task.
include { GIAP_STEP } from './modules/local/giap_step/main'

workflow {
    def inputs = (params.inputs ?: [:]).findAll { k, v -> v }
    def keys = inputs.keySet().sort()
    def files = keys.collect { k -> file(inputs[k], checkIfExists: true) }
    GIAP_STEP(channel.of(['smoke', files]), keys, ['marker=smoke.done'])
}
