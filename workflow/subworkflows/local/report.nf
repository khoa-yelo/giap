// REPORT subworkflow (GIAP_SPEC.md §6.2 rows "plots" and "report", §6.1, §9.1, §10.3; M9-02
// behaviour 8): one channel of [key, path] pairs (every subworkflow's published outputs,
// 'preflight', and `events` emissions under the key 'events') -> plots/ (the `plots` py step)
// -> report/report.html (the `report` py step), both through GIAP_STEP aliases. `preflight` is
// the only required key; a subworkflow that did not run contributes nothing and the plots that
// need it skip themselves. Stub runs name one file inside each output directory (cf. mapping.nf's
// lociSpec()), so -stub creates and publishes plots/ and report/report.html.
include { GIAP_STEP as PLOTS; GIAP_STEP as REPORT_HTML } from '../../modules/local/giap_step/main'
include { stepIn; kv } from './annotate'

workflow REPORT {
    take:
    published   // [key, path] pairs

    main:
    def p = stepIn('plots', [published], ['preflight'])
    PLOTS(p.files, p.keys, workflow.stubRun ? ["plots=plots/stub.done"] : ['plots=plots'])
    def plots = PLOTS.out.outputs.map { s, f -> f }
    def r = stepIn('report', [kv('plots', plots), published], ['plots', 'preflight'])
    REPORT_HTML(r.files, r.keys, workflow.stubRun ? ["report=report/report.html"] : ['report=report'])

    emit:
    plots    = plots
    report   = REPORT_HTML.out.outputs.map { s, f -> f }
    versions = PLOTS.out.versions.mix(REPORT_HTML.out.versions)
}
