// One Python pipeline step through `giap _step` (GIAP_SPEC.md §6.2, §6.3; M2-06).
// inputs[i] is passed as --in in_keys[i]=<staged path>; out_spec items are "key=file name"
// (a glob such as windows_*.fa names a scattered output; the stub creates the listed files; in stub
// runs an item "key=file HEADER" makes the stub write HEADER as the file's only line).
// The optional logs/events.jsonl that `giap _step` writes (§6.7) is published as
// OUT/logs/steps/<step>.events.jsonl; `giap run` appends those to OUT/logs/events.jsonl after
// Nextflow ends, in trace order (giap.runner.nextflow.append_step_events; M3-17 decision 15).

// Stub runs: steps whose real run always logs an event write one §6.7 line of the same shape
// (ts, step, level, event, data). Only what the stub knows goes in: placement's
// mapping.search_db_resolved keeps data.key; the resolved value and query count need the peptide
// table and are left out (data.stub = true marks the line).
def stubEvent(String step) {
    def events = [placement: [event: 'mapping.search_db_resolved', data: [key: 'search_db', stub: true]]]
    def e = events[step]
    if (!e) return ''
    def body = groovy.json.JsonOutput.toJson([data: e.data, event: e.event, level: 'INFO', step: step])
    def json = body.substring(1)   // the ts key goes first: {"ts":"<now>",...}
    return "mkdir -p logs && printf '{\"ts\":\"%s\",%s\\n' \"\$(date -u +%Y-%m-%dT%H:%M:%S+00:00)\" '${json}' > logs/events.jsonl"
}

process GIAP_STEP {
    tag "${step}"
    label 'process_single'
    publishDir "${params.outdir}", mode: 'copy', saveAs: { fn ->
        fn == 'versions.yml' ? null : fn == 'logs/events.jsonl' ? "logs/steps/${step}.events.jsonl" : fn }

    input:
    tuple val(step), path(inputs, stageAs: 'in/*')
    val(in_keys)
    val(out_spec)

    output:
    tuple val(step), path(out_spec.collect { it.tokenize('=')[1].tokenize(' ')[0] }), emit: outputs
    path 'versions.yml', emit: versions
    path 'logs/events.jsonl', optional: true, emit: events

    script:
    def files = inputs instanceof List ? inputs : [inputs]
    def ins = [in_keys, files].transpose().collect { k, f -> "--in ${k}=${f}" }.join(' ')
    def outs = out_spec.collect { "--out '${it}'" }.join(' ')
    def resolved = params.config ?: 'resolved_config.yaml'
    """
    giap _step ${step} --config ${resolved} ${ins} ${outs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        giap: \$(giap version --json | python -c 'import json,sys; print(json.load(sys.stdin)["giap"])')
    END_VERSIONS
    """

    stub:   // .giap_inputs: the --in entries in channel order (staging lines are name-sorted)
    def files = inputs instanceof List ? inputs : [inputs]
    def ins = [in_keys, files].transpose().collect { k, f -> "'${k}=${f}'" }.join(' ')
    def touches = out_spec.collect { o ->
        def (f, hdr) = [o.tokenize('=')[1].tokenize(' ')[0], o.tokenize(' ')[1]]
        "mkdir -p \$(dirname ${f}) && " + (hdr ? "printf '%s\\n' '${hdr}' > ${f}" : "touch ${f}")
    }.join('\n    ')
    """
    printf '%s\\n' ${ins} > .giap_inputs
    ${touches}
    ${stubEvent(step)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        giap: stub
    END_VERSIONS
    """
}
