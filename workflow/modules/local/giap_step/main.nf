// One Python pipeline step through `giap _step` (GIAP_SPEC.md §6.2, §6.3; M2-06).
// inputs[i] is passed as --in in_keys[i]=<staged path>; each input is staged into its own
// directory in1/, in2/, ... (stageAs 'in?/*') because two keys of one step can name the same
// file: context_stats takes both the gene table (genetable/genes.parquet) and the M4-06 C2
// table (genes.parquet), and one `in/` would be a Nextflow input file name collision (MVP-01).
// out_spec items are "key=file name"
// (a glob such as windows_*.fa names a scattered output; the stub creates the listed files; in stub
// runs an item "key=file HEADER" makes the stub write HEADER (everything after the first space;
// printf %b, so a literal \n starts a new line) as the file's content, gzipped for a .gz file).
// The optional logs/events.jsonl that `giap _step` writes (§6.7) is published as
// OUT/logs/steps/<step>.events.jsonl; `giap run` appends those to OUT/logs/events.jsonl after
// Nextflow ends, in trace order (giap.runner.nextflow.append_step_events; M3-17 decision 15).
// Outputs named work/<file> are intermediates (M7-14 variants_qc): passed on, never published.

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
        fn == 'versions.yml' || fn.startsWith('work/') ? null : fn == 'logs/events.jsonl' ? "logs/steps/${step}.events.jsonl" : fn }

    input:
    tuple val(step), path(inputs, stageAs: 'in?/*')
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
        def f = o.tokenize('=')[1].tokenize(' ')[0]
        def hdr = o.contains(' ') ? o.substring(o.indexOf(' ') + 1) : null
        def sink = f.endsWith('.gz') ? "| gzip -n > ${f}" : "> ${f}"
        "mkdir -p \$(dirname ${f}) && " + (hdr ? "printf '%b\\n' '${hdr}' ${sink}" : "touch ${f}")
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
