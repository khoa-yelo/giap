// One Python pipeline step through `gippa _step` (GIPPA_SPEC.md §6.2, §6.3; M2-06).
// inputs[i] is passed as --in in_keys[i]=<staged path>; each input is staged into its own
// directory in1/, in2/, ... (stageAs 'in?/*') because two keys of one step can name the same
// file: context_stats takes both the gene table (genetable/genes.parquet) and the M4-06 C2
// table (genes.parquet), and one `in/` would be a Nextflow input file name collision (MVP-01).
// config_hash is the step's config slice hash (params.config_hashes, gippa.runner.nextflow
// .write_params; stepIn adds it): a task-hash input only, so a changed key re-runs the steps that
// read it (M7-15 decision 7).
// out_spec items are "key=file name"
// (a glob such as windows_*.fa names a scattered output; the stub creates the listed files; in stub
// runs an item "key=file HEADER" makes the stub write HEADER (everything after the first space;
// printf %b, so a literal \n starts a new line) as the file's content, gzipped for a .gz file).
// The optional logs/events.jsonl that `gippa _step` writes (§6.7) is published as
// OUT/logs/steps/<step>.events.jsonl; `gippa run` appends those to OUT/logs/events.jsonl after
// Nextflow ends, in trace order (gippa.runner.nextflow.append_step_events; M3-17 decision 15).
// Outputs named work/<file> are intermediates (M7-14 variants_qc): passed on, never published.
// M9-17: an out spec item '<x>_support=<x>_support' names an optional plot-support directory
// (gippa.io.plot_support.write_support writes <x>_support/support/<name>.tsv[.gz], and only when
// the step has something to write). It is not part of `outputs`; its files are the optional
// `support` emission (a file glob, not the directory: Nextflow publishes a declared directory
// whole and never the files below it), published file by file as OUT/plots/support/<file>
// (never a whole plots/ directory: the plots step publishes that one). supportDir() gives the
// <x>_support directory back, the path the plots step stages as plots_support. Stub runs write
// each step's real support file names (supportStub).

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

// The M9-17 support files each step writes (docs/specs/M9-17-plot-support-analysis.md), by out key.
def supportStub(String key) {
    def names = [function_support     : ['gsea_running.tsv.gz'],
                 determinants_support : ['determinants_shap.tsv.gz', 'determinants_pd.tsv', 'null_determinants.tsv'],
                 context_stats_support: ['null_regions.tsv.gz', 'null_islands.tsv']]
    return names[key] ?: ["${key}.tsv"]
}

def isSupport(String item) { item.tokenize('=')[0].endsWith('_support') }

// The <x>_support directory of a `support` emission (one file or a list of its support/ files).
def supportDir(fs) { (fs instanceof List ? fs[0] : fs).parent.parent }

// OUT-relative publish path of a task output: null for versions.yml, work/<file> and a support
// directory itself; OUT/plots/support/<file> for a support file; logs/steps/<step>.events.jsonl
// for the step's event log; else the output's own path.
def publishAs(String step, String fn) {
    def parts = fn.tokenize('/')
    if (fn == 'versions.yml' || fn.startsWith('work/')) return null
    if (parts[0].endsWith('_support')) return parts.size() == 3 && parts[1] == 'support' ? "plots/support/${parts[2]}" : null
    return fn == 'logs/events.jsonl' ? "logs/steps/${step}.events.jsonl" : fn
}

process GIPPA_STEP {
    tag "${step}"
    label 'process_single'
    cache 'deep'   // content, not path/mtime: a re-run step with unchanged outputs caches its dependants
    publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> publishAs(step, fn) }

    input:
    tuple val(step), path(inputs, stageAs: 'in?/*'), val(config_hash)
    val(in_keys)
    val(out_spec)

    output:
    tuple val(step), path(out_spec.findAll { !isSupport(it) }.collect { it.tokenize('=')[1].tokenize(' ')[0] }), emit: outputs
    path '*_support/support/*', optional: true, emit: support
    path 'versions.yml', emit: versions
    path 'logs/events.jsonl', optional: true, emit: events

    script:
    def files = inputs instanceof List ? inputs : [inputs]
    def ins = [in_keys, files].transpose().collect { k, f -> "--in ${k}=${f}" }.join(' ')
    def outs = out_spec.collect { "--out '${it}'" }.join(' ')
    def resolved = params.config ?: 'resolved_config.yaml'
    """
    gippa _step ${step} --config ${resolved} ${ins} ${outs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gippa: \$(gippa version --json | python -c 'import json,sys; print(json.load(sys.stdin)["gippa"])')
    END_VERSIONS
    """

    stub:   // .gippa_inputs: the --in entries in channel order (staging lines are name-sorted)
    def files = inputs instanceof List ? inputs : [inputs]
    def ins = [in_keys, files].transpose().collect { k, f -> "'${k}=${f}'" }.join(' ')
    def touches = out_spec.findAll { !isSupport(it) }.collect { o ->
        def f = o.tokenize('=')[1].tokenize(' ')[0]
        def hdr = o.contains(' ') ? o.substring(o.indexOf(' ') + 1) : null
        def sink = f.endsWith('.gz') ? "| gzip -n > ${f}" : "> ${f}"
        "mkdir -p \$(dirname ${f}) && " + (hdr ? "printf '%b\\n' '${hdr}' ${sink}" : "touch ${f}")
    }.join('\n    ')
    def supports = out_spec.findAll { isSupport(it) }.collectMany { o ->
        def d = o.tokenize('=')[1].tokenize(' ')[0]
        supportStub(o.tokenize('=')[0]).collect { n ->
            "mkdir -p ${d}/support && " + (n.endsWith('.gz') ? "printf '' | gzip -n > ${d}/support/${n}" : "touch ${d}/support/${n}")
        }
    }.join('\n    ')
    """
    printf '%s\\n' ${ins} > .gippa_inputs
    ${touches}
    ${supports}
    ${stubEvent(step)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gippa: stub
    END_VERSIONS
    """
}
