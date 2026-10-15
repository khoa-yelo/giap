// One Python pipeline step through `giap _step` (GIAP_SPEC.md §6.2, §6.3; M2-06).
// inputs[i] is passed as --in in_keys[i]=<staged path>; out_spec items are "key=file name"
// (a glob such as windows_*.fa names a scattered output; the stub creates the listed files; in stub
// runs an item "key=file HEADER" makes the stub write HEADER as the file's only line).
process GIAP_STEP {
    tag "${step}"
    label 'process_single'
    publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(step), path(inputs, stageAs: 'in/*')
    val(in_keys)
    val(out_spec)

    output:
    tuple val(step), path(out_spec.collect { it.tokenize('=')[1].tokenize(' ')[0] }), emit: outputs
    path 'versions.yml', emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        giap: stub
    END_VERSIONS
    """
}
