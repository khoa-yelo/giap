// Mapping Python steps through `giap _step` (GIAP_SPEC.md §6.2, §6.3; M3-17). One task per step:
// in_keys[i] names inputs[i] (--in), out_spec items are "key=file name" (--out).
// GIAP_MAPPING runs the work-file steps; GIAP_MAPPING_LOCI writes and publishes the seven §7.2.9
// files (params.outdir/mapping/); its stub writes the §9.3 headers (giap.mapping.outputs).
process GIAP_MAPPING {
    tag "${step}"
    label 'process_single'

    input:
    tuple val(step), val(in_keys), path(inputs, stageAs: 'in/*'), val(out_spec)

    output:
    tuple val(step), path(out_spec.collect { it.tokenize('=')[1] }), emit: outputs
    path 'versions.yml', emit: versions

    script:
    def files = inputs instanceof List ? inputs : [inputs]
    def ins = [in_keys, files].transpose().collect { k, f -> "--in ${k}=${f}" }.join(' ')
    def outs = out_spec.collect { "--out ${it}" }.join(' ')
    def resolved = params.config ?: 'resolved_config.yaml'
    """
    giap _step ${step} --config ${resolved} ${ins} ${outs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        giap: \$(giap version --json | python -c 'import json,sys; print(json.load(sys.stdin)["giap"])')
    END_VERSIONS
    """

    stub:
    def touches = out_spec.collect { "touch ${it.tokenize('=')[1]}" }.join('\n    ')
    """
    ${touches}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        giap: stub
    END_VERSIONS
    """
}

// Shell line writing the header (space-separated `cols`, joined by tabs) of mapping/`file`.
def tsvHeader(String cols, String file) {
    return "printf '%s\\n' '${cols.tokenize(' ').join('\t')}' > mapping/${file}"
}

process GIAP_MAPPING_LOCI {
    tag "loci"
    label 'process_single'
    publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(in_keys), path(inputs, stageAs: 'in/*')

    output:
    path 'mapping', emit: mapping
    path 'versions.yml', emit: versions

    script:
    def files = inputs instanceof List ? inputs : [inputs]
    def ins = [in_keys, files].transpose().collect { k, f -> "--in ${k}=${f}" }.join(' ')
    def resolved = params.config ?: 'resolved_config.yaml'
    """
    giap _step loci --config ${resolved} ${ins} --out mapping=mapping

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        giap: \$(giap version --json | python -c 'import json,sys; print(json.load(sys.stdin)["giap"])')
    END_VERSIONS
    """

    stub:
    def peptide = 'peptide peptide_len source best_q_input pep_input score_input contig start end strand frame crosses_origin route n_diff diff_types placement_q length_class locus_class uniqueness gene novel_candidate_id flags'
    def protein = 'protein_id protein_len_aa route contig start end strand rank identity query_coverage reciprocal_overlap class gene flags'
    def unmapped = 'sequence_id sequence kind label reason contaminant_match'
    def fdr = 'length_class n_diff targets decoys pi_wrong fdr_stratum n_accepted decoy_type section score_bin precision ci_low ci_high c_value flags'
    """
    mkdir mapping
    ${tsvHeader(peptide, 'peptide_hits.tsv')}
    ${tsvHeader(protein, 'protein_hits.tsv')}
    ${tsvHeader(unmapped, 'unmapped.tsv')}
    ${tsvHeader(fdr, 'mapping_fdr.tsv')}
    touch mapping/mappings.probed mapping/mappings.bed mapping/mappings.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        giap: stub
    END_VERSIONS
    """
}
