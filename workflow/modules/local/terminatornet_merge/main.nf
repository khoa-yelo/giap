// Sharded TerminatorNet, step 3 (M2-15b): every chunk's calls, validated against plan.tsv, through
// upstream merge_overlapping/output_terms -> terminatornet.tsv (= TerminatorNet.py -s -n >= 2).
process TERMINATORNET_MERGE {
    tag "${meta_id}"
    label 'process_single'
    conda "${projectDir}/envs/terminatornet.yml"

    input:
    tuple val(meta_id), path(genome, stageAs: 'genome.fa'), path(plan)
    path(calls)

    output:
    path 'terminatornet.tsv', emit: tsv
    path 'versions.yml', emit: versions

    script:
    """
    export TERMINATORNET_HOME=\${TERMINATORNET_HOME:-\$CONDA_PREFIX/share/terminatornet}
    gippa_terminatornet_shard merge

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        terminatornet: \$(cat \$TERMINATORNET_HOME/GIPPA_COMMIT)
    END_VERSIONS
    """

    stub:
    """
    printf 'genome\\tNO TERMINATORS FOUND\\n' > terminatornet.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        terminatornet: stub
    END_VERSIONS
    """
}
