// TerminatorNet, genome-wide (GIAP_SPEC.md §7.1.1; M2-15). Its parser matches flags by substring:
// pass exactly -f -s -n -o (always -s); the genome is staged as genome.fa (no spaces).
process TERMINATORNET {
    tag "${meta_id}"
    label 'process_medium'
    conda "${projectDir}/envs/terminatornet.yml"

    input:
    tuple val(meta_id), path(genome, stageAs: 'genome.fa')

    output:
    path 'terminatornet.tsv', emit: tsv
    path 'versions.yml', emit: versions

    script:
    """
    TERMINATORNET_HOME=\${TERMINATORNET_HOME:-\$CONDA_PREFIX/share/terminatornet}
    python \$TERMINATORNET_HOME/TerminatorNet.py -f genome.fa -s -n ${task.cpus} -o terminatornet.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        terminatornet: \$(cat \$TERMINATORNET_HOME/GIAP_COMMIT)
    END_VERSIONS
    """

    stub:
    """
    echo 'NO TERMINATORS FOUND' > terminatornet.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        terminatornet: stub
    END_VERSIONS
    """
}
