// TransTermHP 2.09 cross-check (GIAP_SPEC.md §7.1.1; M2-15); genes.coords from transterm_input.
process TRANSTERMHP {
    tag "${meta_id}"
    label 'process_single'
    conda "${projectDir}/envs/transtermhp.yml"

    input:
    tuple val(meta_id), path(genome, stageAs: 'genome.fa'), path(coords, stageAs: 'genes.coords')

    output:
    path 'transtermhp.tt', emit: tt
    path 'versions.yml', emit: versions

    script:
    """
    EXPTERM=\${EXPTERM:-\$(dirname \$(command -v transterm))/../data/expterm.dat}
    transterm -p \$EXPTERM genome.fa genes.coords > transtermhp.tt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transtermhp: 2.09
    END_VERSIONS
    """

    stub:
    """
    touch transtermhp.tt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transtermhp: 2.09
    END_VERSIONS
    """
}
