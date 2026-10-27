// UniOP operon pair predictions (GIPPA_SPEC.md §7.1.1); input from the uniop_input step (M2-12).
process UNIOP {
    tag "${meta_id}"
    label 'process_low'
    conda "${projectDir}/envs/uniop.yml"

    input:
    tuple val(meta_id), path(gff)

    output:
    path 'uniop.pred', emit: pred
    path 'versions.yml', emit: versions

    script:
    """
    UNIOP_HOME=\${UNIOP_HOME:-\$CONDA_PREFIX/share/uniop}
    mkdir -p out
    python \$UNIOP_HOME/src/UniOP -f ${gff} --gff_source prodigal_gff -t out
    cp out/*.pred uniop.pred

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        uniop: \$(cat \$UNIOP_HOME/GIPPA_COMMIT)
    END_VERSIONS
    """

    stub:
    """
    printf 'operon prediction scores\\nGene A\\tGene B\\tPrediction\\n' > uniop.pred

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        uniop: stub
    END_VERSIONS
    """
}
