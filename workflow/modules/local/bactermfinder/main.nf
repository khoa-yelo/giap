// BacTermFinder, opt-in (annotation.terminators.bactermfinder; GIAP_SPEC.md §7.1.1; M2-15).
// Command per the README at the pinned commit (ADR-0008).
process BACTERMFINDER {
    tag "${meta_id}"
    label 'process_medium'
    conda "${projectDir}/envs/bactermfinder.yml"

    input:
    tuple val(meta_id), path(genome, stageAs: 'genome.fa')

    output:
    path 'bactermfinder', emit: dir
    path 'versions.yml', emit: versions

    script:
    """
    BACTERMFINDER_HOME=\${BACTERMFINDER_HOME:-\$CONDA_PREFIX/share/bactermfinder}
    python \$BACTERMFINDER_HOME/genome_scan.py genome.fa bactermfinder

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bactermfinder: \$(cat \$BACTERMFINDER_HOME/GIAP_COMMIT)
    END_VERSIONS
    """

    stub:
    """
    mkdir bactermfinder

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bactermfinder: stub
    END_VERSIONS
    """
}
