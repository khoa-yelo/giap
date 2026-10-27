// BacTermFinder, opt-in (annotation.terminators.bactermfinder; GIPPA_SPEC.md §7.1.1; M2-15).
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
    # genome_scan.py <fasta> <step> <prefix> <batch> (README example) reads its models and
    # iLearnPlus from, and writes scratch files to, its working directory: run it from a
    # task-local directory of links, never inside the shared (read-only) source directory.
    rm -rf btf_work && mkdir btf_work bactermfinder
    ln -s \$BACTERMFINDER_HOME/genome_scan.py \$BACTERMFINDER_HOME/iLearnPlus \$BACTERMFINDER_HOME/*.h5 btf_work/
    ln -s ../genome.fa btf_work/genome.fa
    cd btf_work
    python genome_scan.py genome.fa 3 out 10000 > log.out
    cd ..
    cp btf_work/outgenome.fa_mean.csv bactermfinder/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bactermfinder: \$(cat \$BACTERMFINDER_HOME/GIPPA_COMMIT)
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
