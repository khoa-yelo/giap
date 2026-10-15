// antiSMASH 8.0.4 [bgc] on the Bakta GBFF, no gene finding (GIAP_SPEC.md §7.1.1; M2-18).
process ANTISMASH {
    tag "${meta_id}"
    label 'process_high'
    memory { 32.GB * task.attempt }
    conda "${projectDir}/envs/antismash.yml"

    input:
    tuple val(meta_id), path(gbff)
    val(db_dir)

    output:
    path 'antismash', emit: dir
    path 'versions.yml', emit: versions

    script:
    """
    antismash --genefinding-tool none --minimal --databases ${db_dir}/antismash --cpus ${task.cpus} --output-dir antismash ${gbff}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(antismash --version 2>&1 | sed 's/^antiSMASH //')
    END_VERSIONS
    """

    stub:
    """
    mkdir antismash
    echo '{"records": []}' > antismash/${meta_id}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: 8.0.4
    END_VERSIONS
    """
}
