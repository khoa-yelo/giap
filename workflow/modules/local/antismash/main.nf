// antiSMASH 8.0.4 [bgc] on the Bakta GBFF, no gene finding (GIPPA_SPEC.md §7.1.1; M2-18).
// The task's shell decides on --databases, so it sees the path as bound in the container (#515).
process ANTISMASH {
    tag "${meta_id}"
    label 'process_high'
    conda "${projectDir}/envs/antismash.yml"

    input:
    tuple val(meta_id), path(gbff)
    val(db_dir)

    output:
    path 'antismash', emit: dir
    path 'versions.yml', emit: versions

    script:
    def gz = gbff.name.toLowerCase().endsWith('.gz')
    def input = gz ? "${meta_id}.gbff" : "${gbff}"
    """
    ${gz ? "gzip -dc ${gbff} > ${input}" : ''}
    dbs=''
    if [ -d "${db_dir}/antismash" ] && [ -n "\$(ls -A "${db_dir}/antismash" 2>/dev/null)" ]; then dbs="--databases ${db_dir}/antismash"; fi
    antismash --genefinding-tool none --minimal \$dbs --cpus ${task.cpus} --output-dir antismash ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(antismash --version 2>/dev/null | sed 's/^antiSMASH //')
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
