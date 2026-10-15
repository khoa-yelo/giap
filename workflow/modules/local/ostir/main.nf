// OSTIR 1.1.3 RBS strength (GIAP_SPEC.md §7.1.1; M2-16); input from the ostir_input step.
process OSTIR {
    tag "${meta_id}"
    label 'process_medium'
    conda "${projectDir}/envs/ostir.yml"

    input:
    tuple val(meta_id), path(csv, stageAs: 'ostir_input.csv')

    output:
    path 'ostir.csv', emit: csv
    path 'versions.yml', emit: versions

    script:
    """
    ostir -i ostir_input.csv -o ostir.csv -j ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ostir: \$(python -c "import importlib.metadata as m; print(m.version('ostir'))")
    END_VERSIONS
    """

    stub:
    """
    echo 'name,start_codon,start_position,expression,RBS_distance_bp,dG_total,dG_rRNA:mRNA,dG_mRNA,dG_spacing,dG_standby,dG_start_codon' > ostir.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ostir: 1.1.3
    END_VERSIONS
    """
}
