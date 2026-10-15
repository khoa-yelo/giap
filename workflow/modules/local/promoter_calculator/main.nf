// Promoter Calculator on one TU upstream region, -j 1 (GIAP_SPEC.md §7.1.1, G11; M2-14).
process PROMOTER_CALCULATOR {
    tag "${meta_id}"
    label 'process_single'
    conda "${projectDir}/envs/promoter_calculator.yml"

    input:
    tuple val(meta_id), path(region)   // meta_id = TU id

    output:
    path "${meta_id}.csv", emit: csv
    path 'versions.yml', emit: versions

    script:
    """
    promoter-calculator -i ${region} -o ${meta_id}.csv -j 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        promoter_calculator: \$(python -c "import importlib.metadata as m; print(m.version('promotercalculator'))")
    END_VERSIONS
    """

    stub:
    """
    printf 'TSS,strand,Tx_rate\\n' > ${meta_id}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        promoter_calculator: 1.2.4
    END_VERSIONS
    """
}
