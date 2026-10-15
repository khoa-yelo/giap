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
    # no promoter in the region: the tool writes no CSV ("No promotors found."); emit its real header
    [ -e ${meta_id}.csv ] || printf '%s\\r\\n' 'promoter_sequence,TSS,UP,hex35,spacer,hex10,disc,ITR,dG_total,dG_10,dG_35,dG_disc,dG_ITR,dG_ext10,dG_spacer,dG_bind,dG_UP,Tx_rate,UP_position,hex35_position,spacer_position,hex10_position,disc_position,fraction,strand,drop,length' > ${meta_id}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        promoter_calculator: \$(python -c "import importlib.metadata as m; print(m.version('promotercalculator'))")
    END_VERSIONS
    """

    stub:
    """
    printf 'promoter_sequence,TSS,UP,hex35,spacer,hex10,disc,ITR,dG_total,dG_10,dG_35,dG_disc,dG_ITR,dG_ext10,dG_spacer,dG_bind,dG_UP,Tx_rate,UP_position,hex35_position,spacer_position,hex10_position,disc_position,fraction,strand,drop,length\\n' > ${meta_id}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        promoter_calculator: 1.2.4
    END_VERSIONS
    """
}
