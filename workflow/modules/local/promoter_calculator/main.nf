// Promoter Calculator on a batch of TU upstream regions, -j 1 (GIPPA_SPEC.md §7.1.1, G11; M2-14, M2-14b).
// One task runs each region of the batch in the given (TU-id) order and writes <TU>.csv per region;
// a failing region fails the batch (Nextflow's default task shell, bash -ue; keep it).
process PROMOTER_CALCULATOR {
    tag "${meta_id}"
    label 'process_single'
    time { 2.h * task.attempt }
    conda "${projectDir}/envs/promoter_calculator.yml"

    input:
    tuple val(meta_id), path(regions)   // meta_id = pc_<first TU>_<last TU>; regions = <TU>.fa in TU-id order

    output:
    path "*.csv", emit: csv
    path 'versions.yml', emit: versions

    script:
    """
    for r in ${regions}; do
        tu="\${r%.*}"
        promoter-calculator -i "\${r}" -o "\${tu}.csv" -j 1
        # no promoter in the region: the tool writes no CSV ("No promotors found."); emit its real header
        [ -e "\${tu}.csv" ] || printf '%s\\r\\n' 'promoter_sequence,TSS,UP,hex35,spacer,hex10,disc,ITR,dG_total,dG_10,dG_35,dG_disc,dG_ITR,dG_ext10,dG_spacer,dG_bind,dG_UP,Tx_rate,UP_position,hex35_position,spacer_position,hex10_position,disc_position,fraction,strand,drop,length' > "\${tu}.csv"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        promoter_calculator: \$(python -c "import importlib.metadata as m; print(m.version('promotercalculator'))")
    END_VERSIONS
    """

    stub:
    """
    for r in ${regions}; do
        tu="\${r%.*}"
        printf 'promoter_sequence,TSS,UP,hex35,spacer,hex10,disc,ITR,dG_total,dG_10,dG_35,dG_disc,dG_ITR,dG_ext10,dG_spacer,dG_bind,dG_UP,Tx_rate,UP_position,hex35_position,spacer_position,hex10_position,disc_position,fraction,strand,drop,length\\n' > "\${tu}.csv"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        promoter_calculator: 1.2.4
    END_VERSIONS
    """
}
