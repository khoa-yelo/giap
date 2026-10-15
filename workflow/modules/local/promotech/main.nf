// Promotech RF-HOT on one chunk of 40-nt windows (GIAP_SPEC.md §7.1.1, §7.1.4; M2-13).
process PROMOTECH {
    tag "${meta_id}"
    label 'process_medium'
    conda "${projectDir}/envs/promotech.yml"

    input:
    tuple val(meta_id), path(chunk)   // meta_id = promotech_NNN for windows_NNN.fa

    output:
    path "${meta_id}", emit: dir
    path 'versions.yml', emit: versions

    script:
    """
    PROMOTECH_HOME=\${PROMOTECH_HOME:-\$CONDA_PREFIX/share/promotech}
    export OMP_NUM_THREADS=${task.cpus}   # no thread flag in -s mode; caps the numeric libraries
    python \$PROMOTECH_HOME/promotech.py -s -m RF-HOT -f ${chunk} -o ${meta_id}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        promotech: \$(cat \$PROMOTECH_HOME/GIAP_COMMIT)
    END_VERSIONS
    """

    stub:
    """
    mkdir ${meta_id}
    printf 'sequence,score\\n' > ${meta_id}/sequence_predictions.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        promotech: stub
    END_VERSIONS
    """
}
