// Bakta 1.12.2 (GIAP_SPEC.md §7.1.1); flags mirror giap.annotate.bakta_cmd.bakta_command.
process BAKTA {
    tag "${meta_id}"
    label 'process_high'
    memory { 48.GB * task.attempt }
    conda "${projectDir}/envs/bakta.yml"

    input:
    tuple val(meta_id), path(genome), path(replicons)
    val(db)
    val(complete)   // true, false or 'auto' (--complete iff every replicon is circular)
    path(regions)   // inputs.regions_gff3 (staged, part of the task hash), or [] for none

    output:
    path "bakta/${meta_id}.json", emit: json
    path "bakta/${meta_id}.gbff", emit: gbff
    path 'versions.yml', emit: versions

    script:
    def reg = regions ? "--regions ${regions}" : ''
    def extra = (params.annotation?.bakta?.extra_args ?: []).join(' ')
    """
    mode="--replicons ${replicons}"
    if [ "${complete}" = true ] || { [ "${complete}" = auto ] && awk -F'\\t' '\$4 != "c" { exit 1 }' ${replicons}; }; then mode=--complete; fi
    bakta --db ${db} --output bakta --prefix ${meta_id} --threads ${task.cpus} --keep-contig-headers \$mode ${reg} ${extra} ${genome}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(bakta --version 2>&1 | sed 's/^bakta //')
    END_VERSIONS
    """

    stub:
    """
    mkdir bakta
    cp ${projectDir}/assets/stub/bakta.json bakta/${meta_id}.json
    touch bakta/${meta_id}.gff3 bakta/${meta_id}.gbff bakta/${meta_id}.faa bakta/${meta_id}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: 1.12.2
    END_VERSIONS
    """
}
