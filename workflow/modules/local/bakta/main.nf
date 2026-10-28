// Bakta 1.12.2 (GIPPA_SPEC.md §7.1.1); flags mirror gippa.annotate.bakta_cmd.bakta_command.
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
    // M11-18: the JSON and GBFF only (hard links; copies where the filesystem has none), passed to CONTEXT_LAYERS as
    // tool_dirs -> OUT/annotation/tools/bakta/. No publishDir: CONTEXT_LAYERS republishes the
    // whole annotation/ directory later and Nextflow's overwrite would delete it.
    path 'tools/bakta', emit: tool_dir
    path 'versions.yml', emit: versions

    script:
    def reg = regions ? "--regions ${regions}" : ''
    def extra = (params.annotation?.bakta?.extra_args ?: []).join(' ')
    def table = params.annotation?.translation_table ?: 11
    def code = table == 11 ? '' : "--translation-table ${table} "  // M11-16 D2: none at 11
    """
    mode="${code}--replicons ${replicons}"
    if [ "${complete}" = true ] || { [ "${complete}" = auto ] && awk -F'\\t' '\$4 != "c" { exit 1 }' ${replicons}; }; then mode="--complete \$mode"; fi
    bakta --db ${db} --output bakta --prefix ${meta_id} --threads ${task.cpus} --keep-contig-headers \$mode ${reg} ${extra} ${genome}
    mkdir -p tools/bakta && { ln bakta/${meta_id}.json bakta/${meta_id}.gbff tools/bakta/ 2>/dev/null || cp bakta/${meta_id}.json bakta/${meta_id}.gbff tools/bakta/; }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(bakta --version 2>/dev/null | sed 's/^bakta //')
    END_VERSIONS
    """

    stub:
    """
    mkdir bakta
    cp ${projectDir}/assets/stub/bakta.json bakta/${meta_id}.json
    touch bakta/${meta_id}.gff3 bakta/${meta_id}.gbff bakta/${meta_id}.faa bakta/${meta_id}.tsv
    mkdir -p tools/bakta && { ln bakta/${meta_id}.json bakta/${meta_id}.gbff tools/bakta/ 2>/dev/null || cp bakta/${meta_id}.json bakta/${meta_id}.gbff tools/bakta/; }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: 1.12.2
    END_VERSIONS
    """
}
