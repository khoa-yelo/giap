// geNOMAD end-to-end [mobile] (GIAP_SPEC.md §7.1.1; M2-18); only genomad/*_summary/ is used.
process GENOMAD {
    tag "${meta_id}"
    label 'process_high'
    memory { 48.GB * task.attempt }
    conda "${projectDir}/envs/genomad.yml"

    input:
    tuple val(meta_id), path(genome, stageAs: 'genome.fa')
    val(db_dir)

    output:
    path 'genomad', emit: dir
    path 'versions.yml', emit: versions

    script:
    """
    genomad end-to-end --cleanup --splits 8 --threads ${task.cpus} genome.fa genomad ${db_dir}/genomad_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(genomad --version 2>/dev/null | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p genomad/genome_summary && cd genomad/genome_summary
    printf 'seq_name\\tlength\\ttopology\\tcoordinates\\tn_genes\\tgenetic_code\\tvirus_score\\tfdr\\tn_hallmarks\\tmarker_enrichment\\ttaxonomy\\n' > genome_virus_summary.tsv
    printf 'seq_name\\tlength\\ttopology\\tn_genes\\tgenetic_code\\tplasmid_score\\tfdr\\tn_hallmarks\\tmarker_enrichment\\tconjugation_genes\\tamr_genes\\n' > genome_plasmid_summary.tsv
    printf 'gene\\tstart\\tend\\tlength\\tstrand\\tgc_content\\tgenetic_code\\trbs_motif\\tmarker\\tevalue\\tbitscore\\tuscg\\tplasmid_hallmark\\tvirus_hallmark\\ttaxid\\ttaxname\\tannotation_conjscan\\tannotation_amr\\tannotation_accessions\\tannotation_description\\n' > genome_virus_genes.tsv
    cd -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: 1.12.0
    END_VERSIONS
    """
}
