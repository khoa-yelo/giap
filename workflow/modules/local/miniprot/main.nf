// miniprot 0.18 (MIT): route 2a protein alignment (GIAP_SPEC.md §7.2.3, §7.1.1; M3-12).
// Mirrors giap.mapping.miniprot.miniprot_index_cmd / miniprot_align_cmd. §7.1.1 gotcha: -T 11 is
// an indexing option; the alignment options (mapping.proteins.miniprot_args) arrive through
// task.ext.args (set in conf/base.config from params.mapping). Container pin: conf/containers.config.
process MINIPROT {
    tag "${meta_id}"
    label 'process_medium'
    conda "${projectDir}/envs/miniprot.yml"

    input:
    val(meta_id)
    path(genome)
    path(proteins)

    output:
    path 'miniprot.gff', emit: gff
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    miniprot -t ${task.cpus} -T 11 -d genome.mpi ${genome}
    miniprot -t ${task.cpus} ${args} --gff genome.mpi ${proteins} > miniprot.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        miniprot: \$(miniprot --version 2>/dev/null)
    END_VERSIONS
    """

    stub:
    """
    printf '##gff-version 3\\n' > miniprot.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        miniprot: stub
    END_VERSIONS
    """
}
