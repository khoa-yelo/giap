// BLAST+ 2.16.0 (public domain): route 2b tblastn for protein leftovers (GIPPA_SPEC.md §7.2.3,
// §6.2 step `tblastn`; M3-13). Mirrors gippa.mapping.tblastn.makeblastdb_cmd / tblastn_cmd. The
// search options (mapping.proteins.tblastn_args) arrive only through task.ext.args (set in
// conf/base.config from params.mapping). Container pin: conf/containers.config.
process TBLASTN {
    tag "${meta_id}"
    label 'process_medium'
    conda "${projectDir}/envs/blast.yml"

    input:
    val(meta_id)
    path(genome)
    path(proteins)

    output:
    path 'tblastn.tsv', emit: hits
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    makeblastdb -in ${genome} -dbtype nucl -out db/genome
    tblastn -db db/genome -query ${proteins} -num_threads ${task.cpus} ${args} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen nident positive' > tblastn.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tblastn: \$(tblastn -version 2>/dev/null | head -n 1 | sed 's/^tblastn: //')
    END_VERSIONS
    """

    stub:
    """
    : > tblastn.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tblastn: stub
    END_VERSIONS
    """
}
