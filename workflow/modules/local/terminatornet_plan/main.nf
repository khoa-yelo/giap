// Sharded TerminatorNet, step 1 (M2-15b): genome.fa -> plan.tsv, one row per 250,000-window core
// of each contig strand (workflow/bin/gippa_terminatornet_shard plan).
process TERMINATORNET_PLAN {
    tag "${meta_id}"
    label 'process_single'
    conda "${projectDir}/envs/terminatornet.yml"

    input:
    tuple val(meta_id), path(genome, stageAs: 'genome.fa')

    output:
    path 'plan.tsv', emit: plan

    script:
    """
    export TERMINATORNET_HOME=\${TERMINATORNET_HOME:-\$CONDA_PREFIX/share/terminatornet}
    gippa_terminatornet_shard plan
    """

    stub:
    def n = params.stub.terminatornet_chunks ?: 1
    def rows = (0..<n).collect { k -> "c00000p${String.format('%04d', k)}\t0\tgenome\t+\t0\t${k}\t${k + 1}\tshort" }
    """
    printf '%s\\n' 'chunk_id\tseq_index\tseq_name\tstrand\tseq_length\tcore_start\tcore_end\tmode' ${rows.collect { "'${it}'" }.join(' ')} > plan.tsv
    """
}
