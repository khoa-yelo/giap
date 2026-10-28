// Sharded TerminatorNet, step 2 (M2-15b): upstream's window scan over one chunk's own windows on
// the whole cleaned strand sequence -> calls_<chunk_id>.pkl. CPUs: params.terminatornet_chunk_cpus.
process TERMINATORNET_CHUNK {
    tag "${chunk_id}"
    conda "${projectDir}/envs/terminatornet.yml"

    input:
    tuple val(chunk_id), path(genome, stageAs: 'genome.fa'), path(plan)

    output:
    path "calls_${chunk_id}.pkl", emit: calls

    script:
    """
    export TERMINATORNET_HOME=\${TERMINATORNET_HOME:-\$CONDA_PREFIX/share/terminatornet}
    ln -sf \$TERMINATORNET_HOME/model.pkl model.pkl   # loaded from the cwd, as the unsharded module
    gippa_terminatornet_shard.py chunk ${chunk_id} --cpus ${task.cpus}
    """

    stub:
    """
    touch calls_${chunk_id}.pkl
    """
}
