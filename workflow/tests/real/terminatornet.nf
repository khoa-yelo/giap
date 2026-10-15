// Real-tool entry for one process (scripts/slurm/nf_process_real.sbatch); inputs from the recordings.
include { TERMINATORNET } from '../../modules/local/terminatornet/main'
workflow {
    TERMINATORNET(channel.of([params.sample, file(params.genome)]))
}
