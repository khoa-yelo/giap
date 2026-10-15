// Real-tool entry for one process (scripts/slurm/nf_process_real.sbatch); inputs from the recordings.
include { TRANSTERMHP } from '../../modules/local/transtermhp/main'
workflow {
    TRANSTERMHP(channel.of([params.sample, file(params.genome), file(params.coords)]))
}
