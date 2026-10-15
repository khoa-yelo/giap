// Real-tool entry for one process (scripts/slurm/nf_process_real.sbatch); inputs from the recordings.
include { UNIOP } from '../../modules/local/uniop/main'
workflow {
    UNIOP(channel.of([params.sample, file(params.gff)]))
}
