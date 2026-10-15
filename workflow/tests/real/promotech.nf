// Real-tool entry for one process (scripts/slurm/nf_process_real.sbatch); inputs from the recordings.
include { PROMOTECH } from '../../modules/local/promotech/main'
workflow {
    PROMOTECH(channel.of(['promotech_000', file(params.windows)]), params.db_dir)
}
