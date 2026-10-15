// Real-tool entry for one process (scripts/slurm/nf_process_real.sbatch); inputs from the recordings.
include { BAKTA } from '../../modules/local/bakta/main'
workflow {
    BAKTA(channel.of([params.sample, file(params.genome), file(params.replicons)]), params.bakta_db, false, [])
}
