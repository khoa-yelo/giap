// Real-tool entry for one process (scripts/slurm/nf_process_real.sbatch); inputs from the recordings.
include { ANTISMASH } from '../../modules/local/antismash/main'
workflow {
    ANTISMASH(channel.of([params.sample, file(params.gbff)]), params.db_dir)
}
