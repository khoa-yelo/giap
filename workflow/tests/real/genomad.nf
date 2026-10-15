// Real-tool entry for one process (scripts/slurm/nf_process_real.sbatch); inputs from the recordings.
include { GENOMAD } from '../../modules/local/genomad/main'
workflow {
    GENOMAD(channel.of([params.sample, file(params.genome)]), params.db_dir)
}
