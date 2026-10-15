// Real-tool entry for one process (scripts/slurm/nf_process_real.sbatch); inputs from the recordings.
include { OSTIR } from '../../modules/local/ostir/main'
workflow {
    OSTIR(channel.of([params.sample, file(params.csv)]))
}
