// Real-tool entry for one process (scripts/slurm/nf_process_real.sbatch); inputs from the recordings.
// --windows is one chunk or a glob of windows_NNN.fa chunks (named as in subworkflows/local/annotate.nf).
include { PROMOTECH } from '../../modules/local/promotech/main'
workflow {
    PROMOTECH(channel.fromPath(params.windows).map { f -> [f.baseName.replace('windows', 'promotech'), f] }, params.db_dir)
}
