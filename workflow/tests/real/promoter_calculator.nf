// Real-tool entry for one process (scripts/slurm/nf_process_real.sbatch); one task per TU region.
include { PROMOTER_CALCULATOR } from '../../modules/local/promoter_calculator/main'
workflow {
    PROMOTER_CALCULATOR(channel.fromPath("${params.regions}/*.fa").map { f -> [f.baseName, f] })
}
