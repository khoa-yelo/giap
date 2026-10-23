// Real-tool entry for one process (scripts/slurm/nf_process_real.sbatch); one task per TU region.
include { PROMOTER_CALCULATOR } from '../../modules/local/promoter_calculator/main'
workflow {
    PROMOTER_CALCULATOR(channel.fromPath("${params.regions}/*.fa").toSortedList { a, b -> a.baseName <=> b.baseName }.flatMap { l -> l.collate(100) }.map { b -> ["pc_${b[0].baseName}_${b[-1].baseName}", b] })
}
