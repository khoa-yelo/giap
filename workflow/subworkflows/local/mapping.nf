// MAPPING subworkflow (GIAP_SPEC.md §6.1, §6.2, §7.2; M3-17): six-frame index -> exact ->
// tolerant (+ decoys) -> placement FDR -> loci; long proteins -> MINIPROT -> leftovers -> TBLASTN.
// py steps run through GIAP_STEP aliases (M2-07; keys as in giap.steps.mapping); the protein route
// runs only when proteins_long emits. Outputs: params.outdir/mapping/ (seven §7.2.9 files).
include { GIAP_STEP as SIXFRAME_INDEX; GIAP_STEP as MAP_EXACT; GIAP_STEP as MAP_TOLERANT;
          GIAP_STEP as PLACEMENT; GIAP_STEP as PROTEIN_LEFTOVERS;
          GIAP_STEP as LOCI } from '../../modules/local/giap_step/main'
include { MINIPROT } from '../../modules/local/miniprot/main'
include { TBLASTN } from '../../modules/local/tblastn/main'
include { stepIn; kv; pick } from './annotate'

// Stub runs: the loci stub writes the seven §7.2.9 files, the TSVs with their §9.3 headers
// (giap.mapping.outputs, giap.mapping.protein_assign); real runs declare the directory itself.
def lociSpec() {
    def h = [
        peptide_hits: 'peptide peptide_len source best_q_input pep_input score_input contig start end strand frame crosses_origin route n_diff diff_types placement_q length_class locus_class uniqueness gene novel_candidate_id flags',
        protein_hits: 'protein_id protein_len_aa route contig start end strand rank identity query_coverage reciprocal_overlap class gene flags',
        unmapped: 'sequence_id sequence kind label reason contaminant_match',
        mapping_fdr: 'length_class n_diff targets decoys pi_wrong fdr_stratum n_accepted decoy_type section score_bin precision ci_low ci_high c_value flags']
    def tsvs = h.collect { f, cols -> "mapping=mapping/${f}.tsv ${cols.tokenize(' ').join('\t')}".toString() }
    return tsvs + ['probed', 'bed', 'gff3'].collect { e -> "mapping=mapping/mappings.${e}".toString() }
}

workflow MAPPING {
    take:
    genome              // genome.norm.fa
    peptides            // peptides.norm.parquet
    proteins_long       // proteins.long.faa, or channel.empty()
    annotation          // [genes.parquet, features.parquet, replicons.parquet]
    reference_peptides  // reference peptide table, or channel.empty()
    contaminants        // contaminant library FASTA, or channel.empty()

    main:
    def g = kv('genome', genome)
    def pep = kv('peptides', peptides)
    def prot = kv('proteins_long', proteins_long)
    def gm = [kv('genes', annotation.map { a -> a[0] }), kv('features', annotation.map { a -> a[1] }), kv('replicons', annotation.map { a -> a[2] })]

    def si = stepIn('sixframe_index', [g, pep, gm[2]])
    SIXFRAME_INDEX(si.files, si.keys, ['sixframe=sixframe.npz'])
    def six = kv('sixframe', pick(SIXFRAME_INDEX.out.outputs, 'sixframe.npz'))

    def me = stepIn('map_exact', [pep, six])
    MAP_EXACT(me.files, me.keys, ['hits_exact=hits_exact.parquet'])
    def hx = kv('hits_exact', pick(MAP_EXACT.out.outputs, 'hits_exact.parquet'))

    def mt = stepIn('map_tolerant', [pep, six, hx])
    MAP_TOLERANT(mt.files, mt.keys, ['hits_tolerant=hits_tolerant.parquet'])
    def ht = kv('hits_tolerant', pick(MAP_TOLERANT.out.outputs, 'hits_tolerant.parquet'))

    def pl = stepIn('placement', [pep, six, hx, ht, kv('reference_peptides', reference_peptides)] + gm)
    // stratum0_classes = the ADR-0019 stratum-0 pi-hat per length class (an input of context_stats)
    PLACEMENT(pl.files, pl.keys, ['placement=placement.parquet', 'mapping_fdr=mapping_fdr.parquet', 'stratum0_classes=stratum0_classes.parquet'])
    def P = PLACEMENT.out.outputs
    def placed = [kv('placement', pick(P, 'placement.parquet')), kv('mapping_fdr', pick(P, 'mapping_fdr.parquet'))]

    // protein route (§7.2: only miniprot failures move on to tblastn)
    MINIPROT(channel.value('genome'), genome, proteins_long)
    def mp = kv('miniprot_gff', MINIPROT.out.gff)
    def lo = stepIn('protein_leftovers', [prot, mp] + gm, ['proteins_long', 'miniprot_gff'])
    PROTEIN_LEFTOVERS(lo.files, lo.keys, ['proteins_leftover=proteins.leftover.faa'])
    TBLASTN(channel.value('genome'), genome, pick(PROTEIN_LEFTOVERS.out.outputs, 'proteins.leftover.faa'))
    def tb = kv('tblastn_tsv', TBLASTN.out.hits)

    def lc = stepIn('loci', [g, pep, hx, ht, prot, mp, tb, kv('contaminants', contaminants)] + gm + placed)
    LOCI(lc.files, lc.keys, workflow.stubRun ? lociSpec() : ['mapping=mapping'])

    def ch_versions = SIXFRAME_INDEX.out.versions.mix(MAP_EXACT.out.versions, MAP_TOLERANT.out.versions,
        PLACEMENT.out.versions, PROTEIN_LEFTOVERS.out.versions, MINIPROT.out.versions,
        TBLASTN.out.versions, LOCI.out.versions)

    emit:
    mapping  = LOCI.out.outputs
    stratum0 = pick(P, 'stratum0_classes.parquet')   // ADR-0019 table, an input of context_stats
    versions = ch_versions
}
