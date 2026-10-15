// MAPPING subworkflow (GIAP_SPEC.md §6.1, §6.2, §7.2; M3-17): six-frame index -> exact ->
// tolerant (+ decoys) -> placement FDR -> loci; long proteins -> MINIPROT -> leftovers -> TBLASTN.
// Each Python step is one `giap _step` task (keys as in giap.steps.mapping); the protein route runs
// only when proteins_long emits. Outputs: params.outdir/mapping/ (seven §7.2.9 files).
include { GIAP_MAPPING as SIXFRAME_INDEX    } from '../../modules/local/giap_mapping/main'
include { GIAP_MAPPING as MAP_EXACT         } from '../../modules/local/giap_mapping/main'
include { GIAP_MAPPING as MAP_TOLERANT      } from '../../modules/local/giap_mapping/main'
include { GIAP_MAPPING as PLACEMENT         } from '../../modules/local/giap_mapping/main'
include { GIAP_MAPPING as PROTEIN_LEFTOVERS } from '../../modules/local/giap_mapping/main'
include { GIAP_MAPPING_LOCI as LOCI         } from '../../modules/local/giap_mapping/main'
include { MINIPROT                          } from '../../modules/local/miniprot/main'
include { TBLASTN                           } from '../../modules/local/tblastn/main'

// [step, keys, files, out_spec] once every channel of `tagged` ([key, file]) has completed; keys
// sorted. Emits nothing unless every key of `required` is present (optional inputs may be absent).
def gather(String step, List required, List tagged, List outSpec) {
    def ch = tagged[0]
    tagged.drop(1).each { c -> ch = ch.mix(c) }
    return ch
        .toSortedList { a, b -> a[0] <=> b[0] }
        .filter { items -> items*.get(0).containsAll(required) }
        .map { items -> [step, items*.get(0), items*.get(1), outSpec] }
}

// The files of a GIAP_MAPPING output tuple, in out_spec order.
def outFiles(ch) {
    return ch.map { s, f -> [f].flatten() }
}

workflow MAPPING {
    take:
    genome              // genome.norm.fa
    peptides            // peptides.norm.parquet
    proteins_long       // proteins.long.faa, or channel.empty()
    annotation          // [genes.parquet, features.parquet, replicons.parquet]
    reference_peptides  // reference peptide table, or channel.empty()

    main:
    def g    = genome.map { ['genome', it] }
    def p    = peptides.map { ['peptides', it] }
    def gene = annotation.map { ['genes', it[0]] }
    def feat = annotation.map { ['features', it[1]] }
    def rep  = annotation.map { ['replicons', it[2]] }
    def ref  = reference_peptides.map { ['reference_peptides', it] }
    def prot = proteins_long.map { ['proteins_long', it] }

    SIXFRAME_INDEX(gather('sixframe_index', ['genome', 'peptides', 'replicons'], [g, rep, p],
        ['sixframe=sixframe.npz']))
    def six = outFiles(SIXFRAME_INDEX.out.outputs).map { ['sixframe', it[0]] }

    MAP_EXACT(gather('map_exact', ['peptides', 'sixframe'], [p, six], ['hits_exact=hits_exact.parquet']))
    def hx = outFiles(MAP_EXACT.out.outputs).map { ['hits_exact', it[0]] }

    MAP_TOLERANT(gather('map_tolerant', ['hits_exact', 'peptides', 'sixframe'], [p, six, hx],
        ['hits_tolerant=hits_tolerant.parquet']))
    def ht = outFiles(MAP_TOLERANT.out.outputs).map { ['hits_tolerant', it[0]] }

    PLACEMENT(gather('placement', ['hits_exact', 'hits_tolerant', 'peptides', 'sixframe'],
        [p, six, hx, ht, ref], ['placement=placement.parquet', 'mapping_fdr=mapping_fdr.parquet']))
    def placed = outFiles(PLACEMENT.out.outputs)
    def pl  = placed.map { ['placement', it[0]] }
    def fdr = placed.map { ['mapping_fdr', it[1]] }

    // protein route (§7.2: only failures move on)
    MINIPROT(channel.value('genome'), genome, proteins_long)
    def mp = MINIPROT.out.gff.map { ['miniprot_gff', it] }
    PROTEIN_LEFTOVERS(gather('protein_leftovers',
        ['features', 'genes', 'miniprot_gff', 'proteins_long', 'replicons'],
        [prot, mp, gene, feat, rep], ['proteins_leftover=proteins.leftover.faa']))
    def left = outFiles(PROTEIN_LEFTOVERS.out.outputs).map { it[0] }
    TBLASTN(channel.value('genome'), genome, left)
    def tb = TBLASTN.out.hits.map { ['tblastn_tsv', it] }

    def loci_in = gather('loci',
        ['features', 'genes', 'genome', 'hits_exact', 'hits_tolerant', 'mapping_fdr', 'peptides',
         'placement', 'replicons'],
        [g, p, gene, feat, rep, hx, ht, pl, fdr, prot, mp, tb], [])
    LOCI(loci_in.map { step, keys, files, outs -> [keys, files] })

    def ch_versions = SIXFRAME_INDEX.out.versions.mix(MAP_EXACT.out.versions, MAP_TOLERANT.out.versions,
        PLACEMENT.out.versions, PROTEIN_LEFTOVERS.out.versions, MINIPROT.out.versions,
        TBLASTN.out.versions, LOCI.out.versions)

    emit:
    mapping  = LOCI.out.mapping
    versions = ch_versions
}
