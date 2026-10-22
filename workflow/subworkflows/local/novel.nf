// NOVEL subworkflow (GIAP_SPEC.md §6.2, §7.6.2, §9.1; M7-08): MAPPING's loci tables -> candidate
// novel ORFs and their second-round windows (novel_candidates), Promotech RF-HOT on each window
// chunk (PROMOTECH_NOVEL, one task per chunk as M2-13) -> novel_promoter_layer, OSTIR on the
// multi-start rows (OSTIR_NOVEL) -> novel_rbs_layer, then novel_orfs -> the four §9.1 files under
// tables/. The py steps run through GIAP_STEP aliases; both tools reuse their modules and envs.
// A non-null annotation.precomputed.promoters / .rbs_novel replaces the tool and its layer step
// (§6.4, ADR-0024), as in annotate.nf. With no candidate start the window directory is empty and
// the OSTIR CSV header-only: neither tool runs, nor its layer step (rule 3). NOVEL runs only with
// a gene table (lead decision 3). Stub runs: params.stub.novel_windows_chunks window files (0 =
// no candidate start).
include { GIAP_STEP as NOVEL_CANDIDATES; GIAP_STEP as NOVEL_PROMOTER_LAYER;
          GIAP_STEP as NOVEL_RBS_LAYER; GIAP_STEP as NOVEL_ORFS } from '../../modules/local/giap_step/main'
include { PROMOTECH as PROMOTECH_NOVEL } from '../../modules/local/promotech/main'
include { OSTIR as OSTIR_NOVEL } from '../../modules/local/ostir/main'
include { stepIn; kv; pick; byName; inDir; scatterSpec } from './annotate'
include { mappingDir } from './genetable'

workflow NOVEL {
    take:
    genome
    gene_model
    contract
    mapping
    peptides
    contaminants
    gene_table
    layers

    main:
    // take: genome FASTA; [genes, features, replicons]; [xrefs, annotation_meta, tus, layers];
    // LOCI outputs; peptides.norm.parquet; contaminant FASTA or empty; genetable/genes.parquet
    // or empty; ANNOTATE's [key, path] layer pairs (promoter_windows, terminators are used)
    def pc = params.annotation.precomputed
    def db = params.db_dir ?: System.getenv('GIAP_DB_DIR') ?: "${launchDir}/db"
    def n_chunks = params.stub?.novel_windows_chunks != null ? params.stub.novel_windows_chunks as int : 1
    def g = kv('genome', genome)
    def genes = kv('genes', gene_model.map { a -> a[0] })
    def reps = kv('replicons', gene_model.map { a -> a[2] })
    def loci = mappingDir(mapping).combine(gene_table).map { m, t -> m }   // no gene table: no NOVEL

    def nc = stepIn('novel_candidates', [kv('mapping', loci), kv('peptides', peptides), genes,
        kv('features', gene_model.map { a -> a[1] }), reps, g,
        kv('annotation_meta', contract.map { c -> c[1] }), kv('contaminants', contaminants)],
        ['mapping', 'peptides', 'genes', 'features', 'replicons', 'genome', 'annotation_meta'])
    NOVEL_CANDIDATES(nc.files, nc.keys, ['novel_candidates=novel_candidates.parquet',
        'novel_windows_map=novel_windows_map.parquet', 'ostir_input_novel=ostir_input_novel.csv'] +
        scatterSpec('novel_windows', 'novel_windows/windows_*.fa', n_chunks, true))
    def C = NOVEL_CANDIDATES.out.outputs
    def cand = kv('novel_candidates', pick(C, 'novel_candidates.parquet'))

    def npw = channel.empty()
    if (pc.promoters) {
        npw = kv('novel_promoter_windows', channel.value(file(pc.promoters, checkIfExists: true)))
    } else {
        PROMOTECH_NOVEL(inDir(C, 'novel_windows').map { f -> [f.baseName.replace('windows', 'promotech'), f] }, db)
        def pl = stepIn('novel_promoter_layer', [kv('novel_windows_map', pick(C, 'novel_windows_map.parquet')),
            kv('promotech_novel', byName(PROMOTECH_NOVEL.out.dir.collect(sort: true)))],
            ['novel_windows_map', 'promotech_novel'])
        NOVEL_PROMOTER_LAYER(pl.files, pl.keys, ['novel_promoter_windows=novel_promoter_windows.parquet'])
        npw = kv('novel_promoter_windows', pick(NOVEL_PROMOTER_LAYER.out.outputs, 'novel_promoter_windows.parquet'))
    }

    def rbn = channel.empty()
    if (pc.rbs_novel) {
        rbn = kv('rbs_novel', channel.value(file(pc.rbs_novel, checkIfExists: true)))
    } else {
        // a header-only CSV (no candidate start) gives no OSTIR task
        def oin = pick(C, 'ostir_input_novel.csv').filter { f ->
            workflow.stubRun ? n_chunks > 0 : f.withReader { r -> r.readLine(); r.readLine() != null }
        }
        OSTIR_NOVEL(oin.map { f -> ['novel', f] })
        def rl = stepIn('novel_rbs_layer', [kv('ostir_novel', OSTIR_NOVEL.out.csv),
            kv('ostir_input_novel', oin), cand], ['ostir_novel', 'ostir_input_novel', 'novel_candidates'])
        NOVEL_RBS_LAYER(rl.files, rl.keys, ['rbs_novel=rbs_novel.parquet'])
        rbn = kv('rbs_novel', pick(NOVEL_RBS_LAYER.out.outputs, 'rbs_novel.parquet'))
    }

    def first = layers.filter { p -> p[0] in ['promoter_windows', 'terminators'] }
    def no = stepIn('novel_orfs', [cand, kv('gene_table', gene_table), genes, reps, g, npw, rbn, first],
        ['novel_candidates', 'gene_table', 'genes', 'replicons', 'genome'])
    NOVEL_ORFS(no.files, no.keys, ['novel_orfs=tables/novel_orfs.tsv',
        'novel_orfs_gff3=tables/novel_orfs.gff3', 'novel_orfs_faa=tables/novel_orfs.faa',
        'annotation_patch=tables/annotation_patch.gff3'])
    def O = NOVEL_ORFS.out.outputs

    emit:
    novel_orfs       = pick(O, 'novel_orfs.tsv')
    novel_orfs_gff3  = pick(O, 'novel_orfs.gff3')
    novel_orfs_faa   = pick(O, 'novel_orfs.faa')
    annotation_patch = pick(O, 'annotation_patch.gff3')
}
