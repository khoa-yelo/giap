// ANNOTATE: the annotation track (GIAP_SPEC.md §6.1, §6.3, §6.4; M2-index step contract; M2-07).
// py steps run through GIAP_STEP aliases; a non-null annotation.precomputed.<key> replaces its
// tool(s) and parse step; geNOMAD / antiSMASH run only with the mobile / bgc module.
include { GIAP_STEP as PARSE_BAKTA; GIAP_STEP as UNIOP_INPUT; GIAP_STEP as OPERON_PAIRS;
          GIAP_STEP as PROMOTER_WINDOWS; GIAP_STEP as PROMOTER_LAYER; GIAP_STEP as TX_LAYER;
          GIAP_STEP as TRANSTERM_INPUT; GIAP_STEP as TERMINATOR_LAYER; GIAP_STEP as OSTIR_INPUT;
          GIAP_STEP as RBS_LAYER; GIAP_STEP as CODON_USAGE; GIAP_STEP as REGION_LAYER;
          GIAP_STEP as CONTEXT_LAYERS } from '../../modules/local/giap_step/main'
include { BAKTA } from '../../modules/local/bakta/main'
include { UNIOP } from '../../modules/local/uniop/main'
include { PROMOTECH } from '../../modules/local/promotech/main'
include { PROMOTER_CALCULATOR } from '../../modules/local/promoter_calculator/main'
include { TERMINATORNET } from '../../modules/local/terminatornet/main'
include { TRANSTERMHP } from '../../modules/local/transtermhp/main'
include { BACTERMFINDER } from '../../modules/local/bactermfinder/main'
include { OSTIR } from '../../modules/local/ostir/main'
include { GENOMAD } from '../../modules/local/genomad/main'
include { ANTISMASH } from '../../modules/local/antismash/main'

// Channels of [key, path | list] → GIAP_STEP inputs: files = [step, paths], keys aligned (key order).
// Emits nothing unless every key of `need` arrived (a step gated on an optional input, e.g. proteins).
def stepIn(String step, List chans, List need = []) {
    chans.inject(channel.empty()) { acc, c -> acc.mix(c) }.toList()
        .filter { kvs -> kvs.collect { kv -> kv[0] }.containsAll(need) }
        .map { kvs -> kvs.sort { kv -> kv[0] }.collectMany { kv -> (kv[1] instanceof List ? kv[1] : [kv[1]]).collect { f -> [kv[0], f] } } }
        .multiMap { kvs ->
            files: [step, kvs.collect { kv -> kv[1] }]
            keys: kvs.collect { kv -> kv[0] }
        }
}

def kv(String key, ch) { ch.map { f -> [key, f] } }

def pick(ch, String name) { ch.map { s, fs -> (fs instanceof List ? fs : [fs]).find { f -> f.name == name } } }

// collect(sort: true) orders by full work path (task hash); re-sort by file name for stable inputs.
def byName(ch) { ch.map { fs -> fs.toSorted { a, b -> a.name <=> b.name } } }

def pre(String key) { kv(key, channel.value(file(params.annotation.precomputed[key], checkIfExists: true))) }

// Stub runs (-stub-run) emit params.stub.<n> chunk / region files; real runs glob the step's output,
// or (dir = true, glob 'D/pattern') declare the output directory D itself as the key's path.
def scatterSpec(String key, String glob, int n, boolean dir) {
    workflow.stubRun ? (0..<n).collect { i -> "${key}=${glob.replace('*', String.format('%03d', i))}" }
        : [dir ? "${key}=${glob.tokenize('/')[0]}" : "${key}=${glob}"]
}

// Files of an output directory (real run) or the directory's files themselves (stub run).
def inDir(ch, String dir) {
    ch.map { s, fs ->
        (fs instanceof List ? fs : [fs]).collectMany { f ->
            f.name == dir && f.isDirectory() ? f.listFiles().toList() : f.parent.name == dir ? [f] : []
        }.toSorted { a, b -> a.name <=> b.name }
    }.flatten()
}

workflow ANNOTATE {
    take:
    genome      // value channel: normalised genome FASTA
    replicons   // value channel: replicons.tsv

    main:
    def ann = params.annotation
    def pc = ann.precomputed
    def db = params.db_dir ?: System.getenv('GIAP_DB_DIR') ?: "${launchDir}/db"
    def gid = genome.map { f -> f.simpleName }
    def gin = genome.map { f -> [f.simpleName, f] }
    def layers = []

    def bjson = channel.empty()
    def gbff = null   // null: no GBFF (precomputed bakta_json alone) -> antiSMASH skipped
    if (pc.bakta_json) {
        bjson = channel.value(file(pc.bakta_json, checkIfExists: true))
        def g = file(pc.bakta_json.toString().replaceFirst(/\.json$/, '.gbff'))
        gbff = g.exists() ? channel.value(g) : null
    } else {
        def light = ann.bakta.db == 'light' || (ann.bakta.db == 'auto' && !file("${db}/bakta/db").exists())
        def bdb = light ? "${db}/bakta/db-light" : ann.bakta.db in ['auto', 'full'] ? "${db}/bakta/db" : ann.bakta.db
        def rg = params.inputs?.regions_gff3 ? file(params.inputs.regions_gff3, checkIfExists: true) : []
        BAKTA(gin.combine(replicons), bdb, ann.bakta.complete, rg)
        bjson = BAKTA.out.json
        gbff = BAKTA.out.gbff
    }
    def pb = stepIn('parse_bakta', [kv('bakta_json', bjson), kv('genome', genome)])
    PARSE_BAKTA(pb.files, pb.keys, ['genes=genes.parquet', 'features=features.parquet',
        'replicons=replicons.parquet', 'rrna16s=rrna16s.fa', 'trna_counts=trna_counts.tsv',
        'xrefs=xrefs.parquet', 'annotation_meta=annotation_meta.json'])
    def P = PARSE_BAKTA.out.outputs
    def genes = kv('genes', pick(P, 'genes.parquet'))
    def reps = kv('replicons', pick(P, 'replicons.parquet'))
    def core = [genes, kv('features', pick(P, 'features.parquet')), reps]

    // M2-12 registry: operon_pairs always gets uniop_cds / uniop_map; a precomputed operon_pairs
    // layer (read by the step from the config) replaces only UniOP itself.
    def ui = stepIn('uniop_input', [genes, reps])
    UNIOP_INPUT(ui.files, ui.keys, ['uniop_cds=uniop_cds.gff', 'uniop_map=uniop_map.tsv'])
    def U = UNIOP_INPUT.out.outputs
    def cds = pick(U, 'uniop_cds.gff')
    def pairs = [kv('uniop_cds', cds), kv('uniop_map', pick(U, 'uniop_map.tsv'))]
    if (!pc.operon_pairs) {
        UNIOP(gid.combine(cds))
        pairs += [kv('uniop_pred', UNIOP.out.pred)]
    }
    def op = stepIn('operon_pairs', [genes, reps, kv('genome', genome)] + pairs)
    OPERON_PAIRS(op.files, op.keys, ['pairs_raw=pairs_raw.parquet', 'tus=tus.parquet', 'operon_meta=operon_meta.json'] +
        scatterSpec('tu_upstream', 'tu_upstream/tu*.fa', params.stub?.tu_regions ?: 1, true))
    def O = OPERON_PAIRS.out.outputs
    def tus = kv('tus', pick(O, 'tus.parquet'))
    layers += [kv('pairs_raw', pick(O, 'pairs_raw.parquet')), tus]

    if (pc.promoters) {
        layers += [pre('promoters')]
    } else {
        def pw = stepIn('promoter_windows', [genes, reps, kv('genome', genome)])
        PROMOTER_WINDOWS(pw.files, pw.keys, ['windows_map=windows_map.parquet'] +
            scatterSpec('windows', 'windows/windows_*.fa', params.stub?.windows_chunks ?: 1, true))
        def W = PROMOTER_WINDOWS.out.outputs
        def chunks = inDir(W, 'windows')
        PROMOTECH(chunks.map { f -> [f.baseName.replace('windows', 'promotech'), f] }, db)
        def pl = stepIn('promoter_layer', [kv('windows_map', pick(W, 'windows_map.parquet')),
            kv('promotech', byName(PROMOTECH.out.dir.collect(sort: true)))])
        PROMOTER_LAYER(pl.files, pl.keys, ['promoter_windows=promoter_windows.parquet'])
        layers += [kv('promoter_windows', pick(PROMOTER_LAYER.out.outputs, 'promoter_windows.parquet'))]
    }

    if (pc.tx_rates) {
        layers += [pre('tx_rates')]
    } else if (ann.promoter_calculator.enabled) {
        def regions = inDir(O, 'tu_upstream')
        PROMOTER_CALCULATOR(regions.map { f -> [f.baseName, f] })
        def tx = stepIn('tx_layer', [kv('pc', byName(PROMOTER_CALCULATOR.out.csv.collect(sort: true))), tus])
        TX_LAYER(tx.files, tx.keys, ['tx_rates=tx_rates.parquet'])
        layers += [kv('tx_rates', pick(TX_LAYER.out.outputs, 'tx_rates.parquet'))]
    }

    if (pc.terminators) {
        layers += [pre('terminators')]
    } else {
        def term = ann.terminators
        def tools = []
        if (term.tool == 'terminatornet') {
            TERMINATORNET(gin)
            tools += [kv('terminatornet', TERMINATORNET.out.tsv)]
        }
        if (term.tool == 'transtermhp' || term.cross_check_transtermhp) {
            def ti = stepIn('transterm_input', [genes])
            TRANSTERM_INPUT(ti.files, ti.keys, ['genes_coords=genes.coords'])
            TRANSTERMHP(gin.combine(pick(TRANSTERM_INPUT.out.outputs, 'genes.coords')))
            tools += [kv('transtermhp', TRANSTERMHP.out.tt)]
        }
        if (term.bactermfinder) {
            BACTERMFINDER(gin)
            tools += [kv('bactermfinder', BACTERMFINDER.out.dir)]
        }
        def tl = stepIn('terminator_layer', [reps] + tools)
        TERMINATOR_LAYER(tl.files, tl.keys, ['terminators=terminators.parquet', 'terminators_crosscheck=terminators_crosscheck.tsv'])
        layers += [kv('terminators', pick(TERMINATOR_LAYER.out.outputs, 'terminators.parquet'))]
    }

    if (pc.rbs) {
        layers += [pre('rbs')]
    } else {
        def oi = stepIn('ostir_input', [genes, kv('annotation_meta', pick(P, 'annotation_meta.json')), reps, kv('genome', genome)])
        OSTIR_INPUT(oi.files, oi.keys, ['ostir_input=ostir_input.csv'])
        def oin = pick(OSTIR_INPUT.out.outputs, 'ostir_input.csv')
        OSTIR(gid.combine(oin))
        def rl = stepIn('rbs_layer', [kv('ostir', OSTIR.out.csv), kv('ostir_input', oin)])
        RBS_LAYER(rl.files, rl.keys, ['rbs=rbs.parquet'])
        layers += [kv('rbs', pick(RBS_LAYER.out.outputs, 'rbs.parquet'))]
    }

    def cu = stepIn('codon_usage', [genes, kv('trna_counts', pick(P, 'trna_counts.tsv'))])
    CODON_USAGE(cu.files, cu.keys, ['codon=codon.parquet'])
    layers += [kv('codon', pick(CODON_USAGE.out.outputs, 'codon.parquet'))]

    if (pc.regions) {
        layers += [pre('regions')]
    } else {
        def optional = []
        if ('mobile' in params.modules) {
            GENOMAD(gin, db)
            optional += [kv('genomad', GENOMAD.out.dir)]
        }
        if ('bgc' in params.modules && gbff != null) {
            ANTISMASH(gid.combine(gbff), db)
            optional += [kv('antismash', ANTISMASH.out.dir)]
        } else if ('bgc' in params.modules) {
            println "WARN: layer_missing:antismash: precomputed bakta_json ${pc.bakta_json} has no GBFF next to it; antiSMASH skipped"
        }
        def rg = stepIn('region_layer', core + optional)
        REGION_LAYER(rg.files, rg.keys, ['regions=regions.parquet', 'region_genes=region_genes.parquet'])
        def R = REGION_LAYER.out.outputs
        layers += [kv('regions', pick(R, 'regions.parquet')), kv('region_genes', pick(R, 'region_genes.parquet'))]
    }

    def cl = stepIn('context_layers', core + layers)
    CONTEXT_LAYERS(cl.files, cl.keys, ['layers=layers.parquet', 'annotation=annotation'])

    emit:
    context    = CONTEXT_LAYERS.out.outputs
    gene_model = pick(P, 'genes.parquet').combine(pick(P, 'features.parquet')).combine(pick(P, 'replicons.parquet'))
}
