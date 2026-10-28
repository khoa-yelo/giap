// Sharded TerminatorNet scan (M2-15b, docs/specs/M2-15b-terminatornet-sharding.md; ADR-0062
// amendment 1): PLAN -> one CHUNK per plan row -> MERGE, byte-identical to the unsharded
// TERMINATORNET module at -n >= 2. annotate.nf runs it when params.terminatornet_sharded is true.
include { TERMINATORNET_PLAN } from '../../modules/local/terminatornet_plan/main'
include { TERMINATORNET_CHUNK } from '../../modules/local/terminatornet_chunk/main'
include { TERMINATORNET_MERGE } from '../../modules/local/terminatornet_merge/main'

workflow TERMINATORNET_SHARDED {
    take:
    gin   // [meta_id, genome.fa]

    main:
    TERMINATORNET_PLAN(gin)
    def plan = TERMINATORNET_PLAN.out.plan
    def ids = plan.splitCsv(header: true, sep: '\t').map { it.chunk_id }
    TERMINATORNET_CHUNK(ids.combine(gin.map { it[1] }).combine(plan))
    // MERGE reads the pickles in plan order, so the collected order does not matter
    TERMINATORNET_MERGE(gin.combine(plan), TERMINATORNET_CHUNK.out.calls.collect())

    emit:
    tsv      = TERMINATORNET_MERGE.out.tsv
    versions = TERMINATORNET_MERGE.out.versions
}
