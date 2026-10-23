// VARIANTS subworkflow (GIAP_SPEC.md §6.2, §7.5, §9.1; M7-04): the `variants` py step (Module 5,
// core) reads the loci directory `mapping` (its peptide_hits.tsv), MAP_TOLERANT's hits_tolerant
// (optional: without it frameshift rows are skipped with a run flag), the C2/C1 gene model, the
// genome and the C9 peptides, and publishes tables/pileup.tsv.gz, tables/variants.tsv and
// tables/variants.vcf (no work Parquet). The required keys go to stepIn's `need` explicitly, so
// the step does not run when one of their channels is empty. peptide_spectra (C17, optional,
// M7-10) carries the paired de novo and reference spectra for the class error rates.
include { GIAP_STEP as VARIANTS_STEP } from '../../modules/local/giap_step/main'
include { stepIn; kv; pick } from './annotate'

// The loci output directory: the `mapping` directory itself in a real run, the parent of its
// files in a stub run (the loci stub emits the files of mapping/).
def lociDir(ch) { ch.map { s, fs -> def f = fs instanceof List ? fs[0] : fs; f.name == 'mapping' ? f : f.parent } }

// Stub runs (rule 6): header-only tables with their §9.3 headers (giap.variants.pileup,
// giap.variants.vcf_records) and a header-only VCF with the §7.5.5 header (giap.io.vcf).
def vcfMeta(String kind, List f) { "##${kind}=<ID=${f[0]},Number=${f[1]},Type=${f[2]},Description=\"${f[3]}\">".toString() }

def variantsSpec(List spec) {
    if (!workflow.stubRun) return spec
    def pile = 'contig codon_start strand gene aa_pos ref_aa n_ref n_il alt_counts call'
    def vars = 'variant_id contig codon_start strand gene aa_pos ref_codon ref_aa alt_aa alt_type call confidence n_ref n_alt nt_changes alt_codon_vcf class_error_rate class_error_n source flags'
    def info = [['GENE', '1', 'String', 'Locus tag of the gene framing the codon'],
                ['AA_POS', '1', 'Integer', '1-based residue position of the codon in the gene'],
                ['AA_REF', '1', 'String', 'Reference residue(s) on the gene strand'],
                ['AA_ALT', 'A', 'String', 'Alternative residue(s) observed in peptides (J = I or L)'],
                ['CALL', '1', 'String', 'Codon call: ALT or mixed'],
                ['CONF', 'A', 'String', 'Confidence of each ALT allele: high or candidate'],
                ['NT_CHANGES', 'A', 'Integer', 'Nucleotide changes of the minimal-change ALT codon(s)'],
                ['ALT_TYPE', 'A', 'String', 'ALT type: substitution, swap or block'],
                ['FLAGS', '.', 'String', 'Variant flags of the alleles, deduplicated, in GIAP §9.4 order'],
                ['CLASS_ERR', 'A', 'String', 'Class-specific error rate of each allele, or unestimated']]
    def fmt = [['GT', '1', 'String', 'Haploid genotype: the best-supported ALT allele; . for mixed or tie'],
               ['AD', 'R', 'Integer', 'Distinct peptides supporting REF, then each ALT allele'],
               ['DP', '1', 'Integer', 'Distinct eligible peptides covering the codon']]
    def vcf = ['##fileformat=VCFv4.3', '##source=giap-stub', '##reference=stub'] +
        info.collect { f -> vcfMeta('INFO', f) } + fmt.collect { f -> vcfMeta('FORMAT', f) } +
        ['#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample'.tokenize(' ').join('\\t')]
    def hdr = [pileup: pile.tokenize(' ').join('\t'), variants: vars.tokenize(' ').join('\t'),
               variants_vcf: vcf.join('\\n')]
    return spec.collect { o -> "${o} ${hdr[o.tokenize('=')[0]]}".toString() }
}

workflow VARIANTS {
    take:
    genome
    gene_model
    mapping
    hits_tolerant
    peptides
    reference_peptides
    peptide_spectra

    main:
    def chans = [kv('mapping', lociDir(mapping)), kv('genome', genome),
                 kv('genes', gene_model.map { a -> a[0] }), kv('replicons', gene_model.map { a -> a[2] }),
                 kv('peptides', peptides), kv('hits_tolerant', hits_tolerant),
                 kv('reference_peptides', reference_peptides), kv('peptide_spectra', peptide_spectra)]
    def needed = ['mapping', 'genome', 'genes', 'replicons', 'peptides']
    def v = stepIn('variants', chans, needed)
    VARIANTS_STEP(v.files, v.keys, variantsSpec(['pileup=tables/pileup.tsv.gz',
        'variants=tables/variants.tsv', 'variants_vcf=tables/variants.vcf']))
    def T = VARIANTS_STEP.out.outputs

    emit:
    pileup       = pick(T, 'pileup.tsv.gz')
    variants     = pick(T, 'variants.tsv')
    variants_vcf = pick(T, 'variants.vcf')
    versions     = VARIANTS_STEP.out.versions
}
