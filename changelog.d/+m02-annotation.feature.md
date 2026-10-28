<!-- m02 -->
Genome annotation: `gippa setup` downloads and checks the Bakta (light or full), geNOMAD, antiSMASH,
Promotech, GO, KEGG, COG and contaminant databases and records them in `db_manifest.json`. A run
annotates the genome with Bakta and adds operon pairs and transcription units (UniOP), promoters,
terminators, ribosome binding sites, codon usage and regions (prophage, plasmid, biosynthetic gene
cluster, CRISPR, AMR). A Bakta JSON or your own GFF3/GBFF annotation can be supplied instead.
The pipeline runs through Nextflow (`-profile farmshare,apptainer`) or a local executor.
