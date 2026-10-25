<!-- m07 -->
Selenocysteine and pyrrolysine codons: annotated Sec (`U` over TGA) and Pyl (`O` over TAG) codons
are now consumable sites for tolerant peptide placement, with decoys mirrored as for annotated
stops. Rows over them are flagged `sec_site` / `pyl_site` in `peptide_hits.tsv` (GFF3
`recoded_codon`), `variants.tsv` reports the `U`/`O` residue as an ALT site, and
`detection_summary.tsv` counts the two new `placement_flag` groups.
