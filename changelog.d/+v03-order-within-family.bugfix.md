V03 locus plots are ordered within each test family instead of on one mixed scale. Loci
significant at their own family's level come first, then non-significant loci, then loci with
no estimate, and the families (prophage, BGC, island, operon discordance, missed operon) take
turns within each group, so under the `report.max_locus_plots` cap one family's significant loci
can no longer crowd out another's.
Operon flags from the two flag families are no longer merged into one locus. The index figure
names each estimate (`q`, `FWER p`, `decoy-FDR q`) and labels the groups.
