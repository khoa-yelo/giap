<!-- m10 -->
Pilot validation: GIAP was run on public E. coli Nissle 1917 proteomics data (iProX PXD078459)
searched with DIA-NN against a generic K-12 reference proteome and against the Nissle database
from `giap makedb`, and the two runs were compared with `giap compare`. `giap convert --from
cascadia` reads de novo peptides. A reproducibility runbook lists every command and input hash,
and runtimes on 16 threads were measured on FarmShare.
