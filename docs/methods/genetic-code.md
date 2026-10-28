# Genetic codes 4, 11 and 25: the genetic-code check and its limits

GIPPA translates with one NCBI genetic code per run, set by `annotation.translation_table`
(4, 11 or 25; [Inputs › Genetic code](../user-guide/inputs.md#genetic-code)). The three tables
differ only at TGA: a stop under 11, W (Trp) under 4 and G (Gly) under 25. This page describes
what the preflight check `genome.genetic_code` can and cannot tell from the genome, and the
evidence behind its cut.

## The statistic

The check scans all six frames for open reading frames of at least 300 codons twice: once with
TGA as a stop (`bp11`, the bases they cover) and once with TGA as a sense codon (`bp4`, with `n4`
such frames). The ratio is `bp11 / bp4`, compared on exact integers as `20 * bp11 < 11 * bp4`
(ratio < 0.55). The check is evaluated only when `n4 >= 20`; otherwise it is `info`.

- At table 11 it warns when the ratio is below 0.55: the genome probably uses code 4 or 25.
- At tables 4 and 25 it warns when the ratio is 0.55 or more: the genome probably uses code 11.

A ratio of exactly 0.55 is `ok` at 11 and `warn` at 4 and 25, so a genome is "code-11-like" if
and only if its ratio is at least 0.55, and the two decisions never disagree. The warning never
changes a result; `--strict` turns it into exit code 3. The printed ratio is truncated to two
decimals, never rounded up to the cut.

## Code 4 versus code 25

Codes 4 and 25 both read TGA as a sense codon, so this statistic cannot separate them, and no
other genome-only rule in GIPPA does. The table at 4 or 25 must come from the organism's
taxonomy. A code-25 genome run at table 4 (or the reverse) gives no warning. Its peptides with G
(or W) at TGA codons map only through the tolerant route with one difference
(`strain_difference` under `related`).

## Evidence for the 0.55 cut

The statistician's review (2026-10-28) recomputed the ratio on the 61-assembly panel of official
NCBI assemblies used for the code-11 warning, and confirmed that the same cut serves the reverse
warning at tables 4 and 25.

- **Whole assemblies.** The 22 code-4/25 assemblies score 0.006–0.470; the highest are Nasuia
  0.470 (n4 = 41), an M. pneumoniae draft 0.384 and M. genitalium 0.370. The 39 code-11
  assemblies score 0.679–0.999; the lowest are Sodalis 0.679 and M. leprae 0.717. False warnings
  at table 4/25: 0/22 (95 % Clopper–Pearson upper bound 15 %). A code-11 genome run at table 4
  is always warned: 39/39, missed 0/39 (upper bound 9 %). These bounds reflect a panel chosen to
  stress the rule, not a random sample of genomes; the 0.21 gap between the classes is what makes
  the result robust.
- **Drafts** (exponential contig lengths with means of 1, 2 and 5 kb, contigs under 200 bp
  dropped, 5 seeds each). Nasuia at a mean of 1 kb reaches 0.520 (n4 = 22–29), still 0.03 below
  the cut. M. pneumoniae and M. genitalium stay at or below 0.39. Sodalis goes down to 0.635 and
  M. leprae to 0.662.
- **One cut, not two.** A separate cut of 0.60 for tables 4/25 would lower false warnings on
  partial inputs from 0.52 % to 0.05 % but double the missed warnings from 0.75 % to 1.50 %, and
  shrink the margin to the Sodalis drafts from 0.09 to 0.04. A missed warning at table 4/25 (a
  code-11 genome read with TGA as sense) silently gives read-through fusions, while a false one
  only prints a message. 0.55 sits near the centre of the gap and is roughly balanced on partial
  inputs.

## Limit 1: small inputs (20–40 long open reading frames)

On whole genomes the check made no error on the panel. On partial inputs it can be wrong. The
review drew 60 random single windows of 20–100 kb per assembly and kept those with n4 ≥ 20:

| Error | Windows | Rate | 95 % upper bound |
|---|---|---|---|
| false warning on a code-4/25 window at table 4/25 (ratio ≥ 0.55) | 11/2109 | 0.52 % | 0.93 % |
| code-11 window below the cut (a missed warning at table 4/25, or a false warning at table 11) | 59/7912 | 0.75 % | 0.96 % |

The highest ratio of a code-4/25 window was 0.611. False warnings at tables 4/25 are
concentrated at small n4: 10/1056 windows at n4 20–29, 1/447 at 30–39 and 0/606 at 40 or more.

**So, on inputs with only 20–40 long open reading frames, such as a plasmid alone or a partial
assembly, either warning can be wrong about 1 % of the time.** On such an input, treat the
warning as a prompt to check the organism's genetic code, not as a finding. Raising the gate
from 20 to 40 frames was considered and rejected: it would put Nasuia (n4 = 41) and its drafts
(n4 = 22–40) at `info`.

## Limit 2: code 25 rests on one assembly

The panel holds one code-25 genome (Ca. Absconditicoccus praedator, ratio 0.006). The code-25
side of the cut therefore rests on a single assembly. Any official code-25 assembly (SR1,
Gracilibacteria) used in later validation is added to this evidence; if one scores 0.40 or more,
the statistician reviews the cut again.

## Chance hits in six-frame mapping

Under codes 4 and 25 there are fewer stops per frame, so stop-free segments are longer and there
are more chance hits; the chance model takes its residue frequencies from the six-frame index, so
it follows the table. The statistician's note on the chance-hit calibration: the residual six-frame chance-model ratio of
about +3 % at codes 11, 4 and 25 comes from codon dependence within a window and does not depend
on the table.
