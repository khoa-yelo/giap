"""giapsim six-frame target rows name only real makedb entries, ``nterm_extension`` included
(docs/specs/M7-11b-class-fdr-wiring.md lead decisions 9-10; the M7-13b fast-tier run, Slurm
1873888, still had 22 ``nterm_extension`` ``no_makedb_pair`` ORFs over R = 10 at 1M).

``giap makedb`` keeps an ``nterm_extension`` entry only when the upstream part is >= 7 aa
(``giap.makedb.build.NTERM_MIN_UPSTREAM_AA``). A planted N-terminal extension of 5-6 codons with
a stop right before its true start therefore has no entry, and a search against that FASTA cannot
report a peptide reaching upstream of the gene's annotated start. In ``novel_orf_power(30,
"100k")`` one such gene has a detected upstream peptide (giapsim before this fix named it
``GIAPORF_nterm_extension_chr1_18104_18208_m``, an entry makedb does not build). Checked against
the database GIAP itself builds (``m7_round3_oracle``).
"""

from __future__ import annotations

from pathlib import Path

from helpers import m7_round3_oracle as Q
from helpers import m7_sim_oracle as O

import giapsim
import giapsim.scenarios as S

SEED = 30  # 100k: a detected upstream peptide of an N-terminal extension without an entry


def test_target_ids_are_makedb_entries_and_entryless_extensions_are_unsearched(
    tmp_path: Path,
) -> None:
    out = tmp_path / "sim"
    giapsim.simulate(S.novel_orf_power(SEED, "100k"), SEED, out)
    entries = Q.makedb_entries(out, tmp_path / "db")
    tags = {e.locus_tag for e in entries.values() if e.cls == "nterm_extension"}
    feats = {r["feature_id"]: r for r in O.rows(out / "truth" / "features.tsv")}
    miss = {k for k, r in feats.items() if r["true_start"] and k not in tags}
    assert miss, "the seed no longer has an N-terminal extension without a makedb entry"

    rows = O.rows(out / "peptides.tsv")
    orf_ids = {i.removeprefix("rev_") for r in rows for i in Q.ids(r["protein_ids"])}
    bad = sorted(i for i in orf_ids if i.startswith("GIAPORF_") and i not in entries)
    assert not bad, f"{len(bad)} protein IDs name no makedb entry: {bad[:5]}"

    targets = {r["peptide"] for r in rows if r["is_decoy"] != "true"}
    occ: dict[str, list[bool]] = {}  # peptide -> per truth occurrence: upstream of a missing entry
    for t in O.rows(out / "truth" / "peptides.tsv"):
        if t["origin"] != "gene":
            continue
        f = feats.get(t["locus_tag"])
        up = (f is not None and t["locus_tag"] in miss
              and (int(t["start"]) < int(f["start"]) if f["strand"] == "+"
                   else int(t["end"]) > int(f["end"])))  # fmt: skip
        occ.setdefault(t["peptide"], []).append(up)
    unsearched = sorted(p for p in targets if occ.get(p) and all(occ[p]))
    assert not unsearched, f"target rows with no makedb entry to come from: {unsearched}"
