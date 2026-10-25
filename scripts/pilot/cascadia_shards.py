"""Cascadia 0.0.7 ``sequence`` split into prep / predict-shard / merge, for a Slurm array on CPU.

    python -I cascadia_shards.py UNIMOD_XML prep    MZML MODEL WORKDIR [-b 32] [-w 2] [-c 4]
    python -I cascadia_shards.py UNIMOD_XML predict MODEL WORKDIR SHARD N_SHARDS [--limit B]
    python -I cascadia_shards.py UNIMOD_XML merge   WORKDIR OUTPREFIX [-t 0.8]

``prep`` runs Cascadia's own augmentation and builds the spectrum index once (WORKDIR/index.hdf5,
meta.json). ``predict`` runs one contiguous block of whole batches: the same batches, in the same
order, as the unsplit run (greedy decoding, one spectrum per row, so a batch's predictions do not
depend on any other batch). ``merge`` concatenates the blocks in shard order and calls Cascadia's
``write_results`` once, so OUTPREFIX.ssl has the same rows and row numbers as
``cascadia sequence``. The model hyper-parameters and preprocessing are copied from
``cascadia.cascadia.sequence`` (0.0.7). Every step installs ``cascadia_offline``'s socket guard and
offline Unimod table first (R6).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cascadia_offline

META = "meta.json"
INDEX = "index.hdf5"


def shard_batches(n_batches: int, shard: int, n_shards: int) -> range:
    """Contiguous batch block of ``shard`` (0-based); the blocks partition range(n_batches)."""
    if not 0 <= shard < n_shards:
        raise ValueError(f"shard {shard} not in [0, {n_shards})")
    size, extra = divmod(n_batches, n_shards)
    start = shard * size + min(shard, extra)
    return range(start, start + size + (shard < extra))


def _dataset(tokenizer: Any, asf: str | None, index: Path) -> Any:
    from cascadia.depthcharge.data.preprocessing import (  # type: ignore[import-not-found]
        scale_intensity,
        scale_to_unit_norm,
    )
    from cascadia.depthcharge.data.spectrum_datasets import (  # type: ignore[import-not-found]
        AnnotatedSpectrumDataset,
    )

    return AnnotatedSpectrumDataset(
        tokenizer, asf, index_path=str(index),
        preprocessing_fn=[scale_intensity(scaling="root"), scale_to_unit_norm],
    )  # fmt: skip


def _tokenizer() -> Any:
    from cascadia.depthcharge.tokenizers import PeptideTokenizer  # type: ignore[import-not-found]

    return PeptideTokenizer.from_massivekb(reverse=False, replace_isoleucine_with_leucine=True)


def prep(mzml: str, workdir: Path, batch_size: int, width: int, max_charge: int) -> int:
    from cascadia.augment import augment_spectra  # type: ignore[import-not-found]

    workdir.mkdir(parents=True, exist_ok=True)
    if (workdir / INDEX).exists() or (workdir / META).exists():
        print(f"cascadia_shards: {workdir} already prepared", file=sys.stderr)
        return 2
    asf, iso, cycle = augment_spectra(mzml, str(workdir), max_charge=max_charge)
    ds = _dataset(_tokenizer(), asf, workdir / INDEX)
    n = len(ds)
    os.remove(asf)
    meta = {"mzml": mzml, "n_spectra": n, "batch_size": batch_size,
            "n_batches": -(-n // batch_size), "isolation_window_size": float(iso),
            "cycle_time": float(cycle), "width": width, "max_charge": max_charge}  # fmt: skip
    (workdir / META).write_text(json.dumps(meta, indent=1) + "\n")
    print(json.dumps(meta))
    return 0


def predict(model_ckpt: str, workdir: Path, shard: int, n_shards: int, limit: int | None) -> int:
    import pytorch_lightning as pl  # type: ignore[import-not-found]
    import torch  # type: ignore[import-not-found]
    from cascadia.model import AugmentedSpec2Pep  # type: ignore[import-not-found]
    from torch.utils.data import DataLoader, Subset  # type: ignore[import-not-found]

    meta = json.loads((workdir / META).read_text())
    bs, n = meta["batch_size"], meta["n_spectra"]
    block = shard_batches(meta["n_batches"], shard, n_shards)
    if limit is not None:
        block = block[:limit]
    rows = range(block.start * bs, min(block.stop * bs, n)) if len(block) else range(0)
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "0")) or None
    if threads:
        torch.set_num_threads(threads)
    tokenizer = _tokenizer()
    ds = _dataset(tokenizer, None, workdir / INDEX)
    loader = DataLoader(Subset(ds, list(rows)), batch_size=bs, collate_fn=ds.collate_fn,
                        num_workers=2, pin_memory=False)  # fmt: skip
    model = AugmentedSpec2Pep.load_from_checkpoint(
        model_ckpt, d_model=512, n_layers=9, n_head=8, dim_feedforward=1024, dropout=0,
        rt_width=2, tokenizer=tokenizer, max_charge=10,
    )  # fmt: skip
    print(f"shard {shard}/{n_shards}: batches {block.start}..{block.stop} rows {len(rows)} "
          f"threads {threads}")  # fmt: skip
    trainer = pl.Trainer(accelerator="cpu", devices=1, logger=False, enable_checkpointing=False)
    preds = trainer.predict(model, dataloaders=loader) if len(rows) else []
    out = workdir / f"shard_{shard:03d}_of_{n_shards:03d}.pt"
    torch.save({"block": [block.start, block.stop], "rows": len(rows), "preds": preds},
               str(out) + ".tmp")  # fmt: skip
    os.replace(str(out) + ".tmp", out)
    print(f"wrote {out}")
    return 0


def merge(workdir: Path, outprefix: str, threshold: float) -> int:
    import torch  # type: ignore[import-not-found]
    from cascadia.utils import write_results  # type: ignore[import-not-found]

    meta = json.loads((workdir / META).read_text())
    shards = sorted(workdir.glob("shard_*_of_*.pt"))
    counts = {p.name.rsplit("_of_", 1)[1] for p in shards}
    if len(counts) != 1:
        print(f"cascadia_shards: mixed shard counts {sorted(counts)}", file=sys.stderr)
        return 2
    preds: list[Any] = []
    stop = 0
    for p in shards:
        rec = torch.load(p, weights_only=False)
        if rec["block"][0] != stop:
            print(f"cascadia_shards: gap or overlap at {p.name} ({rec['block']}, want {stop})",
                  file=sys.stderr)  # fmt: skip
            return 2
        stop = rec["block"][1]
        preds += rec["preds"]
    if stop != meta["n_batches"] or len(preds) != meta["n_batches"]:
        print(f"cascadia_shards: {len(preds)} of {meta['n_batches']} batches", file=sys.stderr)
        return 2
    write_results(preds, outprefix, meta["mzml"], meta["isolation_window_size"], threshold,
                  meta["width"] * meta["cycle_time"])  # fmt: skip
    print(f"wrote {outprefix}.ssl from {len(shards)} shards")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        print(__doc__, file=sys.stderr)
        return 2
    table, rest = args[0], args[1:]
    if "://" in table or not os.path.isfile(table):
        print(f"cascadia_shards: not an existing local file: {table!r}", file=sys.stderr)
        return 2
    cascadia_offline.install(os.path.abspath(table))
    ap = argparse.ArgumentParser(prog="cascadia_shards.py UNIMOD_XML")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prep")
    p.add_argument("mzml"), p.add_argument("model"), p.add_argument("workdir", type=Path)
    p.add_argument("-b", type=int, default=32), p.add_argument("-w", type=int, default=2)
    p.add_argument("-c", type=int, default=4)
    q = sub.add_parser("predict")
    q.add_argument("model"), q.add_argument("workdir", type=Path)
    q.add_argument("shard", type=int), q.add_argument("n_shards", type=int)
    q.add_argument("--limit", type=int, default=None)
    m = sub.add_parser("merge")
    m.add_argument("workdir", type=Path), m.add_argument("outprefix")
    m.add_argument("-t", type=float, default=0.8)
    a = ap.parse_args(rest)
    if a.cmd == "prep":
        return prep(a.mzml, a.workdir, a.b, a.w, a.c)
    if a.cmd == "predict":
        return predict(a.model, a.workdir, a.shard, a.n_shards, a.limit)
    return merge(a.workdir, a.outprefix, a.t)


if __name__ == "__main__":
    sys.exit(main())
