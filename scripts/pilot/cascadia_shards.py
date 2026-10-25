"""Cascadia 0.0.7 ``sequence`` split into prep / predict-shard / merge, for a Slurm array on CPU.

    python -I cascadia_shards.py UNIMOD_XML prep    MZML WORKDIR [-b 32] [-w 2] [-c 4]
    python -I cascadia_shards.py UNIMOD_XML predict MODEL WORKDIR SHARD N_SHARDS [--limit B]
    python -I cascadia_shards.py UNIMOD_XML work    MODEL GEN DEADLINE_S CHUNK WORKDIR [WORKDIR ...]
    python -I cascadia_shards.py UNIMOD_XML merge   WORKDIR OUTPREFIX [-t 0.8]

``prep`` runs Cascadia's own augmentation and builds the spectrum index once (WORKDIR/index.hdf5,
meta.json). ``predict`` runs one contiguous block of whole batches: the same batches, in the same
order, as the unsplit run (greedy decoding, one spectrum per row, so a batch's predictions do not
depend on any other batch). ``merge`` concatenates the blocks in shard order and calls Cascadia's
``write_results`` once, so OUTPREFIX.ssl has the same rows and row numbers as
``cascadia sequence``. ``work`` is the load-balanced form: batch cost grows several-fold along the
gradient, so equal blocks were up to 6x apart. Each worker claims CHUNK-batch chunks (an atomic
mkdir under claims/g<GEN>/) of the work dirs in order, saves each chunk as it finishes, and claims
no new chunk after DEADLINE_S. A later generation picks up whatever is left. Chunks lie on the
global grid [k*CHUNK, (k+1)*CHUNK), so the grid never shifts with the markers. Batches listed in a
reserved_<start>_<stop> file (a half-open block an earlier ``predict`` task is still writing) or
already held by a shard_*_of_*.pt / chunk_*_*.pt file are skipped (``pending_chunks``, which the
watcher calls too; it needs no torch). Files are written to a per-job temporary name next to the
target and renamed into place.
The model hyper-parameters and preprocessing are copied from
``cascadia.cascadia.sequence`` (0.0.7). Every step installs ``cascadia_offline``'s socket guard and
offline Unimod table first (R6).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cascadia_offline

META = "meta.json"
INDEX = "index.hdf5"
_RESERVED = re.compile(r"reserved_(\d+)_(\d+)")
_SHARD = re.compile(r"shard_(\d+)_of_(\d+)\.pt")
_CHUNK = re.compile(r"chunk_(\d+)_(\d+)\.pt")


def shard_batches(n_batches: int, shard: int, n_shards: int) -> range:
    """Contiguous batch block of ``shard`` (0-based); the blocks partition range(n_batches)."""
    if not 0 <= shard < n_shards:
        raise ValueError(f"shard {shard} not in [0, {n_shards})")
    size, extra = divmod(n_batches, n_shards)
    start = shard * size + min(shard, extra)
    return range(start, start + size + (shard < extra))


def _save_atomic(obj: Any, out: Path) -> None:
    """``torch.save`` to a per-job temporary file next to ``out``, then rename it into place."""
    import torch  # type: ignore[import-not-found]

    job = os.environ.get("SLURM_JOB_ID", "nojob")
    tmp = out.with_name(f"{out.name}.{job}.{os.getpid()}.tmp")
    try:
        torch.save(obj, str(tmp))
        os.replace(tmp, out)
    finally:
        if tmp.exists():
            tmp.unlink()


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
    _save_atomic({"block": [block.start, block.stop], "rows": len(rows), "preds": preds}, out)
    print(f"wrote {out}")
    return 0


def _reserved(workdir: Path) -> list[tuple[int, int]]:
    """Half-open blocks named by ``reserved_<a>_<b>`` markers; any other file name is ignored."""
    out = []
    for p in workdir.glob("reserved_*"):
        m = _RESERVED.fullmatch(p.name)
        if m:
            out.append((int(m[1]), int(m[2])))
    return sorted(out)


def _covered(workdir: Path, n_batches: int) -> list[tuple[int, int]]:
    """Blocks held by finished ``shard_*_of_*.pt`` and ``chunk_*_*.pt`` files, from their names.

    A shard file's block is ``shard_batches(n_batches, shard, n_shards)``; a ``predict --limit``
    shard holds less than that, which ``merge`` then reports as a gap.
    """
    out = []
    for p in workdir.iterdir():
        if m := _CHUNK.fullmatch(p.name):
            out.append((int(m[1]), int(m[2])))
        elif (m := _SHARD.fullmatch(p.name)) and int(m[1]) < int(m[2]):
            blk = shard_batches(n_batches, int(m[1]), int(m[2]))
            out.append((blk.start, blk.stop))
    return sorted(out)


def chunk_plan(n_batches: int, reserved: list[tuple[int, int]], size: int) -> list[range]:
    """The batches of range(n_batches) outside ``reserved``, cut on the global grid of ``size``.

    Each chunk is a maximal run of free batches inside one cell ``[k*size, (k+1)*size)``, so a
    reserved block only trims the cells it touches. Blocks may be unsorted, overlapping or reach
    past ``n_batches`` (they are clamped).
    """
    if size < 1:
        raise ValueError(f"chunk size {size} < 1")
    free: list[tuple[int, int]] = []
    pos = 0
    for a, b in sorted((max(0, a), min(b, n_batches)) for a, b in reserved):
        if a >= b:
            continue
        if a > pos:
            free.append((pos, a))
        pos = max(pos, b)
    if pos < n_batches:
        free.append((pos, n_batches))
    plan: list[range] = []
    for a, b in free:
        s = a
        while s < b:
            e = min((s // size + 1) * size, b)
            plan.append(range(s, e))
            s = e
    return plan


def pending_chunks(workdir: Path, size: int) -> list[range]:
    """The chunks of ``workdir`` still to compute: not reserved and not held by a finished file."""
    n_batches = int(json.loads((workdir / META).read_text())["n_batches"])
    return chunk_plan(n_batches, _reserved(workdir) + _covered(workdir, n_batches), size)


def work(model_ckpt: str, gen: int, deadline_s: float, size: int, workdirs: list[Path]) -> int:
    import time

    import pytorch_lightning as pl  # type: ignore[import-not-found]
    import torch  # type: ignore[import-not-found]
    from cascadia.model import AugmentedSpec2Pep  # type: ignore[import-not-found]
    from torch.utils.data import DataLoader, Subset  # type: ignore[import-not-found]

    t0 = time.monotonic()
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "0")) or None
    if threads:
        torch.set_num_threads(threads)
    tokenizer = _tokenizer()
    model = AugmentedSpec2Pep.load_from_checkpoint(
        model_ckpt, d_model=512, n_layers=9, n_head=8, dim_feedforward=1024, dropout=0,
        rt_width=2, tokenizer=tokenizer, max_charge=10,
    )  # fmt: skip
    done = 0
    for workdir in workdirs:
        meta = json.loads((workdir / META).read_text())
        bs, n = meta["batch_size"], meta["n_spectra"]
        ds = _dataset(tokenizer, None, workdir / INDEX)
        claims = workdir / "claims" / f"g{gen}"
        claims.mkdir(parents=True, exist_ok=True)
        for block in pending_chunks(workdir, size):
            out = workdir / f"chunk_{block.start:05d}_{block.stop:05d}.pt"
            if out.exists():
                continue
            if time.monotonic() - t0 > deadline_s:
                print(f"deadline: stopped after {done} chunks")
                return 0
            try:
                (claims / out.stem).mkdir()
            except FileExistsError:
                continue
            rows = list(range(block.start * bs, min(block.stop * bs, n)))
            loader = DataLoader(Subset(ds, rows), batch_size=bs, collate_fn=ds.collate_fn,
                                num_workers=2, pin_memory=False)  # fmt: skip
            trainer = pl.Trainer(accelerator="cpu", devices=1, logger=False,
                                 enable_checkpointing=False, enable_progress_bar=False)  # fmt: skip
            t = time.monotonic()
            preds = trainer.predict(model, dataloaders=loader)
            _save_atomic({"block": [block.start, block.stop], "rows": len(rows), "preds": preds},
                         out)  # fmt: skip
            done += 1
            print(f"{workdir.parent.name} {out.name} {time.monotonic() - t:.0f}s", flush=True)
    print(f"no chunks left: {done} done")
    return 0


def merge(workdir: Path, outprefix: str, threshold: float) -> int:
    import torch  # type: ignore[import-not-found]
    from cascadia.utils import write_results  # type: ignore[import-not-found]

    meta = json.loads((workdir / META).read_text())
    bs, n = meta["batch_size"], meta["n_spectra"]
    files = sorted(p for p in workdir.iterdir()
                   if _SHARD.fullmatch(p.name) or _CHUNK.fullmatch(p.name))  # fmt: skip
    recs = [torch.load(p, weights_only=False) for p in files]
    recs.sort(key=lambda r: r["block"][0])
    preds: list[Any] = []
    stop = 0
    for rec in recs:
        a, b = rec["block"]
        if a != stop or b <= a:
            print(f"cascadia_shards: gap or overlap at block {rec['block']} (want start {stop})",
                  file=sys.stderr)  # fmt: skip
            return 2
        want_rows = max(0, min(b * bs, n) - a * bs)
        if len(rec["preds"]) != b - a or rec["rows"] != want_rows:
            print(f"cascadia_shards: block {rec['block']} has {len(rec['preds'])} predictions "
                  f"and {rec['rows']} rows (want {b - a} and {want_rows})",
                  file=sys.stderr)  # fmt: skip
            return 2
        stop = b
        preds += rec["preds"]
    if stop != meta["n_batches"] or len(preds) != meta["n_batches"]:
        print(f"cascadia_shards: {len(preds)} of {meta['n_batches']} batches", file=sys.stderr)
        return 2
    write_results(preds, outprefix, meta["mzml"], meta["isolation_window_size"], threshold,
                  meta["width"] * meta["cycle_time"])  # fmt: skip
    print(f"wrote {outprefix}.ssl from {len(recs)} blocks")
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
    p.add_argument("mzml"), p.add_argument("workdir", type=Path)
    p.add_argument("-b", type=int, default=32), p.add_argument("-w", type=int, default=2)
    p.add_argument("-c", type=int, default=4)
    q = sub.add_parser("predict")
    q.add_argument("model"), q.add_argument("workdir", type=Path)
    q.add_argument("shard", type=int), q.add_argument("n_shards", type=int)
    q.add_argument("--limit", type=int, default=None)
    w = sub.add_parser("work")
    w.add_argument("model"), w.add_argument("gen", type=int)
    w.add_argument("deadline_s", type=float), w.add_argument("size", type=int)
    w.add_argument("workdirs", type=Path, nargs="+")
    m = sub.add_parser("merge")
    m.add_argument("workdir", type=Path), m.add_argument("outprefix")
    m.add_argument("-t", type=float, default=0.8)
    a = ap.parse_args(rest)
    if a.cmd == "prep":
        return prep(a.mzml, a.workdir, a.b, a.w, a.c)
    if a.cmd == "predict":
        return predict(a.model, a.workdir, a.shard, a.n_shards, a.limit)
    if a.cmd == "work":
        return work(a.model, a.gen, a.deadline_s, a.size, a.workdirs)
    return merge(a.workdir, a.outprefix, a.t)


if __name__ == "__main__":
    sys.exit(main())
