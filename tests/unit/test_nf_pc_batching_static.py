"""M2-14b static checks: annotate.nf batches the Promoter Calculator regions and the module runs a
batch (docs/specs/M2-14b-promoter-calculator-batching.md items 1, 2, 4, 6, 7). No Nextflow needed;
the behaviour itself is checked by workflow/tests/annotate_pc_batch.nf.test and
workflow/tests/promoter_calculator_batch.nf.test.

Pinned by the test author: the batch id and task tag is ``pc_<first TU>_<last TU>`` (item 4, first
form), built in annotate.nf and passed as the module's ``meta_id``.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SUBWF = REPO / "workflow" / "subworkflows" / "local" / "annotate.nf"
MODULE = REPO / "workflow" / "modules" / "local" / "promoter_calculator" / "main.nf"
DEFAULTS = REPO / "src" / "gippa" / "config" / "defaults.yaml"


def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(?m)//.*$", "", text)


def _sub() -> str:
    return _strip_comments(SUBWF.read_text(encoding="utf-8"))


def _module() -> str:
    return _strip_comments(MODULE.read_text(encoding="utf-8"))


def _pc_wiring() -> str:
    """annotate.nf from ``inDir(O, 'tu_upstream')`` to the end of the PROMOTER_CALCULATOR(...) call
    (balanced parentheses, so a call spanning several lines is included)."""
    sub = _sub()
    m = re.search(r"inDir\s*\(\s*O\s*,\s*['\"]tu_upstream['\"]\s*\)", sub)
    assert m, "annotate.nf: no inDir(O, 'tu_upstream') regions"
    call = re.compile(r"\bPROMOTER_CALCULATOR\s*\(").search(sub, m.start())
    assert call, "annotate.nf: no PROMOTER_CALCULATOR(...) call after the tu_upstream regions"
    depth, i = 1, call.end()
    while depth and i < len(sub):
        depth += {"(": 1, ")": -1}.get(sub[i], 0)
        i += 1
    return sub[m.start() : i]


def test_pc_batch_is_a_constant_of_100_in_annotate_nf() -> None:
    """Item 7: PC_BATCH = 100 is a constant in annotate.nf."""
    assert re.search(r"(?m)^\s*(?:[\w ]+\s)?PC_BATCH\s*=\s*100\b",
                     _sub()), "annotate.nf must define PC_BATCH = 100"  # fmt: skip


def test_pc_batch_is_not_a_config_key() -> None:
    """Item 7: not a config key (defaults.yaml is protected; a key needs an ADR)."""
    assert "pc_batch" not in DEFAULTS.read_text(encoding="utf-8").lower()


def test_regions_are_sorted_by_tu_id_then_collated_in_batches_of_pc_batch() -> None:
    """Item 1: the regions, sorted by TU id (file base name), are collated in batches of PC_BATCH
    keeping the remainder (Groovy List.collate / channel collate keep it by default). The sort
    is either explicit in the wiring or inDir's own by-name sort of the step's region files."""
    wiring = _pc_wiring()
    coll = re.search(r"\bcollate\s*\(\s*PC_BATCH\b([^)]*)\)", wiring)
    assert coll, f"no collate(PC_BATCH ...) between the regions and PROMOTER_CALCULATOR: {wiring}"
    rest = coll.group(1)
    assert not re.search(r"\bfalse\b", rest), f"remainder dropped: collate(PC_BATCH{rest})"
    by_name = r"(?:toSortedList|toSorted|sort)\s*(?:\(\s*)?\{[^}]*\b(?:baseName|name|simpleName)\b"
    explicit = re.search(by_name, wiring[: coll.start()])
    indir = re.search(r"def\s+inDir\b.*?\n\}", _sub(), flags=re.S)
    assert explicit or (indir and re.search(by_name, indir.group(0))), (
        "regions must be sorted by TU id (file name) before collate(PC_BATCH)"
    )


def test_one_promoter_calculator_call_in_annotate() -> None:
    assert len(re.findall(r"\bPROMOTER_CALCULATOR\s*\(", _sub())) == 1


def test_batch_id_is_pc_first_last_and_the_module_tags_with_it() -> None:
    """Item 4: each batch is ``[pc_<first TU>_<last TU>, regions]``; the module keeps
    ``tuple val(meta_id), path(...)`` and tag "${meta_id}" (test_nf_modules_lint.py)."""
    wiring = _pc_wiring()
    assert re.search(r"pc_\$\{[^}]*\bbaseName\b[^}]*\}_\$\{[^}]*\bbaseName\b[^}]*\}", wiring), (
        f'annotate.nf must build the batch id "pc_${{<first TU>}}_${{<last TU>}}": {wiring}'
    )
    inp = re.search(r"\binput\s*:(.*?)\boutput\s*:", _module(), flags=re.S)
    assert inp, "no input: block"
    assert re.search(r"tuple\s+val\s*\(\s*meta_id\s*\)\s*,\s*path\s*\(", inp.group(1)), (
        f"input must be tuple val(meta_id), path(<regions>): {inp.group(1).strip()}"
    )


def test_module_emits_every_csv_of_the_batch() -> None:
    """Item 2: output ``path "*.csv"`` (one CSV per TU, the same names) emitted as csv."""
    out = re.search(r"\boutput\s*:(.*?)\b(?:script|when)\s*:", _module(), flags=re.S)
    assert out, "no output: block"
    assert re.search(
        r"path\s*\(?\s*['\"]\*\.csv['\"]\s*\)?\s*,\s*emit\s*:\s*csv\b", out.group(1)
    ), f'output must be path "*.csv", emit: csv: {out.group(1).strip()}'


def test_module_keeps_process_single_and_two_hours_per_attempt() -> None:
    """Item 6: label process_single; time { 2.h * task.attempt }."""
    head = re.search(r"\bprocess\s+PROMOTER_CALCULATOR\s*\{(.*?)\binput\s*:", _module(), re.S)
    assert head, "process block with an input: section not found"
    d = head.group(1)
    assert re.search(r"label\s+['\"]process_single['\"]", d)
    assert re.search(r"\btime\s*\{\s*2\.h\s*\*\s*task\.attempt\s*\}", d), (
        f"time must be {{ 2.h * task.attempt }}: {d.strip()}"
    )


def test_script_runs_each_region_with_the_unchanged_command_and_fallback() -> None:
    """Item 2: per region `promoter-calculator -i <region> -o <TU>.csv -j 1`, then the header
    fallback for that TU, inside a loop over the batch."""
    m = re.search(r"\bscript\s*:(.*?)\bstub\s*:", _module(), flags=re.S)
    assert m, "no script: block"
    script = m.group(1)
    assert re.search(r"\bfor\b|\.each\b|\.collect\b", script), "no per-region loop in the script"
    calls = re.findall(r"promoter-calculator\b[^\n]*", script)
    assert calls, "no promoter-calculator call"
    for c in calls:
        assert re.search(r"(?<![\w-])-i\s+\S+", c) and re.search(r"-o\s+\S*\.csv", c), c
        assert re.search(r"(?<![\w-])-j\s+1(?!\d)", c), c
    assert re.search(r"\[\s*-e\s+[^\]]*\.csv['\"]?\s*\]\s*\|\|\s*printf", script), (
        "the no-promoter header fallback must stay per region"
    )
