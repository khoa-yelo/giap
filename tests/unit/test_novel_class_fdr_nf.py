"""Static checks of the M7-11b Nextflow wiring (docs/specs/M7-11b-class-fdr-wiring.md "Nextflow
wiring"; no Nextflow needed, as tests/unit/test_novel_nf_registry.py does for M7-08):

* ``mapping.nf``: PLACEMENT's out spec gains ``placement_costs=placement_costs.parquet``;
  ``emit:`` gains ``hits_exact`` and ``placement_costs`` (``hits_tolerant`` stays);
* ``novel.nf``: ``take:`` gains ``hits`` (read here: appended after ``layers``); the
  ``novel_orfs`` stepIn passes ``peptides``, ``features``, ``hits_exact``, ``hits_tolerant`` and
  ``placement_costs``, none of them in its required-keys list; its out spec gains
  ``novel_class_fdr=novel/novel_class_fdr.parquet`` and ``emit:`` gains ``novel_class_fdr``;
* ``main.nf``: NOVEL gets ``MAPPING.out.hits_exact``, ``.hits_tolerant`` and
  ``.placement_costs``; ``kv('novel_class_fdr', NOVEL.out.novel_class_fdr)`` is published;
* ``runner/local.py``: ``PUBLISHED_FILES["novel_class_fdr"] = "novel/novel_class_fdr.parquet"``.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
NOVEL_NF = WF / "subworkflows" / "local" / "novel.nf"
MAPPING_NF = WF / "subworkflows" / "local" / "mapping.nf"
MAIN_NF = WF / "main.nf"
NEW_KEYS = ("peptides", "features", "hits_exact", "hits_tolerant", "placement_costs")


def _strip(text: str) -> str:
    return re.sub(r"(?m)//.*$", "", text)


def _read(p: Path) -> str:
    return _strip(p.read_text(encoding="utf-8"))


def _call_body(text: str, start: int) -> str:
    depth, i = 1, start
    while depth:
        depth += {"(": 1, ")": -1, "[": 0, "]": 0}.get(text[i], 0)
        i += 1
    return text[start : i - 1]


def _emit(text: str, workflow: str) -> dict[str, str]:
    m = re.search(rf"workflow\s+{workflow}\s*\{{.*?emit:(.*)\}}\s*$", text, flags=re.S)
    assert m, f"no emit: block in workflow {workflow}"
    return dict(re.findall(r"(?m)^\s*(\w+)\s*=\s*(.+?)\s*$", m.group(1)))


def _stepin_body(text: str, step: str) -> str:
    m = re.search(rf"stepIn\(\s*'{step}'\s*,", text)
    assert m, f"no stepIn('{step}', ...)"
    return _call_body(text, m.end())


def _keys_in(text: str, body: str) -> set[str]:
    keys = set(re.findall(r"kv\(\s*'(\w+)'", body))
    var_keys = dict(re.findall(r"def\s+(\w+)\s*=\s*kv\(\s*'(\w+)'", text))
    for tok in re.findall(r"\b(\w+)\b", body):
        if tok in var_keys:
            keys.add(var_keys[tok])
    return keys


# ---------------------------------------------------------------- mapping.nf


def test_placement_writes_placement_costs() -> None:
    text = _read(MAPPING_NF)
    m = re.search(r"PLACEMENT\(\s*pl\.files\s*,\s*pl\.keys\s*,", text)
    assert m
    spec = _call_body(text, m.end())
    assert "'placement_costs=placement_costs.parquet'" in spec, spec


def test_mapping_emits_the_hit_tables_and_costs() -> None:
    emit = _emit(_read(MAPPING_NF), "MAPPING")
    assert "hits_exact.parquet" in emit.get("hits_exact", ""), emit
    assert "hits_tolerant.parquet" in emit.get("hits_tolerant", ""), emit
    assert "placement_costs.parquet" in emit.get("placement_costs", ""), emit


# ---------------------------------------------------------------- novel.nf


def test_novel_takes_hits_last() -> None:
    m = re.search(r"workflow\s+NOVEL\s*\{\s*take:(.*?)main:", _read(NOVEL_NF), flags=re.S)
    assert m
    names = re.findall(r"(?m)^\s*(\w+)\s*$", m.group(1))
    assert names[-1] == "hits", names
    assert names[:-1] == ["genome", "gene_model", "contract", "mapping", "peptides",
                          "contaminants", "gene_table", "layers"]  # fmt: skip


def test_novel_orfs_gets_the_new_optional_keys() -> None:
    text = _read(NOVEL_NF)
    body = _stepin_body(text, "novel_orfs")
    keys = _keys_in(text, body)
    assert set(NEW_KEYS) <= keys, sorted(keys)
    req = re.search(r"\[\s*((?:'\w+'\s*,?\s*)+)\]\s*$", body)
    assert req, body
    required = set(re.findall(r"'(\w+)'", req.group(1)))
    assert required == {"novel_candidates", "gene_table", "genes", "replicons", "genome"}
    assert not required & set(NEW_KEYS)


def test_novel_orfs_publishes_and_emits_the_side_table() -> None:
    text = _read(NOVEL_NF)
    m = re.search(r"NOVEL_ORFS\(\s*no\.files\s*,\s*no\.keys\s*,", text)
    assert m
    spec = _call_body(text, m.end())
    assert "'novel_class_fdr=novel/novel_class_fdr.parquet'" in spec, spec
    emit = _emit(text, "NOVEL")
    assert "novel_class_fdr.parquet" in emit.get("novel_class_fdr", ""), emit


# ---------------------------------------------------------------- main.nf + local runner


def test_main_passes_mapping_hits_to_novel_and_publishes_the_side_table() -> None:
    text = _read(MAIN_NF)
    n = text.find("NOVEL(")
    assert n >= 0
    body = _call_body(text, n + len("NOVEL("))
    for out in ("hits_exact", "hits_tolerant", "placement_costs"):
        assert f"MAPPING.out.{out}" in body, (out, body)
    assert re.search(r"kv\(\s*'novel_class_fdr'\s*,\s*NOVEL\.out\.novel_class_fdr\s*\)", text)


def test_local_runner_publishes_the_side_table() -> None:
    from giap.runner.local import PUBLISHED_FILES

    assert PUBLISHED_FILES.get("novel_class_fdr") == "novel/novel_class_fdr.parquet"
