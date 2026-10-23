"""PR #205 review round B4: under Nextflow, REPORT receives every published key that the local
executor hands to the ``plots`` and ``report`` steps (GIAP_SPEC.md §6.4 executor equality; M9-02
behaviour 8). Static parse of ``workflow/main.nf`` and ``workflow/subworkflows/local/*.nf``
against the step registry, in the style of test_report_wiring.py / test_nf_operons_registry.py;
no Nextflow needed.

Keys under test, generated from the registry: the inputs and optional inputs of ``plots`` and
``report`` that are published (``PUBLISHED`` / ``PUBLISHED_FILES``) and produced by a registered
step: the local runner stages each of them from ``work/`` whenever its producer ran. For each:

* main.nf adds ``kv('<key>', ...)`` or ``outs('<key>', ...)`` to the ``pub`` list given to
  ``REPORT``;
* when that channel is ``SUB.out.<name>`` (optionally inside ``pick(...)``), the subworkflow
  ``SUB`` emits ``<name>``; if a subworkflow's GIAP_STEP out spec publishes ``'<key>=...'``,
  ``SUB`` is that subworkflow; and the published file's name is picked, in the emit line or in
  main.nf (so the channel carries that file, not another output of the same task).

``gene_context``, ``genes_pre`` and ``island_diagnostics`` are the keys the review found
missing.
"""

from __future__ import annotations

import re
from pathlib import Path, PurePosixPath

import pytest

from giap.runner.local import PUBLISHED, PUBLISHED_FILES
from giap.steps import load_step_modules
from giap.steps import registry as reg

REPO = Path(__file__).resolve().parents[2]
MAIN = REPO / "workflow" / "main.nf"
LOCAL = REPO / "workflow" / "subworkflows" / "local"


def _strip(text: str) -> str:
    return re.sub(r"(?m)//.*$", "", text)


def _report_keys() -> list[str]:
    load_step_modules()
    produced = {k for s in reg.all_steps() for k in s.outputs}
    published = set(PUBLISHED) | set(PUBLISHED_FILES)
    keys: set[str] = set()
    for name in ("plots", "report"):
        spec = reg.get(name)
        keys |= {*spec.inputs, *spec.optional_inputs}
    return sorted((keys & published & produced) - {"plots", "report"})


KEYS = _report_keys()


def _main() -> str:
    return _strip(MAIN.read_text(encoding="utf-8"))


def _pub_text() -> str:
    """main.nf from the genome branch up to the REPORT call (where ``pub`` is built)."""
    text = _main()
    end = text.index("REPORT(pub")
    return text[:end]


def _includes() -> dict[str, Path]:
    """Subworkflow name -> its .nf file, from main.nf's ``include { A; b } from '...'``."""
    out: dict[str, Path] = {}
    for m in re.finditer(r"include\s*\{([^}]*)\}\s*from\s*'\./subworkflows/local/(\w+)'", _main()):
        for name in (n.strip() for n in m.group(1).split(";")):
            if re.fullmatch(r"[A-Z][A-Z0-9_]*", name):
                out[name] = LOCAL / f"{m.group(2)}.nf"
    return out


def _emits(path: Path) -> dict[str, str]:
    """``emit:`` name -> its right-hand side."""
    text = _strip(path.read_text(encoding="utf-8"))
    m = re.search(r"\bemit:\s*\n(.*?)\n\}", text, re.S)
    assert m, f"{path.name}: no emit block"
    out = {}
    for line in m.group(1).splitlines():
        lm = re.match(r"\s*(\w+)\s*=\s*(.+)$", line)
        if lm:
            out[lm.group(1)] = lm.group(2)
    return out


def _out_spec_owner(key: str) -> list[str]:
    """Subworkflows whose text publishes ``'<key>=`` (a GIAP_STEP out spec)."""
    found = []
    for name, path in _includes().items():
        if re.search(rf"['\"]{key}=", _strip(path.read_text(encoding="utf-8"))):
            found.append(name)
    return sorted(set(found))


def _call(text: str, start: int) -> str:
    """The text of a call from just after its ``(`` to the matching ``)``."""
    depth, i = 1, start
    while depth:
        depth += {"(": 1, ")": -1}.get(text[i], 0)
        i += 1
    return text[start : i - 1]


def test_key_list_is_generated_and_covers_the_review_findings() -> None:
    assert {"gene_context", "genes_pre", "island_diagnostics"} <= set(KEYS), KEYS
    assert {"islands", "genes_tsv", "operon_flags", "novel_orfs"} <= set(KEYS), KEYS


@pytest.mark.parametrize("key", KEYS)
def test_main_nf_passes_every_report_key_to_report(key: str) -> None:
    pub = _pub_text()
    assert re.search(rf"\b(?:kv|outs)\(\s*'{key}'\s*,", pub), (
        f"main.nf: no kv('{key}', ...) / outs('{key}', ...) in the channels given to REPORT, "
        f"though the local executor stages '{key}' for plots/report"
    )


@pytest.mark.parametrize("key", KEYS)
def test_report_key_channel_is_emitted_by_its_producing_subworkflow(key: str) -> None:
    pub = _pub_text()
    m = re.search(rf"\b(kv|outs)\(\s*'{key}'\s*,", pub)
    if m is None:
        pytest.fail(f"main.nf does not pass '{key}' to REPORT")
    arg = _call(pub, m.start() + len(m.group(1)) + 1)
    ref = re.search(r"\b([A-Z][A-Z0-9_]*)\.out\.(\w+)", arg)
    if ref is None:
        return  # a root (opt(params...)) or ANNOTATE's gene model; no subworkflow emit to check
    sub, name = ref.group(1), ref.group(2)
    includes = _includes()
    assert sub in includes, f"'{key}': {sub} is not an included subworkflow"
    emits = _emits(includes[sub])
    assert name in emits, f"'{key}': {sub} ({includes[sub].name}) does not emit '{name}'"
    owners = _out_spec_owner(key)
    if owners:
        assert sub in owners, f"'{key}' is published by {owners}, but main.nf takes it from {sub}"
    if key in PUBLISHED_FILES:
        base = PurePosixPath(PUBLISHED_FILES[key]).name
        assert base in emits[name] or base in arg, (
            f"'{key}': neither {sub}.{name} ({emits[name]!r}) nor main.nf picks {base}"
        )
