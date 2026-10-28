"""Per-step config slices and their task hashes (M7-15 lead decision 7; review #278 B2, lead
decision 11): ``config_slice`` and ``config_hashes`` in ``gippa.steps.registry``.

* ``threads`` changes no hash (resources only, §6.6);
* ``variants.*`` and ``sample_name`` change only the ``variants`` hash;
* ``novel_orfs.*`` changes only the ``novel_orfs`` hash;
* any other section changes every hash;
* the hash is the same across two dumps of one config, and does not depend on key order.
"""

from __future__ import annotations

import copy
import hashlib
import json
from typing import Any

import pytest
from hypothesis import given
from hypothesis import strategies as st

from gippa.config import resolve_config
from gippa.steps.registry import config_hashes, config_slice

STEPS = ["variants", "novel_orfs", "gene_table", "detection", "mapping", "report"]


def _dumped(*overrides: str) -> dict[str, Any]:
    dumped: dict[str, Any] = resolve_config(overrides=list(overrides)).model_dump(mode="json")
    return dumped


def _changed(before: dict[str, Any], after: dict[str, Any]) -> set[str]:
    a, b = config_hashes(before, STEPS), config_hashes(after, STEPS)
    assert set(a) == set(b) == set(STEPS)
    return {n for n in STEPS if a[n] != b[n]}


def test_hash_is_sha256_hex_per_step() -> None:
    hashes = config_hashes(_dumped(), STEPS)
    assert list(hashes) == STEPS
    for h in hashes.values():
        assert len(h) == 64 and int(h, 16) >= 0


def test_hash_is_the_canonical_json_of_the_slice() -> None:
    d = _dumped()
    want = hashlib.sha256(
        json.dumps(config_slice(d, "variants"), sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    assert config_hashes(d, ["variants"]) == {"variants": want}


def test_two_dumps_of_one_config_hash_identically() -> None:
    assert config_hashes(_dumped(), STEPS) == config_hashes(_dumped(), STEPS)


def test_key_order_does_not_change_the_hash() -> None:
    d = _dumped()
    rev = {k: d[k] for k in reversed(list(d))}
    rev["variants"] = {k: d["variants"][k] for k in reversed(list(d["variants"]))}
    assert config_hashes(rev, STEPS) == config_hashes(d, STEPS)


def test_threads_changes_no_hash() -> None:
    d = _dumped()
    e = copy.deepcopy(d)
    e["threads"] = d["threads"] + 7
    assert _changed(d, e) == set()
    for n in STEPS:
        assert "threads" not in config_slice(d, n)


@pytest.mark.parametrize(
    ("section", "key"), [("variants", "high_min_peptides"), ("novel_orfs", "min_codons")]
)
def test_a_step_section_changes_only_that_step(section: str, key: str) -> None:
    d = _dumped()
    e = copy.deepcopy(d)
    e[section][key] = d[section][key] + 1
    assert _changed(d, e) == {section}


def test_variants_change_via_override_reruns_only_variants() -> None:
    """Example 3 of the spec (item 3): ``variants.high_min_peptides`` re-runs only VARIANTS."""
    base = _dumped()
    n = base["variants"]["high_min_peptides"] + 1
    assert _changed(base, _dumped(f"variants.high_min_peptides={n}")) == {"variants"}


def test_sample_name_changes_only_variants() -> None:
    d = _dumped()
    e = copy.deepcopy(d)
    e["sample_name"] = "other_sample"
    assert _changed(d, e) == {"variants"}


@pytest.mark.parametrize("section", ["seed", "mapping", "detection", "annotation", "stats"])
def test_any_other_section_changes_every_step(section: str) -> None:
    d = _dumped()
    e = copy.deepcopy(d)
    e[section] = {"changed": True} if isinstance(d[section], dict) else 12345
    assert _changed(d, e) == set(STEPS)


def test_slice_contents() -> None:
    d = _dumped()
    v, n, g = (config_slice(d, s) for s in ("variants", "novel_orfs", "gene_table"))
    assert "variants" in v and "sample_name" in v and "novel_orfs" not in v
    assert "novel_orfs" in n and "variants" not in n and "sample_name" not in n
    assert not {"variants", "novel_orfs", "sample_name", "threads"} & set(g)
    shared = set(d) - {"variants", "novel_orfs", "sample_name", "threads"}
    for s in (v, n, g):
        assert shared <= set(s)
        assert all(s[k] == d[k] for k in s)


def test_slice_does_not_mutate_its_input() -> None:
    d = _dumped()
    before = copy.deepcopy(d)
    config_slice(d, "variants")
    config_hashes(d, STEPS)
    assert d == before


def test_unknown_top_level_key_is_in_every_slice() -> None:
    """Over-including is safe (a change re-runs too much, never too little)."""
    d = {**_dumped(), "future_section": {"x": 1}}
    for n in STEPS:
        assert config_slice(d, n)["future_section"] == {"x": 1}


@given(
    st.dictionaries(
        st.sampled_from(["seed", "mapping", "variants", "novel_orfs", "sample_name", "threads"]),
        st.integers(-5, 5),
        min_size=1,
    ),
    st.sampled_from(STEPS),
)
def test_property_slice_is_a_sub_dict_without_threads(d: dict[str, int], step: str) -> None:
    s = config_slice(d, step)
    assert "threads" not in s
    assert all(d[k] == s[k] for k in s)
    for k in d:
        mine = {"variants": {"variants"}, "novel_orfs": {"novel_orfs"}}.get(k)
        if k == "sample_name":
            mine = {"variants"}
        expected = k != "threads" and (mine is None or step in mine)
        assert (k in s) == expected, (k, step)
