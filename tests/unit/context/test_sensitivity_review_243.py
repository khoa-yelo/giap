"""M6-10 review #243 B2 at unit level: ``Variant.families`` and ``Variant.patch_only``
(docs/specs/M6-10-sensitivity-operon-flags.md, lead decision 8; GIAP_SPEC.md §8.6).

Interface (lead decision 8, B2 option a): ``Variant`` gains ``families: tuple[str, ...] | None``
(``None`` = every family, the default); the ``pair_threshold=*`` variants carry
``("operon_flags",)``. ``Variant.patch_only`` is true exactly when the variant has a config patch
and changes nothing else (no overrides, no single-peptide widening). Step-level checks are in
tests/integration/test_sensitivity_review_243_step.py.
"""

from __future__ import annotations

import dataclasses
from typing import Any

import pytest
from helpers import sensitivity_operons_data as so

PATCH = {so.PATCH_KEY: 0.3}


def _S() -> Any:
    from giap.context import sensitivity as S

    return S


def test_variant_has_a_families_field_defaulting_to_none() -> None:
    S = _S()
    names = [f.name for f in dataclasses.fields(S.Variant)]
    assert "families" in names
    assert S.Variant("missing=min", {"abundance.missing": "min"}).families is None


@pytest.mark.parametrize("tier", ["mayu", "two_peptide"])
def test_pair_threshold_variants_carry_operon_flags_only(tier: str) -> None:
    S = _S()
    vs = S.variants(so.config(output={"sensitivity": "full"}), tier)
    pt = [v for v in vs if v.name.startswith("pair_threshold=")]
    assert [v.name for v in pt] == [so.variant_name(t) for t in so.UNIOP_THRESHOLDS]
    for v in pt:
        assert v.families == ("operon_flags",), v.name
        assert isinstance(v.families, tuple)


@pytest.mark.parametrize("level", ["basic", "full"])
@pytest.mark.parametrize("tier", ["mayu", "two_peptide"])
def test_other_variants_run_for_every_family(level: str, tier: str) -> None:
    """Every variant that changes the abundances (``missing=*``, ``single_peptide_detected``,
    ``top3``) keeps ``families = None``."""
    S = _S()
    vs = S.variants(so.config(output={"sensitivity": level}), tier)
    others = [v for v in vs if not v.name.startswith("pair_threshold=")]
    assert others
    for v in others:
        assert v.families is None, v.name


@pytest.mark.parametrize(
    ("kw", "want"),
    [
        ({"overrides": {}, "config_patch": PATCH}, True),
        ({"overrides": {"abundance.missing": "exclude"}, "config_patch": PATCH}, False),
        ({"overrides": {}, "config_patch": PATCH, "single_peptide_detected": True}, False),
        ({"overrides": {"abundance.missing": "exclude"}}, False),
        ({"overrides": {}, "single_peptide_detected": True}, False),
        ({"overrides": {}}, False),
    ],
    ids=["patch", "patch+override", "patch+single", "override", "single", "nothing"],
)
def test_patch_only_property(kw: dict[str, Any], want: bool) -> None:
    S = _S()
    assert isinstance(S.Variant.__dict__.get("patch_only"), property)
    v = S.Variant("v", **kw)
    assert v.patch_only is want


def test_patch_only_matches_the_shipped_variants() -> None:
    S = _S()
    for v in S.variants(so.config(output={"sensitivity": "full"}), "two_peptide"):
        assert v.patch_only is v.name.startswith("pair_threshold="), v.name


def test_patch_only_variant_genes_is_the_primary_table() -> None:
    """N6: ``variant_genes`` and ``patch_only`` agree: a patch-only variant gives the primary
    table back, even with imputed genes (a ``missing`` rerun would change them)."""
    import pandas as pd

    S = _S()
    cfg = so.config()
    c2 = pd.DataFrame({
        "locus_tag": [f"G{i}" for i in range(12)], "gene": None, "product": "p",
        "feature_type": "CDS", "contig": "c1", "start": [i * 100 for i in range(12)],
        "end": [i * 100 + 90 for i in range(12)], "strand": "+", "partial": False,
        "crosses_origin": False, "hypothetical": False,
    })  # fmt: skip
    gt = so.gene_table(c2, imputed=(1, 2))
    v = S.Variant("p", {}, config_patch=PATCH, families=("operon_flags",))
    assert v.patch_only
    pd.testing.assert_frame_equal(S.variant_genes(gt, v, cfg), gt)
