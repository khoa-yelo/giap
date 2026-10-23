"""Lead decision 10 (M9-24): a title that starts with a declared phrase passes the "title states
the finding" item without a number; V08's no-abundance title is the one declared phrase."""

from __future__ import annotations

from types import SimpleNamespace

from giap.viz import checklist, detection


def _inst(title: str) -> SimpleNamespace:
    return SimpleNamespace(title=title)


def test_v08_no_abundance_title_is_declared() -> None:
    assert detection.NO_ABUNDANCE in checklist.DECLARED_TITLES
    t = f"{detection.NO_ABUNDANCE}; no COG annotation (layer_missing:cog)"
    ok, _ = checklist._title(_inst(t), empty=False)  # type: ignore[arg-type]
    assert ok


def test_other_numberless_titles_still_fail() -> None:
    ok, _ = checklist._title(_inst("Detected fraction per COG category"), empty=False)  # type: ignore[arg-type]
    assert not ok
