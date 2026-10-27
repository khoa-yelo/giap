"""M9-16b behaviour 8 and lead decision 4: the two #178 guards in ``gippa.report`` go once #178
(M8-10) is in the same tree as #228 (M9-16) — whichever merges second removes them.

``"determinants_model" in PATHS`` (build.py, content.py) and ``getattr(gippa.determinants,
"CAVEATS", ())`` (content.py) exist only so the report runs without M8-10. While
``gippa.determinants.CAVEATS`` and the ``determinants_model`` RunView key are both absent, the
guards must stay (the report must still build); this test then skips. Once both exist, the guards
must be gone. Reads the source files only.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

GUARDS = re.compile(r'"determinants_model" in PATHS|getattr\(gippa\.determinants')


def _report_sources() -> dict[str, str]:
    import gippa.report

    root = Path(gippa.report.__file__).parent
    return {p.name: p.read_text(encoding="utf-8") for p in sorted(root.glob("*.py"))}


def test_the_178_guards_are_removed_once_m8_10_is_in_the_tree() -> None:
    import gippa.determinants
    from gippa.io.run_view import PATHS

    has_178 = hasattr(gippa.determinants, "CAVEATS") and "determinants_model" in PATHS
    if not has_178:
        pytest.skip("#178 (M8-10) is not in this tree: the guards stay (lead decision 4)")
    found = {name: GUARDS.findall(text) for name, text in _report_sources().items()}
    assert not any(found.values()), f"#178 guards still present: {found}"


def test_partial_178_is_not_a_silent_state() -> None:
    """Either both M8-10 names exist or neither does; one without the other means a half merge."""
    import gippa.determinants
    from gippa.io.run_view import PATHS

    assert hasattr(gippa.determinants, "CAVEATS") == ("determinants_model" in PATHS)
