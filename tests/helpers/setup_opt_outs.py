"""Explicit ``gippa setup`` opt-outs for tests that want "nothing but X" (M11-15c2).

M11-15c2 (docs/specs/M11-15-portable-clusters.md item 6) turns antiSMASH, GO, KEGG, COG,
contaminants and the images on by default and adds ``--no-go``, ``--no-kegg``, ``--no-cog``,
``--no-contaminants`` and ``--no-images``. Tests that relied on "off by default" state the items
they do not want instead. Before 15c2 those ``--no-*`` options do not exist (and the items are off
anyway), so :func:`opt_outs` keeps only the forms the current ``setup`` command accepts: the same
argv then means the same set of items on either side of the change.
"""

from __future__ import annotations


def _setup_options() -> set[str]:
    from gippa.cli import cli

    cmd = cli.commands["setup"]  # type: ignore[attr-defined]
    return {o for p in cmd.params for o in (*p.opts, *p.secondary_opts)}


def opt_outs(*names: str) -> list[str]:
    """``--no-<name>`` for each name the ``setup`` command accepts (``_`` becomes ``-``)."""
    known = _setup_options()
    flags = ["--no-" + n.replace("_", "-") for n in names]
    return [f for f in flags if f in known]


#: no resource download and no antiSMASH, whatever the defaults
NO_RESOURCE_FLAGS: tuple[str, ...] = ("antismash", "go", "kegg", "cog", "contaminants")
