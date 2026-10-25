"""Run Cascadia 0.0.7 with every Unimod lookup pointed at a local table (M10-08, deliverable 4).

Cascadia 0.0.7 reaches ``http://www.unimod.org`` twice per run; R6 forbids both:
at the start, ``PeptideTokenizer.from_massivekb()`` -> ``pyteomics.proforma`` -> psims
``obo_cache.resolve(UNIMOD_OBO_URL)`` -> ``Unimod()`` -> ``lxml.etree.parse(<http URL>)``; at the
end, pyteomics ``mass.unimod.Unimod()`` with no path. Inside ``giap-cascadia``:

    micromamba run -n giap-cascadia python -I scripts/pilot/cascadia_offline.py \\
        UNIMOD_XML sequence MZML MODEL -o OUTPREFIX -b 32

(a) refuses a ``UNIMOD_XML`` that is not an existing local file (exit 2, Cascadia not imported);
(b) refuses every Python-level IP socket connection (AF_UNIX stays open for torch's DataLoader
    workers; defence in depth; libxml2's C I/O is not covered);
(c) registers a psims resolver for the Unimod OBO URL that reads the local table;
(d) sets ``pyteomics.mass.unimod._unimod_xml_download_url`` to its absolute path;
(e) hands the rest of ``argv`` unchanged to ``cascadia.cascadia.main()``.
The only GIAP file that imports Cascadia (Apache-2.0); never imported from ``src/`` or
``validation/``.
"""

from __future__ import annotations

import os
import socket
import sys
from collections.abc import Sequence
from typing import NoReturn


def _refuse(*_args: object, **_kwargs: object) -> NoReturn:
    raise OSError("cascadia_offline: network access is disabled (R6)")


def _guard_sockets() -> None:
    """Refuse IP connections; keep AF_UNIX, which torch's DataLoader workers use to share tensors
    (job 1889042 failed in ``rebuild_storage_fd`` with every connect refused)."""
    connect, connect_ex = socket.socket.connect, socket.socket.connect_ex

    def _connect(self: socket.socket, address: object) -> None:
        if getattr(socket, "AF_UNIX", None) is not None and self.family == socket.AF_UNIX:
            return connect(self, address)  # type: ignore[arg-type]
        _refuse()

    def _connect_ex(self: socket.socket, address: object) -> int:
        if getattr(socket, "AF_UNIX", None) is not None and self.family == socket.AF_UNIX:
            return connect_ex(self, address)  # type: ignore[arg-type]
        _refuse()

    socket.create_connection = _refuse
    socket.socket.connect = _connect  # type: ignore[method-assign]
    socket.socket.connect_ex = _connect_ex  # type: ignore[method-assign]


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        print("usage: cascadia_offline.py UNIMOD_XML CASCADIA_ARGS...", file=sys.stderr)
        return 2
    table, rest = args[0], args[1:]
    if "://" in table or not os.path.isfile(table):
        print(f"cascadia_offline: not an existing local file: {table!r}", file=sys.stderr)
        return 2
    table_abs = os.path.abspath(table)
    _guard_sockets()
    import pyteomics.mass.unimod as unimod  # type: ignore[import-not-found]
    from psims.controlled_vocabulary import unimod as psims_unimod  # type: ignore[import-not-found]
    from psims.controlled_vocabulary.controlled_vocabulary import (  # type: ignore[import-not-found]
        obo_cache,
    )

    obo_cache.set_resolver(
        psims_unimod.UNIMOD_OBO_URL, lambda _cache: psims_unimod.Unimod(None, table_abs)
    )
    unimod._unimod_xml_download_url = table_abs
    # After both are set: Cascadia may read them at import.
    from cascadia import cascadia  # type: ignore[import-not-found]

    sys.argv = [cascadia.__file__ or "cascadia", *rest]
    cascadia.main()
    return 0


if __name__ == "__main__":
    sys.exit(main())
