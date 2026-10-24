"""Run Cascadia 0.0.7 with pyteomics pointed at a local Unimod table (M10-08, deliverable 4).

Cascadia 0.0.7 downloads ``http://www.unimod.org/xml/unimod_tables.xml`` at the end of every run
(pyteomics ``mass.unimod.Unimod()`` with no path); R6 forbids that. Inside ``giap-cascadia``:

    micromamba run -n giap-cascadia python -I scripts/pilot/cascadia_offline.py \\
        UNIMOD_XML sequence MZML MODEL -o OUTPREFIX -b 32

(a) refuses a ``UNIMOD_XML`` that is not an existing local file (exit 2, Cascadia not imported);
(b) sets ``pyteomics.mass.unimod._unimod_xml_download_url`` to its absolute path;
(c) hands the rest of ``argv`` unchanged to ``cascadia.cascadia.main()``.
The only GIAP file that imports Cascadia (Apache-2.0); never imported from ``src/`` or
``validation/``.
"""

from __future__ import annotations

import os
import sys
from collections.abc import Sequence


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        print("usage: cascadia_offline.py UNIMOD_XML CASCADIA_ARGS...", file=sys.stderr)
        return 2
    table, rest = args[0], args[1:]
    if "://" in table or not os.path.isfile(table):
        print(f"cascadia_offline: not an existing local file: {table!r}", file=sys.stderr)
        return 2
    import pyteomics.mass.unimod as unimod  # type: ignore[import-not-found]

    unimod._unimod_xml_download_url = os.path.abspath(table)
    # After the URL is set: Cascadia may read it at import.
    from cascadia import cascadia  # type: ignore[import-not-found]

    sys.argv = [cascadia.__file__ or "cascadia", *rest]
    cascadia.main()
    return 0


if __name__ == "__main__":
    sys.exit(main())
