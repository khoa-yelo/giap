"""Regenerate src/gippa/schemas/config.schema.json from the pydantic model.

Usage: uv run python scripts/write_config_schema.py [--check]
--check exits 1 when the committed file is out of date.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from gippa.config import json_schema

TARGET = Path(__file__).resolve().parents[1] / "src" / "gippa" / "schemas" / "config.schema.json"


def render() -> str:
    return json.dumps(json_schema(), indent=2, sort_keys=True, allow_nan=False) + "\n"


def main(argv: list[str]) -> int:
    text = render()
    if "--check" in argv:
        current = TARGET.read_text(encoding="utf-8") if TARGET.exists() else ""
        if current != text:
            print(f"{TARGET} is out of date; run scripts/write_config_schema.py", file=sys.stderr)
            return 1
        return 0
    TARGET.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
