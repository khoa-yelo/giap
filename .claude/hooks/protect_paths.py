#!/usr/bin/env python3
"""Block edits to protected files listed in ci/protected.txt (GIAP_SPEC.md R3, §14.6)."""
import fnmatch, json, os, sys

event = json.load(sys.stdin)
tool_input = event.get("tool_input") or {}
path = tool_input.get("file_path") or tool_input.get("notebook_path") or ""
root = os.environ.get("CLAUDE_PROJECT_DIR", os.getcwd())
if path:
    rel = os.path.relpath(os.path.abspath(path), root)
    try:
        with open(os.path.join(root, "ci", "protected.txt")) as fh:
            globs = [l.split("#")[0].strip() for l in fh if l.split("#")[0].strip()]
    except FileNotFoundError:
        globs = []
    # New files may be created; existing protected files may not be changed.
    if os.path.exists(os.path.join(root, rel)) and any(fnmatch.fnmatch(rel, g) for g in globs):
        print(f"Blocked: {rel} is protected (GIAP_SPEC.md R3). Draft an ADR and ask Khoa.", file=sys.stderr)
        sys.exit(2)
sys.exit(0)
