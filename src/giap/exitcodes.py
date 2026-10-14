"""Process exit codes of the ``giap`` command line (GIAP_SPEC.md §5.1)."""

from __future__ import annotations

from enum import IntEnum


class ExitCode(IntEnum):
    SUCCESS = 0
    INVALID_INPUT = 1  # invalid input (preflight errors, bad configuration)
    ENVIRONMENT = 2  # a tool, database or model is missing
    WARNINGS_STRICT = 3  # warnings were raised and --strict was given
    STEP_FAILED = 4  # a pipeline step failed
    INTERNAL = 5  # internal error (a bug; the message asks to report it)
