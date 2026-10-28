#!/usr/bin/env python3
"""The GIAP -> GIPPA rename (ADR-0064 decision 2; inventory docs/specs/rename-gippa.md).

Deterministic and idempotent: it reads the tracked files of a git checkout (``git ls-files``)
and can be re-run on any later main. Three modes::

    python scripts/dev/rename_gippa.py --plan  [--repo DIR] [-v]   # counts, moves, .github list
    python scripts/dev/rename_gippa.py --apply [--repo DIR]        # content map, then git mv
    python scripts/dev/rename_gippa.py --check [--repo DIR] [-v]   # exit 1 on any leftover

Content map: plain case-exact substrings ``giap``->``gippa``, ``GIAP``->``GIPPA``,
``Giap``->``Gippa`` (no ``\\b``: ``"\\tgiap"`` and ``b"giap"`` must change too), outside the keep
list and the protected spans below. Paths: every tracked path with ``giap`` in a component,
outside the path keep list, is moved with ``git mv`` (shallowest component first, so
``src/giap`` and ``tests/sim/giapsim`` move as whole directories).

Q-K32 (ADR-0064 amendment): every document is rewritten, the spec moves to ``GIPPA_SPEC.md``
and the keep list holds data only, plus the rename record (ADR-0064, docs/specs/rename-gippa.md).
The old expansion becomes "genome integrated prokaryotic proteomics analysis" everywhere.

Never edited: the keep paths (``KEEP_PATHS``, every ``PROVENANCE.md`` / ``*.PROVENANCE.md``),
``.github/`` (listed for Khoa only), binary files (a NUL byte, not UTF-8, or a ``BINARY_EXT``
suffix such as ``.gz``/``.parquet``/``.png``/``.pdf``/``.sif``) and sequence files
(``SEQUENCE_EXT``: FASTA and GenBank content is never touched, because GIAP is a tetrapeptide).

Protected spans inside renamed files (kept as written; ``--check`` reports them by category):
the checkout path ``bios/270/khoa/giap`` (it
does not move, decision 5), the recording marker ``.giap-complete`` (on-disk under ``$GIAP_ROOT``,
whose contents do not move), the setup lock ``.giap_setup.lock`` (old and new ``setup`` must lock
each other out on a shared ``$GIAP_DB_DIR``), the SVG salt ``"giap"`` (decision 7), uppercase
``GIAP`` inside a longer uppercase run (a possible peptide; ``GIAPORF``/``GIAPNOVEL``/``GIAPSIM``
are identifiers and are renamed), sequence-like cells of ``seq``/``peptide`` columns in TSV/CSV,
and the ADR-0064 transition aliases that ``TRANSITIONS`` adds (old ``GIAP_*`` exports in
farmshare_env.sh, ``giap-*`` job counting in the R14 checks and probe_env.py, ``.giap/`` in
.gitignore).

Every transition asserts its end state and raises ``TransitionError`` when its anchor no longer
matches (a later main changed the line), so ``--apply`` refuses and ``--check`` reports it.
``--check`` also fails on any case-insensitive ``giap`` left in the documentation (``is_doc``:
the root docs, docs/ and every tracked ``*.md`` / ``*.rst``) outside the rename record, the data
keep paths and the kept data tokens (``DOC_ALLOWED``).
"""

from __future__ import annotations

import argparse
import re
import subprocess
from collections import Counter
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath

REPO = Path(__file__).resolve().parents[2]
ANY = re.compile(r"giap", re.I)

# ---------------------------------------------------------------- keep list (ADR-0064 decision 3)

#: (pattern, category). A pattern ending in "/" is a directory prefix; else an exact path or regex.
#: Q-K32: data only, plus the rename record. Documentation (spec, ADRs, docs/specs, PROGRESS.md,
#: backlog, docs/validation, validation/reports) is rewritten. validation/release/ holds only
#: .gitkeep: a file added there is rewritten unless it is sha256-pinned (then list it here).
KEEP_PATHS: tuple[tuple[str, str], ...] = (
    ("docs/adr/ADR-0064-rename-giap-to-gippa.md", "rename-record"),  # maps the old name
    ("docs/specs/rename-gippa.md", "rename-record"),  # its inventory: giap->gippa tables
    ("tests/data/recorded/", "recorded-tool-data"),
    ("tests/data/uniop/", "recorded-tool-data"),  # verbatim UniOP input/output excerpts (spot 3)
    ("tests/data/golden/r/VERSIONS.txt", "recorded-tool-data"),  # the R env as recorded
    ("validation/pilot/REPORT.md", "pilot-artefact"),
    ("validation/pilot/metrics/", "pilot-artefact"),
    ("validation/pilot/db_a2/", "pilot-artefact"),
    ("validation/pilot/inputs/", "pilot-artefact"),
    ("validation/pilot/diann_A1.cfg", "pilot-artefact"),  # sha256-pinned in reproduce.md
    ("validation/pilot/diann_A2.cfg", "pilot-artefact"),
    (".github/", "github-khoa"),
    ("scripts/dev/rename_gippa.py", "rename-tooling"),  # the map itself
    ("scripts/dev/rename_compare.py", "rename-tooling"),
    ("tests/unit/test_rename_compare.py", "rename-tooling"),
    ("tests/unit/test_rename_gippa_tooling.py", "rename-tooling"),
    ("tests/unit/test_rename_gippa_docs.py", "rename-tooling"),
)
#: keep paths whose *path* also stays (all other keep categories only keep their content)
PATH_KEEP_CATEGORIES = frozenset(
    {"rename-record", "recorded-tool-data", "github-khoa", "pilot-artefact", "rename-tooling"}
)  # fmt: skip
BINARY_EXT = frozenset(
    {".gz", ".bgz", ".bz2", ".xz", ".zst", ".zip", ".tar", ".parquet", ".png", ".pdf", ".sif",
     ".npz", ".npy", ".pkl", ".pickle", ".jpg", ".jpeg", ".gif", ".ico", ".woff", ".woff2",
     ".ttf", ".otf", ".bam", ".cram", ".bai", ".xlsx", ".h5", ".hdf5", ".sqlite", ".db"}
)  # fmt: skip
SEQUENCE_EXT = frozenset(
    {".fa", ".faa", ".fasta", ".fas", ".fna", ".ffn", ".frn", ".gbff", ".gbk", ".gb", ".embl"}
)
#: never sniffed as FASTA: a page may open with a blockquote (``> **Note** ...``)
DOC_SUFFIXES = frozenset({".md", ".yml", ".yaml", ".cff", ".txt", ".rst", ".html"})

# ---------------------------------------------------------------- protected spans


@dataclass(frozen=True)
class Protect:
    """A regex whose matches stay as written, in every file or in the files ``paths`` names
    (new, post-rename paths)."""

    pattern: re.Pattern[str]
    category: str
    paths: tuple[str, ...] = ()

    def applies(self, path: str) -> bool:
        return not self.paths or path in self.paths


def _p(rx: str, flags: int = 0) -> re.Pattern[str]:
    return re.compile(rx, flags)


#: an uppercase GIAP that is not a whole word or a known identifier prefix: a possible peptide
PEPTIDE = _p(r"[A-Z]+GIAP[A-Z]*|GIAP(?!ORF|NOVEL|SIM)[A-Z]+")
ADR64_FILE = "ADR-0064-rename-giap-to-gippa"  # the kept rename record's file name
GLOBAL_PROTECT: tuple[Protect, ...] = (
    Protect(_p(re.escape(ADR64_FILE)), "link to the kept ADR-0064 file"),
    Protect(_p(r"bios([/-])270\1khoa\1giap(?![A-Za-z0-9_])"), "checkout path (does not move)"),
    Protect(_p(r"\.giap-complete"), "recording marker under $GIAP_ROOT (does not move)"),
    Protect(_p(r"\.giap_setup\.lock"), "setup lock (old and new setup exclude each other)"),
    Protect(PEPTIDE, "possible peptide (review)"),
)

PLOT_FILES = "src/gippa/io/plot_files.py"
STYLE = "src/gippa/viz/style.py"
SALT_COMMENT = "#: the pre-rename salt stays: a new salt changes every SVG id (ADR-0064 item 7)"
FARMSHARE = "scripts/farmshare_env.sh"
PROBE = "scripts/probe_env.py"
ALIAS_MARK = "# ADR-0064 decision 4: the old names, same values, for jobs in flight, old worktrees"
ALIAS_MARK2 = "# and the installed scrontab; removed after the first publish."
R14_NOTE = "(R14 counts `giap-*` jobs too until none is left, ADR-0064)"

FILE_PROTECT: tuple[Protect, ...] = (
    Protect(_p(r'HASHSALT = "giap"'), "svg hashsalt (decision 7)", (PLOT_FILES,)),
    Protect(_p(r'"svg\.hashsalt": "giap"'), "svg hashsalt (decision 7)", (STYLE,)),
    Protect(_p(r'\["svg\.hashsalt"\] == "giap"'), "svg hashsalt (decision 7)",
            ("tests/unit/viz/test_style.py",)),
    Protect(_p(r"\$\{GIPPA_[A-Z_]+:-\$\{GIAP_[A-Z_]+:-"), "transition: old override honoured",
            (FARMSHARE,)),
    Protect(_p(r"envs/gippa-wf/bin:\$MAMBA_ROOT_PREFIX/envs/giap-wf/bin"),
            "transition: old wf env on PATH", (FARMSHARE,)),
    Protect(_p(r'^export GIAP_[A-Z_]+="\$GIPPA_[A-Z_]+"$', re.M), "transition: old GIAP_* export",
            (FARMSHARE,)),
    Protect(_p(r"\(giap\|gippa\)-"), "transition: R14 counts giap-* too",
            ("scripts/resume_status.sh", "scripts/session_start.sh",
             "scripts/slurm/benchmark_standard.sbatch")),
    Protect(_p(r'"\$GIPPA_LOG_DIR"/giap-\*\.out'), "transition: R14 counts giap-* too",
            ("scripts/resume_status.sh",)),
    Protect(_p(r"giap-\*\|gippa-\*\)"), "transition: R14 counts giap-* too",
            ("scripts/slurm/calibration_full_launch.sh", "scripts/slurm/calibration_sharded.sh")),
    Protect(_p(re.escape(R14_NOTE)), "transition: R14 counts giap-* too", ("AGENTS.md",)),
    Protect(_p(r'\(" giap-", " gippa-"\)'), "transition: R14 counts giap-* too", (PROBE,)),
    Protect(_p(r"^\.giap/$", re.M), "transition: old probe dir stays ignored", (".gitignore",)),
    Protect(_p(r'old="GIAP_\$\{name#GIPPA_\}"'), "transition: kept pilot cfgs name old vars",
            ("scripts/slurm/diann_pilot.sbatch",)),
    Protect(_p(r'key\.replace\("GIPPA_", "GIAP_", 1\)'),
            "transition: kept pilot cfgs name old vars", ("tests/helpers/pilot_kit.py",)),
    Protect(_p(r"^- \[ADR-0064:.*$", re.M), "ADR-0064 index line (the rename record)",
            ("docs/adr/index.md",)),
)  # fmt: skip

# ---------------------------------------------------------------- transitions (idempotent edits)


class TransitionError(ValueError):
    """A transition's anchor no longer matches and its end state is absent."""


def _require(text: str, *posts: str) -> str:
    for post in posts:
        if post not in text:
            raise TransitionError(f"end state missing (anchor changed? update TRANSITIONS): {post}")
    return text


def _salt_comment(text: str) -> str:
    line = 'HASHSALT = "giap"'
    if line in text and SALT_COMMENT not in text:
        text = text.replace(line, f"{SALT_COMMENT}\n{line}", 1)
    return _require(text, f"{SALT_COMMENT}\n{line}")


def _style_comment(text: str) -> str:
    line = '"svg.hashsalt": "giap",'
    note = "  # the pre-rename salt stays (ADR-0064 item 7)"
    if line in text and line + note not in text:
        text = text.replace(line, line + note, 1)
    return _require(text, line + note)


_DEFAULT = re.compile(r'^([ \t]*export GIPPA_([A-Z_]+)=")\$\{GIPPA_\2:-(?!\$\{GIAP_)(.*)\}"', re.M)
_EXPORT = re.compile(r"^[ \t]*export (GIPPA_[A-Z_]+)=", re.M)
_ANY_DEFAULT = re.compile(r"\$\{GIPPA_([A-Z_]+):-(?!\$\{GIAP_\1:-)")


def _farmshare(text: str) -> str:
    # an old GIAP_X set by the caller (old worktree, job in flight) is still honoured
    text = _DEFAULT.sub(lambda m: f'{m[1]}${{GIPPA_{m[2]}:-${{GIAP_{m[2]}:-{m[3]}}}}}"', text)
    wf = "$MAMBA_ROOT_PREFIX/envs/gippa-wf/bin:"
    old = "$MAMBA_ROOT_PREFIX/envs/giap-wf/bin:"
    if wf in text and wf + old not in text:
        text = text.replace(wf, wf + old, 1)  # until build_envs.sbatch has built gippa-wf
    lines = text.split("\n")
    if ALIAS_MARK in lines:  # drop the old block; it is rebuilt from the current exports
        i = lines.index(ALIAS_MARK)
        j = i + 1
        while j < len(lines) and (lines[j] == ALIAS_MARK2 or lines[j].startswith("export GIAP_")):
            j += 1
        del lines[i:j]
    names = list(dict.fromkeys(_EXPORT.findall("\n".join(lines))))
    block = [ALIAS_MARK, ALIAS_MARK2]
    block += [f'export {n.replace("GIPPA_", "GIAP_", 1)}="${n}"' for n in names]
    end = len(lines) - 1 if lines and lines[-1] == "" else len(lines)
    lines[end:end] = block
    text = "\n".join(lines)
    if m := _ANY_DEFAULT.search(text):
        raise TransitionError(f"GIPPA_{m[1]} default without the old GIAP_{m[1]} fallback")
    return _require(text, wf + old) if wf in text else text


def _ensure(old: str, new: str) -> Callable[[str], str]:
    """Replace ``old`` by ``new`` unless ``new`` is there; raise when neither is."""

    def edit(text: str) -> str:
        return _require(text.replace(old, new) if new not in text else text, new)

    return edit


def _gitignore(text: str) -> str:
    lines = text.split("\n")
    if ".gippa/" not in lines and ".giap/" in lines:
        lines.insert(lines.index(".giap/"), ".gippa/")
    if ".gippa/" not in lines or ".giap/" not in lines:
        raise TransitionError("end state missing: .gippa/ and .giap/ both ignored")
    return "\n".join(lines)


#: the old expansion in every variant (Q-K32); the loose form is the end-state assertion and is
#: anchored to "genome-integrated", so other prose about proteomics data never trips it
NEW_EXPANSION = "genome integrated prokaryotic proteomics analysis"
_OLD_EXPANSION = re.compile(
    r"genome[-\s]+integrated\s+(?:analysis\s+of\s+proteomics\s+data(?:\s+for\s+bacterial"
    r"\s+genomes)?|annotation\s+of\s+proteomics)",
    re.I,
)
_OLD_LOOSE = re.compile(r"genome[-\s]+integrated(?!\s+prokaryotic\s+proteomics\s+analysis)", re.I)


def _expansion(text: str) -> str:
    def new(m: re.Match[str]) -> str:
        return NEW_EXPANSION.capitalize() if m[0][0].isupper() else NEW_EXPANSION

    text = _OLD_EXPANSION.sub(new, text)
    if m := _OLD_LOOSE.search(text):
        raise TransitionError(f"old expansion left (update _OLD_EXPANSION): {m[0]!r}")
    # "the GIAP → GIPPA rename" reads "GIPPA → GIPPA" after the map
    text = _SELF_ARROW.sub(lambda m: "rename to GIPPA" if m[1] else "GIPPA rename", text)
    # "renames the image prefix (`giap-` → `gippa-`)" -> "renames the image prefix to `gippa-`"
    text = _SELF_PAREN.sub(lambda m: f" to {m[1]}{m[3]}{m[1]}", text)
    for m in _SELF_LOOSE.finditer(text):
        if m[1].lower() == m[2].lower():
            raise TransitionError(f"'GIPPA → GIPPA' left (update _SELF_ARROW): {m[0]!r}")
    return text


_SELF_ARROW = re.compile(r"(rename\s+)GIPPA\s+→\s+GIPPA|GIPPA\s+→\s+GIPPA\s+rename")
_TOKEN = r"\w*gippa[\w*-]*"
_SELF_PAREN = re.compile(
    rf"\s*\(\s*([`'\"]*)({_TOKEN})\1\s*(?:→|->)\s*\1(\2)\1\s*\)", re.I
)  # fmt: skip
#: any "X → X" left after the map (any case, code spans, quotes, -/_ suffixes), compared by token;
#: the right side is a lookahead so chains ("a → b → c") are checked pair by pair
_SELF_LOOSE = re.compile(
    rf"({_TOKEN})[`'\"*,;:]*\s*(?:→|->|=>|\bto\b)\s*[`'\"*]*(?=({_TOKEN}))", re.I
)  # fmt: skip


#: README "Status": the rename is the PR itself, so its pending item goes (post-map wording)
_CLAUSE = (
    ", and the code rename: until it lands, the package, the command and the environment "
    "variables are `gippa` (`gippa run`, `GIPPA_DB_DIR`);"
)
_README_CLAUSE = re.compile(r"\s+".join(re.escape(w) for w in _CLAUSE.split()))


def _readme(text: str) -> str:
    text = _README_CLAUSE.sub(";", text)
    if m := re.search(r"until it lands|the code rename", text):
        raise TransitionError(f"README rename clause left (update _README_CLAUSE): {m[0]!r}")
    return text


_R14_GREP = _ensure("grep ' gippa-'", "grep -E ' (giap|gippa)-'")
_R14_CASE = _ensure('case "$name" in gippa-*)', 'case "$name" in giap-*|gippa-*)')
_PROBE_OLD = (
    """    gippa_jobs = [line for line in queue.splitlines() """
    """if " gippa-" in f" {line.split(' ', 1)[-1]}"]"""
)
_PROBE_NEW = """\
    gippa_jobs = [  # ADR-0064: R14 counts the old job prefix too until none is left
        line
        for line in queue.splitlines()
        if any(p in f" {line.split(' ', 1)[-1]}" for p in (" giap-", " gippa-"))
    ]"""

#: new path -> idempotent edit applied after the substring map (ADR-0064 decisions 4, 5, 7);
#: each one raises ``TransitionError`` when its end state cannot be reached
TRANSITIONS: dict[str, Callable[[str], str]] = {
    PLOT_FILES: _salt_comment,
    STYLE: _style_comment,
    FARMSHARE: _farmshare,
    "scripts/resume_status.sh": lambda t: _ensure(
        '"$GIPPA_LOG_DIR"/gippa-*.out ', '"$GIPPA_LOG_DIR"/gippa-*.out "$GIPPA_LOG_DIR"/giap-*.out '
    )(_R14_GREP(t)),
    "scripts/session_start.sh": _R14_GREP,
    "scripts/slurm/benchmark_standard.sbatch": _ensure("$2 ~ /^gippa-/", "$2 ~ /^(giap|gippa)-/"),
    "scripts/slurm/calibration_full_launch.sh": _R14_CASE,
    "scripts/slurm/calibration_sharded.sh": _R14_CASE,
    PROBE: _ensure(_PROBE_OLD, _PROBE_NEW),
    "AGENTS.md": _ensure("jobs named `gippa-*`,", f"jobs named `gippa-*` {R14_NOTE},"),
    ".gitignore": _gitignore,
    "README.md": _readme,
}  # fmt: skip

# ---------------------------------------------------------------- the map


def map_text(s: str) -> str:
    """The case-exact substring map (no word boundaries)."""
    return s.replace("giap", "gippa").replace("GIAP", "GIPPA").replace("Giap", "Gippa")


def map_path(path: str) -> str:
    return "/".join(map_text(c) for c in path.split("/"))


def keep_category(path: str) -> str | None:
    """The keep category of a tracked path (content stays as written), else None."""
    name = PurePosixPath(path).name
    if name == "PROVENANCE.md" or name.endswith(".PROVENANCE.md"):
        return "provenance"
    for pat, cat in KEEP_PATHS:
        if pat.startswith("re:"):
            if re.match(pat[3:], path):
                return cat
        elif (pat.endswith("/") and path.startswith(pat)) or path == pat:
            return cat
    return None


def path_stays(path: str) -> bool:
    return keep_category(path) in PATH_KEEP_CATEGORIES


def file_kind(path: str, data: bytes) -> str:
    """``binary``, ``sequence`` or ``text``."""
    suffix = PurePosixPath(path).suffix.lower()
    if suffix in BINARY_EXT or b"\0" in data:
        return "binary"
    try:
        data.decode("utf-8")
    except UnicodeDecodeError:
        return "binary"
    sniff = suffix not in DOC_SUFFIXES and not is_doc(path)
    if suffix in SEQUENCE_EXT or (sniff and data.lstrip()[:1] == b">"):
        return "sequence"
    return "text"


_SEQ_COL = re.compile(r"seq|peptide", re.I)
_MODS = re.compile(r"\[[^\]]*\]|\([^)]*\)")
_AA = re.compile(r"[A-Z*]+")


def _table_spans(path: str, text: str) -> list[tuple[int, int, str]]:
    """Sequence-like cells (amino-acid letters, modifications removed) of seq/peptide columns."""
    sep = {".tsv": "\t", ".csv": ","}.get(PurePosixPath(path).suffix.lower())
    if sep is None:
        return []
    spans: list[tuple[int, int, str]] = []
    cols: set[int] | None = None
    pos = 0
    for line in text.splitlines(keepends=True):
        body = line.rstrip("\r\n")
        if body.startswith("#") or not body:
            pos += len(line)
            continue
        cells = body.split(sep)
        if cols is None:
            cols = {i for i, c in enumerate(cells) if _SEQ_COL.search(c)}
        else:
            off = pos
            for i, c in enumerate(cells):
                if i in cols and "GIAP" in c and _AA.fullmatch(_MODS.sub("", c).strip("_.-")):
                    spans.append((off, off + len(c), "sequence cell (never renamed)"))
                off += len(c) + len(sep)
        pos += len(line)
    return spans


def protected_spans(path: str, text: str) -> list[tuple[int, int, str]]:
    """Non-overlapping (start, end, category) spans that stay as written, in text order."""
    found: list[tuple[int, int, str]] = list(_table_spans(path, text))
    for p in (*FILE_PROTECT, *GLOBAL_PROTECT):
        if p.applies(path):
            found += [(m.start(), m.end(), p.category) for m in p.pattern.finditer(text)]
    found.sort(key=lambda s: (s[0], -s[1]))
    out: list[tuple[int, int, str]] = []
    for s in found:
        if s[1] > s[0] and (not out or s[0] >= out[-1][1]):
            out.append(s)
    return out


def transform(path: str, text: str) -> str:
    """The renamed content of a text file whose post-rename path is ``path``."""
    pieces: list[str] = []
    pos = 0
    for start, end, _ in protected_spans(path, text):
        pieces += [_expansion(map_text(text[pos:start])), text[start:end]]
        pos = end
    pieces.append(_expansion(map_text(text[pos:])))  # protected spans never pass _expansion
    out = "".join(pieces)
    edit = TRANSITIONS.get(path)
    return edit(out) if edit else out


# ---------------------------------------------------------------- the repository


def git(repo: Path, *args: str) -> str:
    res = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)
    if res.returncode:
        raise SystemExit(f"git {' '.join(args)}: {res.stderr.strip()}")
    return res.stdout


def tracked(repo: Path) -> list[str]:
    """Tracked regular files (symlinks and submodules excluded), sorted."""
    out = []
    for rec in git(repo, "ls-files", "-s", "-z").split("\0"):
        if not rec:
            continue
        meta, path = rec.split("\t", 1)
        if meta.split()[0] in {"100644", "100755"}:
            out.append(path)
    return sorted(out)


def area(path: str) -> str:
    parts = path.split("/")
    if len(parts) == 1:
        return "(root)"
    if parts[0] in {"src", "tests", "docs", "workflow", "scripts", "validation"} and len(parts) > 2:
        return f"{parts[0]}/{map_path(parts[1])}"
    return parts[0]


@dataclass
class Report:
    changes: dict[str, tuple[str, str]] = field(default_factory=dict)  # path -> (old, new) text
    moves: list[tuple[str, str]] = field(default_factory=list)  # file-level, sorted
    kept: Counter[str] = field(default_factory=Counter)  # category -> occurrences left
    kept_files: dict[str, list[tuple[str, int, str]]] = field(default_factory=dict)
    skipped_binary: list[str] = field(default_factory=list)  # binary files with giap inside
    github: list[tuple[str, int, str, str]] = field(default_factory=list)
    broken: list[tuple[str, str]] = field(default_factory=list)  # (new path, TransitionError)


def _line_hits(text: str) -> list[tuple[int, str]]:
    return [(i + 1, ln) for i, ln in enumerate(text.splitlines()) if ANY.search(ln)]


def survey(repo: Path) -> Report:
    """Everything ``--apply`` would do, and every ``giap`` it leaves, by category."""
    rep = Report()
    for path in tracked(repo):
        new_path = path if path_stays(path) else map_path(path)
        if new_path != path:
            rep.moves.append((path, new_path))
        data = (repo / path).read_bytes()
        quick = data.decode("latin-1")
        if not ANY.search(quick) and not _OLD_LOOSE.search(quick) and new_path not in TRANSITIONS:
            continue
        kind = file_kind(path, data)
        cat = keep_category(path)
        if kind != "text" or cat is not None:
            label = cat or ("binary (never edited)" if kind == "binary" else "sequence file")
            text = data.decode("utf-8", "replace")
            n = len(ANY.findall(text))
            rep.kept[label] += n
            rep.kept_files.setdefault(label, []).append((path, n, ""))
            if kind == "binary":
                rep.skipped_binary.append(path)
            if cat == "github-khoa":
                for no, ln in _line_hits(text):
                    if (mapped := transform(new_path, ln)) != ln:
                        rep.github.append((path, no, ln, mapped))
            continue
        text = data.decode("utf-8")
        try:
            new = transform(new_path, text)
        except TransitionError as e:
            rep.broken.append((new_path, str(e)))
            continue
        if new != text:
            rep.changes[path] = (text, new)
        for start, end, label in protected_spans(new_path, new):
            n = len(ANY.findall(new[start:end]))
            if n:
                rep.kept[label] += n
                line = new.count("\n", 0, start) + 1
                rep.kept_files.setdefault(label, []).append((new_path, line, new[start:end]))
    return rep


def apply(repo: Path) -> int:
    """Content first (at the old paths), then ``git mv``. Returns the number of files written."""
    rep = survey(repo)
    if rep.broken:
        raise SystemExit("\n".join(f"rename_gippa: {p}: {e}" for p, e in rep.broken))
    collisions = [(p, d) for p, d in rep.moves if (repo / d).is_file()]
    if collisions:
        raise SystemExit(
            "\n".join(f"rename_gippa: {d} exists, cannot move {p}" for p, d in collisions)
        )
    for path, (_, new) in sorted(rep.changes.items()):
        (repo / path).write_bytes(new.encode("utf-8"))
    move_paths(repo)
    return len(rep.changes)


def move_paths(repo: Path) -> list[tuple[str, str]]:
    """``git mv`` every renamable path, the shallowest renamable component first."""
    done: list[tuple[str, str]] = []
    while True:
        todo = [p for p in tracked(repo) if not path_stays(p) and map_path(p) != p]
        if not todo:
            return done
        prefixes: dict[str, str] = {}
        for p in todo:
            parts = p.split("/")
            i = next(k for k, c in enumerate(parts) if map_text(c) != c)
            src = "/".join(parts[: i + 1])
            prefixes[src] = "/".join([*parts[:i], map_text(parts[i])])
        for src, dst in sorted(prefixes.items()):
            if (repo / dst).exists():  # never nest a directory into an existing one
                for p in (q for q in todo if q == src or q.startswith(src + "/")):
                    d = dst + p[len(src) :]
                    if (repo / d).exists():  # git mv -k would skip it and loop forever
                        raise SystemExit(f"rename_gippa: {d} exists, cannot move {p}")
                    (repo / d).parent.mkdir(parents=True, exist_ok=True)
                    git(repo, "mv", p, d)
                    done.append((p, d))
            else:
                (repo / dst).parent.mkdir(parents=True, exist_ok=True)
                git(repo, "mv", src, dst)
                done.append((src, dst))


_HEX64 = re.compile(r"\b[0-9a-f]{64}\b")
_SIMREF = re.compile(r"(giap|gippa)sim|\bsim\d+k\b|SIM\d+K", re.I)


def repin_candidates(repo: Path) -> list[tuple[str, int]]:
    """Test and benchmark files outside the keep list that pin 64-hex digests and use the
    simulator: their pins may change where the simulator writes its own name (risky spot 8)."""
    out = []
    for path in tracked(repo):
        if keep_category(path) or not path.endswith(".py"):
            continue
        if not path.startswith(("tests/", "benchmarks/", "validation/")):
            continue
        text = (repo / path).read_text(encoding="utf-8", errors="replace")
        n = len(_HEX64.findall(text))
        if n and _SIMREF.search(text):
            out.append((path, n))
    return out


# ---------------------------------------------------------------- reports


def _print_kept(rep: Report, verbose: bool) -> None:
    for label in ("possible peptide (review)", "sequence cell (never renamed)"):
        for path, no, span in rep.kept_files.get(label, []):
            print(f"REVIEW {label}: {path}:{no}: {span}")
    if verbose:
        print("\n## Kept spans and files (category, path, line or count, span)")
        for cat in sorted(rep.kept_files):
            for path, no, span in rep.kept_files[cat]:
                print(f"{cat}\t{path}\t{no}\t{span}")


def print_plan(rep: Report, repo: Path, verbose: bool) -> None:
    files: Counter[str] = Counter()
    occ: Counter[str] = Counter()
    for path, (old, _) in rep.changes.items():
        files[area(path)] += 1
        occ[area(path)] += len(ANY.findall(old))
    print("## Content changes by area (files / occurrences before the map)")
    for a in sorted(files):
        print(f"{a:40s} {files[a]:6d} {occ[a]:7d}")
    print(f"{'TOTAL':40s} {sum(files.values()):6d} {sum(occ.values()):7d}")
    print("\n## Kept as written (occurrences)")
    for cat, n in sorted(rep.kept.items()):
        print(f"{cat:55s} {n:6d}  ({len(rep.kept_files[cat])} files/spans)")
    print(f"\n## Moves ({len(rep.moves)} files)")
    dirs: Counter[tuple[str, str]] = Counter()
    for src, dst in rep.moves:
        sp, dp = src.split("/"), dst.split("/")
        i = next(k for k, (x, y) in enumerate(zip(sp, dp, strict=True)) if x != y)
        dirs["/".join(sp[: i + 1]), "/".join(dp[: i + 1])] += 1
    for (s_dir, d_dir), n in sorted(dirs.items()):
        print(f"git mv {s_dir} {d_dir}" + (f"   ({n} files)" if n > 1 else ""))
    for path, err in rep.broken:
        print(f"TRANSITION BROKEN {path}: {err}")
    print("\n## .github/ lines for Khoa (agents never edit .github/)")
    for path, no, old, new in rep.github:
        print(f"{path}:{no}:\n  - {old.strip()}\n  + {new.strip()}")
    _print_kept(rep, verbose)
    print("\n## Re-pin candidates (pin 64-hex digests and use the simulator; spot 8)")
    for path, n in repin_candidates(repo):
        print(f"{path} ({n} digests)")


DOC_FILES = frozenset(
    {"README.md", "CITATION.cff", "mkdocs.yml", "AGENTS.md", "CLAUDE.md", "GIPPA_SPEC.md",
     "PROGRESS.md"}
)  # fmt: skip
#: the only giap a document may keep (Q-K32): data tokens, the rename record's file name and the
#: AGENTS.md R14 note (decision 5: agents count old `giap-*` jobs until none is left)
DOC_ALLOWED = _p(
    r"\.giap-complete|\.giap_setup\.lock|bios([/-])270\1khoa\1giap(?![A-Za-z0-9_])|"
    + re.escape(ADR64_FILE)
    + "|"
    + re.escape(R14_NOTE)
)
_ADR64_LINE = _p(r"^- \[ADR-0064:")


def is_doc(path: str) -> bool:
    """Documentation (Q-K32): the root docs, docs/ and every ``*.md`` / ``*.rst`` (data keep
    paths among them are reported as kept by ``doc_leftovers``)."""
    suffix = PurePosixPath(path).suffix.lower()
    return path in DOC_FILES or path.startswith("docs/") or suffix in {".md", ".rst"}


def doc_leftovers(repo: Path) -> tuple[list[tuple[str, int, str]], list[str]]:
    """``giap`` left in the documentation outside the allowed tokens, and the documentation
    files kept as data (``keep_category``: PROVENANCE notes, the rename record)."""
    hits: list[tuple[str, int, str]] = []
    kept: list[str] = []
    for path in tracked(repo):
        if not is_doc(path):
            continue
        data = (repo / path).read_bytes()
        if file_kind(path, data) != "text" or not ANY.search(text := data.decode("utf-8")):
            continue
        if keep_category(path):
            kept.append(path)
            continue
        for no, line in _line_hits(text):
            if path == "docs/adr/index.md" and _ADR64_LINE.match(line):
                continue
            if ANY.search(DOC_ALLOWED.sub("", line)):
                hits.append((path, no, line.strip()))
    return hits, kept


def check(repo: Path, verbose: bool) -> int:
    """Exit status 1 when anything that should be renamed is left (content or path)."""
    rep = survey(repo)
    bad = 0
    for path, (old, new) in sorted(rep.changes.items()):
        a, b = old.splitlines(), new.splitlines()
        diff = [(i + 1, x) for i, (x, y) in enumerate(zip(a, b, strict=False)) if x != y]
        if len(a) != len(b) and not diff:
            diff = [(min(len(a), len(b)) + 1, "(line count differs)")]
        for no, line in diff:
            what = line.strip() if ANY.search(line) else f"transition edit missing: {line.strip()}"
            print(f"LEFTOVER {path}:{no}: {what}")
            bad += 1
    for src, dst in rep.moves:
        print(f"LEFTOVER path {src} (should be {dst})")
        bad += 1
    for path, err in rep.broken:
        print(f"LEFTOVER transition {path}: {err}")
        bad += 1
    docs, docs_kept = doc_leftovers(repo)
    for path, no, line in docs:
        print(f"LEFTOVER doc {path}:{no}: {line}")
    bad += len(docs)
    for path in docs_kept:
        print(f"kept in docs ({keep_category(path)}): {path}")
    print(f"docs leftovers: {len(docs)}")
    print("## giap left after the rename, by keep category (occurrences)")
    for cat, n in sorted(rep.kept.items()):
        print(f"{cat:55s} {n:6d}")
    _print_kept(rep, verbose)
    print(f"uncategorised: {bad}")
    return 1 if bad else 0


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--plan", action="store_true")
    mode.add_argument("--apply", action="store_true")
    mode.add_argument("--check", action="store_true")
    ap.add_argument("--repo", type=Path, default=REPO)
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    repo = a.repo.resolve()
    if a.plan:
        print_plan(survey(repo), repo, a.verbose)
        return 0
    if a.apply:
        n = apply(repo)
        print(f"rename_gippa: {n} files rewritten; paths moved; now run --check")
        return 0
    return check(repo, a.verbose)


if __name__ == "__main__":
    raise SystemExit(main())
