# M9-01: The style module and the plot framework

Milestone M9 · **pass 1** · depends on: none · estimate ≈ 360 code lines (+ ≈ 400 test lines) ·
not statistical

## The gap this closes

Sixteen plots written independently would drift in colour, size, file layout and empty-state
handling. This issue fixes all of that once, before any plot is written: one style module, one
registry, one reader of a run's tables, one writer of the three image formats plus the data TSV.
Every later plot issue is "a pure function from frames to a figure and a data frame" on top of it.

## Spec sections implemented

- §10.1 (files `plots/Vnn_<name>.{svg,png,pdf}` + `plots/data/Vnn_<name>.tsv`; multi-instance
  suffix `V03_locus_<id>`; "One style module, `gippa.viz.style` … for the CLI and the report alike";
  widths 3.5 in / 7.2 in; Okabe–Ito categories, viridis/cividis continuous; the colour table; COG
  super-categories with letter labels; "Titles state the finding … generated from the tables";
  "Plots aggregate into bins when a track would draw more than 20,000 elements").
- §5.2 `report.formats` ([svg, png, pdf]), `report.dpi` (300).
- §9.2 (the data TSVs follow the table conventions), §6.5 (coordinates), §6.6 (number formats,
  no timestamps), §6.7 (events), §1 rule 6 (wording), conventions "file I/O only in io/ and steps/".

## Public interface

`src/gippa/viz/style.py` (pure constants and functions; no I/O)
```python
SINGLE_IN: float = 3.5
DOUBLE_IN: float = 7.2
OKABE_ITO: Mapping[str, str]       # black, orange, sky_blue, bluish_green, yellow, blue, vermillion, reddish_purple, grey (#999999)
DETECTION: Mapping[str, str]       # detected #0072B2, not_detected #BBBBBB, unestimated #E69F00 (§10.1)
STATUS_COLOUR: Mapping[str, str]   # detection_status enum -> colour (index decision 9)
REGION: Mapping[str, str]          # prophage #D55E00, plasmid #009E73, bgc #CC79A7, crispr #56B4E9, amr #F0E442
ISLAND: Mapping[str, tuple[str, str]]  # high ("#0072B2", "///"), low ("#D55E00", "\\\\\\")
COG_SUPER: Mapping[str, str]       # letter -> super-category key (index decision 8)
COG_SUPER_COLOUR: Mapping[str, str]
CONTINUOUS: str = "viridis"; LOCUS_ABUNDANCE: str = "cividis"
def apply() -> None                # rcParams: Agg, DejaVu Sans, sizes, svg.hashsalt, no metadata dates
def new_figure(width: Literal["single", "double"], height_in: float, **kw) -> Figure
def cog_colour(letter: str) -> str
```

`src/gippa/viz/base.py` (registry and result types; pure)
```python
@dataclass(frozen=True)
class Instance:                    # one output file set
    suffix: str | None             # None for single plots; "<id>" for V03_locus_<id>
    figure: Figure
    data: pd.DataFrame             # exactly what is drawn
    title: str

@dataclass(frozen=True)
class Skip:
    reason: str                    # "module_off:<m>" | "layer_missing:<l>" | "input_absent:<file>"

PlotFn = Callable[[RunView, GippaConfig], list[Instance] | Skip]

def plot(pid: str, name: str, *, needs: Sequence[str], optional: Sequence[str] = (),
         module: str = "core", layers: Sequence[str] = ()) -> Callable[[PlotFn], PlotFn]
def registered() -> Mapping[str, PlotSpec]       # "V01" -> spec, sorted by ID
def normalise_ids(text: str) -> list[str]        # "V1,v02" -> ["V01", "V02"]; ValueError on unknown
def empty_state(pid: str, name: str, message: str, columns: Sequence[str], width: ...) -> Instance
def binned(starts: np.ndarray, ends: np.ndarray, length: int, max_elements: int = 20_000) -> np.ndarray | None
def title_pct(n: int, d: int) -> str              # "12 %"; d == 0 -> "NA"
```

`src/gippa/io/run_view.py` (all file reads)
```python
class RunView:
    @classmethod
    def from_outdir(cls, out: Path) -> RunView
    @classmethod
    def from_inputs(cls, inputs: Mapping[str, Path]) -> RunView   # step keys (PUBLISHED / PUBLISHED_FILES)
    def table(self, key: str) -> pd.DataFrame | None    # None = file absent; header-only -> 0 rows
    def preflight(self) -> Mapping[str, Any] | None
    def support(self, name: str) -> pd.DataFrame | None # plots/support/<name>.tsv[.gz] (M9-03)
    def run_flags(self) -> frozenset[str]               # from events: data.flag and data.flags
    def modules_off(self) -> frozenset[str]             # from resolved_config.yaml when readable
```

`src/gippa/io/plot_files.py`
```python
def write_instance(inst: Instance, pid: str, name: str, plots_dir: Path,
                   formats: Sequence[str], dpi: int) -> list[Path]
def remove_plot_files(plots_dir: Path, pid: str) -> None   # before a rebuild (M9-02)
```

## Behaviour

1. **Style.** `apply()` sets `backend Agg`, `font.family = DejaVu Sans` (bundled with matplotlib),
   base 7 pt, titles 8 pt, tick labels 6 pt, `axes.prop_cycle` = the Okabe–Ito order, `svg.fonttype
   = none` (text stays text in the SVG, so the report can be searched), `svg.hashsalt =
   "gippa"`, `pdf.fonttype = 42`. It is called once by the `plots` step and by every test fixture.
2. **Colours.** Exactly the §10.1 table. COG letters map to four super-categories (index decision
   8); a gene with several letters (e.g. `KL`) is drawn in each letter, the convention of
   `gippa.context.summary._cog_letters` (M4-16). Letter descriptions come from the bundled
   `src/gippa/data/cog_categories.tsv` (26 letters; NCBI COG, public domain).
3. **RunView keys.** One table maps a key to its `OUT/` relative path, built from
   `gippa.runner.local.PUBLISHED` and `PUBLISHED_FILES` plus the entries that module does not list
   yet (`preflight/preflight.json`, `mapping/*.tsv`, `annotation/*.{tsv,bed}`, `plots/support/*`,
   `logs/events.jsonl`, `logs/steps/*.events.jsonl`, `resolved_config.yaml` or
   `work/resolved_config.yaml`). `from_outdir` and `from_inputs` therefore see the same keys, and
   a test pins them against `PUBLISHED_FILES` (the `test_local_publish_names.py` pattern).
4. **Reading tables.** TSVs are read with their Frictionless schema's types (`NA` → missing,
   `true`/`false` → bool, lists split on `;`, ADR-0006); a missing required column raises
   `ValueError` naming the table and column (a schema drift must fail loudly). Header-only files
   return zero-row frames with the schema's columns.
5. **Run flags.** `run_flags()` reads every events file and collects `data["flag"]` (string) and
   `data["flags"]` (list) from events, sorted and de-duplicated. JSON is parsed with Python's
   `json` (it accepts the `-Infinity` that `operons.missed_operons` writes; backlog M6-06).
6. **Skip logic** (index decision 4), in the registry wrapper, before the plot function runs: if
   the spec's `module` is in `modules_off()` → `Skip("module_off:<m>")`; if any of the spec's
   `layers` has `layer_missing:<layer>` in `run_flags()` → `Skip("layer_missing:<layer>")`; if a
   `needs` table is absent → `Skip("input_absent:<path>")`. Optional tables absent are passed as
   `None` and the plot omits that panel or track, noting it in the data TSV (`panel` column) and
   in the figure as a small grey "not available" text.
7. **Writing.** `write_instance` writes `Vnn_<name>[_<suffix>].<fmt>` for each configured format
   and `data/Vnn_<name>[_<suffix>].tsv` through `gippa.io.tsv` (§9.2: `%.6g`, `NA`, `true`/`false`,
   lists `;`-joined), sorted as the plot spec says. Suffixes are sanitised to `[A-Za-z0-9._-]`.
   PNG at `report.dpi`. Metadata: SVG `{"Date": None}`, PDF `{"CreationDate": None, "ModDate":
   None}`, PNG `{"Software": None}`, so two renders of the same data are byte-identical on one
   platform (index decision 6). Figures are closed after writing.
8. **Binning helper.** `binned()` returns `None` when `len(starts) ≤ max_elements`, otherwise the edges
   of `max_elements` equal-width bins over `[0, length)`;
   callers draw counts or means per bin and write the bins, not the elements, to the data TSV
   (index decision 12).
9. **Empty state.** `empty_state` draws an axes-free figure with the message centred and returns a
   header-only data frame with the plot's columns; the title is the message ("No islands at
   FWER ≤ 0.05 in this run").

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `OUT/` has no `tables/` and no `preflight/` | `RunView.from_outdir` raises `NotAGippaRun` (the CLI maps it to exit 1, M9-02) |
| a table header-only | zero-row frame; plots draw their empty state |
| a table missing a schema column | `ValueError` naming table and column |
| events files absent (a hand-built `OUT/`) | `run_flags()` = ∅ |
| `resolved_config.yaml` absent | `modules_off()` = ∅; absence of a module's tables still skips via `input_absent` |
| a COG letter not in the bundled table | drawn in the "poorly characterised" grey, data TSV keeps the letter |
| a suffix with `/` or spaces | sanitised to `_` |
| `report.formats = [svg]` | only SVG + data TSV written |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `style.STATUS_COLOUR["single_peptide_unestimated"]` | `"#E69F00"` |
| 2 | `style.cog_colour("J")`, `("C")`, `("S")`, `("D")` | `#0072B2`, `#E69F00`, `#999999`, `#009E73` |
| 3 | `normalise_ids("V1, v18,V03")` with V01, V03, V18 registered | `["V01", "V03", "V18"]` |
| 4 | `normalise_ids("V04")` (V4 is deferred, §10) | `ValueError: unknown plot ID V04` |
| 5 | `binned(starts of 4,500 genes, ends, 5_000_000)` | `None` (under 20,000) |
| 6 | `binned(25,000 intervals, …, 1_000_000)` | 20,000 equal bins of 50 bp each |
| 7 | an events file with `{"event": "context_layers.layer_missing", "data": {"flag": "layer_missing:genomad"}}` and `{"event": "function.run_flags", "data": {"flags": ["layer_missing:kegg"]}}` | `run_flags() == {"layer_missing:genomad", "layer_missing:kegg"}` |
| 8 | a plot spec `layers=["genomad"]` and example 7's flags | `Skip("layer_missing:genomad")`; the plot function is not called |
| 9 | the same `Instance` written twice into two directories | identical bytes for `.svg`, `.pdf`, `.png` and `.tsv` |
| 10 | `RunView.from_outdir(local/runs/val-mvp-2026-10-22/OUT)`, `table("islands")` | zero rows, the 24 schema columns |

## Acceptance criteria

- [ ] `tests/unit/viz/test_style.py`: examples 1–2; every §10.1 colour present with its exact hex;
  widths 3.5 and 7.2; `apply()` sets Agg and DejaVu Sans.
- [ ] `tests/unit/viz/test_registry.py`: examples 3, 4, 8; IDs unique; skip order (module, layer,
  input).
- [ ] `tests/unit/io/test_run_view.py`: examples 7, 10; keys pinned against `PUBLISHED_FILES`;
  header-only and missing-column cases; `from_inputs` and `from_outdir` give equal frames on one
  tree.
- [ ] `tests/unit/io/test_plot_files.py`: example 9 (byte identity), names and suffix sanitising,
  formats honoured, data TSV encoding (ADR-0006).
- [ ] `tests/unit/viz/test_binning.py`: examples 5, 6; property test: bins ≤ 20,000 and cover
  `[0, length)`.
- [ ] `src/gippa/data/cog_categories.tsv` with provenance line; `mypy --strict` clean.

## Test ideas

- Hypothesis: any list of `(start, end)` on a 10 Mb contig → `binned` bins cover every element.
- A fake plot registered in a test module proves the registry wrapper never calls a skipped plot.

## Hardening owed (pass 2)

- M9-24 turns the style rules into checks over every registered plot.
- M9-25 freezes the rcParams; after Khoa's approval a change to `style.py` is a reviewed change.

## Questions for the lead

1. File I/O placement: the convention says I/O only in `io/` and `steps/`. **Proposed default:**
   `RunView` and the writers live in `gippa/io/`; `gippa/viz/` holds pure figure code.
2. A plot that raises: fail the step (exit 4) or skip it with `reason = error`? **Proposed
   default:** fail (no silent fallback, the M4-17 rule); declared skips cover the expected cases.
3. Font: DejaVu Sans is the only font guaranteed on every CI leg and on FarmShare without a system
   package. **Proposed default:** DejaVu Sans; Khoa may ask for another at the M9-25 review.
4. COG super-category colours (index decision 8) reuse hues the detection and region tables use.
   §10.1 has only eight Okabe–Ito colours and assigns seven. **Proposed default:** reuse is
   acceptable because COG colours never appear in the same plot as detection colours (V8 and V14
   only); M9-24 checks that no plot mixes the two palettes.
5. The 20,000-element rule (index decision 12): per track, as §10.1 words it. **Proposed default:**
   per track.
