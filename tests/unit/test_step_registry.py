"""Step registry, topological order and ``run_step`` (M2-02; GIAP_SPEC.md §6.2, §6.6)."""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.config import GiapConfig, resolve_config
from giap.rng import for_step


@contextmanager
def fresh_registry() -> Iterator[ModuleType]:
    """An empty registry for use inside Hypothesis examples (same as the ``registry`` fixture)."""
    from giap.steps import registry as reg

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(reg, "_REGISTRY", {})
        yield reg


def _noop(ctx: Any) -> None:
    return None


def _copy(ctx: Any) -> None:
    ctx.outputs["b"].write_bytes(ctx.inputs["a"].read_bytes())


# ---------------------------------------------------------------- decorator (example 1)


def test_step_decorator_registers_a_py_spec(registry: ModuleType) -> None:
    def toy(ctx: Any) -> None:
        return None

    returned = registry.step("toy", inputs=["a"], outputs=["b"])(toy)
    assert returned is toy
    spec = registry.get("toy")
    assert spec == registry.StepSpec(
        name="toy",
        kind="py",
        inputs=("a",),
        outputs=("b",),
        optional_inputs=(),
        module="core",
        precomputed=None,
        func=toy,
    )


def test_step_decorator_keeps_module_optional_inputs_and_precomputed(registry: ModuleType) -> None:
    registry.step(
        "toy",
        inputs=["a"],
        outputs=["b", "c"],
        optional_inputs=["o"],
        module="bgc",
        precomputed="regions",
    )(_noop)
    spec = registry.get("toy")
    assert spec.kind == "py"
    assert spec.outputs == ("b", "c")
    assert spec.optional_inputs == ("o",)
    assert spec.module == "bgc"
    assert spec.precomputed == "regions"


def test_external_registers_an_ext_spec_without_func(registry: ModuleType) -> None:
    assert (
        registry.external(
            "tool", inputs=["genome.norm.fa"], outputs=["tool.tsv"], precomputed="terminators"
        )
        is None
    )
    spec = registry.get("tool")
    assert spec.kind == "ext"
    assert spec.func is None
    assert spec.inputs == ("genome.norm.fa",)
    assert spec.outputs == ("tool.tsv",)
    assert spec.precomputed == "terminators"
    assert spec.module == "core"


def test_step_spec_is_frozen(registry: ModuleType) -> None:
    registry.step("toy", inputs=["a"], outputs=["b"])(_noop)
    with pytest.raises(AttributeError):
        registry.get("toy").name = "other"


def test_all_steps_sorted_by_name(registry: ModuleType) -> None:
    for name in ["zeta", "alpha", "mid"]:
        registry.step(name, inputs=[], outputs=[f"{name}.out"])(_noop)
    registry.external("beta", inputs=[], outputs=["beta.out"])
    assert [s.name for s in registry.all_steps()] == ["alpha", "beta", "mid", "zeta"]


# ---------------------------------------------------------------- duplicates (example 2)


def test_duplicate_step_name_rejected(registry: ModuleType) -> None:
    registry.step("toy", inputs=["a"], outputs=["b"])(_noop)
    with pytest.raises(ValueError, match=r"step 'toy' already registered"):
        registry.step("toy", inputs=["x"], outputs=["y"])(_noop)


def test_duplicate_between_py_and_ext_rejected(registry: ModuleType) -> None:
    registry.external("toy", inputs=[], outputs=["b"])
    with pytest.raises(ValueError, match=r"step 'toy' already registered"):
        registry.step("toy", inputs=[], outputs=["c"])(_noop)


# ---------------------------------------------------------------- topological order (3, 4)


def _register_xyz(reg: ModuleType, order: list[str]) -> None:
    edges = {"x": (["a"], ["b"]), "y": (["b"], ["c"]), "z": (["a"], ["d"])}
    for name in order:
        ins, outs = edges[name]
        reg.step(name, inputs=ins, outputs=outs)(_noop)


def test_topo_order_example(registry: ModuleType) -> None:
    _register_xyz(registry, ["x", "y", "z"])
    assert registry.topo_order(["y", "z", "x"]) == ["x", "y", "z"]


def test_topo_order_producer_before_consumer_beats_name_order(registry: ModuleType) -> None:
    registry.step("a_consumer", inputs=["k"], outputs=["k2"])(_noop)
    registry.step("z_producer", inputs=[], outputs=["k"])(_noop)
    assert registry.topo_order(["a_consumer", "z_producer"]) == ["z_producer", "a_consumer"]


def test_topo_order_optional_inputs_are_edges(registry: ModuleType) -> None:
    registry.step("a_use", inputs=[], outputs=["u"], optional_inputs=["opt"])(_noop)
    registry.external("z_make", inputs=[], outputs=["opt"])
    assert registry.topo_order(["a_use", "z_make"]) == ["z_make", "a_use"]


def test_topo_order_ignores_inputs_without_producer(registry: ModuleType) -> None:
    registry.step("only", inputs=["genome.norm.fa"], outputs=["x"])(_noop)
    assert registry.topo_order(["only"]) == ["only"]


@settings(max_examples=30, deadline=None)
@given(st.permutations(["x", "y", "z"]), st.permutations(["x", "y", "z"]))
def test_topo_order_independent_of_registration_and_argument_order(
    reg_order: list[str], arg_order: list[str]
) -> None:
    with fresh_registry() as reg:
        _register_xyz(reg, list(reg_order))
        assert reg.topo_order(list(arg_order)) == ["x", "y", "z"]


def test_topo_order_cycle_rejected_naming_both_steps(registry: ModuleType) -> None:
    registry.step("x", inputs=["a"], outputs=["b"])(_noop)
    registry.step("y", inputs=["b"], outputs=["a"])(_noop)
    with pytest.raises(ValueError) as info:
        registry.topo_order(["x", "y"])
    msg = str(info.value)
    assert "x" in msg and "y" in msg
    assert "cycle" in msg.lower()


# ---------------------------------------------------------------- annotation contract table


# (name, kind, inputs, optional_inputs, outputs, module, precomputed) from M2-index.md
CONTRACT: list[tuple[str, str, list[str], list[str], list[str], str, str | None]] = [
    ("bakta", "ext", ["genome.norm.fa", "replicons.tsv"], ["regions.gff3"],
     ["bakta.json", "bakta.gff3", "bakta.gbff", "bakta.faa", "bakta.tsv"], "core", "bakta_json"),
    ("parse_bakta", "py", ["bakta.json", "genome.norm.fa"], [],
     ["genes.parquet", "features.parquet", "replicons.parquet", "rrna16s.fa", "trna_counts.tsv",
      "xrefs.parquet", "annotation_meta.json"], "core", None),
    ("uniop_input", "py", ["genes.parquet"], [], ["uniop_cds.gff", "uniop_map.tsv"], "core",
     "operon_pairs"),
    ("uniop", "ext", ["uniop_cds.gff"], [], ["uniop.pred"], "core", "operon_pairs"),
    ("operon_pairs", "py", ["genes.parquet", "uniop_map.tsv"], ["uniop.pred", "uniop_cds.gff"],
     ["pairs_raw.parquet", "tus.parquet", "tu_upstream/", "operon_meta.json"], "core", None),
    ("promoter_windows", "py", ["genes.parquet", "genome.norm.fa"], [],
     ["windows.fa", "windows_map.parquet"], "core", "promoters"),
    ("promotech", "ext", ["windows.fa"], [], ["promotech/"], "core", "promoters"),
    ("promoter_layer", "py", ["windows_map.parquet", "promotech/"], [],
     ["promoter_windows.parquet"], "core", "promoters"),
    ("promoter_calculator", "ext", ["tu_upstream/"], [], ["pc/"], "core", "tx_rates"),
    ("tx_layer", "py", ["pc/", "tus.parquet"], [], ["tx_rates.parquet"], "core", "tx_rates"),
    ("terminatornet", "ext", ["genome.norm.fa"], [], ["terminatornet.tsv"], "core", "terminators"),
    ("transterm_input", "py", ["genes.parquet"], [], ["genes.coords"], "core", "terminators"),
    ("transtermhp", "ext", ["genome.norm.fa", "genes.coords"], [], ["transtermhp.tt"], "core",
     "terminators"),
    ("bactermfinder", "ext", ["genome.norm.fa"], [], ["bactermfinder/"], "opt_in", "terminators"),
    ("terminator_layer", "py", ["terminatornet.tsv", "transtermhp.tt"], ["bactermfinder/"],
     ["terminators.parquet", "terminators_crosscheck.tsv"], "core", "terminators"),
    ("ostir_input", "py", ["genes.parquet", "annotation_meta.json"], [], ["ostir_input.csv"],
     "core", "rbs"),
    ("ostir", "ext", ["ostir_input.csv"], [], ["ostir.csv"], "core", "rbs"),
    ("rbs_layer", "py", ["ostir.csv", "ostir_input.csv"], [], ["rbs.parquet"], "core", "rbs"),
    ("codon_usage", "py", ["genes.parquet", "trna_counts.tsv"], [], ["codon.parquet"], "core",
     None),
    ("genomad", "ext", ["genome.norm.fa"], [], ["genomad/"], "mobile", "regions"),
    ("antismash", "ext", ["bakta.gbff"], [], ["antismash/"], "bgc", "regions"),
    ("region_layer", "py", ["genes.parquet", "features.parquet", "replicons.parquet"],
     ["genomad/", "antismash/"], ["regions.parquet", "region_genes.parquet"], "core", "regions"),
    ("context_layers", "py",
     ["genes.parquet", "tus.parquet", "promoter_windows.parquet", "tx_rates.parquet",
      "terminators.parquet", "rbs.parquet", "codon.parquet", "regions.parquet"], [],
     ["layers.parquet"], "core", None),
]  # fmt: skip


def _register_contract(reg: ModuleType, rows: list[Any] = CONTRACT) -> None:
    for name, kind, ins, opt, outs, module, pre in rows:
        kwargs = dict(inputs=ins, outputs=outs, optional_inputs=opt, module=module, precomputed=pre)
        if kind == "py":
            reg.step(name, **kwargs)(_noop)
        else:
            reg.external(name, **kwargs)


def test_every_contract_step_registers_without_conflict(registry: ModuleType) -> None:
    _register_contract(registry)
    names = [row[0] for row in CONTRACT]
    assert [s.name for s in registry.all_steps()] == sorted(names)
    order = registry.topo_order(names)
    assert sorted(order) == sorted(names)
    pos = {n: i for i, n in enumerate(order)}
    for before, after in [
        ("bakta", "parse_bakta"),
        ("parse_bakta", "uniop_input"),
        ("uniop", "operon_pairs"),
        ("operon_pairs", "promoter_calculator"),
        ("promotech", "promoter_layer"),
        ("transterm_input", "transtermhp"),
        ("bactermfinder", "terminator_layer"),
        ("ostir", "rbs_layer"),
        ("genomad", "region_layer"),
        ("antismash", "region_layer"),
        ("region_layer", "context_layers"),
        ("codon_usage", "context_layers"),
    ]:
        assert pos[before] < pos[after], (before, after)


def test_contract_topo_order_is_deterministic(registry: ModuleType) -> None:
    _register_contract(registry)
    names = [row[0] for row in CONTRACT]
    first = registry.topo_order(names)
    assert registry.topo_order(list(reversed(names))) == first
    with fresh_registry() as other:
        _register_contract(other, list(reversed(CONTRACT)))
        assert other.topo_order(sorted(names)) == first


def _register_bakta_like(reg: ModuleType, with_producer: bool) -> None:
    if with_producer:
        reg.step("make_regions", inputs=[], outputs=["regions.gff3"])(_noop)
    reg.external(
        "bakta",
        inputs=["genome.norm.fa"],
        outputs=["bakta.json"],
        optional_inputs=["regions.gff3"],
        precomputed="bakta_json",
    )
    reg.step("parse", inputs=["bakta.json"], outputs=["genes.parquet"])(_noop)


@pytest.mark.parametrize("with_producer", [True, False])
def test_external_optional_input_plans_with_and_without_producer(
    with_producer: bool, tmp_path: Path
) -> None:
    from giap.runner.local import LayerMissing, plan_local

    with fresh_registry() as reg:
        _register_bakta_like(reg, with_producer)
        spec = reg.get("bakta")
        assert spec.kind == "ext" and spec.func is None
        assert spec.inputs == ("genome.norm.fa",)
        assert spec.optional_inputs == ("regions.gff3",)
        names = ["parse", "bakta", *(["make_regions"] if with_producer else [])]
        expected = [*(["make_regions"] if with_producer else []), "bakta", "parse"]
        assert reg.topo_order(names) == expected
        with pytest.raises(LayerMissing):
            plan_local(GiapConfig(), ["parse"])
        layer = tmp_path / "bakta.json"
        layer.write_text("{}")
        cfg = resolve_config(None, [f"annotation.precomputed.bakta_json={layer}"])
        assert plan_local(cfg, ["parse"]) == ["parse"]


# ---------------------------------------------------------------- run_step


@pytest.fixture
def cfg() -> GiapConfig:
    return GiapConfig(seed=7)


def test_run_step_copies_and_creates_output_dirs(
    registry: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    registry.step("toy", inputs=["a"], outputs=["b"])(_copy)
    src = tmp_path / "in.txt"
    src.write_text("hello\n")
    out = tmp_path / "deep" / "er" / "out.txt"
    registry.run_step("toy", cfg, {"a": src}, {"b": out})
    assert out.read_text() == "hello\n"


def test_run_step_context_fields(registry: ModuleType, cfg: GiapConfig, tmp_path: Path) -> None:
    seen: dict[str, Any] = {}

    def toy(ctx: Any) -> None:
        seen["ctx"] = ctx
        seen["draw"] = ctx.rng.random()
        ctx.outputs["b"].write_text("x")

    registry.step("toy", inputs=["a"], outputs=["b"])(toy)
    src = tmp_path / "in.txt"
    src.write_text("1")
    registry.run_step("toy", cfg, {"a": src}, {"b": tmp_path / "out.txt"})
    ctx = seen["ctx"]
    assert isinstance(ctx, registry.StepContext)
    assert ctx.cfg == cfg
    assert dict(ctx.inputs) == {"a": src}
    # outputs are temporary paths renamed into place on success (review decision 3)
    assert set(ctx.outputs) == {"b"}
    assert (tmp_path / "out.txt").read_text() == "x"
    assert seen["draw"] == for_step("toy", cfg.seed).random()
    assert hasattr(ctx.log, "info")


def test_run_step_optional_input_absent(
    registry: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    seen: dict[str, Any] = {}

    def toy(ctx: Any) -> None:
        seen["keys"] = set(ctx.inputs)
        ctx.outputs["b"].write_text("x")

    registry.step("toy", inputs=["a"], outputs=["b"], optional_inputs=["opt"])(toy)
    src = tmp_path / "in.txt"
    src.write_text("1")
    registry.run_step("toy", cfg, {"a": src}, {"b": tmp_path / "out.txt"})
    assert seen["keys"] == {"a"}


def test_run_step_directory_input_allowed(
    registry: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    def toy(ctx: Any) -> None:
        names = sorted(p.name for p in ctx.inputs["genomad/"].iterdir())
        ctx.outputs["b"].write_text(",".join(names))

    registry.step("toy", inputs=["genomad/"], outputs=["b"])(toy)
    d = tmp_path / "genomad"
    d.mkdir()
    (d / "x.tsv").write_text("1")
    registry.run_step("toy", cfg, {"genomad/": d}, {"b": tmp_path / "out.txt"})
    assert (tmp_path / "out.txt").read_text() == "x.tsv"


def test_run_step_unknown_name_is_wiring_error(registry: ModuleType, cfg: GiapConfig) -> None:
    with pytest.raises(registry.WiringError, match="nope"):
        registry.run_step("nope", cfg, {}, {})


def test_run_step_undeclared_input_key_is_wiring_error(
    registry: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    registry.step("toy", inputs=["a"], outputs=["b"])(_copy)
    src = tmp_path / "in.txt"
    src.write_text("1")
    with pytest.raises(registry.WiringError) as info:
        registry.run_step("toy", cfg, {"a": src, "zz": src}, {"b": tmp_path / "o"})
    assert "toy" in str(info.value) and "zz" in str(info.value)


def test_run_step_undeclared_output_key_is_wiring_error(
    registry: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    registry.step("toy", inputs=["a"], outputs=["b"])(_copy)
    src = tmp_path / "in.txt"
    src.write_text("1")
    with pytest.raises(registry.WiringError) as info:
        registry.run_step("toy", cfg, {"a": src}, {"b": tmp_path / "o", "qq": tmp_path / "q"})
    assert "toy" in str(info.value) and "qq" in str(info.value)


def test_run_step_missing_required_input_is_wiring_error(
    registry: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    registry.step("toy", inputs=["a"], outputs=["b"])(_copy)
    with pytest.raises(registry.WiringError) as info:
        registry.run_step("toy", cfg, {}, {"b": tmp_path / "o"})
    assert "toy" in str(info.value) and "a" in str(info.value)


def test_run_step_nonexistent_input_path_is_wiring_error(
    registry: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    registry.step("toy", inputs=["a"], outputs=["b"])(_copy)
    with pytest.raises(registry.WiringError, match="toy"):
        registry.run_step("toy", cfg, {"a": tmp_path / "missing.txt"}, {"b": tmp_path / "o"})


def test_run_step_wraps_exception_in_step_error(
    registry: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    def toy(ctx: Any) -> None:
        raise RuntimeError("boom")

    registry.step("toy", inputs=[], outputs=["b"])(toy)
    with pytest.raises(registry.StepError) as info:
        registry.run_step("toy", cfg, {}, {"b": tmp_path / "o"})
    assert "toy" in str(info.value)
    assert "boom" in str(info.value)
    assert isinstance(info.value.__cause__, RuntimeError) or any(
        isinstance(a, RuntimeError) for a in info.value.args
    )


def test_run_step_missing_declared_output_is_wiring_error(
    registry: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    def toy(ctx: Any) -> None:
        ctx.outputs["b"].write_text("x")  # forgets "c"

    registry.step("toy", inputs=[], outputs=["b", "c"])(toy)
    with pytest.raises(registry.WiringError) as info:
        registry.run_step("toy", cfg, {}, {"b": tmp_path / "b", "c": tmp_path / "c"})
    assert "toy" in str(info.value) and "c" in str(info.value)


def test_run_step_ext_step_has_no_function(
    registry: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    registry.external("tool", inputs=[], outputs=["b"])
    with pytest.raises(registry.WiringError, match="tool"):
        registry.run_step("tool", cfg, {}, {"b": tmp_path / "b"})


# ---------------------------------------------------------------- determinism (§6.6)


def test_run_step_outputs_byte_identical(
    registry: ModuleType, cfg: GiapConfig, tmp_path: Path
) -> None:
    def toy(ctx: Any) -> None:
        data = ctx.inputs["a"].read_bytes()
        draws = ctx.rng.integers(0, 1_000_000, size=20)
        ctx.outputs["b"].write_bytes(data + ",".join(map(str, draws)).encode())

    registry.step("toy", inputs=["a"], outputs=["b"])(toy)
    src = tmp_path / "in.txt"
    src.write_text("seeded\n")
    registry.run_step("toy", cfg, {"a": src}, {"b": tmp_path / "r1" / "b"})
    registry.run_step("toy", cfg, {"a": src}, {"b": tmp_path / "r2" / "b"})
    assert (tmp_path / "r1" / "b").read_bytes() == (tmp_path / "r2" / "b").read_bytes()
    registry.run_step("toy", GiapConfig(seed=8), {"a": src}, {"b": tmp_path / "r3" / "b"})
    assert (tmp_path / "r3" / "b").read_bytes() != (tmp_path / "r1" / "b").read_bytes()
