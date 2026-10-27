"""The image-cache helper (M11-14 item 3, Q-L11; shared by preflight 14f1, fgsea 14c, setup 14f2).

``gippa.io.image_key.image_cache(environ, db_dir) -> Path`` is ``NXF_APPTAINER_CACHEDIR``, else
``$GIPPA_ROOT/containers``, else ``<db_dir>/containers``. Reading chosen here: an empty value
counts as unset (as ``os.environ.get(...) or ...`` does in preflight today); the function is pure:
it reads only ``environ`` (never ``os.environ``) and creates nothing.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st


def _image_cache():  # type: ignore[no-untyped-def]
    from gippa.io.image_key import image_cache

    return image_cache


DB = Path("/fake/db")


def test_nxf_apptainer_cachedir_wins() -> None:
    env = {"NXF_APPTAINER_CACHEDIR": "/fake/nxf", "GIPPA_ROOT": "/fake/root"}
    assert _image_cache()(env, DB) == Path("/fake/nxf")


def test_gippa_root_containers_second() -> None:
    assert _image_cache()({"GIPPA_ROOT": "/fake/root"}, DB) == Path("/fake/root/containers")


def test_db_dir_containers_last() -> None:
    assert _image_cache()({}, DB) == Path("/fake/db/containers")


def test_empty_values_count_as_unset() -> None:
    env = {"NXF_APPTAINER_CACHEDIR": "", "GIPPA_ROOT": ""}
    assert _image_cache()(env, DB) == Path("/fake/db/containers")


def test_returns_a_path() -> None:
    assert isinstance(_image_cache()({}, DB), Path)


def test_reads_only_the_given_environ(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("NXF_APPTAINER_CACHEDIR", "/fake/from-os")
    monkeypatch.setenv("GIPPA_ROOT", "/fake/os-root")
    assert _image_cache()({}, DB) == Path("/fake/db/containers")


def test_creates_nothing(tmp_path: Path) -> None:
    db = tmp_path / "db"
    got = _image_cache()({"GIPPA_ROOT": str(tmp_path / "root")}, db)
    assert got == tmp_path / "root" / "containers"
    assert list(tmp_path.iterdir()) == []


_val = st.one_of(st.none(), st.just(""), st.sampled_from(["/a/nxf", "/b/root", "rel/dir"]))


@given(nxf=_val, root=_val)
def test_property_precedence(nxf: str | None, root: str | None) -> None:
    env = {}
    if nxf is not None:
        env["NXF_APPTAINER_CACHEDIR"] = nxf
    if root is not None:
        env["GIPPA_ROOT"] = root
    got = _image_cache()(env, DB)
    if nxf:
        assert got == Path(nxf)
    elif root:
        assert got == Path(root) / "containers"
    else:
        assert got == DB / "containers"
