from pathlib import Path
from b410._b49 import B49_DIR, verify_manifest   # shim puts the vendored B4.9 first on sys.path
import proca_dispersion


def test_vendored_b49_is_byte_identical_to_handoff_manifest():
    res = verify_manifest()
    assert len(res) == 16
    bad = [rel for rel, ok in res if not ok]
    assert bad == []


def test_b49_is_imported_from_vendored_copy():
    assert Path(proca_dispersion.__file__).resolve().parent.parent == B49_DIR.resolve()
