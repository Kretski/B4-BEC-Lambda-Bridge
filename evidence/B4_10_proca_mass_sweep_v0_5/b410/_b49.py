"""Import shim: always use the vendored, hash-verified B4.9 package.

B4.9 is not modified.  The vendored copy is put first on sys.path so that a
different locally installed ``proca_dispersion`` cannot be picked up silently.
"""
from __future__ import annotations
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
B49_DIR = ROOT / "vendor" / "B4_9_proca_eigenvector_tracking_v0_1"
MANIFEST = ROOT / "vendor" / "SHA256SUMS_B4_9.txt"

if str(B49_DIR) not in sys.path:
    sys.path.insert(0, str(B49_DIR))

import proca_dispersion  # noqa: E402
from proca_dispersion.model import (  # noqa: E402,F401
    Background, local_roots, full_roots, full_polynomial, local_polynomial)
from proca_dispersion.eigenmodes import (  # noqa: E402,F401
    coupled_matrix, null_vector, mode_data, track_eigenmodes)

if Path(proca_dispersion.__file__).resolve().parent.parent != B49_DIR.resolve():
    raise ImportError("proca_dispersion was not loaded from the vendored B4.9 copy: "
                      f"{proca_dispersion.__file__}")


def verify_manifest():
    """Return list of (relative_path, ok) for every file in the B4.9 manifest."""
    out = []
    for line in MANIFEST.read_text().splitlines():
        if not line.strip():
            continue
        digest, rel = line.split(None, 1)
        rel = rel.strip()
        path = ROOT / "vendor" / rel
        ok = path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() == digest
        out.append((rel, ok))
    return out
