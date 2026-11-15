"""B4.12 common-M check (specification: B4_12_COMMON_M_TEST_PROPOSAL_v0.1.2.md, commit 3cba04e).

Usage:
    python -I b412_common_m_check.py <evidence_dir> <outdir>

<evidence_dir> is evidence/B4_10_proca_mass_sweep_v0_5 (commit 7ae9c5c, PROVENANCE fixed in 511f6fe).

Standard library only. Does not import the b410 package, numpy or scipy.

For each of the 41 historical B4.12 backgrounds and each set A, B, C:
  * Original criterion (historical, reproduced):  sqrtX_max <= g_max * rho_max            (b410/uv.py, controlled)
  * Common-M criterion (new):                     c_max * G <= g_max * rho_max
        G     = max(|q|, |w|, |r|, |v|)
        c_max = max(sqrt(p), |kappa|/sqrt(p))       (code coefficients: g_phi = M sqrt(p), g_psi = -M kappa/sqrt(p))
    classified PASS / FAIL / BORDERLINE with margin m = g_max*rho_max - c_max*G and tolerance 1e-12*g_max*rho_max.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

SPEC = "B4_12_COMMON_M_TEST_PROPOSAL_v0.1.2.md"

# Thresholds copied from b410/uv.py, line 17 (SETS). Same names, same order, same values.
SETS = {
    "A_g1_rho0.3": (1.0, 0.3),
    "B_g1_rho0.1": (1.0, 0.1),
    "C_gStrong_rho0.3": (math.sqrt(4 * math.pi), 0.3),
}
BORDER_REL_TOL = 1e-12   # spec §3
K3_TOL = 1e-12           # spec §5, K3
K3_M = 7.0               # spec §5, K3 (as run_b412.py, line 20)
K4_REPORT_THRESHOLD = 1e-9  # spec §5, K4 (report only)
# Expected counts for K1 (spec §2 and §5 K1). Checked, not only reported.
K1_EXPECTED = dict(manifest_selected=40, total=41, negative_controls=12, noncontrol=29)

# Files read. INPUTS feed the criteria; REFERENCE is read only for K1 (historical counts).
INPUTS = ("data/B3_max_speed_7_315096.json", "data/e3_manifest.json")
REFERENCE = "results_b412/b412_summary.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_sums(evidence: Path) -> dict:
    sums = {}
    for line in (evidence / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
        h, rel = line.split("  ", 1)
        sums[rel] = h
    return sums


def load_cases(evidence: Path) -> tuple:
    """The 41 backgrounds, in the order of run_b412.py lines 15-17: anchor first, then manifest 'selected'.
    Also returns the number of manifest 'selected' records and the manifest's own 'n_selected' (K1)."""
    b3 = json.loads((evidence / INPUTS[0]).read_text(encoding="utf-8"))
    man = json.loads((evidence / INPUTS[1]).read_text(encoding="utf-8"))
    keys = ("p", "s", "kappa", "w", "q", "v", "r")
    cases = [dict(tag="ANCHOR_B3", cls="anchor", mode=b3["mode"], **{k: float(b3[k]) for k in keys})]
    for row in man["selected"]:
        cases.append(dict(tag=f"{row['mode']}_{row['cls']}_{row['index']}", cls=row["cls"], mode=row["mode"],
                          **{k: float(row[k]) for k in keys}))
    return cases, len(man["selected"]), man.get("n_selected")


def original_sqrtX_max(c: dict) -> float:
    # identical expression to b410/uv.py, controlled(), line 31 (X uses p for phi and s for psi)
    return max(math.sqrt(c["p"] * c["q"] ** 2), math.sqrt(c["s"] * c["r"] ** 2),
               math.sqrt(c["p"] * c["w"] ** 2), math.sqrt(c["s"] * c["v"] ** 2))


def common_m_quantities(c: dict) -> tuple:
    G = max(abs(c["q"]), abs(c["w"]), abs(c["r"]), abs(c["v"]))
    coef_phi = math.sqrt(c["p"])                     # g_phi / M
    coef_psi = abs(c["kappa"]) / math.sqrt(c["p"])   # |g_psi| / M, code coefficient (always used, spec §5 K4)
    return G, max(coef_phi, coef_psi), coef_phi, coef_psi


def classify(margin: float, scale: float) -> str:
    tol = BORDER_REL_TOL * scale
    if margin > tol:
        return "PASS"
    if margin < -tol:
        return "FAIL"
    return "BORDERLINE"


def k3_identity_dev(c: dict) -> float:
    """|g_j rho_j - sqrt(X_j)| at M = 7, as run_b412.py (uv.invariants)."""
    M = K3_M
    g_phi = M * math.sqrt(c["p"])
    g_psi = abs(-M * c["kappa"] / math.sqrt(c["p"]))
    prods = dict(Xq=g_phi * abs(c["q"]) / M, Xw=g_phi * abs(c["w"]) / M,
                 Xr=g_psi * abs(c["r"]) / M, Xv=g_psi * abs(c["v"]) / M)
    sx = dict(Xq=math.sqrt(c["p"] * c["q"] ** 2), Xw=math.sqrt(c["p"] * c["w"] ** 2),
              Xr=math.sqrt(c["s"] * c["r"] ** 2), Xv=math.sqrt(c["s"] * c["v"] ** 2))
    return max(abs(prods[k] - sx[k]) for k in sx)


def main(evidence_dir: str, outdir: str) -> int:
    evidence = Path(evidence_dir)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    checks = {}

    # K5 — input (and K1 reference) hashes against the package manifest; stop on mismatch
    sums = read_sums(evidence)
    k5 = {rel: dict(expected=sums.get(rel), actual=sha256(evidence / rel)) for rel in INPUTS + (REFERENCE,)}
    k5_ok = all(v["expected"] is not None and v["expected"] == v["actual"] for v in k5.values())
    checks["K5_inputs"] = dict(pass_=k5_ok, files=k5)
    if not k5_ok:
        print("K5 FAIL: input hash mismatch; stopping")
        (out / "b412_common_m.json").write_text(json.dumps(dict(spec=SPEC, verdict="K5_FAIL", checks=checks), indent=1) + "\n",
                                                encoding="utf-8", newline="\n")
        return 2

    cases, n_manifest, n_manifest_declared = load_cases(evidence)
    ref = json.loads((evidence / REFERENCE).read_text(encoding="utf-8"))

    rows = []
    for c in cases:
        sxm = original_sqrtX_max(c)
        G, cmax, cphi, cpsi = common_m_quantities(c)
        row = dict(tag=c["tag"], cls=c["cls"], mode=c["mode"], p=c["p"], s=c["s"], kappa=c["kappa"],
                   sqrtX_max=sxm, G=G, coef_phi=cphi, coef_psi=cpsi, c_max=cmax,
                   delta_rank1=abs(c["kappa"] ** 2 / (c["p"] * c["s"]) - 1.0),
                   k3_identity_dev=k3_identity_dev(c))
        for name, (gmax, rmax) in SETS.items():
            scale = gmax * rmax
            row[f"original_{name}"] = bool(sxm <= scale)                 # Original criterion (independent of Common-M)
            m = scale - cmax * G                                         # Common-M criterion
            row[f"commonM_{name}"] = classify(m, scale)
            row[f"margin_{name}"] = m
            row[f"M_low_{name}"] = G / rmax
            row[f"M_high_{name}"] = gmax / cmax
        rows.append(row)

    n = len(rows)
    negctrl = [r for r in rows if r["cls"] == "negative_control"]
    nonctrl = [r for r in rows if r["cls"] != "negative_control"]

    # K1 — reproduction of the historical Original counts, with explicit denominators
    # Denominators are checked against the expected counts of the specification, not only reported.
    k1_counts = dict(manifest_selected=n_manifest, manifest_n_selected_field=n_manifest_declared, total=n,
                     negative_controls=len(negctrl), noncontrol=len(nonctrl),
                     anchor_in_noncontrol=any(r["tag"] == "ANCHOR_B3" for r in nonctrl))
    k1_counts_ok = (n_manifest == K1_EXPECTED["manifest_selected"]
                    and n_manifest_declared == K1_EXPECTED["manifest_selected"]
                    and n == K1_EXPECTED["total"]
                    and len(negctrl) == K1_EXPECTED["negative_controls"]
                    and len(nonctrl) == K1_EXPECTED["noncontrol"]
                    and k1_counts["anchor_in_noncontrol"])
    k1 = {}
    for name in SETS:
        pass_all = sum(r[f"original_{name}"] for r in rows)
        pass_non = sum(r[f"original_{name}"] for r in nonctrl)
        hist_all = ref["controlled_counts"][name]
        hist_non = ref["controlled_counts_noncontrol"][name]
        k1[name] = dict(denominator_all=n, passing_all=pass_all, negative_controls=len(negctrl),
                        denominator_noncontrol=len(nonctrl), passing_noncontrol=pass_non,
                        historical_n=ref["n"], historical_all=hist_all, historical_noncontrol=hist_non,
                        agreement=(n == ref["n"] and pass_all == hist_all and pass_non == hist_non))
    checks["K1_original_reproduction"] = dict(pass_=k1_counts_ok and all(v["agreement"] for v in k1.values()),
                                              counts=k1_counts, expected=K1_EXPECTED, counts_ok=k1_counts_ok, sets=k1)

    # Common-M table, same layout, three categories
    cm = {}
    for name in SETS:
        cats = {}
        for label, sub in (("all", rows), ("noncontrol", nonctrl)):
            cats[label] = dict(denominator=len(sub),
                               PASS=sum(r[f"commonM_{name}"] == "PASS" for r in sub),
                               FAIL=sum(r[f"commonM_{name}"] == "FAIL" for r in sub),
                               BORDERLINE=sum(r[f"commonM_{name}"] == "BORDERLINE" for r in sub))
        cm[name] = dict(**cats, borderline_tags=[r["tag"] for r in rows if r[f"commonM_{name}"] == "BORDERLINE"])

    # K2 — Common-M PASS implies Original true
    viol = [(r["tag"], name) for r in rows for name in SETS
            if r[f"commonM_{name}"] == "PASS" and not r[f"original_{name}"]]
    checks["K2_implication"] = dict(pass_=not viol, violations=viol)

    # K3 — identity
    k3max = max(r["k3_identity_dev"] for r in rows)
    checks["K3_identity"] = dict(pass_=k3max <= K3_TOL, max_dev=k3max, tolerance=K3_TOL, M=K3_M)

    # K4 — rank-one deviation (report only; formula never switched)
    d = max(r["delta_rank1"] for r in rows)
    checks["K4_rank1_report"] = dict(
        max_delta_rank1=d, threshold=K4_REPORT_THRESHOLD,
        statement=("|kappa|/sqrt(p) agrees numerically with sqrt(s)" if d <= K4_REPORT_THRESHOLD
                   else "deviation above threshold; code coefficient used unchanged"),
        above_threshold=[r["tag"] for r in rows if r["delta_rank1"] > K4_REPORT_THRESHOLD])

    # K6 — anchor values for a hand check
    a = rows[0]
    checks["K6_anchor_hand_check"] = {k: a[k] for k in ("tag", "p", "s", "kappa", "G", "coef_phi", "coef_psi", "c_max",
                                                         "sqrtX_max", "margin_C_gStrong_rho0.3",
                                                         "M_low_C_gStrong_rho0.3", "M_high_C_gStrong_rho0.3",
                                                         "commonM_C_gStrong_rho0.3")}

    verdict = "ALL_CHECKS_PASS" if all(checks[k]["pass_"] for k in
                                       ("K1_original_reproduction", "K2_implication", "K3_identity", "K5_inputs")) \
        else "CHECKS_FAILED"
    summary = dict(spec=SPEC, n=n, thresholds={k: dict(g_max=v[0], rho_max=v[1], product=v[0] * v[1]) for k, v in SETS.items()},
                   original=k1, common_m=cm, checks=checks, verdict=verdict)
    (out / "b412_common_m.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8", newline="\n")

    cols = ["tag", "cls", "mode", "p", "s", "kappa", "sqrtX_max", "G", "c_max", "delta_rank1"]
    for name in SETS:
        cols += [f"original_{name}", f"commonM_{name}"]
    cols += ["margin_C_gStrong_rho0.3", "M_low_C_gStrong_rho0.3", "M_high_C_gStrong_rho0.3"]
    with open(out / "b412_common_m.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    for name in SETS:
        o, c_ = k1[name], cm[name]["all"]
        print(f"{name}: Original {o['passing_all']}/{o['denominator_all']} (noncontrol {o['passing_noncontrol']}/"
              f"{o['denominator_noncontrol']}; historical {o['historical_all']}, {o['historical_noncontrol']}) | "
              f"Common-M PASS {c_['PASS']} FAIL {c_['FAIL']} BORDERLINE {c_['BORDERLINE']}")
    for k in ("K1_original_reproduction", "K2_implication", "K3_identity", "K5_inputs"):
        print(k, "PASS" if checks[k]["pass_"] else "FAIL")
    print("K4", checks["K4_rank1_report"]["statement"], f"(max {d:.3e})")
    print("K6", json.dumps(checks["K6_anchor_hand_check"]))
    print(verdict)
    return 0 if verdict == "ALL_CHECKS_PASS" else 1


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(64)
    sys.exit(main(sys.argv[1], sys.argv[2]))
