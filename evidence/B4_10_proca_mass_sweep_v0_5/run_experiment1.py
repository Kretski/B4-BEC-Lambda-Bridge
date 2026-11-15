"""Run B4.10 experiment 1 and write results/*.csv and results/exp1_summary.json.

Usage:  python run_experiment1.py
"""
import csv
import json
from pathlib import Path
import numpy as np
from b410 import experiment1 as e1
from b410._b49 import track_eigenmodes, verify_manifest
from b410.exact import discriminant_report
from b410.fixed_eft import (load_b3, background, domain_report, leading_shift_coefficients,
                            contact_coefficients)
from b410.spectrum import universal_c, universal_omega_over_M

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)


def write_csv(name, rows):
    with (OUT / name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


b3 = load_b3()
masses = [1e-3, 0.1, 1.0, 7.0, 1e3, 1e6]

local = e1.local_limit_table(b3, masses)
recon = e1.coupling_reconstruction(b3, masses)
collapse = e1.scaling_collapse(b3, [0.01, 0.05, 0.14, 0.3, 0.5, 1.0, 2.0], [0.3, 1.0, 7.0])
conv = e1.shift_convergence(b3, [0.0025, 0.005, 0.01, 0.02, 0.04])
pts = [(1.0, 0.01), (1.0, 0.14), (1.0, 0.2077), (1.0, 0.5), (3.0, 1.0), (0.2, 0.05)]
verif = e1.verification_table(b3, pts)
vmap = e1.b49_validity_map(b3, [0.3, 1.0, 7.0, 40.0, 1e3, 1e5], [0.01, 0.14, 0.5, 2.0, 50.0])
for name, rows in (("exp1_local_limit.csv", local), ("exp1_coupling_reconstruction.csv", recon),
                   ("exp1_scaling_collapse.csv", collapse), ("exp1_shift_convergence.csv", conv),
                   ("exp1_verification.csv", verif), ("exp1_b49_validity_map.csv", vmap)):
    write_csv(name, rows)

# Exploratory (no claims beyond the numbers): the universal x-sweep of all six roots.
xs = np.logspace(-3, np.log10(50), 400)
sweep = []
for x in xs:
    om = universal_omega_over_M(b3, x)
    row = {"x": float(x)}
    for i, z in enumerate(om):
        row[f"omega_over_M_{i}_re"] = float(z.real); row[f"omega_over_M_{i}_im"] = float(z.imag)
    sweep.append(row)
write_csv("exp1b_universal_x_sweep.csv", sweep)

res = track_eigenmodes(background(b3, 1.0), np.linspace(0.02, 3.0, 300))
comp = []
for ik, k in enumerate(res["k"]):
    row = {"k_over_M": float(k), "mean_overlap": float(res["mean_overlap"][ik])}
    for j in range(6):
        row[f"mode{j}_omega_over_M"] = float(res["omega"][ik, j].real)
        row[f"mode{j}_scalar_fraction"] = float(res["scalar_fraction"][ik, j])
        row[f"mode{j}_null_residual"] = float(res["null_residual"][ik, j])
    comp.append(row)
write_csv("exp1b_tracked_composition.csv", comp)

c0, A = leading_shift_coefficients(b3)
summary = {
    "prescription": "fixed p,s,kappa,Lambda,(w,q,v,r); g_phi=M sqrt(p), g_psi=-M kappa/sqrt(p)",
    "b49_manifest_all_ok": all(ok for _, ok in verify_manifest()),
    "domain_invariants": domain_report(b3),
    "contact_coefficients_M_independent": contact_coefficients(b3),
    "local_limit_max_abs_diff_vs_B3": max(r["max_abs_diff_vs_B3"] for r in local),
    "local_limit_max_abs_diff_vs_independent_T17": max(r["max_abs_diff_vs_independent_T17"] for r in local),
    "coupling_reconstruction_max_rel_err": max(max(r["rel_err_p"], r["rel_err_s"], r["rel_err_kappa"]) for r in recon),
    "scaling_collapse_max_abs_dc": max(r["max_abs_dc_over_masses"] for r in collapse),
    "leading_shift": {"c_B3": c0.tolist(), "A_i_in_c_i=c_B3+A_i*x^2": A.tolist()},
    "dressed_proca_frequency_k_to_0_over_M": float(np.sqrt(1 + b3.p*b3.q**2 + b3.s*b3.r**2)),
    "verification_max_abs_diff_qz": max(r["abs_diff_qz"] for r in verif),
    "verification_max_abs_diff_mpmath": max(r["abs_diff_mpmath"] for r in verif),
    "verification_max_residual": max(r["residual"] for r in verif),
    "verification_min_simple_root_ratio": min(r["simple_root_ratio"] for r in verif),
    "verification_max_root_condition": max(r["root_condition"] for r in verif),
    "b49_dropped_root_cases": [r for r in vmap if r["b49_n_roots"] != 6][:6],
    "b49_dropped_root_case_count": sum(r["b49_n_roots"] != 6 for r in vmap),
    "universal_vs_mpmath_max_rel_err_over_map": max(r["universal_vs_mpmath_rel_err"] for r in vmap),
    "exact_discriminant": discriminant_report(b3),
    "exploratory_note": "exp1b_* files are exploratory outputs; no tested claims rest on them",
}
(OUT / "exp1_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps({k: v for k, v in summary.items() if k not in ("domain_invariants", "b49_dropped_root_cases")}, indent=2))
