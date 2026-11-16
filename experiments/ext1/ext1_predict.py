"""EXT-1 prediction script (registration step only; NO comparison).

Computes the mean-field prediction of the spin/density sound-speed ratio
    gamma_pred = c_s / c_d
for a symmetric two-component 23Na mixture (equal masses, equal densities,
V1 = V2 = 0, Omega_R = 0) from the scattering lengths a and a_updn in
ext1_data_v3.json, using theory/bec_dispersion.py of the repository.

It does NOT read gamma_exp or any measured sound speed. It makes NO comparison.

Usage:  python3 -I ext1_predict.py <data_json> <repo_root> <out_json>
"""
import hashlib
import json
import math
import platform
import sys
from fractions import Fraction

DATA_SHA = "e7752cf0e949c5abe613298eafd5b8bac78f2787d6f8cdcc8e00f47d117b19ad"
MODULE_SHA = "228086f4a1a3ab292c1798910dffcc4ba69cec9e17c88dd7046369f52b878865"
SEED = 20261116
N_MC = 10_000
RHOS = [0.0, -0.9, -0.5, 0.5, 0.9]

# CODATA 2018
HBAR = 1.054571817e-34      # J s
A0 = 5.29177210903e-11      # m
U = 1.66053906660e-27       # kg
M_NA = 22.98976928 * U      # kg


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def get_entry(data, entry_id):
    """Return the unique entry with the given id; abort if missing or duplicated."""
    hits = [e for e in data.get("entries", []) if isinstance(e, dict) and e.get("id") == entry_id]
    if len(hits) != 1:
        sys.exit(f"ABORT: entry '{entry_id}' found {len(hits)} times (expected 1)")
    return hits[0]


def read_length(data, entry_id):
    """Scattering length and 1-sigma uncertainty (Bohr radii) from an entry."""
    e = get_entry(data, entry_id)
    v, u, unit = e.get("value"), e.get("uncertainty_as_printed"), e.get("unit")
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v <= 0:
        sys.exit(f"ABORT: invalid value in '{entry_id}': {v!r}")
    if not isinstance(u, str):
        sys.exit(f"ABORT: missing uncertainty in '{entry_id}'")
    try:
        u = float(u)
    except ValueError:
        sys.exit(f"ABORT: non-numeric uncertainty in '{entry_id}': {u!r}")
    if not math.isfinite(u) or u <= 0:
        sys.exit(f"ABORT: invalid uncertainty in '{entry_id}': {u!r}")
    if not isinstance(unit, str) or not unit.startswith("a0"):
        sys.exit(f"ABORT: unexpected unit in '{entry_id}': {unit!r}")
    return float(v), u


def read_positive_number(data, entry_id):
    e = get_entry(data, entry_id)
    v = e.get("value")
    if isinstance(v, str):
        try:
            v = float(Fraction(v))
        except (ValueError, ZeroDivisionError):
            sys.exit(f"ABORT: non-numeric value in '{entry_id}': {v!r}")
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v <= 0:
        sys.exit(f"ABORT: invalid value in '{entry_id}': {v!r}")
    return float(v)


def main():
    data_path, repo, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    if sha(data_path) != DATA_SHA:
        sys.exit("ABORT: data sha256 mismatch")
    mod_path = repo.rstrip("/") + "/theory/bec_dispersion.py"
    if sha(mod_path) != MODULE_SHA:
        sys.exit("ABORT: module sha256 mismatch")
    sys.path.insert(0, repo.rstrip("/") + "/theory")
    import numpy as np
    import bec_dispersion as bd

    data = json.load(open(data_path, encoding="utf-8"))
    return data, np, bd, data_path, mod_path, out_path


def run():
    data, np, bd, data_path, mod_path, out_path = main()

    # ---- inputs: read from the data file by entry id (gamma_exp NOT read) ----
    a_c, s_a = read_length(data, "a")
    b_c, s_b = read_length(data, "a_updn")
    # cross-check against the primary-source (Knoop) records
    ka, ks_a = read_length(data, "knoop_a")
    kb, ks_b = read_length(data, "knoop_a_updn")
    if (ka, ks_a, kb, ks_b) != (a_c, s_a, b_c, s_b):
        sys.exit("ABORT: D1 inputs differ from the Knoop records")
    if not a_c > b_c:
        sys.exit("ABORT: a <= a_updn (no real spin sound speed)")
    n0a3 = read_positive_number(data, "n0_a3")
    neff = read_positive_number(data, "n_eff_factor")
    input_ids = ["a", "a_updn", "knoop_a", "knoop_a_updn", "n0_a3", "n_eff_factor"]

    def g_of(a_bohr):
        return 4.0 * math.pi * HBAR**2 * (a_bohr * A0) / M_NA

    # reference total density (illustrative only, enters no ratio):
    # n0 a^3 ~ 1.4e-6 (no uncertainty published); n_eff = (2/3) n0
    n0 = n0a3 / (a_c * A0) ** 3
    n_tot = neff * n0
    n_i = n_tot / 2.0

    def speeds(a_bohr, b_bohr, ni):
        """c_d, c_s in mm/s via the repository module (hydrodynamic roots).
        Scaled units: m1=m2=1, n1=n2=1, g_ii = g*n_i/m * 1e6 (mm^2/s^2)."""
        k = ni / M_NA * 1e6
        p = bd.BECParams(m1=1.0, m2=1.0, n1=1.0, n2=1.0,
                         g11=g_of(a_bohr) * k, g22=g_of(a_bohr) * k,
                         g12=g_of(b_bohr) * k)
        r = bd.hydrodynamic_roots(p)        # sorted [-c_d, -c_s, c_s, c_d]
        return float(r[3]), float(r[2])

    def gamma_closed(a, b):
        return math.sqrt((a - b) / (a + b))

    # ---- central prediction ----
    c_d, c_s = speeds(a_c, b_c, n_i)
    gam = c_s / c_d
    gam_cf = gamma_closed(a_c, b_c)

    # ---- internal checks ----
    cd_exp = math.sqrt((g_of(a_c) + g_of(b_c)) * n_i / M_NA) * 1e3
    cs_exp = math.sqrt((g_of(a_c) - g_of(b_c)) * n_i / M_NA) * 1e3
    E1 = abs(c_d / cd_exp - 1) < 1e-10 and abs(c_s / cs_exp - 1) < 1e-10
    E2 = abs(gam / gam_cf - 1) < 1e-10
    cd2, cs2 = speeds(a_c, b_c, 10.0 * n_i)
    g_ratio = g_of(a_c) / g_of(b_c)
    E3 = (abs((cs2 / cd2) / gam - 1) < 1e-10
          and abs(g_ratio / (a_c / b_c) - 1) < 1e-14
          and abs(cd2 / c_d / math.sqrt(10.0) - 1) < 1e-10)

    # ---- analytic linear propagation ----
    dA = b_c / (gam_cf * (a_c + b_c) ** 2)
    dB = -a_c / (gam_cf * (a_c + b_c) ** 2)
    analytic = {}
    for rho in RHOS:
        var = (dA * s_a) ** 2 + (dB * s_b) ** 2 + 2 * rho * dA * dB * s_a * s_b
        analytic[f"{rho:+.1f}"] = {"rho": rho, "sigma_gamma": math.sqrt(var)}

    # ---- Monte Carlo ----
    mc = {}
    for rho in RHOS:
        rng = np.random.default_rng(SEED)
        cov = [[s_a**2, rho * s_a * s_b], [rho * s_a * s_b, s_b**2]]
        smp = rng.multivariate_normal([a_c, b_c], cov, size=N_MC)
        vals, fail = [], 0
        for a, b in smp:
            try:
                if not (a > b > 0):
                    raise ValueError("outside a > a_updn > 0")
                cd, cs = speeds(a, b, n_i)
                vals.append(cs / cd)
            except ValueError:
                fail += 1
        v = np.array(vals)
        mc[f"{rho:+.1f}"] = {
            "rho": rho, "n_samples": N_MC, "n_failed": fail,
            "mean": float(v.mean()), "std": float(v.std(ddof=1)),
            "median": float(np.median(v)),
            "p16": float(np.percentile(v, 16)), "p84": float(np.percentile(v, 84)),
        }

    out = {
        "file": "ext1_predictions.json",
        "status": "PREDICTION COMPUTED — AWAITING REVIEW AND REGISTRATION. NO COMPARISON MADE. gamma_exp was not read. Evidence level 3 at most; not evidence for B4.",
        "model": "Symmetric two-component mean-field BEC, equal masses (23Na), equal densities, V1=V2=0, Omega_R=0, no LHY, hydrodynamic (k->0) limit; c^2 = (g +- g12) n_i / m.",
        "observable": "gamma_pred = c_s / c_d (spin / density sound-speed ratio)",
        "inputs": {
            "source_entry_ids": input_ids,
            "a_bohr": a_c, "sigma_a_bohr": s_a,
            "a_updn_bohr": b_c, "sigma_a_updn_bohr": s_b,
            "uncertainties_treated_as_1sigma": True,
            "primary_rho": 0.0, "sensitivity_rhos": [-0.9, -0.5, 0.5, 0.9],
            "constants_CODATA2018": {"hbar": HBAR, "a0": A0, "u": U, "m_Na_u": 22.98976928},
            "n0_a3": n0a3, "n_eff_factor": neff,
            "reference_density_note": "n_i = (1/2) n_eff_factor n0 with n0 = n0_a3 / a^3 (no uncertainty published); illustrative only, ratio is density-independent (E3).",
        },
        "central": {
            "gamma_pred_module": gam, "gamma_pred_closed_form": gam_cf,
            "c_d_mm_s_reference_density": c_d, "c_s_mm_s_reference_density": c_s,
            "c_absolute_not_a_prediction": True,
        },
        "analytic_linear": analytic,
        "monte_carlo": {"seed": SEED, "generator": "numpy.random.default_rng / multivariate_normal, same seed per rho", "per_rho": mc},
        "internal_checks": {"E1_module_vs_sqrt_formula": E1, "E2_ratio_vs_closed_form": E2,
                            "E3_g_conversion_and_density_independence": E3},
        "provenance": {
            "script_sha256": sha(__file__),
            "data_file_sha256": sha(data_path), "module_sha256": sha(mod_path),
            "repo_head_at_run": "a9b950d1d83d7a27dc434dc9471164af48293b32",
            "python": platform.python_version(), "numpy": np.__version__,
        },
        "open_items": ["Journal versions of Kim et al. (PRA 101, 061601(R)) and Knoop et al. (PRA 83, 042704) not yet checked; required before comparison.",
                       "Comparison test NOT authorized."],
    }
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print("E1", E1, "E2", E2, "E3", E3)
    print("gamma_pred", gam, "closed", gam_cf)
    for k in analytic:
        print(k, "analytic", analytic[k]["sigma_gamma"], "MC", mc[k]["mean"], mc[k]["std"], "fail", mc[k]["n_failed"])
    if not (E1 and E2 and E3):
        sys.exit("INTERNAL CHECK FAILED")


if __name__ == "__main__":
    run()
