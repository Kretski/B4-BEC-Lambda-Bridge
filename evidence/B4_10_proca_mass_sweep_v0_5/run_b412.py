"""B4.12: control parameters for the 41 backgrounds + controlled-population front velocities."""
import csv, json, sys
from pathlib import Path
import numpy as np
import b410._b49 as _b49
from b410 import uv, causality as ca, experiment3 as e3
from b410.fixed_eft import load_b3
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/"vendor_t173"))
from robustness_scan import generate_samples  # noqa: E402


def main(out="results_b412", limit=None, population=True, n=10000):
    out = Path(out); out.mkdir(parents=True, exist_ok=True); _b49.verify_manifest()
    man = json.loads((ROOT/"data"/"e3_manifest.json").read_text())
    sel = man["selected"][:limit] if limit else man["selected"]
    cases = [("ANCHOR_B3", "anchor", load_b3())] + [(f"{r['mode']}_{r['cls']}_{r['index']}", r["cls"], e3.make_point(r)) for r in sel]
    rows = []
    for tag, cls, b3 in cases:
        c = uv.controlled(b3); inv = uv.invariants(b3, 7.0)
        allr, real = ca.principal_roots(b3)
        row = dict(tag=tag, cls=cls, p=b3.p, s=b3.s, Lambda_star_phi=inv["Lambda_star_phi"], Lambda_star_psi=inv["Lambda_star_psi"], **c,
                   identity_max_dev=max(abs(inv["g_times_rho"][k] - inv["sqrtX"][k]) for k in inv["sqrtX"]),
                   front_velocity=float(np.abs(allr.real).max()) if real else None)
        if real:
            xw = uv.x_window(b3); row["x_dev1pct"] = xw[0.01]; row["x_dev10pct"] = xw[0.1]
            row["front_minus_1_over_sqrtX"] = (row["front_velocity"] - 1)/c["sqrtX_max"]
        rows.append(row)
    keys = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(out/"b412_control_parameters.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    summ = dict(n=len(rows), thresholds={k: dict(g_max=v[0], rho_max=v[1], sqrtX_max=v[0]*v[1]) for k, v in uv.SETS.items()},
                controlled_counts={k: sum(r[k] for r in rows) for k in uv.SETS},
                controlled_counts_noncontrol={k: sum(r[k] for r in rows if r["cls"] != "negative_control") for k in uv.SETS},
                anchor=next(r for r in rows if r["tag"] == "ANCHOR_B3"),
                identity_max_dev=max(r["identity_max_dev"] for r in rows))
    if population:
        pop = {}
        for cut in (0.01, 0.04, 0.09, 0.25, 0.95):
            for mode, sg in (("rank1_minus", -1.0), ("rank1_plus", 1.0)):
                fr = []; sx = []; real = 0
                for bg in generate_samples(n, (0.05, 1.0), (0.05, 1.0), (0.0, 1.0), cut, 172026):
                    b3 = e3.make_point(dict(p=bg["p"], s=bg["s"], kappa=float(sg*np.sqrt(bg["p"]*bg["s"])), Lambda=bg["Lambda"], w=float(bg["w"]), q=float(bg["q"]), v=float(bg["v"]), r=float(bg["r"])))
                    allr, ok = ca.principal_roots(b3)
                    if ok:
                        real += 1; fr.append(float(np.abs(allr.real).max())); sx.append(max(np.sqrt(bg["Xq"]), np.sqrt(bg["Xr"]), np.sqrt(bg["Xw"]), np.sqrt(bg["Xv"])))
                fr = np.array(fr); sx = np.array(sx)
                pop[f"cut{cut}/{mode}"] = dict(n=n, principal_real=real, subluminal=int((fr <= 1).sum()), front_min=float(fr.min()),
                                              front_median=float(np.median(fr)), excess_over_sqrtX_median=float(np.median((fr - 1)/sx)),
                                              excess_over_sqrtX_min=float((( fr - 1)/sx).min()), excess_over_sqrtX_max=float(((fr - 1)/sx).max()))
        summ["population_front_velocity"] = pop
    (out/"b412_summary.json").write_text(json.dumps(summ, indent=1, default=float))
    print(json.dumps({k: summ[k] for k in ("controlled_counts", "controlled_counts_noncontrol", "identity_max_dev")}, indent=1))
    if population:
        for k, v in summ["population_front_velocity"].items(): print(k, {a: round(b, 4) if isinstance(b, float) else b for a, b in v.items()})
    print(summ["anchor"])


if __name__ == "__main__":
    main(*(sys.argv[1:2] or ["results_b412"]))
