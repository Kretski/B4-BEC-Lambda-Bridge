"""B4.11 front-velocity (causality) analysis for the 41 backgrounds: python run_b411.py [outdir]"""
import csv, json, sys
from pathlib import Path
import numpy as np
import b410._b49 as _b49
from b410 import causality as ca, experiment3 as e3
from b410.fixed_eft import load_b3

ROOT = Path(__file__).resolve().parent


def main(out="results_b411", limit=None, time_domain=True):
    out = Path(out); out.mkdir(parents=True, exist_ok=True); _b49.verify_manifest()
    man = json.loads((ROOT/"data"/"e3_manifest.json").read_text())
    sel = man["selected"][:limit] if limit else man["selected"]
    cases = [("ANCHOR_B3", "rank1_minus", "anchor", load_b3())] + [(f"{r['mode']}_{r['cls']}_{r['index']}", r["mode"], r["cls"], e3.make_point(r)) for r in sel]
    rows = []
    for tag, mode, cls, b3 in cases:
        allr, real = ca.principal_roots(b3); lr, lreal = ca.local_roots(b3)
        row = dict(tag=tag, mode=mode, cls=cls, principal_roots_real=real, local_roots_real=lreal,
                   principal_c_min=float(allr.real.min()) if real else None, principal_c_max=float(allr.real.max()) if real else None,
                   front_velocity=float(np.abs(allr.real).max()) if real else None,
                   n_principal_speeds_above_light=int(np.sum(np.abs(allr.real) > 1)) if real else None,
                   local_front=float(np.abs(lr.real).max()) if lreal else None,
                   max_imag_c_principal=float(np.abs(allr.imag).max()))
        row["front_over_local"] = (row["front_velocity"]/row["local_front"]) if real and lreal else None
        row["subluminal_in_minkowski_cone"] = bool(real and row["front_velocity"] <= 1.0)
        if time_domain and (cls != "negative_control"):
            td, maxre = ca.evolve_front(b3)
            row["td_measured_speed_right"] = (td[8.0]["right"] - td[4.0]["right"])/4.0
            row["td_measured_speed_left"] = (td[8.0]["left"] - td[4.0]["left"])/4.0
            row["td_max_growth_rate_at_kmax"] = maxre
            row["td_blowup"] = bool(td[8.0]["max_amp"] > 1e6)
        rows.append(row)
    keys = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(out/"b411_front_velocity.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    ok = [r for r in rows if r["principal_roots_real"]]
    summ = dict(n=len(rows), principal_real=len(ok),
                principal_complex=[r["tag"] for r in rows if not r["principal_roots_real"] and r["cls"] != "negative_control"],
                negative_controls_principal_complex=sum((not r["principal_roots_real"]) for r in rows if r["cls"] == "negative_control"),
                front_velocity_min=min(r["front_velocity"] for r in ok), front_velocity_max=max(r["front_velocity"] for r in ok),
                front_velocity_median=float(np.median([r["front_velocity"] for r in ok])),
                n_subluminal=sum(r["subluminal_in_minkowski_cone"] for r in ok),
                front_over_local_min=min(r["front_over_local"] for r in ok if r["front_over_local"]),
                front_over_local_max=max(r["front_over_local"] for r in ok if r["front_over_local"]),
                td_max_rel_dev_right=max((abs(r["td_measured_speed_right"] - r["principal_c_max"])/abs(r["principal_c_max"]) for r in ok if "td_measured_speed_right" in r), default=None),
                td_max_rel_dev_left=max((abs(r["td_measured_speed_left"] - r["principal_c_min"])/abs(r["principal_c_min"]) for r in ok if "td_measured_speed_left" in r), default=None))
    (out/"b411_summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main(*(sys.argv[1:2] or ["results_b411"]))
