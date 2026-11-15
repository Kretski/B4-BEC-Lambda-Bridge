"""Three separate results per background: lab energy, shift-symmetric (fixed-charge / fixed-mu) energy, dynamical reality.
python run_energy_triplet.py [outdir]"""
import csv, json, sys
from pathlib import Path
import numpy as np
import b410._b49 as _b49
from b410 import experiment3 as e3, energy as en, charges as ch
from b410.fixed_eft import load_b3, couplings

ROOT = Path(__file__).resolve().parent


def main(out="results_e3", limit=None):
    out = Path(out); out.mkdir(parents=True, exist_ok=True); _b49.verify_manifest()
    man = json.loads((ROOT/"data"/"e3_manifest.json").read_text())
    sel = man["selected"][:limit] if limit else man["selected"]
    cases = [("ANCHOR_B3", "rank1_minus", "anchor", load_b3())] + [(f"{r['mode']}_{r['cls']}_{r['index']}", r["mode"], r["cls"], e3.make_point(r)) for r in sel]
    rows = []; maxdev = 0.0
    for tag, mode, cls, b3 in cases:
        res = e3.run_background(b3, n0=121, mp_every=10**6)
        g = couplings(b3, 1.0); P = (1.0, *g, b3.w, b3.q, b3.v, b3.r, b3.Lam)
        n_real = 0; neg_lab = neg_sh = 0; samples = 0; min_lab_rel = np.inf; min_sh_rel = np.inf
        for n in res["nodes"]:
            if not n["real"]:
                continue
            n_real += 1
            for b in n["branches"]:
                e = ch.energies(P, en.u_matrix() @ b["vec"], b["omega"], n["x"])
                maxdev = max(maxdev, abs(e["E_shift"] - b["E"])/abs(b["E"]) if abs(b["E"]) > 0 else 0)
                samples += 1; neg_lab += e["E_lab"] <= 0; neg_sh += e["E_shift"] <= 0
                min_lab_rel = min(min_lab_rel, e["E_lab"]/abs(b["E"])); min_sh_rel = min(min_sh_rel, e["E_shift"]/abs(b["E"]))
        rows.append(dict(tag=tag, mode=mode, cls=cls, dynamic_real_fraction=n_real/len(res["nodes"]),
                         dynamic_complex_somewhere=n_real < len(res["nodes"]),
                         samples=samples, lab_energy_negative_samples=int(neg_lab), shift_energy_negative_samples=int(neg_sh),
                         min_Elab_over_absEshift=(None if samples == 0 else float(min_lab_rel)),
                         h_positive_all_nodes=all(n["h_positive"] for n in res["nodes"])))
    keys = list(rows[0].keys())
    with open(out/"e3_energy_triplet.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    summ = dict(max_rel_dev_Eshift_from_nonlinear_expansion_vs_Whitham=maxdev, n=len(rows))
    (out/"e3_energy_triplet_summary.json").write_text(json.dumps(summ, indent=1))
    return rows, summ


if __name__ == "__main__":
    rows, s = main(*(sys.argv[1:2] or ["results_e3"]))
    print(s)
    for r in rows:
        if r["shift_energy_negative_samples"] or r["lab_energy_negative_samples"] or r["dynamic_complex_somewhere"]:
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
