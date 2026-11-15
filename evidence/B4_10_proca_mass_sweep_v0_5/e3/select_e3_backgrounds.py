"""E3 background selection (pre-registered rules; NO Proca sweep is run here).

Population: the T1.7.2 scan, generate_samples(10000, p,s in [0.05,1], Lambda in [0,1], cut 0.95, seed 172026),
both rank-one modes.  Classification uses only the LOCAL T1.7 quartic (analyze_state).
Rules (fixed before looking at results):
  negative_control : max|Im c| > 1e-6                       (complex local roots)
  boundary         : four real roots and rel_gap < 0.10     (rel_gap = min root separation / max|root|)
  stable           : four real roots and rel_gap >= 0.30
  (0.10 <= rel_gap < 0.30 is a counted but unselected buffer)
Selection: within a (mode,class) list sorted by Xw+Xv (and index for ties), take evenly spaced
quantile positions (np.linspace) -> deterministic, covers the range of p w^2 + s v^2.
Counts per mode: 8 stable, 6 boundary, 6 negative_control.
"""
import json, sys, hashlib
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent/"vendor_t173"))
from robustness_scan import generate_samples, analyze_state

N, SEED, CUT = 10000, 172026, 0.95
QUOTA = {"stable": 8, "boundary": 6, "negative_control": 6}
TH_BOUNDARY, TH_STABLE = 0.10, 0.30


def classify(res):
    if res["max_abs_imag_c"] > 1e-6:
        return "negative_control"
    roots = np.array([z[0] for z in res["roots_real"]])
    rel = res["min_root_separation"]/max(np.abs(roots).max(), 1e-300)
    if rel < TH_BOUNDARY:
        return "boundary"
    return "stable" if rel >= TH_STABLE else "buffer"


def main(out="."):
    pop = {"rank1_minus": {}, "rank1_plus": {}}
    stats = {m: {} for m in pop}
    for bg in generate_samples(N, (0.05, 1.0), (0.05, 1.0), (0.0, 1.0), CUT, SEED):
        for mode, sg in (("rank1_minus", -1.0), ("rank1_plus", 1.0)):
            k = sg*np.sqrt(bg["p"]*bg["s"])
            res = analyze_state(bg["w"], bg["q"], bg["v"], bg["r"], bg["p"], bg["s"], k, bg["Lambda"])
            cls = classify(res)
            roots = np.array([z[0] for z in res["roots_real"]])
            row = dict(index=bg["index"], mode=mode, cls=cls, p=bg["p"], s=bg["s"], kappa=float(k), Lambda=bg["Lambda"],
                       w=float(bg["w"]), q=float(bg["q"]), v=float(bg["v"]), r=float(bg["r"]),
                       Xq=float(bg["Xq"]), Xr=float(bg["Xr"]), Xw=float(bg["Xw"]), Xv=float(bg["Xv"]), temporal_sum=float(bg["Xw"] + bg["Xv"]),
                       local_roots_real_part=sorted(roots.tolist()), max_abs_imag_c=res["max_abs_imag_c"],
                       rel_gap=res["min_root_separation"]/max(np.abs(roots).max(), 1e-300))
            pop[mode].setdefault(cls, []).append(row)
            stats[mode][cls] = stats[mode].get(cls, 0) + 1
    chosen = []
    for mode in pop:
        for cls, n in QUOTA.items():
            lst = sorted(pop[mode].get(cls, []), key=lambda r: (r["temporal_sum"], r["index"]))
            if cls == "negative_control" or cls == "boundary":
                pass
            pos = sorted(set(np.round(np.linspace(0, len(lst) - 1, n)).astype(int).tolist()))
            for pz in pos:
                r = dict(lst[pz]); r["selection_rank_in_class"] = pz; r["class_size"] = len(lst)
                chosen.append(r)
    for r in chosen:
        mx = max(r["Xq"], r["Xr"], r["Xw"], r["Xv"])
        r["eft_control_cut0.5"] = bool(mx < 0.5); r["eft_control_cut0.7"] = bool(mx < 0.7); r["eft_control_cut0.95"] = bool(mx < 0.95)
    manifest = dict(population=dict(n=N, seed=SEED, cut=CUT), rules=dict(TH_BOUNDARY=TH_BOUNDARY, TH_STABLE=TH_STABLE, QUOTA=QUOTA,
                    negative_control="max|Im c|>1e-6"), class_counts=stats, n_selected=len(chosen), selected=chosen,
                    anchor_not_in_population="B3_max_speed_7_315096 (row 8899): not reproducible from seed 172026; reported separately")
    Path(out, "e3_manifest.json").write_text(json.dumps(manifest, indent=1))
    return manifest


if __name__ == "__main__":
    m = main(*(sys.argv[1:2] or ["."]))
    print(json.dumps(m["class_counts"], indent=1), m["n_selected"])
    import collections
    print(collections.Counter((r["mode"], r["cls"]) for r in m["selected"]))
