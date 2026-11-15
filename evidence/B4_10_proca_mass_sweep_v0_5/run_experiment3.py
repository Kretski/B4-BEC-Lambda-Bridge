"""B4.10 Experiment 3 driver: python run_experiment3.py [outdir]  (uses data/e3_manifest.json unchanged)"""
import csv, json, sys, collections
from pathlib import Path
import numpy as np
import b410._b49 as _b49
from b410 import experiment3 as e3
from b410.fixed_eft import load_b3

ROOT = Path(__file__).resolve().parent


def main(out="results_e3", limit=None):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    _b49.verify_manifest()
    man = json.loads((ROOT/"data"/"e3_manifest.json").read_text())
    sel = man["selected"][:limit] if limit else man["selected"]
    anchor = load_b3()
    cases = [dict(tag="ANCHOR_B3", mode="rank1_minus", cls="anchor", index=8899, point=anchor, row=None)]
    for r in sel:
        cases.append(dict(tag=f"{r['mode']}_{r['cls']}_{r['index']}", mode=r["mode"], cls=r["cls"], index=r["index"],
                          point=e3.make_point(r), row=r))
    f_nodes = open(out/"e3_nodes.csv", "w", newline=""); wn = csv.writer(f_nodes)
    wn.writerow("tag x branch omega_over_M E K sign_omega sign_K proca_energy_share group_velocity".split())
    summaries = []; cross_rows = []
    for c in cases:
        res = e3.run_background(c["point"])
        s = e3.summarize_background(c["point"], res)
        for sd in res["segments"]:
            if sd["tbl"] is None:
                continue
            for n, rows in zip(sd["nodes"], sd["tbl"]):
                for br, b in enumerate(rows):
                    wn.writerow([c["tag"], f"{n['x']:.10g}", br, f"{b['omega']:.10g}", f"{b['E']:.8g}", f"{b['K']:.8g}",
                                 int(np.sign(b['omega'])), int(np.sign(b['K'])), f"{e3.e2.shares(b)['f_proca']:.8g}", f"{b['vg']:.8g}"])
        rowinfo = dict(tag=c["tag"], mode=c["mode"], cls=c["cls"], bg_index=c["index"], **{k: v for k, v in s.items() if k not in ("crossings", "transitions")})
        if c["row"]:
            rowinfo.update(local_rel_gap=c["row"]["rel_gap"], local_max_imag_c=c["row"]["max_abs_imag_c"],
                           local_real=bool(c["row"]["max_abs_imag_c"] <= 1e-6), Xw_plus_Xv=c["row"]["temporal_sum"])
        else:
            rowinfo.update(local_real=True, Xw_plus_Xv=anchor.p*anchor.w**2 + anchor.s*anchor.v**2)
        rowinfo["transition_x"] = ";".join(f"{t['x']:.4g}->{'real' if t['to_real'] else 'complex'}" for t in s["transitions"])
        summaries.append(rowinfo)
        for cr in s["crossings"]:
            cross_rows.append(dict(tag=c["tag"], mode=c["mode"], cls=c["cls"], **{k: v for k, v in cr.items() if k != "in_window"},
                                   in_window_eta0_1=cr["in_window"]["0.1"], in_window_eta0_3=cr["in_window"]["0.3"]))
    f_nodes.close()
    for name, rows in (("e3_background_summary.csv", summaries), ("e3_exchange_crossings.csv", cross_rows)):
        if rows:
            keys = list(rows[0].keys())
            for r in rows:
                for k in r:
                    if k not in keys: keys.append(k)
            with open(out/name, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    # pre-registered outcome table by (mode, class) -- anchor kept separate
    table = {}
    for (mode, cls), grp in collections.defaultdict(list, {}).items():
        pass
    groups = collections.defaultdict(list)
    for r in summaries:
        groups[(r["mode"], r["cls"])].append(r)
    for (mode, cls), g in sorted(groups.items()):
        table[f"{mode}/{cls}"] = dict(
            n=len(g), local_real=sum(r["local_real"] for r in g), full_real_all_x=sum(r["full_real_all_x"] for r in g),
            full_real_never=sum(r["frac_nodes_real"] == 0 for r in g), full_partly_real=sum(0 < r["frac_nodes_real"] < 1 for r in g),
            h_positive_all_nodes=sum(r["h_positive_all_nodes"] for r in g),
            h_positive_at_complex_node=sum(r["h_positive_at_any_complex_node"] for r in g),
            h_indefinite_at_real_node=sum(r["h_indefinite_at_real_nodes"] > 0 for r in g),
            inertia_mismatch=sum(r['inertia_mismatch_nodes'] > 0 for r in g),
            energy_negative_somewhere=sum(r["n_negative_energy_samples"] > 0 for r in g),
            krein_sign_change=sum(r["krein_sign_changes"] > 0 for r in g),
            exchange_present=sum(r["exchange_present"] for r in g),
            exchange_in_window_eta0_1=sum(r["exchange_inside_window"]["0.1"] for r in g),
            exchange_in_window_eta0_3=sum(r["exchange_inside_window"]["0.3"] for r in g),
            group_velocity_gt1_somewhere=sum(r["max_abs_group_velocity"] > 1 for r in g if r["max_abs_group_velocity"] is not None),
            mp_reality_agrees=sum(r["mp_reality_agrees"] for r in g),
            unresolved_tracking=sum(r["unresolved_overlap_intervals"] > 0 for r in g))
    summary = dict(version="B4.10 Experiment 3 (package v0.3.0)", manifest_rules=man["rules"], n_backgrounds=len(sel),
                   anchor_reported_separately=True, outcome_table=table,
                   consistency=dict(
                       h_positive_but_complex_roots=sum(r["h_positive_at_any_complex_node"] for r in summaries),
                       h_indefinite_but_real_roots=sum(r["h_indefinite_at_real_nodes"] > 0 for r in summaries),
                       note="h positive definite must imply real roots and E>0; the converse (real roots, indefinite h) can occur only with mixed-sign energies"),
                   caveats=["stratified diagnostic sample, not a random sample of the background population",
                            "one grid (241 log-spaced x plus adaptive refinement), x in [1e-3,1e3]",
                            "group velocities are not front velocities",
                            "EFT validity is a separate flag (|omega|/M<=eta and x<=eta; protocol cuts 0.5/0.7/0.95), never used to exclude a background",
                            "g ~ M_A: perturbative UV completeness not established (B4.12)"])
    (out/"e3_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(table, indent=1))


if __name__ == "__main__":
    main(*(sys.argv[1:2] or ["results_e3"]))
