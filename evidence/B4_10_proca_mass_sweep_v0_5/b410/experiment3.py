"""Experiment 3: multi-background analysis (pre-registered manifest, no criteria changed).

Per background the SAME audited machinery as Experiment 2 is used: QZ roots (+ mpmath), reduced
Hamiltonian h(k) from the original action (positivity independent of Krein), canonical energy,
Krein value, Proca energy share and composition exchange, branch gaps, group velocities, EFT flags.
Reality of the LOCAL quartic and of the FULL six-root system are recorded separately; EFT validity
is reported but never used to exclude a background.
"""
from __future__ import annotations
import numpy as np
from . import experiment2 as e2, energy as en, reference as ref
from . import audit_hamiltonian as ah
from .fixed_eft import B3Point, couplings

REAL_TOL = 1e-7          # max|Im omega| / max|omega|
ETAS = (0.1, 0.3)


def make_point(row) -> B3Point:
    return B3Point(p=row["p"], s=row["s"], kappa=row["kappa"], Lam=row["Lambda"], w=row["w"], q=row["q"],
                   v=row["v"], r=row["r"], roots=(), record={})


def light_node(b3: B3Point, x: float):
    g = couplings(b3, 1.0)
    a = (1.0, *g, b3.w, b3.q, b3.v, b3.r, b3.Lam)
    bl = ref.H_blocks(1.0, x, *g, b3.w, b3.q, b3.v, b3.r, b3.Lam)
    roots = ref.qz_roots(bl)
    scale = float(np.abs(roots).max())
    im = float(np.abs(roots.imag).max()/scale)
    h = ah.hamiltonian_matrix(x, *a)
    ev = np.linalg.eigvalsh(h)
    node = dict(x=float(x), real=bool(im <= REAL_TOL), max_rel_imag=im, h_min=float(ev[0]),
                h_norm=float(ev[-1]), h_neg=int(np.sum(ev < -1e-13*max(1.0, ev[-1]))), h_positive=bool(ev[0] > 1e-13*max(1.0, ev[-1])))
    if node["real"]:
        node.update(e2.evaluate(b3, x))
        node["real"] = True
    return node


def run_background(b3: B3Point, n0=241, x_lo=1e-3, x_hi=1e3, mp_every=20):
    xs = list(np.geomspace(x_lo, x_hi, n0))
    cache = {}
    nd = lambda x: cache.setdefault(x, light_node(b3, x))
    transitions = []
    unresolved = []
    i = 0
    while i < len(xs) - 1 and len(xs) < 20000:
        A, B = nd(xs[i]), nd(xs[i+1])
        dlog = np.log(xs[i+1]/xs[i])
        if A["real"] != B["real"]:
            if dlog > 1e-7:
                xs.insert(i+1, float(np.sqrt(xs[i]*xs[i+1]))); continue
            transitions.append(dict(x=xs[i], to_real=B["real"]))
        elif A["real"] and B["real"]:
            Pa = np.column_stack([b["vec"] for b in A["branches"]]); Pb = np.column_stack([b["vec"] for b in B["branches"]])
            perm, asg, mar = e2.match(Pa, Pb)
            if asg.min() < 0.98:
                if dlog > 1e-5:
                    xs.insert(i+1, float(np.sqrt(xs[i]*xs[i+1]))); continue
                unresolved.append(dict(x_left=xs[i], x_right=xs[i+1], min_assigned=float(asg.min())))
        i += 1
    nodes = [nd(x) for x in xs]
    # contiguous real segments
    segs, cur = [], []
    for n in nodes:
        if n["real"]:
            cur.append(n)
        elif cur:
            segs.append(cur); cur = []
    if cur:
        segs.append(cur)
    seg_data = []
    for sg in segs:
        if len(sg) < 3:
            seg_data.append(dict(nodes=sg, tbl=None, lab=None)); continue
        lab, st = e2.track(sg)
        seg_data.append(dict(nodes=sg, tbl=e2.tracked_table(sg, lab), lab=lab, min_assigned=min(s[0] for s in st)))
    # mpmath stability check (also on complex nodes: reality classification must agree)
    mp = []
    for n in nodes[::mp_every]:
        m = e2.mp_check(b3, n["x"])
        mp.append(dict(x=n["x"], diff=m["max_diff"], mp_complex=bool(m["mp_max_abs_imag"] > REAL_TOL),
                       qz_complex=not n["real"]))
    return dict(nodes=nodes, segments=seg_data, transitions=transitions, unresolved=unresolved, mp=mp)


def summarize_background(b3: B3Point, res, row=None):
    nodes, segs = res["nodes"], res["segments"]
    real = [n for n in nodes if n["real"]]
    cplx = [n for n in nodes if not n["real"]]
    out = dict(n_nodes=len(nodes), frac_nodes_real=len(real)/len(nodes),
               full_real_all_x=bool(len(cplx) == 0),
               full_real_small_x=bool(nodes[0]["real"]), full_real_large_x=bool(nodes[-1]["real"]),
               n_real_segments=len(segs), transitions=res["transitions"],
               max_rel_imag=float(max((n["max_rel_imag"] for n in nodes), default=0.0)),
               # Hamiltonian positivity vs reality (must be consistent)
               h_positive_all_nodes=bool(all(n["h_positive"] for n in nodes)),
               h_positive_at_all_real_nodes=bool(all(n["h_positive"] for n in real)) if real else None,
               h_positive_at_any_complex_node=bool(any(n["h_positive"] for n in cplx)),
               h_indefinite_at_real_nodes=int(sum(1 for n in real if not n["h_positive"])),
               mp_max_diff_real=float(max((m["diff"] for m in res["mp"] if not m["qz_complex"]), default=np.nan)),
               mp_reality_agrees=bool(all(m["mp_complex"] == m["qz_complex"] for m in res["mp"])),
               unresolved_overlap_intervals=len(res["unresolved"]))
    # energy / Krein over real nodes
    Es = [b["E"] for n in real for b in n["branches"]]
    out["n_real_branch_samples"] = len(Es)
    out["n_negative_energy_samples"] = int(sum(1 for E in Es if E <= 0))
    out["inertia_mismatch_nodes"] = int(sum(1 for n in real if n["h_neg"] != sum(1 for b in n["branches"] if b["E"] < 0)))
    out["krein_equals_sign_omega_all"] = bool(all(np.sign(b["K"]) == np.sign(b["omega"]) for n in real for b in n["branches"]))
    ksc = 0; ecr = []; ex_branches = set(); gaps = []; vg_max = 0.0; c_max = 0.0
    window = {eta: 0 for eta in ETAS}
    for sd in segs:
        for n in sd["nodes"]:
            om = np.array([b["omega"] for b in n["branches"]]); d = np.diff(om)
            gaps.append((n["x"], float(d.min()), float((d/np.maximum(np.abs(om[:-1]), np.abs(om[1:]))).min())))
            vg_max = max(vg_max, max(abs(b["vg"]) for b in n["branches"]))
            c_max = max(c_max, max(abs(b["c"]) for b in n["branches"]))
            for eta in ETAS:
                if n["x"] <= eta and sum(abs(b["omega"]) <= eta for b in n["branches"]) >= 4:
                    window[eta] = max(window[eta], n["x"])
        if sd["tbl"] is None:
            continue
        tbl = sd["tbl"]; xs = np.array([n["x"] for n in sd["nodes"]])
        for br in range(6):
            K = np.array([r[br]["K"] for r in tbl]); ksc += int(np.sum(np.diff(np.sign(K)) != 0))
            fe = np.array([e2.shares(r[br])["f_proca"] for r in tbl])
            for cx, dirn in e2.crossings(xs, fe):
                j = int(np.argmin(np.abs(np.log(xs/cx))))
                om_c = abs(tbl[j][br]["omega"])
                ecr.append(dict(branch=br, x=cx, direction=dirn, abs_omega_over_M=float(om_c),
                                in_window={str(eta): bool(cx <= eta and om_c <= eta) for eta in ETAS}))
            if (fe[0] - 0.5)*(fe[-1] - 0.5) < 0:
                ex_branches.add(br)
    out.update(krein_sign_changes=ksc, n_crossings=len(ecr), crossings=ecr,
               exchanged_branches_segment_ends=sorted(ex_branches),
               exchange_present=bool(len(ecr) > 0),
               exchange_inside_window={str(eta): bool(any(c["in_window"][str(eta)] for c in ecr)) for eta in ETAS},
               min_gap_abs=min((g[1] for g in gaps), default=None), min_gap_rel=min((g[2] for g in gaps), default=None),
               max_abs_group_velocity=vg_max, max_abs_phase_c=c_max,
               frac_real_samples_with_group_velocity_gt_1=float(np.mean([abs(b["vg"]) > 1 for n in real for b in n["branches"]])) if real else None,
               eft_window_x_max={str(eta): float(window[eta]) for eta in ETAS})
    mx = max(abs(b3.p*b3.q**2), abs(b3.s*b3.r**2), b3.p*b3.w**2, b3.s*b3.v**2)
    out["eft_control_max_X"] = float(mx)
    out["eft_control_cut0.5"] = bool(mx < 0.5); out["eft_control_cut0.7"] = bool(mx < 0.7); out["eft_control_cut0.95"] = bool(mx < 0.95)
    try:
        cinf, im = e2.c_infinity(b3)
        out["c_infinity_real"] = bool(im < 1e-9); out["c_infinity_max_abs"] = float(np.abs(cinf).max())
    except Exception:
        out["c_infinity_real"] = None
    return out
