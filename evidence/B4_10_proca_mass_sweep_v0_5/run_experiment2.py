"""B4.10 Experiment 2 driver:  python run_experiment2.py [outdir]"""
import csv, json, sys
from pathlib import Path
import numpy as np
import b410._b49 as _b49
from b410 import experiment2 as e2, energy as en, reference as ref
from b410.fixed_eft import load_b3, couplings

ETAS = (0.1, 0.3)


def main(out="results", n0=241):
    out = Path(out); out.mkdir(exist_ok=True, parents=True)
    _b49.verify_manifest()
    b3 = load_b3()
    nodes, unresolved = e2.sweep(b3, n0=n0)
    lab, st = e2.track(nodes)
    tbl = e2.tracked_table(nodes, lab)
    xs = np.array([n["x"] for n in nodes])
    mark = e2.eft_mark(b3, tbl, ETAS)
    dom, X = e2.eft_domain_flags(b3); dom_strict, _ = e2.eft_domain_flags(b3, e2.STRICT_CUTS)
    a = e2._args(b3)

    # ---- energy tracking / mixing / krein / canonical energy tables -------------------
    f_trk = open(out/"exp2_energy_tracking.csv", "w", newline="")
    f_krn = open(out/"exp2_krein_signature.csv", "w", newline="")
    f_can = open(out/"exp2_canonical_energy.csv", "w", newline="")
    f_mix = open(out/"exp2_mixing_diagnostics.csv", "w", newline="")
    W = {k: None for k in "trk krn can mix".split()}
    w_trk = csv.writer(f_trk); w_krn = csv.writer(f_krn); w_can = csv.writer(f_can); w_mix = csv.writer(f_mix)
    w_trk.writerow("x branch omega_over_M c=omega/k group_velocity overlap_with_previous margin_to_next_best".split())
    w_krn.writerow("x branch omega_over_M krein_value sign_krein sign_omega sign_canonical_energy "
                   "krein_equals_sign_omega sigma_min sigma_2nd_smallest cond_sigma2_over_sigma1 "
                   "relative_root_condition residual".split())
    w_can.writerow("x branch omega_over_M E_whitham E_noether_terms E_time_domain E_phi E_psi E_A E_phiA E_psiA".split())
    w_mix.writerow("x branch f_scalar_energy f_proca_energy f_interaction_energy amp_proca_euclid a0_constraint_weight "
                   "gap_to_lower_branch gap_to_upper_branch background_in_protocol_domain status "
                   f"in_window_eta{ETAS[0]} in_window_eta{ETAS[1]}".split())
    tl = en.terms(*a)
    max_E_mis = 0.0
    for i, (nd, rows) in enumerate(zip(nodes, tbl)):
        for br, r in enumerate(rows):
            if i == 0:
                ov, mg = "", ""
            else:
                ov, mg = st[i-1][0], st[i-1][1]
            w_trk.writerow([f"{nd['x']:.12g}", br, f"{r['omega']:.12g}", f"{r['c']:.12g}", f"{r['vg']:.9g}", ov, mg])
            sK, sO, sE = int(np.sign(r["K"])), int(np.sign(r["omega"])), int(np.sign(r["E"]))
            w_krn.writerow([f"{nd['x']:.12g}", br, f"{r['omega']:.12g}", f"{r['K']:.12g}", sK, sO, sE,
                            int(sK == sO), f"{r['sigma_min']:.4g}", f"{r['sigma_2nd']:.6g}", f"{r['cond']:.4g}",
                            f"{r['root_cond']:.6g}", f"{r['residual']:.3g}"])
            z = en.u_matrix() @ r["vec"]
            Etd = en.energy_density_time_domain(tl, z, r["omega"], nd["x"], n=1024)
            g = r["groups"]
            max_E_mis = max(max_E_mis, abs(r["E"] - r["E_noether"])/abs(r["E"]), abs(r["E"] - Etd)/abs(r["E"]))
            w_can.writerow([f"{nd['x']:.12g}", br, f"{r['omega']:.12g}", f"{r['E']:.12g}", f"{r['E_noether']:.12g}",
                            f"{Etd:.12g}", f"{g['phi']:.9g}", f"{g['psi']:.9g}", f"{g['A']:.9g}",
                            f"{g['phiA']:.9g}", f"{g['psiA']:.9g}"])
            sh = e2.shares(r)
            lo = r["omega"] - rows[br-1]["omega"] if br > 0 else ""
            up = rows[br+1]["omega"] - r["omega"] if br < 5 else ""
            w_mix.writerow([f"{nd['x']:.12g}", br, f"{sh['f_scalar']:.9g}", f"{sh['f_proca']:.9g}", f"{sh['f_int']:.9g}",
                            f"{sh['amp_proca']:.9g}", f"{sh['a0_weight']:.9g}", lo, up, int(dom), mark[i]["status"],
                            int(mark[i][f"per_branch_in_window_eta{ETAS[0]}"][br]),
                            int(mark[i][f"per_branch_in_window_eta{ETAS[1]}"][br])])
    for f in (f_trk, f_krn, f_can, f_mix):
        f.close()

    # ---- stability checks: QZ vs mpmath at every node ----------------------------------
    mpc = [e2.mp_check(b3, float(x)) for x in xs]
    # ---- degeneracy report -------------------------------------------------------------
    pairs, clustered = e2.gap_report(nodes)
    sub = []
    for j in range(5):
        vals = []
        for i in range(len(nodes) - 1):
            pa = [tbl[i][j]["vec"], tbl[i][j+1]["vec"]]; pb = [tbl[i+1][j]["vec"], tbl[i+1][j+1]["vec"]]
            vals.append(e2.subspace_overlap(np.column_stack(pa), np.column_stack(pb)))
        sub.append(dict(pair=[j, j+1], min_subspace_overlap_between_nodes=float(min(vals))))
    cinf, cinf_im = e2.c_infinity(b3)
    deg = dict(
        n_nodes=len(nodes), x_range=[float(xs[0]), float(xs[-1])],
        exact_degeneracies_possible=False,
        reason="discriminant of the universal polynomial (degree 9 in eps=x^2/(1+x^2)) has only the real root eps=0 (Exp. 1); "
               "all six roots are simple for every x>0",
        near_degeneracy_threshold_relative=1e-3, near_degenerate_points=clustered,
        min_adjacent_gaps=pairs, subspace_overlap_tracking_adjacent_pairs=sub,
        unresolved_overlap_intervals=unresolved,
        note="the subspace tracker is exercised on a synthetic exact crossing in tests/test_experiment2.py",
    )
    (out/"exp2_degeneracy_report.json").write_text(json.dumps(deg, indent=1))

    # ---- summary -------------------------------------------------------------------------
    ex = e2.exchange_summary(tbl, xs)
    exch_branches = [e["branch"] for e in ex if e["exchanged_energy_metric"]]
    # eft windows
    win = {}
    for eta in ETAS:
        omax = max(abs(r["omega"]) for rows in tbl for r in rows)
        win[f"eta{eta}"] = dict(
            x_max_all_four_scalar_modes_inside=float(max([m["x"] for m in mark if m[f"modes_in_window_eta{eta}"] >= 4], default=0.0)),
            note="x such that |omega|/M<=eta for >=4 branches and x<=eta (operational criterion)")
    cross_x = sorted(c for e in ex for c in e["crossing_x_energy_metric"])
    omega_at_cross = []
    for e in ex:
        for cx in e["crossing_x_energy_metric"]:
            j = int(np.argmin(np.abs(np.log(xs/cx))))
            omega_at_cross.append(dict(branch=e["branch"], x=cx, abs_omega_over_M=abs(tbl[j][e["branch"]]["omega"])))
    summary = dict(
        version="B4.10 Experiment 2 (package v0.2.0)",
        background="B3_max_speed_7_315096 (row 8899, rank1_minus), p,s,kappa fixed, g_phi,g_psi ~ M",
        independent_variable="x = k/M_A (M=1), 1e-3..1e3, adaptive log grid",
        n_nodes=len(nodes),
        max_abs_imag_root_over_scale_qz=max(n["max_imag"] for n in nodes),
        mpmath_vs_qz_max_rel_diff=max(m["max_diff"] for m in mpc),
        mpmath_vs_qz_by_x=[dict(x=m['x'], rel_diff=m['max_diff']) for m in mpc[::20]],
        mpmath_max_rel_imag=max(m["mp_max_abs_imag"] for m in mpc),
        complex_roots_found=bool(max(m["mp_max_abs_imag"] for m in mpc) > 1e-9 or max(n["max_imag"] for n in nodes) > 1e-9),
        tracking=dict(label_permutation_is_omega_order=bool((lab == lab[0]).all()),
                      min_assigned_overlap=min(s[0] for s in st), min_margin_euclid=min(s[1] for s in st),
                      unresolved_intervals=len(unresolved),
                      note="Euclidean margins are small because null vectors of a polynomial pencil are not orthogonal; "
                           "tracking is confirmed instead by assigned overlap >= 0.98 on refined intervals and by agreement "
                           "with omega-ordering (roots are real and simple, so order is invariant)"),
        krein=dict(all_sign_krein_equals_sign_omega=bool(all(e["krein_sign"] == e["omega_sign"] for e in ex)),
                   sign_changes_total=int(sum(e["krein_sign_changes"] for e in ex)),
                   energy_sign_changes_total=int(sum(e["energy_sign_changes"] for e in ex)),
                   min_abs_krein_unit_norm=min(e["min_abs_krein_unit_norm"] for e in ex)),
        canonical_energy=dict(all_positive=bool(all(r["E"] > 0 for rows in tbl for r in rows)),
                              max_rel_mismatch_whitham_vs_noether_vs_timedomain=max_E_mis),
        conditioning=dict(max_sigma_min=max(r["sigma_min"] for rows in tbl for r in rows),
                          min_sigma2_over_sigma1=min(r["cond"] for rows in tbl for r in rows),
                          max_residual=max(r["residual"] for rows in tbl for r in rows),
                          max_relative_root_condition=max(r["root_cond"] for rows in tbl for r in rows)),
        asymptotics=dict(c_large_x=[float(r['c']) for r in tbl[-1]], c_infinity_exact=[float(c) for c in cinf],
                         max_abs_diff_c_large_x_vs_exact=float(np.abs(np.array([r['c'] for r in tbl[-1]]) - cinf).max()),
                         c_small_x_scalar_branches=[float(tbl[0][j]['c']) for j in (1,2,3,4)], c_B3=[float(c) for c in b3.roots]),
        exchange=ex, exchanged_branches_energy_metric=exch_branches,
        exchange_crossing_x_range=[cross_x[0], cross_x[-1]] if cross_x else None,
        abs_omega_over_M_at_exchange=omega_at_cross,
        eft=dict(protocol_domain_recovered=dict(definition="p q^2, s r^2, p w^2, s v^2 all < 0.95 (T1.7.1 rank_one_grid defaults)",
                                                satisfied_by_B3_row=dom, Xvalues=X),
                 strict_cut_0p5_satisfied=dom_strict,
                 windows=win,
                 uv_completeness="NOT established: g_phi, g_psi ~ M_A so M_A -> infinity at fixed p,s,kappa is a formal limit; "
                                 "perturbativity of the UV completion is not shown",
                 verdict="every point is marked 'unestablished' for physical EFT validity"),
        caveats=["L2 was first reconstructed from the B4.9 symbol; the audit (AUDIT_E2.md, tests/test_audit.py) re-derives it from the T1.7 action by sympy and builds the reduced Hamiltonian independently. The action itself is taken from the T1.7 text (J_phi^1=phi_t phi_x)",
                 "one background point (B3 row), not a global statement",
                 "composition metrics (energy shares, Euclid amplitude) are decompositions, not unique mixing angles; "
                 "reported both and compared",
                 "positive energy for every sampled branch means no negative-energy hybridisation: exchange is "
                 "composition exchange between positive-energy branches, not an instability"],
    )
    (out/"exp2_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: summary[k] for k in ("n_nodes", "complex_roots_found", "mpmath_vs_qz_max_rel_diff", "krein", "canonical_energy",
                                              "exchanged_branches_energy_metric", "exchange_crossing_x_range")}, indent=1))


if __name__ == "__main__":
    main(*(sys.argv[1:2]))
