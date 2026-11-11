"""HeavyScalar_TwoField: population test of T2, T3, T4-proxy, T5-static, T5-Hamiltonian (inertia) and T5-Krein.

Model (M = 1): L = 1/2 X_phi + 1/2 X_psi + 1/2 (d chi)^2 - 1/2 M^2 chi^2 + chi (a X_phi + b X_psi + e Y),
background phi0 = w t + q x, psi0 = v t + r x, chi_bar = O/M^2.

Criteria and tolerances (fixed before the runs; the official T1-T6 in docs/UV_ACCEPTANCE_PROTOCOL.md are not modified):
  T2          real roots of the quartic factor of det K on X = k/M in [1e-3, 1e3] (61 log points), |Im c| <= 1e-7 max(1,|c|).
              det K = (c^2 - 1) * quartic exactly; the luminal pair is factored out symbolically (remainder asserted zero).
  T3          principal symbol diag(Z,1)(omega^2 - k^2): front = 1 iff min eig Z > 1e-12.
  T4-proxy    max |eig(Z) - 1| <= 0.3.  NOT the official T4 (g <= 1, rho <= 0.3), which stays unchecked for this model.
  T5-static   Z > 0 and min eig (Z - u u^T / M^2) >= -1e-12, u = (2a q + e r, 2b r + e q).
  T5-Hamiltonian  inertia of the 6x6 Hermitian quadratic Hamiltonian H2(k) in (q, pi) at k in {1e-2, 0.1, 1, 10}: number of negative eigenvalues.
  T5-Krein    for each normal mode at the same k: sign of v^dag H2 v (mode energy); number of negative-energy modes.
Outputs: one CSV row per background and a JSON summary (seed, tolerances, counts).
"""
import argparse, csv, json, os
import numpy as np, sympy as sp

KS_H = (1e-2, 0.1, 1.0, 10.0)
XGRID = np.geomspace(1e-3, 1e3, 61)
TOL_IM, TOL_EIG, PROXY = 1e-7, 1e-12, 0.3

w, q, v, r, a, b, e, M, k, c = sp.symbols('w q v r a b e M k c', real=True)
ft, fx, gt, gx, Ct, Cx, Cc = sp.symbols('ft fx gt gx Ct Cx Cc', real=True)
_L = (sp.Rational(1, 2)*(ft**2-fx**2) + sp.Rational(1, 2)*(gt**2-gx**2) + sp.Rational(1, 2)*(Ct**2-Cx**2)
      - sp.Rational(1, 2)*M**2*Cc**2 + Cc*(a*(ft**2-fx**2) + b*(gt**2-gx**2) + e*(ft*gt-fx*gx)))
_jets = [ft, fx, gt, gx, Ct, Cx, Cc]
_chi = (a*(w**2-q**2) + b*(v**2-r**2) + e*(w*v-q*r))/M**2
_H = sp.hessian(_L, _jets).subs({ft: w, fx: q, gt: v, gx: r, Ct: 0, Cx: 0, Cc: _chi})
_T = sp.zeros(7, 3)
_T[0, 0] = -sp.I*c*k; _T[1, 0] = sp.I*k; _T[2, 1] = -sp.I*c*k; _T[3, 1] = sp.I*k; _T[4, 2] = -sp.I*c*k; _T[5, 2] = sp.I*k; _T[6, 2] = 1
_D = sp.expand(sp.simplify((_T.H*_H*_T).applyfunc(sp.expand).det()))
_Q4, _rem = sp.div(sp.Poly(_D, c), sp.Poly(c**2-1, c))
assert sp.simplify(_rem.as_expr()) == 0, "det K is not divisible by (c^2 - 1)"
_ARGS = (w, q, v, r, a, b, e, M, k)
_coef4 = [sp.lambdify(_ARGS, co, 'numpy') for co in _Q4.all_coeffs()]
_Hnum = sp.lambdify((w, q, v, r, a, b, e, M), _H, 'numpy')


def quartic_roots(bg, th, kv, Mv=1.0):
    return np.roots([complex(f(*bg, *th, Mv, kv)) for f in _coef4])


def hamiltonian(bg, th, kv, Mv=1.0):
    """6x6 Hermitian H2(k) in (q, pi), q = (phi, psi, chi) Fourier amplitudes."""
    H = np.array(_Hnum(*bg, *th, Mv), dtype=complex)
    Tt = np.zeros((7, 3), dtype=complex); Ty = np.zeros((7, 3), dtype=complex)
    Tt[0, 0] = Tt[2, 1] = Tt[4, 2] = 1
    Ty[1, 0] = Ty[3, 1] = Ty[5, 2] = 1j*kv; Ty[6, 2] = 1
    A = Tt.conj().T@H@Tt; B = Tt.conj().T@H@Ty; C = -(Ty.conj().T@H@Ty)
    Ai = np.linalg.inv(A)
    Hm = np.block([[C + B.conj().T@Ai@B, -B.conj().T@Ai], [-Ai@B, Ai]])
    return 0.5*(Hm + Hm.conj().T), A


def evaluate(bg, th, Mv=1.0):
    w_, q_, v_, r_ = bg; a_, b_, e_ = th
    ch = (a_*(w_**2-q_**2) + b_*(v_**2-r_**2) + e_*(w_*v_-q_*r_))/Mv**2
    Z = np.array([[1+2*a_*ch, e_*ch], [e_*ch, 1+2*b_*ch]]); ez = np.linalg.eigvalsh(Z)
    u = np.array([2*a_*q_+e_*r_, 2*b_*r_+e_*q_])
    t3 = bool(ez.min() > TOL_EIG)
    t4p = bool(np.abs(ez-1).max() <= PROXY)
    es = float(np.linalg.eigvalsh(Z - np.outer(u, u)/Mv**2).min())
    t5s = bool(t3 and es >= -TOL_EIG)
    worst = 0.0
    for X in XGRID:
        rt = quartic_roots(bg, th, X*Mv, Mv)
        worst = max(worst, float(np.abs(rt.imag).max()/max(1.0, np.abs(rt).max())))
    t2 = bool(worst <= TOL_IM)
    neg_inertia = 0; neg_krein = 0; complex_modes = 0
    J = np.block([[np.zeros((3, 3)), np.eye(3)], [-np.eye(3), np.zeros((3, 3))]])
    for kv in KS_H:
        Hm, A = hamiltonian(bg, th, kv*Mv, Mv)
        neg_inertia = max(neg_inertia, int((np.linalg.eigvalsh(Hm) < -1e-10).sum()))
        om, V = np.linalg.eig(1j*J@Hm)          # v ~ exp(-i omega t)
        for j in range(6):
            if abs(om[j].imag) > 1e-7*max(1.0, abs(om[j])):
                complex_modes += 1; continue
            en = float(np.real(V[:, j].conj()@Hm@V[:, j]))
            if en < -1e-10*np.linalg.norm(Hm): neg_krein += 1
    return dict(w=w_, q=q_, v=v_, r=r_, chi_bar=ch, eigZ_min=float(ez.min()), eigZ_max=float(ez.max()), static_min_eig=es,
                T2=t2, T3=t3, T4_proxy=t4p, T5_static=t5s, max_rel_imag=worst,
                H_negative_eigs=neg_inertia, krein_negative_modes=neg_krein, complex_modes_at_Hk=complex_modes)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=21); ap.add_argument('--n', type=int, default=500)
    ap.add_argument('--out', default='results/heavy_scalar_population')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    sets = [(0.2, 0.2, 0.4), (0.1, 0.1, 0.2), (0.3, -0.2, 0.1)]
    rng = np.random.default_rng(args.seed)      # one generator, coupling sets in this order (reproduces the first run)
    summary = dict(seed=args.seed, n=args.n, M=1.0, background_range=[-1, 1], couplings=sets, tolerances=dict(
        T2_rel_imag=TOL_IM, eig=TOL_EIG, T4_proxy=PROXY, X_grid=[1e-3, 1e3, 61], k_for_H=list(KS_H)), results={})
    for th in sets:
        rows = [evaluate(tuple(rng.uniform(-1, 1, 4)), th) for _ in range(args.n)]
        tag = "a%g_b%g_e%g" % th
        with open(os.path.join(args.out, f"seed{args.seed}_{tag}.csv"), 'w', newline='') as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
        B = lambda key: np.array([x[key] for x in rows], dtype=bool)
        t2, t3, t4, t5 = B('T2'), B('T3'), B('T4_proxy'), B('T5_static')
        hin = np.array([x['H_negative_eigs'] for x in rows]); kr = np.array([x['krein_negative_modes'] for x in rows])
        summary['results'][str(th)] = dict(T2=int(t2.sum()), T3=int(t3.sum()), T4_proxy=int(t4.sum()), T5_static=int(t5.sum()),
            all_four=int((t2 & t3 & t4 & t5).sum()), T5static_but_not_T2=int((t5 & ~t2).sum()), T2_but_not_T5static=int((t2 & ~t5).sum()),
            H_positive_definite=int((hin == 0).sum()), H_posdef_equals_T5static=bool(((hin == 0) == t5).all()),
            backgrounds_with_negative_krein=int((kr > 0).sum()), negative_krein_within_T5static=int(((kr > 0) & t5).sum()),
            max_rel_imag_when_T5static=float(max([x['max_rel_imag'] for x in rows if x['T5_static']] or [0.0])))
    with open(os.path.join(args.out, f"seed{args.seed}_summary.json"), 'w') as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary['results'], indent=1))


if __name__ == "__main__":
    main()
