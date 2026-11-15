"""Experiment 2: energy-weighted mode mixing at fixed EFT coefficients (B3 background).

Universal variable x = k/M_A at M = 1 (exact scaling omega/M = f(x), see Exp. 1).
Three diagnostics are computed and stored SEPARATELY:
  * eigenvector overlap  (tracking only)
  * Krein value K_n = z^dag H'(omega) z   (and its sign)
  * canonical energy E_n = (1/4) omega z^dag H' z  (= Noether T^00 average)
"""
from __future__ import annotations
import numpy as np
import scipy.linalg as sla
from scipy.optimize import linear_sum_assignment
from . import reference as ref
from . import energy as en
from .fixed_eft import B3Point, couplings, _polys

# --- protocol EFT domain (recovered, see AUDIT_NOTE) ------------------------------
PROTOCOL_CUTS = {"max_normalized": 0.95, "max_w2": 0.95, "max_v2": 0.95}
STRICT_CUTS = {"max_normalized": 0.5, "max_w2": 0.5, "max_v2": 0.5}


def eft_domain_flags(b3: B3Point, cuts=PROTOCOL_CUTS):
    X = dict(Xq=b3.p*b3.q**2, Xr=b3.s*b3.r**2, Xw=b3.p*b3.w**2, Xv=b3.s*b3.v**2)
    ok = (X["Xq"] < cuts["max_normalized"] and X["Xr"] < cuts["max_normalized"]
          and X["Xw"] < cuts["max_w2"] and X["Xv"] < cuts["max_v2"])
    return bool(ok), X


# --- single-x evaluation ---------------------------------------------------------
def _args(b3):
    gphi, gpsi = couplings(b3, 1.0)
    return (1.0, gphi, gpsi, b3.w, b3.q, b3.v, b3.r, b3.Lam)


def evaluate(b3: B3Point, x: float):
    """All six branches at k/M = x (M=1), ordered by omega."""
    a = _args(b3)
    M, gphi, gpsi, w, q, v, r, Lam = a
    bl = ref.H_blocks(M, x, gphi, gpsi, w, q, v, r, Lam)
    roots = ref.qz_roots(bl)
    if len(roots) != 6:
        raise RuntimeError(f"expected 6 finite roots, got {len(roots)} at x={x}")
    max_im = float(np.abs(roots.imag).max())
    om = roots.real
    rows = []
    for o in om:
        H = ref.H_of(o, bl)
        xv, s1, s2 = en.null_vector_b49(H)
        xv = xv/np.linalg.norm(xv)
        Hp = 2*o*bl[2] + bl[1]
        K = en.krein_value(Hp, xv)
        E = 0.25*o*K
        ge = en.group_energies(*a, o, x, xv)
        diag = ref.root_diagnostics(o, bl)
        rows.append(dict(omega=float(o), c=float(o/x), vec=xv, K=K, E=E, groups=ge,
                         E_noether=sum(ge.values()), sigma_min=s1, sigma_2nd=s2,
                         cond=s2/max(s1, 1e-300), root_cond=diag["root_condition"],
                         residual=diag["residual"],
                         vg=float(-(xv@_Hk(a, o, x)@xv)/(xv@Hp@xv))))
    return dict(x=x, max_imag=max_im, branches=rows)


def _Hk(a, om, k, h=1e-6):
    M, gphi, gpsi, w, q, v, r, Lam = a
    f = lambda kk: ref.H_of(om, ref.H_blocks(M, kk, gphi, gpsi, w, q, v, r, Lam))
    hh = h*max(1.0, k)
    return (f(k + hh) - f(k - hh))/(2*hh)


def shares(row):
    """Composition metrics of one branch.

    Energy shares (signed, sum to 1): scalar = (phi+psi groups), proca = A group
    (a0,a1 with A0 slaved by its constraint), interaction = phiA+psiA.
    Amplitude share: Euclidean |a1|^2/(|phi|^2+|psi|^2+|a1|^2) of the unit null
    vector (A0 excluded: it is constrained), a cross-check metric only.
    """
    g, E = row["groups"], row["E_noether"]
    v = np.abs(row["vec"])**2
    amp = v[3]/(v[0] + v[1] + v[3])
    return dict(f_scalar=(g["phi"] + g["psi"])/E, f_proca=g["A"]/E,
                f_int=(g["phiA"] + g["psiA"])/E, f_phi=g["phi"]/E, f_psi=g["psi"]/E,
                amp_proca=float(amp), a0_weight=float(v[2]))


# --- tracking --------------------------------------------------------------------
def overlap_matrix(Va, Vb):
    return np.abs(Va.conj().T @ Vb)


def subspace_overlap(Va, Vb):
    """Overlap of the spans of orthonormalized columns: ||Qa^† Qb||_F^2 / dim."""
    Qa, _ = np.linalg.qr(Va); Qb, _ = np.linalg.qr(Vb)
    return float(np.linalg.norm(Qa.conj().T @ Qb, 'fro')**2/Qa.shape[1])


def match(Pa, Pb):
    """Hungarian assignment on overlap matrix; returns (perm, assigned, margin)."""
    Ov = overlap_matrix(Pa, Pb)
    ri, ci = linear_sum_assignment(-Ov)
    assigned = Ov[ri, ci]
    margin = np.empty(len(ri))
    for n, (i, j) in enumerate(zip(ri, ci)):
        other = np.delete(Ov[i], j)
        margin[n] = Ov[i, j] - other.max()
    return ci, assigned, margin


def sweep(b3: B3Point, x_lo=1e-3, x_hi=1e3, n0=121, min_assigned=0.98, min_margin=0.0,
          min_dlog=1e-5, max_nodes=20000):
    """Adaptive log-x sweep; refines an interval until overlap tracking is unambiguous."""
    xs = list(np.geomspace(x_lo, x_hi, n0))
    cache = {}
    ev = lambda x: cache.setdefault(x, evaluate(b3, x))
    unresolved = []
    i = 0
    while i < len(xs) - 1 and len(xs) < max_nodes:
        A, B = ev(xs[i]), ev(xs[i+1])
        Pa = np.column_stack([b["vec"] for b in A["branches"]])
        Pb = np.column_stack([b["vec"] for b in B["branches"]])
        perm, asg, mar = match(Pa, Pb)
        bad = asg.min() < min_assigned or mar.min() < min_margin
        dlog = np.log(xs[i+1]/xs[i])
        if bad and dlog > min_dlog:
            xs.insert(i+1, float(np.sqrt(xs[i]*xs[i+1])))
            continue
        if bad:
            unresolved.append(dict(x_left=xs[i], x_right=xs[i+1],
                                   min_assigned=float(asg.min()), min_margin=float(mar.min())))
        i += 1
    nodes = [ev(x) for x in xs]
    return nodes, unresolved


def track(nodes):
    """Branch labels by overlap (Hungarian); returns label[node][position]."""
    n = len(nodes)
    lab = np.zeros((n, 6), int)
    lab[0] = np.arange(6)
    stats = []
    for i in range(1, n):
        Pa = np.column_stack([b["vec"] for b in nodes[i-1]["branches"]])
        Pb = np.column_stack([b["vec"] for b in nodes[i]["branches"]])
        perm, asg, mar = match(Pa, Pb)       # position j in node i-1 -> perm[j] in node i
        new = np.empty(6, int)
        for j in range(6):
            new[perm[j]] = lab[i-1][j]
        lab[i] = new
        stats.append((float(asg.min()), float(mar.min())))
    return lab, stats


def tracked_table(nodes, lab):
    """rows[node][branch] with the branch at its tracked label."""
    out = []
    for nd, lb in zip(nodes, lab):
        by = {int(lb[pos]): nd["branches"][pos] for pos in range(6)}
        out.append([by[b] for b in range(6)])
    return out


# --- degeneracy / near-degeneracy ---------------------------------------------------
def gap_report(nodes, rel_tol=1e-3):
    """Minimum adjacent-root gaps (absolute in omega/M and relative to scale) per sweep."""
    res = []
    for nd in nodes:
        om = np.array([b["omega"] for b in nd["branches"]])
        d = np.diff(om)
        rel = d/np.maximum(np.abs(om[:-1]), np.abs(om[1:]))
        res.append((nd["x"], d, rel))
    xs = np.array([r[0] for r in res])
    D = np.array([r[1] for r in res]); R = np.array([r[2] for r in res])
    pairs = []
    for j in range(5):
        i = int(np.argmin(D[:, j]))
        pairs.append(dict(pair=[j, j+1], min_gap_omega_over_M=float(D[i, j]), at_x=float(xs[i]),
                          min_rel_gap=float(R[:, j].min()), at_x_rel=float(xs[int(np.argmin(R[:, j]))])))
    clustered = [dict(x=float(x), pair=j) for x, r in zip(xs, R) for j in range(5) if r[j] < rel_tol]
    return pairs, clustered


# --- stability: QZ vs mpmath --------------------------------------------------------
def mp_check(b3: B3Point, x: float, dps=60):
    a = _args(b3)
    M, gphi, gpsi, w, q, v, r, Lam = a
    mr = ref.mp_roots(M, x, gphi, gpsi, w, q, v, r, Lam, dps=dps)
    mr = np.array(sorted([complex(z) for z in mr], key=lambda z: (z.real, z.imag)))
    qz = np.array(sorted(ref.qz_roots(ref.H_blocks(M, x, gphi, gpsi, w, q, v, r, Lam)),
                         key=lambda z: (z.real, z.imag)))
    scale = np.abs(qz).max()
    return dict(x=x, max_diff=float(np.abs(mr - qz).max()/scale),
                mp_max_abs_imag=float(np.abs(mr.imag).max()/scale),
                qz_max_abs_imag=float(np.abs(qz.imag).max()/scale))


# --- exact asymptotics ------------------------------------------------------------
def c_infinity(b3: B3Point):
    """x -> infinity roots of Q_inf(c) = F(c) - c^2 d1 d2 (eps -> 1)."""
    P = np.polynomial.polynomial
    d1, d2, F = _polys(b3)
    Q = P.polysub(F, P.polymul(P.polymul([0, 0, 1.0], d1), d2))
    return np.sort(P.polyroots(Q).real), np.abs(P.polyroots(Q).imag).max()


# --- EFT marking -------------------------------------------------------------------
def eft_mark(b3: B3Point, nodes_tbl, etas=(0.1, 0.3)):
    """Per-node flags.  Operational only: a mode is 'in window' if |omega|/M<=eta and
    x<=eta (derivative-expansion parameters).  Domain flag from the recovered protocol
    cuts is a property of the background (x independent).  UV completeness is never
    established (g ~ M)."""
    dom, X = eft_domain_flags(b3)
    rows = []
    for tb in nodes_tbl:
        x = tb[0]["omega"]/tb[0]["c"]
        r = dict(x=x, background_in_protocol_domain=dom)
        for eta in etas:
            ok = [abs(b["omega"]) <= eta and x <= eta for b in tb]
            r[f"modes_in_window_eta{eta}"] = int(sum(ok))
            r[f"all_in_window_eta{eta}"] = bool(all(ok))
            r[f"per_branch_in_window_eta{eta}"] = ok
        r["status"] = "unestablished"
        rows.append(r)
    return rows


# --- composition exchange summary ---------------------------------------------------
def crossings(xs, f, level=0.5):
    """Log-x linear-interpolated crossings of f through `level` (list of (x, direction))."""
    out = []
    lx = np.log(xs)
    for i in range(len(xs) - 1):
        a, b = f[i] - level, f[i+1] - level
        if a == 0:
            out.append((float(xs[i]), 0))
        elif a*b < 0:
            t = a/(a - b)
            out.append((float(np.exp(lx[i] + t*(lx[i+1] - lx[i]))), 1 if b > a else -1))
    return out


def transition_width(xs, f, lo=0.25, hi=0.75):
    """log10-width between first crossings of the `lo` and `hi` levels (None if absent)."""
    a = crossings(xs, f, lo); b = crossings(xs, f, hi)
    if not a or not b:
        return None
    return float(abs(np.log10(a[0][0]) - np.log10(b[0][0])))


def exchange_summary(tbl, xs):
    """Per branch: Proca-share at both ends, crossings, whether composition is exchanged."""
    res = []
    for br in range(6):
        fe = np.array([shares(r[br])["f_proca"] for r in tbl])
        fa = np.array([shares(r[br])["amp_proca"] for r in tbl])
        K = np.array([r[br]["K"] for r in tbl]); om = np.array([r[br]["omega"] for r in tbl])
        E = np.array([r[br]["E"] for r in tbl])
        ce, ca = crossings(xs, fe), crossings(xs, fa)
        res.append(dict(
            branch=br,
            f_proca_small_x=float(fe[0]), f_proca_large_x=float(fe[-1]),
            amp_proca_small_x=float(fa[0]), amp_proca_large_x=float(fa[-1]),
            exchanged_energy_metric=bool((fe[0] - 0.5)*(fe[-1] - 0.5) < 0),
            exchanged_amplitude_metric=bool((fa[0] - 0.5)*(fa[-1] - 0.5) < 0),
            crossing_x_energy_metric=[c[0] for c in ce],
            crossing_x_amplitude_metric=[c[0] for c in ca],
            width_log10_energy_metric_25_75=transition_width(xs, fe),
            krein_sign_changes=int(np.sum(np.diff(np.sign(K)) != 0)),
            energy_sign_changes=int(np.sum(np.diff(np.sign(E)) != 0)),
            krein_sign=int(np.sign(K[0])), omega_sign=int(np.sign(om[0])),
            min_abs_krein_unit_norm=float(np.abs(K).min()),
            min_energy_over_omega2=float(np.min(E/om**2)),
        ))
    return res
