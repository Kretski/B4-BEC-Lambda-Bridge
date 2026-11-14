"""T7a step 4b (DRAFT v0.9, unregistered): independent IBP-matrix implementation.

Specification: T7A_STEP4B_IBP_MATRIX_PROPOSAL.md v0.9 (commit 4074199), canonical rule v0.7.1, option A rules (v0.6).
Pure Python with fractions.Fraction: no SymPy, no Euler operator, no import from t7a_step4_classify.py, no floating point.
The step 4 result is read ONLY for the comparison K7.
Status: DRAFT - unregistered, not a T7a result. No numerical verification, no M^-6, no B4.12 conclusions.

Usage: python -I t7a_step4b_ibp_matrix.py <step4_json> <output_dir>
"""
import hashlib
import itertools
import json
import platform
import sys
from fractions import Fraction

SPEC_VERSION = "v0.9"
STEP4_JSON_SHA256 = "aaf07cc060d5508101bf66b711a1a5ba184b34bd0ddcb548a32c4f6e07f57269"  # file as committed in 98780c9
FIELDS = ("f", "g")
CONTENT = {"f4": ("f", "f", "f", "f"), "f2g2": ("f", "f", "g", "g"), "g4": ("g", "g", "g", "g")}
SECTORS = ("f4", "f2g2", "g4")


# ---------------------------------------------------------------- section 1: combinatorial monomials
# A factor is a triple (field, a, b) = d_t^a d_x^b field, a, b >= 0. A monomial is a sorted tuple of factors
# (repeats kept); a polynomial is a dict {monomial: Fraction}.

def add_to(poly, mono, coeff):
    v = poly.get(mono, Fraction(0)) + coeff
    if v == 0:
        poly.pop(mono, None)
    else:
        poly[mono] = v


def D(poly, direction):
    """Total derivative on combinatorial monomials: sum over positions, raise a (t) or b (x) by one, resort."""
    out = {}
    for mono, c in poly.items():
        for i, (p, a, b) in enumerate(mono):
            new = list(mono)
            new[i] = (p, a + 1, b) if direction == "t" else (p, a, b + 1)
            add_to(out, tuple(sorted(new)), c)
    return out


def mul(P, Q):
    out = {}
    for m1, c1 in P.items():
        for m2, c2 in Q.items():
            add_to(out, tuple(sorted(m1 + m2)), c1 * c2)
    return out


def mono(*factors):
    return {tuple(sorted(factors)): Fraction(1)}


def label(m):
    return " ".join(f"{p}_{'t' * a}{'x' * b}" if a + b else p for p, a, b in m)


def key(m):  # v0.7.1: k1, k2, k3
    orders = [a + b for _, a, b in m]
    return (max(orders), sum(o * o for o in orders), tuple(sorted(m)))


# ---------------------------------------------------------------- section 2: spaces (two independent enumerations)
def factor_list(nt, nx, underived):
    return sorted((p, a, b) for p in FIELDS for a in range(nt + 1) for b in range(nx + 1) if underived or a + b >= 1)


def enum_cwr(nt, nx, sector, underived):
    want = tuple(sorted(CONTENT[sector]))
    out = set()
    for combo in itertools.combinations_with_replacement(factor_list(nt, nx, underived), 4):
        if sum(t[1] for t in combo) == nt and sum(t[2] for t in combo) == nx and tuple(sorted(t[0] for t in combo)) == want:
            out.add(tuple(combo))
    return out


def enum_product(nt, nx, sector, underived):
    want = tuple(sorted(CONTENT[sector]))
    out = set()
    for combo in itertools.product(factor_list(nt, nx, underived), repeat=4):
        if sum(t[1] for t in combo) == nt and sum(t[2] for t in combo) == nx and tuple(sorted(t[0] for t in combo)) == want:
            out.add(tuple(sorted(combo)))
    return out


def has_underived(m):
    return any(a == 0 and b == 0 for _, a, b in m)


# ---------------------------------------------------------------- exact linear algebra over Q
class Echelon:
    """Fully reduced echelon form over Q with sparse rows; each row carries a tag vector (dict) for coordinates.
    Pivot = smallest column index of the row in the given column order."""

    def __init__(self, col_index):
        self.ci = col_index
        self.rows = {}  # pivot column index -> (row dict colidx->Fraction, tag dict)

    def _vec(self, poly):
        return {self.ci[m]: c for m, c in poly.items()}

    def reduce(self, vec, tag=None):
        vec = dict(vec)
        tag = dict(tag or {})
        for pc in sorted(self.rows):
            if pc in vec:
                f = vec[pc]
                row, rtag = self.rows[pc]
                for c, v in row.items():
                    add_to(vec, c, -f * v)
                for c, v in rtag.items():
                    add_to(tag, c, -f * v)
        return vec, tag

    def insert(self, vec, tag=None):
        vec, tag = self.reduce(vec, tag)
        if not vec:
            return False
        pc = min(vec)
        inv = 1 / vec[pc]
        vec = {c: v * inv for c, v in vec.items()}
        tag = {c: v * inv for c, v in tag.items()}
        for k in list(self.rows):
            row, rtag = self.rows[k]
            if pc in row:
                f = row[pc]
                row, rtag = dict(row), dict(rtag)
                for c, v in vec.items():
                    add_to(row, c, -f * v)
                for c, v in tag.items():
                    add_to(rtag, c, -f * v)
                self.rows[k] = (row, rtag)
        self.rows[pc] = (vec, tag)
        return True


def frac_str(q):
    return f"{q.numerator}/{q.denominator}" if q.denominator != 1 else f"{q.numerator}"


# ---------------------------------------------------------------- per-sector computation
def build_WS(sector, VS, OS, preimages_t, preimages_x):
    """Section 4: matrix T (rows O_S then V_S, columns relations); RREF of T^T with O_S columns first."""
    US = OS + VS  # fixed order: O_S block first, then V_S block
    ci = {m: i for i, m in enumerate(US)}
    nO = len(OS)
    relations = [D({Q: Fraction(1)}, "t") for Q in preimages_t] + [D({Q: Fraction(1)}, "x") for Q in preimages_x]
    for r in relations:
        for m in r:
            if m not in ci:
                raise RuntimeError(f"relation monomial outside U_S: {m}")
    raw_in_VS = sum(1 for r in relations if r and all(ci[m] >= nO for m in r))
    raw_zero = sum(1 for r in relations if not r)
    ech = Echelon(ci)
    for r in relations:
        ech.insert({ci[m]: c for m, c in r.items()})
    WS = [row for pc, (row, _) in sorted(ech.rows.items()) if pc >= nO]  # pivot in V_S block
    k5 = all(all(c >= nO for c in row) for row in WS)  # zero O_S part
    # convert W_S basis to V_S-only polynomials
    WS_polys = [{US[c]: v for c, v in row.items()} for row in WS]
    return dict(relations=relations, n_relations=len(relations), raw_in_VS=raw_in_VS, raw_zero=raw_zero,
                raw_leave_VS=len(relations) - raw_in_VS - raw_zero, WS=WS_polys, k5=k5, rank_T=len(ech.rows))


def classify(sector, VS, WS_polys):
    """Section 5: greedy v0.7.1 basis modulo W_S, coordinates via tags."""
    order = sorted(VS, key=key)
    ci = {m: i for i, m in enumerate(order)}
    ech = Echelon(ci)
    for w in WS_polys:
        if not ech.insert({ci[m]: c for m, c in w.items()}):
            raise RuntimeError("W_S basis not independent")
    B = []
    for m in order:
        if ech.insert({ci[m]: Fraction(1)}, {len(B): Fraction(1)}):
            B.append(m)
    WS_only = Echelon(ci)
    for w in WS_polys:
        WS_only.insert({ci[m]: c for m, c in w.items()})
    return dict(order=order, ci=ci, ech=ech, B=B, WS_only=WS_only)


def coords(cls, poly):
    """X = sum c_k b_k + w, w in W_S. Returns c (list) or raises if poly not in span."""
    ci = cls["ci"]
    for m in poly:
        if m not in ci:
            raise RuntimeError(f"monomial outside V_S: {m}")
    vec, tag = cls["ech"].reduce({ci[m]: c for m, c in poly.items()})
    if vec:
        raise RuntimeError("not in span(W_S + B_S)")
    # reduction subtracted f*rows; tag accumulated -(sum f*tag_row); coordinates are the negative
    return [-tag.get(k, Fraction(0)) for k in range(len(cls["B"]))]


def residual_zero(cls, poly, c):
    """Independent re-reduction of X - sum c_k b_k modulo W_S only (no B rows)."""
    R = dict(poly)
    for ck, m in zip(c, cls["B"]):
        add_to(R, m, -ck)
    vec, _ = cls["WS_only"].reduce({cls["ci"][m]: v for m, v in R.items()})
    return not vec


def X_sector(sector):
    A = D(mono(("f", 1, 0), ("f", 0, 1)), "t")
    Bp = D(mono(("g", 1, 0), ("g", 0, 1)), "t")
    return {"f4": mul(A, A), "f2g2": mul(A, Bp), "g4": mul(Bp, Bp)}[sector]


# ---------------------------------------------------------------- main
def main():
    step4_path, outdir = sys.argv[1], sys.argv[2]
    with open(__file__, "rb") as fh:
        script_sha = hashlib.sha256(fh.read()).hexdigest()
    with open(step4_path, "rb") as fh:
        raw = fh.read()
    step4_sha = hashlib.sha256(raw).hexdigest()
    step4 = json.loads(raw.decode("utf-8"))

    ok_all, failures, controls = True, [], {}

    def ctrl(name, passed, desc, extra=None):
        nonlocal ok_all
        controls[name] = {"result": "PASS" if passed else "FAIL", "description": desc}
        if extra:
            controls[name].update(extra)
        if not passed:
            ok_all = False
            failures.append(name)

    res, k6a_ok, k6a_detail = {}, True, {}
    for S in SECTORS:
        # two independent enumerations (K6a) of U_S and preimages
        U1, U2 = enum_cwr(4, 2, S, True), enum_product(4, 2, S, True)
        Pt1, Pt2 = enum_cwr(3, 2, S, True), enum_product(3, 2, S, True)
        Px1, Px2 = enum_cwr(4, 1, S, True), enum_product(4, 1, S, True)
        V1 = enum_cwr(4, 2, S, False)
        same = (U1 == U2) and (Pt1 == Pt2) and (Px1 == Px2) and (V1 == {m for m in U1 if not has_underived(m)})
        n_und = sum(1 for Q in Pt1 if has_underived(Q)) + sum(1 for Q in Px1 if has_underived(Q))
        k6a_ok = k6a_ok and same and n_und > 0
        k6a_detail[S] = {"enumerations_agree": same, "n_preimages_with_underived": n_und}

        VS = sorted((m for m in U1 if not has_underived(m)), key=key)
        OS = sorted((m for m in U1 if has_underived(m)), key=key)
        Pt, Px = sorted(Pt1, key=key), sorted(Px1, key=key)
        W = build_WS(S, VS, OS, Pt, Px)
        cls = classify(S, VS, W["WS"])
        dim = len(VS) - len(W["WS"])
        X = X_sector(S)
        c = coords(cls, X)
        rz = residual_zero(cls, X, c)
        # K6b ablation: preimages without underived factors only (report only)
        Wab = build_WS(S, VS, OS, [Q for Q in Pt if not has_underived(Q)], [Q for Q in Px if not has_underived(Q)])
        res[S] = dict(VS=VS, OS=OS, Pt=Pt, Px=Px, W=W, cls=cls, dim=dim, X=X, c=c, rz=rz,
                      ablation_dim=len(VS) - len(Wab["WS"]), ablation_dimWS=len(Wab["WS"]))
        if not rz or len(cls["B"]) != dim:
            ok_all = False
            failures.append(f"sector {S}: residual_zero={rz}, |B_S|={len(cls['B'])}, dim={dim}")

    # K1-K3 total derivatives
    Y1 = D(mono(("f", 1, 0), ("f", 0, 1), ("f", 1, 0), ("f", 1, 1)), "t")
    Y2 = D(D(mono(("f", 1, 0), ("f", 1, 0), ("f", 1, 0), ("f", 1, 0)), "x"), "x")
    Y2_expected = {tuple(sorted([("f", 1, 0), ("f", 1, 0), ("f", 1, 1), ("f", 1, 1)])): Fraction(12),
                   tuple(sorted([("f", 1, 0), ("f", 1, 0), ("f", 1, 0), ("f", 1, 2)])): Fraction(4)}
    Y3 = D(mono(("f", 1, 0), ("f", 0, 1), ("g", 1, 0), ("g", 1, 1)), "t")

    def zero_mod_W(S, Y):
        cl = res[S]["cls"]
        vec, _ = cl["WS_only"].reduce({cl["ci"][m]: v for m, v in Y.items()})
        return not vec and all(v == 0 for v in coords(cl, Y))

    ctrl("K1", zero_mod_W("f4", Y1), "D_t(f_t f_x f_t f_tx) reduces to 0 mod W_f4, c = 0")
    ctrl("K2", Y2 == Y2_expected and zero_mod_W("f4", Y2), "D_x^2(f_t^4) = 12 f_t^2 f_tx^2 + 4 f_t^3 f_txx; reduces to 0 mod W_f4")
    ctrl("K3", zero_mod_W("f2g2", Y3), "D_t(f_t f_x g_t g_tx) reduces to 0 mod W_f2g2, c = 0")
    k4 = True
    for S in SECTORS:
        cl = res[S]["cls"]
        for k, m in enumerate(cl["B"]):
            ck = coords(cl, {m: Fraction(1)})
            k4 = k4 and all((v == 1) if i == k else (v == 0) for i, v in enumerate(ck))
    ctrl("K4", k4, "every basis monomial b_k gives the unit vector e_k (all sectors)")
    ctrl("K5", all(res[S]["W"]["k5"] for S in SECTORS), "every W_S basis vector has zero O_S part; raw relation counts reported per sector",
         {"raw_relations": {S: {"total": res[S]["W"]["n_relations"], "zero": res[S]["W"]["raw_zero"],
                                "wholly_in_VS": res[S]["W"]["raw_in_VS"], "leave_VS": res[S]["W"]["raw_leave_VS"]} for S in SECTORS}})
    ctrl("K6a", k6a_ok, "two independent enumerations of U_S, P_t, P_x agree; preimages with underived factors present (>0)",
         {"per_sector": k6a_detail})
    controls["K6b_ablation_report_only"] = {"description": "dim(V_S/~) using only preimages without underived factors; report only",
                                            "per_sector": {S: {"ablation_dim": res[S]["ablation_dim"], "full_dim": res[S]["dim"]} for S in SECTORS}}
    # K7 comparison with step 4 (98780c9)
    sha_ok = step4_sha == STEP4_JSON_SHA256
    k7_detail, k7 = {}, sha_ok
    for S in SECTORS:
        s4 = step4["sectors"][S]
        dims = (s4["dim"] == res[S]["dim"])
        basis_same = [tuple(tuple(t) for t in b["triples"]) for b in s4["basis"]] == list(res[S]["cls"]["B"])
        c4 = [Fraction(x) for x in s4["coords"]]
        coords_same = c4 == res[S]["c"]
        k7 = k7 and dims and basis_same and coords_same
        k7_detail[S] = {"dim_step4": s4["dim"], "dim_ibp": res[S]["dim"], "dim_equal": dims, "basis_identical": basis_same,
                        "coords_identical": coords_same}
    ctrl("K7", k7, "step 4 JSON sha256 verified; dims, bases (same monomials, same order) and exact coords identical, all sectors",
         {"step4_json_sha256_ok": sha_ok, "per_sector": k7_detail})

    out = {"meta": {"spec_version": SPEC_VERSION, "python": platform.python_version(), "script_sha256": script_sha,
                    "step4_json_sha256": step4_sha, "status": "DRAFT - unregistered, not a T7a result",
                    "spec": "T7A_STEP4B_IBP_MATRIX_PROPOSAL.md v0.9 (4074199); canonical rule v0.7.1; option A (v0.6)",
                    "independence": "pure Python Fraction; no SymPy, no Euler operator, no import from step 4",
                    "scope_exclusions": "no registration, no numerical verification, no M^-6, no B4.12 conclusions"},
           "sectors": {}, "controls": controls}
    for S in SECTORS:
        r = res[S]
        out["sectors"][S] = {
            "n_VS": len(r["VS"]), "n_OS": len(r["OS"]),
            "n_preimages_t": len(r["Pt"]), "n_preimages_t_with_underived": sum(1 for Q in r["Pt"] if has_underived(Q)),
            "n_preimages_x": len(r["Px"]), "n_preimages_x_with_underived": sum(1 for Q in r["Px"] if has_underived(Q)),
            "n_relations": r["W"]["n_relations"], "rank_T": r["W"]["rank_T"], "dim_WS": len(r["W"]["WS"]), "dim": r["dim"],
            "basis": [{"index": k, "triples": [list(t) for t in m], "label": label(m)} for k, m in enumerate(r["cls"]["B"])],
            "coords": [frac_str(v) for v in r["c"]], "residual_zero": r["rz"], "ablation_dim": r["ablation_dim"]}
    out["verdict"] = "ALL_CHECKS_PASS" if ok_all else "CHECK_FAILED: " + ", ".join(failures)

    with open(f"{outdir}/t7a_step4b_ibp.json", "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=1, sort_keys=True, ensure_ascii=False)
        fh.write("\n")
    L = [f"T7a step 4b independent IBP-matrix implementation - {out['meta']['status']}",
         f"spec {SPEC_VERSION}; python {out['meta']['python']}; pure Python Fraction (no SymPy, no Euler operator)",
         f"script sha256 {script_sha}", f"step4 json sha256 {step4_sha} (expected {STEP4_JSON_SHA256})", ""]
    for S in SECTORS:
        s = out["sectors"][S]
        L += [f"[sector {S}] |V_S| = {s['n_VS']}, |O_S| = {s['n_OS']}",
              f"  preimages: D_t on (3,2): {s['n_preimages_t']} ({s['n_preimages_t_with_underived']} with underived); "
              f"D_x on (4,1): {s['n_preimages_x']} ({s['n_preimages_x_with_underived']} with underived)",
              f"  relations: {s['n_relations']}, rank T = {s['rank_T']}, dim W_S = {s['dim_WS']}, dim(V_S/~) = {s['dim']}, "
              f"ablation dim = {s['ablation_dim']}",
              "  basis (greedy order, modulo W_S):"]
        L += [f"    b{b['index']}: {b['label']}" for b in s["basis"]]
        L += [f"  coords c = [{', '.join(s['coords'])}], residual_zero = {s['residual_zero']}", ""]
    L.append("controls:")
    for k in sorted(controls):
        if "result" in controls[k]:
            L.append(f"  {k}: {controls[k]['result']} - {controls[k]['description']}")
    L.append("  K6b (report only): " + ", ".join(f"{S}: ablation {res[S]['ablation_dim']} vs full {res[S]['dim']}" for S in SECTORS))
    L.append("  K5 raw relations: " + ", ".join(
        f"{S}: {res[S]['W']['raw_in_VS']} wholly in V_S, {res[S]['W']['raw_leave_VS']} leave V_S, {res[S]['W']['raw_zero']} zero" for S in SECTORS))
    L += ["", out["verdict"]]
    with open(f"{outdir}/t7a_step4b_output.txt", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
