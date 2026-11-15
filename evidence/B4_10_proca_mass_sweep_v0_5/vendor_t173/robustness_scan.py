"""Reproducible randomized EFT robustness scan for the rank-one mediator model."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np

def matrices(w,q,v,r,p,s,k,Lam):
    A=np.array([[1+p*q*q, -k*q*r],[-k*q*r,1+s*r*r]],float)
    B=np.array([[-k*v*r+2*p*q*w, -k*q*v],
                [-k*w*r, -k*q*w+2*s*r*v]],float)
    C=np.array([[-1-3*Lam*q*q+p*w*w, -k*w*v],
                [-k*w*v, -1+s*v*v]],float)
    return A,B,C

def characteristic_coeffs(A,B,C):
    D=-(B+B.T)
    # ascending powers of c for each matrix entry
    def entry(i,j): return np.array([C[i,j],D[i,j],A[i,j]],float)
    m00,m11,m01=entry(0,0),entry(1,1),entry(0,1)
    poly=np.polynomial.polynomial.polymul(m00,m11)-np.polynomial.polynomial.polymul(m01,m01)
    # return descending order, padded if leading cancellation occurs
    poly=np.pad(poly,(0,max(0,5-len(poly))))
    return poly[:5][::-1]

def analyze_state(w,q,v,r,p,s,k,Lam,real_tol=1e-6):
    A,B,C=matrices(w,q,v,r,p,s,k,Lam)
    eigA=np.linalg.eigvalsh(A)
    coeff=characteristic_coeffs(A,B,C)
    roots=np.roots(coeff)
    max_imag=float(np.max(np.abs(roots.imag)))
    # A double real root is returned by np.roots with an imaginary error of
    # about sqrt(machine eps) ~ 1e-8, so the real-root tolerance must be larger
    # than that (default 1e-6); otherwise exactly degenerate backgrounds
    # (e.g. p=s with equal amplitudes) are misclassified as non-hyperbolic.
    hyper=bool(max_imag <= real_tol)
    sep=float(min(abs(roots[i]-roots[j]) for i in range(4) for j in range(i+1,4)))
    # Companion matrix eigenvector conditioning is a diagnostic, not a proof.
    try:
        # numpy eig returns right eigenvectors
        _, V=np.linalg.eig(np.polynomial.polynomial.polycompanion(coeff[::-1]))
        cond=float(np.linalg.cond(V))
    except Exception:
        cond=float("inf")
    return {
        "kinetic_positive": bool(eigA[0] > 0),
        "min_eig_A": float(eigA[0]),
        "hyperbolic_at_sample": hyper,
        "max_abs_imag_c": max_imag,
        "roots_real": [[float(z.real),float(z.imag)] for z in roots],
        "eigenvector_condition": cond,
        "min_root_separation": sep,
        "near_degenerate_roots": bool(sep <= 1e-3*max(1.0,float(np.max(np.abs(roots))))),
        "coefficients_desc": [float(x) for x in coeff],
    }

def generate_samples(n,p_range,s_range,Lam_range,cut,seed):
    rng=np.random.default_rng(seed)
    for i in range(n):
        # Positive mediator self-couplings, sampled log-uniformly.
        p=float(np.exp(rng.uniform(np.log(p_range[0]),np.log(p_range[1]))))
        s=float(np.exp(rng.uniform(np.log(s_range[0]),np.log(s_range[1]))))
        Lam=float(rng.uniform(*Lam_range))
        # Independently sample all four dimensionless EFT amplitudes.
        Xq,Xr,Xw,Xv=rng.uniform(0,cut,4)
        signs=rng.choice([-1.0,1.0],4)
        q=signs[0]*np.sqrt(Xq/p); r=signs[1]*np.sqrt(Xr/s)
        w=signs[2]*np.sqrt(Xw/p); v=signs[3]*np.sqrt(Xv/s)
        yield dict(index=i,p=p,s=s,Lambda=Lam,q=q,r=r,w=w,v=v,
                   Xq=Xq,Xr=Xr,Xw=Xw,Xv=Xv)

def run_scan(n=10000,cut=0.5,seed=172026,p_range=(0.05,1.0),
             s_range=(0.05,1.0),Lam_range=(0.0,1.0),real_tol=1e-6,
             k_modes=("rank1_plus","rank1_minus","decoupled")):
    rows=[]
    for bg in generate_samples(n,p_range,s_range,Lam_range,cut,seed):
        for mode in k_modes:
            p,s=bg["p"],bg["s"]
            if mode=="rank1_plus": k=np.sqrt(p*s)
            elif mode=="rank1_minus": k=-np.sqrt(p*s)
            elif mode=="decoupled": k=0.0
            else: raise ValueError(mode)
            result=analyze_state(bg["w"],bg["q"],bg["v"],bg["r"],p,s,k,bg["Lambda"],real_tol)
            row={**bg,"mode":mode,"kappa":float(k),**result}
            row["eft_control"]=bool(max(bg["Xq"],bg["Xr"],bg["Xw"],bg["Xv"]) < cut)
            rows.append(row)
    return rows

def _with_ext(path, ext):
    """Append an extension without touching dots already in the name.

    Path.with_suffix would turn 'cut_0.7' into 'cut_0.json' and make runs
    with different cuts overwrite each other.
    """
    return Path(str(path) + ext)


def sign_table(rows):
    """Failure counts for rank-one rows by (sign k, sigma=sign(w q v r))."""
    table = {}
    for r in rows:
        if r["mode"] == "decoupled":
            continue
        sg = lambda x: 1 if x > 0 else -1
        key = (sg(r["kappa"]), sg(r["w"]*r["q"]*r["v"]*r["r"]))
        cell = table.setdefault(key, [0, 0])
        cell[0] += 1
        cell[1] += (not r["hyperbolic_at_sample"])
    return {f"k{'+' if k > 0 else '-'},sigma{'+' if s > 0 else '-'}":
            {"samples": n, "nonhyperbolic": b, "fraction": b/n}
            for (k, s), (n, b) in sorted(table.items())}


def summarize(rows):
    out={}
    for mode in sorted({r["mode"] for r in rows}):
        rr=[r for r in rows if r["mode"]==mode]
        bad=[r for r in rr if not r["hyperbolic_at_sample"]]
        out[mode]={
            "samples":len(rr),
            "kinetic_positive":sum(r["kinetic_positive"] for r in rr),
            "hyperbolic":len(rr)-len(bad),
            "nonhyperbolic":len(bad),
            "nonhyperbolic_fraction":len(bad)/len(rr) if rr else None,
            "max_imaginary_part":max((r["max_abs_imag_c"] for r in rr),default=0.0),
            "max_eigenvector_condition":max((r["eigenvector_condition"] for r in rr),default=0.0),
        }
    return out

def main(argv=None):
    ap=argparse.ArgumentParser(description="T1.7.2 randomized rank-one EFT robustness scan")
    ap.add_argument("--samples",type=int,default=10000)
    ap.add_argument("--cut",type=float,default=0.5)
    ap.add_argument("--seed",type=int,default=172026)
    ap.add_argument("--p-min",type=float,default=0.05); ap.add_argument("--p-max",type=float,default=1.0)
    ap.add_argument("--s-min",type=float,default=0.05); ap.add_argument("--s-max",type=float,default=1.0)
    ap.add_argument("--lambda-min",type=float,default=0.0); ap.add_argument("--lambda-max",type=float,default=1.0)
    ap.add_argument("--real-tol",type=float,default=1e-6)
    ap.add_argument("--out",default="results/t172_robustness")
    ap.add_argument("--no-rows-json",action="store_true",help="write only the summary to JSON (rows stay in the CSV)")
    args=ap.parse_args(argv)
    if args.samples<1 or not (0<args.cut<1): ap.error("samples must be positive and cut in (0,1)")
    if min(args.p_min,args.s_min)<=0 or args.p_max<args.p_min or args.s_max<args.s_min:
        ap.error("p and s ranges must be positive and ordered")
    rows=run_scan(args.samples,args.cut,args.seed,(args.p_min,args.p_max),
                  (args.s_min,args.s_max),(args.lambda_min,args.lambda_max),args.real_tol)
    summary={"config":vars(args),"summary":summarize(rows)}
    out=Path(args.out); out.parent.mkdir(parents=True,exist_ok=True)
    summary["sign_table"]=sign_table(rows)
    payload={"summary":summary} if args.no_rows_json else {"summary":summary,"rows":rows}
    _with_ext(out,".json").write_text(json.dumps(payload,indent=2))
    fields=[k for k in rows[0] if k!="roots_real" and k!="coefficients_desc"]
    with _with_ext(out,".csv").open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=fields); wr.writeheader()
        for row in rows: wr.writerow({k:row[k] for k in fields})
    print(json.dumps(summary,indent=2))
    print(f"Wrote {_with_ext(out,'.json')} and {_with_ext(out,'.csv')}")

if __name__=="__main__": main()
