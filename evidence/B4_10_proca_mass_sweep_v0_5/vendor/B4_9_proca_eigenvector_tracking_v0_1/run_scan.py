import csv, json
from pathlib import Path
import numpy as np
from proca_dispersion import (Background, scan_spectrum, avoided_crossing_diagnostic,
    track_eigenmodes, summarize_eigenmode_scan, local_roots)

OUT=Path("results")
OUT.mkdir(exist_ok=True)

base=dict(g_phi=np.sqrt(0.05),g_psi=0.0,w=1.0,q=1.0,v=0.0,r=0.0,Lambda=0.0)
cases={
 "B3_illustrative": Background(M=1.0,**base),
 "B1_boundary": Background(M=1.0,g_phi=1.0,g_psi=1.0,
     w=np.sqrt(0.5),q=0.0,v=np.sqrt(0.5),r=0.0,Lambda=0.0),
 "complex_frequency_stress_outside_B1_domain": Background(M=1.0,
     g_phi=0.8108811308204898,g_psi=1.3198751273700329,
     w=-0.5549437639433696,q=0.39273626882885226,
     v=-1.7629934306179855,r=-0.4494727955570852,Lambda=0.16151817312910333),
}
summaries={}
for name,bg in cases.items():
    ratios=np.geomspace(1e-3,2.0,100) if name.startswith("complex_") else np.geomspace(1e-3,10.0,240)
    ks=ratios*bg.M
    tracks,summary=scan_spectrum(bg,ks)
    if name=="B3_illustrative":
        summary["candidate_avoided_crossing"]=avoided_crossing_diagnostic(bg,2.7,3.6,901)
        # Eigenvector-overlap tracking through the avoided crossing.
        k_local=np.linspace(2.7,3.6,181)*bg.M
        modes=track_eigenmodes(bg,k_local)
        eig_summary=summarize_eigenmode_scan(modes)
        upper=np.argsort(modes["omega"][0].real)[-2:]
        eig_summary["upper_pair_scalar_fraction_start"]=modes["scalar_fraction"][0,upper].tolist()
        eig_summary["upper_pair_scalar_fraction_end"]=modes["scalar_fraction"][-1,upper].tolist()
        summary["eigenvector_tracking_across_avoided_crossing"]=eig_summary
        with (OUT/"B3_avoided_crossing_eigenvectors.csv").open("w",newline="",encoding="utf-8") as f:
            w=csv.writer(f); w.writerow(["k_over_M","mode_id","omega_real","omega_imag","scalar_fraction","proca_fraction","null_residual","mean_overlap_from_previous"])
            for i,k in enumerate(k_local):
                for j,om in enumerate(modes["omega"][i]):
                    w.writerow([k/bg.M,j,float(om.real),float(om.imag),
                        float(modes["scalar_fraction"][i,j]),float(modes["proca_fraction"][i,j]),
                        float(modes["null_residual"][i,j]),float(modes["mean_overlap"][i])])
    summaries[name]=summary
    with (OUT/f"{name}_tracked_roots.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["k_over_M","mode_id","omega_real","omega_imag","phase_speed_real","phase_speed_imag"])
        for i,k in enumerate(ks):
            for j,om in enumerate(tracks[i]):
                w.writerow([k/bg.M,j,float(om.real),float(om.imag),
                    float((om/k).real),float((om/k).imag)])

# Hold all non-mediator background quantities fixed; vary only M_A.
mass_sweep=[]
for M in [0.5,0.75,1.0,1.5,2.0,3.0]:
    bg=Background(M=M,**base)
    local=local_roots(bg)
    superluminal=[float(z.real) for z in local if abs(z.imag)<1e-7 and z.real>1.000001]
    cstar=max(superluminal) if superluminal else None
    kpred=M/np.sqrt(cstar*cstar-1.0) if cstar is not None else None
    if kpred is None:
        continue
    kvals=np.linspace(max(0.02,kpred*0.65),kpred*1.35,221)
    modes=track_eigenmodes(bg,kvals)
    gaps=[]
    for i,row in enumerate(modes["omega"]):
        pos=np.sort(row.real[(row.real>0)&(np.abs(row.imag)<1e-8)])
        if len(pos)>=2: gaps.append((float(pos[-1]-pos[-2]),float(kvals[i])))
    gap,kmin=min(gaps,key=lambda z:z[0])
    ss=summarize_eigenmode_scan(modes)
    mass_sweep.append({"M_A":M,"local_superluminal_c":cstar,
        "uncoupled_crossing_k":kpred,"minimum_gap_omega":gap,
        "k_at_minimum_gap":kmin,"k_over_M_at_minimum_gap":kmin/M,
        "all_real_within_tolerance":ss["all_real_within_tolerance"],
        "max_null_residual":ss["max_null_residual"]})
with (OUT/"mediator_mass_sweep.csv").open("w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=list(mass_sweep[0])); w.writeheader(); w.writerows(mass_sweep)
(OUT/"scan_summary.json").write_text(json.dumps({"scans":summaries,"fixed_background_mass_sweep":mass_sweep},indent=2),encoding="utf-8")
print(json.dumps({"fixed_background_mass_sweep":mass_sweep},indent=2))
