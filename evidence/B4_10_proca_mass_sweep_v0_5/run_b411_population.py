"""Population scan: how the UV front velocity depends on the EFT amplitude cut (T1.7.2 generator, seed 172026)."""
import json, sys
from pathlib import Path
import numpy as np
import b410._b49  # noqa
from b410 import causality as ca, experiment3 as e3
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/"vendor_t173"))
from robustness_scan import generate_samples  # noqa: E402


def main(out="results_b411", n=10000, cuts=(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.95)):
    res = {}
    for cut in cuts:
        for mode, sg in (("rank1_minus", -1.0), ("rank1_plus", 1.0)):
            fr = []; real = 0; sub = 0
            for bg in generate_samples(n, (0.05, 1.0), (0.05, 1.0), (0.0, 1.0), cut, 172026):
                k = sg*np.sqrt(bg["p"]*bg["s"])
                b3 = e3.make_point(dict(p=bg["p"], s=bg["s"], kappa=float(k), Lambda=bg["Lambda"], w=float(bg["w"]), q=float(bg["q"]), v=float(bg["v"]), r=float(bg["r"])))
                allr, ok = ca.principal_roots(b3)
                if ok:
                    real += 1; f = float(np.abs(allr.real).max()); fr.append(f); sub += f <= 1.0
            res[f"cut{cut}/{mode}"] = dict(n=n, principal_real=real, subluminal=int(sub),
                                           front_median=float(np.median(fr)) if fr else None,
                                           front_p05=float(np.percentile(fr, 5)) if fr else None,
                                           front_min=float(min(fr)) if fr else None)
    Path(out).mkdir(exist_ok=True, parents=True)
    Path(out, "b411_population_front_velocity.json").write_text(json.dumps(res, indent=1))
    for k, v in res.items(): print(k, v)


if __name__ == "__main__":
    main(*(sys.argv[1:2] or ["results_b411"]))
