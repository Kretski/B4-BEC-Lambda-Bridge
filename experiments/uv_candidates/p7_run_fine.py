# Driver for the third calibration resolution (32768, 0.0125); uses the frozen functions of p7_free_calibration.py unchanged.
import sys, json, numpy as np
sys.path.insert(0, "/mnt/user-data/outputs")
import p7_free_calibration as P
r, ts, nrm, dr = P.analyse(32768, 0.0125, False)
np.savetxt("p7_protocol_N32768.csv", np.c_[ts, nrm, dr], header="t norm_fluct drift_total", delimiter=",")
json.dump(r, open("p7_protocol_N32768.json", "w"), indent=1); print(json.dumps(r), flush=True)
