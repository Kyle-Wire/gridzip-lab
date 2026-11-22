"""Post-hoc checks, not in the protocol (run after reading aij_c.json): LBM D/32 on its own, RANS D/20."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import aij_c as A
from codelen import lawsearch as L

v = Path(__file__).resolve().parents[2] / "WindSim" / "docs" / "validation"
out = {}
for lab, f, col in (("lbm_d32_z0", v / "lbm_case_c_1D_resolution/d32_z0tunnel/points_1D.csv", "vs_pred"),
                    ("lbm_d10", v / "lbm_case_c_d10_blind/points_1D.csv", "vs_pred"),
                    ("rans_d20", v / "rans_case_c_d20/results_points.csv", "pred_vector"),
                    ("rans_d10_1D", v / "rans_case_c_d10/points_1D.csv", "pred_vector")):
    ds = A.dataset(f.parent, {"1D": f.name}, col)
    r = np.log(ds["m"] / ds["p"]); lib, ex, cls = A.features(ds["x"], ds["y"], ds["hc"], ds["p"], ds["pk"])
    lib = {k: val for k, val in lib.items() if "hc" not in k}        # one layout: hc is the intercept
    top, n = L.search(lib, r, A.G, 3, exclusive=ex, keep=5); b0 = L.score(lib, [], r, A.G)
    st = A.spatial_test(ds, r, 199)
    lp_only = L.score(lib, ["lp"], r, A.G)
    out[lab] = {"n": len(r), "bias_only_bpr": b0["bits_per_reading"], "bias": b0["fit_coef"]["c"], "sd_after_bias": b0["sd"],
                "top": [{"terms": t["terms"], "bpr": t["bits_per_reading"], "fit_coef": t["fit_coef"], "fit_alpha": t["fit_alpha"],
                         "sd": t["sd"]} for t in top],
                "lp_only": {"bpr": lp_only["bits_per_reading"], "fit_coef": lp_only["fit_coef"]},
                "spatial": {k: val for k, val in st.items() if k != "null"}}
    print(lab, json.dumps(out[lab], default=float)[:900], flush=True)
Path(__file__).resolve().parents[1].joinpath("results/physics/aij_c_posthoc.json").write_text(json.dumps(out, indent=1, default=float))
