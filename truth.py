"""Score grown surfaces against a reference (human-made) mesh: per surface, the share of its vertices near the reference
(<= 4 vx) among those the reference covers (<= 60 vx), and by distance from the surface's own seed.
  truth.py REF_DIR OUT.json SURF_DIR [SURF_DIR ...]"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seedcoh import load, dense, mesh_seed, RINGS
from scipy.spatial import cKDTree
ref, out, dirs = sys.argv[1], sys.argv[2], sys.argv[3:]
T = cKDTree(dense(*load(ref)))
res = {}
for d in dirs:
    P, v = load(d); p = P[v]; dist, _ = T.query(p, distance_upper_bound=60.0); cov = np.isfinite(dist)
    s = mesh_seed(d); r = np.linalg.norm(p - s, axis=1)
    rings = {}
    for lo, hi in zip(RINGS[:-1], RINGS[1:]):
        m = cov & (r >= lo) & (r < hi)
        if m.sum() >= 20: rings[f"{int(lo)}-{'inf' if hi > 1e8 else int(hi)}"] = round(float(np.mean(dist[m] <= 4)), 3)
    res[d] = dict(on_ref_le4=round(float(np.mean(dist[cov] <= 4)), 3) if cov.any() else None,
                  one_wrap_off=round(float(np.mean((dist[cov] > 11) & (dist[cov] <= 19))), 3) if cov.any() else None,
                  covered=round(float(cov.mean()), 3), seed_on_ref=bool(np.min(np.where(cov, dist, 99)[r < 100], initial=99) <= 4), rings=rings)
    print(d, res[d], flush=True)
json.dump(res, open(out, "w"), indent=1)
