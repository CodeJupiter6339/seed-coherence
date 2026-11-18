#!/usr/bin/env python3
"""seed-coherence: label-free measure of sheet switches in auto-grown tifxyz surfaces.

Grow several surfaces from seeds placed on the same sheet, then measure how far each
surface's vertices sit from the other surfaces. Where they agree, they followed the same
sheet; a cluster of distances near one sheet spacing means one of them switched wraps.

  seedcoh.py seeds  MESH [--offset 3]            -> 4 seed points on MESH around its meta.json seed
  seedcoh.py score  DIR DIR [DIR ...] [--json F] -> pairwise agreement, distance histogram, agreement vs. distance from seeds
  seedcoh.py select DIR DIR [DIR ...] [--json F] -> medoid surface (best median agreement, one-wrap mass penalised), or "none"
  seedcoh.py trim   DIR OTHER [OTHER ...] --out M -> mask (M.npy + M.json) of the seed-connected part supported by >= 2/3 of the others
"""
import argparse, itertools, json, os, sys
import numpy as np, tifffile
from scipy.ndimage import zoom, label, binary_propagation
from scipy.spatial import cKDTree

BINS = [0, 2, 4, 8, 11, 19, 26, 34, 41, 60]
RINGS = [0, 100, 200, 400, 800, 1600, 1e9]
NEAR = "0-100"


def load(d):
    P = np.stack([tifffile.imread(os.path.join(d, c + ".tif")).astype(np.float32) for c in "xyz"], -1)
    v = (P[..., 0] > 0) & (P[..., 2] > 0)
    return P, v


def dense(P, v, f=8):
    """Vertices plus bilinear points inside each quad, so distances are point-to-surface, not point-to-vertex."""
    Pz = np.stack([zoom(np.where(v, P[..., k], np.nan), f, order=1) for k in range(3)], -1)
    vz = zoom(v.astype(np.float32), f, order=0) > 0.5
    return Pz[vz & np.all(np.isfinite(Pz), -1)]


def mesh_seed(d):
    m = json.load(open(os.path.join(d, "meta.json")))
    s = np.asarray(m.get("seed", [0, 0, 0]), float)
    g = os.path.join(d, "generations.tif")
    if not np.any(s) and os.path.exists(g):  # a resumed run may write seed [0,0,0]: use the oldest vertex instead
        P, v = load(d); G = tifffile.imread(g).astype(float); G[~v] = np.inf
        s = P[np.unravel_index(np.argmin(G), G.shape)]
    return s if np.any(s) else np.full(3, np.nan)


def cmd_seeds(a):
    P, v = load(a.mesh)
    if a.center:  # fractional grid position, e.g. 0.5,0.5; snapped to the nearest valid vertex
        fi, fj = map(float, a.center.split(","))
        ii, jj = np.nonzero(v); t = np.argmin((ii - fi * v.shape[0]) ** 2 + (jj - fj * v.shape[1]) ** 2); i, j = ii[t], jj[t]
    else:
        s = mesh_seed(a.mesh); D = np.linalg.norm(P - s, axis=-1); D[~v] = np.inf
        i, j = np.unravel_index(np.argmin(D), D.shape)
    k = a.offset
    for n, (di, dj) in zip("ABCD", ((k, 0), (-k, 0), (0, k), (0, -k))):
        ii, jj = i + di, j + dj
        if not (0 <= ii < v.shape[0] and 0 <= jj < v.shape[1] and v[ii, jj]):
            sys.exit(f"seed {n} falls outside the mesh; use a smaller --offset")
        print(n, *map(int, P[ii, jj]))


def cmd_score(a):
    S = {d: load(d) for d in a.dirs}; T = {d: cKDTree(dense(*S[d])) for d in S}
    seeds = {d: mesh_seed(d) for d in S}
    pairs, hist, rings = [], np.zeros(len(BINS)), {}
    for x, y in itertools.permutations(S, 2):
        P, v = S[x]; p = P[v]
        d, _ = T[y].query(p, distance_upper_bound=BINS[-1])
        f = np.isfinite(d)
        h, _ = np.histogram(d[f], bins=BINS); hist[:-1] += h / len(d); hist[-1] += np.mean(~f)
        pairs.append(dict(a=x, b=y, le2=float(np.mean(d <= 2)), le4=float(np.mean(d <= 4)),
                          one_wrap_11_19=float(np.mean((d > 11) & (d <= 19))), none_within_60=float(np.mean(~f))))
        c = (seeds[x] + seeds[y]) / 2; pr = pairs[-1]["rings"] = {}
        if np.all(np.isfinite(c)):
            r = np.linalg.norm(p - c, axis=1)
            for lo, hi in zip(RINGS[:-1], RINGS[1:]):
                m = (r >= lo) & (r < hi)
                if m.sum() >= 20: pr[f"{int(lo)}-{'inf' if hi > 1e8 else int(hi)}"] = float(np.mean(d[m] <= 4))
    # a pair "starts on the same sheet" when they agree near their seeds; only those pairs measure divergence during growth
    same = [p for p in pairs if p["rings"].get(NEAR, 0) >= 0.9]
    for p in same:
        for k, v in p["rings"].items(): rings.setdefault(k, []).append(v)
    n = len(pairs); hist /= n
    out = dict(surfaces=list(S), n_directed_pairs=n,
               mean={k: round(float(np.mean([p[k] for p in pairs])), 3) for k in ("le2", "le4", "one_wrap_11_19", "none_within_60")},
               histogram_vx={f"{BINS[i]}-{BINS[i+1]}": round(float(hist[i]), 3) for i in range(len(BINS) - 1)} | {">60": round(float(hist[-1]), 3)},
               pairs_starting_on_same_sheet=f"{len(same)}/{n}",
               mean_same_start={k: round(float(np.mean([p[k] for p in same])), 3) for k in ("le2", "le4", "one_wrap_11_19", "none_within_60")} if same else None,
               agreement_le4_vs_distance_from_seeds_vx_same_start={k: round(float(np.mean(v)), 3) for k, v in rings.items()},
               pairs=pairs)
    print(json.dumps({k: out[k] for k in out if k != "pairs"}, indent=1))
    if a.json: json.dump(out, open(a.json, "w"), indent=1)

def dists(S, T):
    """D[x][y]: distance of each vertex of surface x to surface y (60 when farther)."""
    return {(x, y): np.minimum(T[y].query(S[x][0][S[x][1]], distance_upper_bound=60.0)[0], 60.0) for x in S for y in S if x != y}


def cmd_select(a):
    S = {d: load(d) for d in a.dirs}; T = {d: cKDTree(dense(*S[d])) for d in S}; D = dists(S, T)
    zone = {x: np.linalg.norm(S[x][0][S[x][1]] - mesh_seed(x), axis=1) < 100 for x in S}  # ring 0-100 around its own seed
    res = {}
    for x in S:
        o = [y for y in S if y != x]
        le4 = [np.mean(D[x, y] <= 4) for y in o]; wrap = np.mean([np.mean((D[x, y] > 11) & (D[x, y] <= 19)) for y in o])
        near = [np.mean(D[x, y][zone[x]] <= 4) >= 0.5 if zone[x].any() else False for y in o]  # y agrees on x's seed zone
        res[x] = dict(median_le4=float(np.median(le4)), one_wrap_11_19=float(wrap), score=float(np.median(le4) - wrap), support=float(np.mean(near)))
    best = max(res, key=lambda x: res[x]["score"]); sup = max(r["support"] for r in res.values())
    out = dict(choice=best if sup >= 2 / 3 - 1e-9 else None, max_support=sup, surfaces=res)
    print(out["choice"] or "none", json.dumps({k: round(v, 3) for k, v in res[best].items()}), f"max_support={sup:.2f}")
    if a.json: json.dump(out, open(a.json, "w"), indent=1)


def trim_mask(S, T, x, others, low=2 / 3):
    """Support of each vertex of x = share of `others` within 4 vx; keep the seed-connected part (>= 2/3, grown down to `low`)."""
    P, v = S[x]; n = np.zeros(v.shape); d = np.zeros(v.shape)
    for y in others:
        dy = np.minimum(T[y].query(P[v], distance_upper_bound=60.0)[0], 60.0); n[v] += dy <= 4; d[v] += dy / len(others)
    sup = n / len(others); R = np.linalg.norm(P - mesh_seed(x), axis=-1)
    strong = v & (sup >= 2 / 3 - 1e-9); lab, _ = label(strong)  # the seed vertex itself may sit on a ragged edge: use the nearest supported vertex within 20 vx
    R[~strong] = np.inf; seed = np.unravel_index(np.argmin(R), R.shape); keep = lab == lab[seed] if R[seed] < 20 else np.zeros(v.shape, bool)
    if low < 2 / 3 - 1e-9 and keep.any(): keep = binary_propagation(keep, mask=v & (sup >= low - 1e-9))
    return keep, sup, d


def cmd_trim(a):
    S = {d: load(d) for d in [a.surf] + a.others}; T = {d: cKDTree(dense(*S[d])) for d in a.others}
    keep, sup, _ = trim_mask(S, T, a.surf, a.others, a.low); v = S[a.surf][1]
    np.save(a.out + ".npy", keep); out = dict(surface=a.surf, others=len(a.others), kept_vertices=int(keep.sum()), valid_vertices=int(v.sum()), kept_share=float(keep.sum() / v.sum()))
    json.dump(out, open(a.out + ".json", "w"), indent=1); print(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("seeds"); s.add_argument("mesh"); s.add_argument("--offset", type=int, default=3, help="grid cells from the seed vertex")
    s.add_argument("--center", help="use this fractional grid position (i,j) of MESH instead of its meta.json seed")
    s = sp.add_parser("score"); s.add_argument("dirs", nargs="+"); s.add_argument("--json")
    s = sp.add_parser("select"); s.add_argument("dirs", nargs="+"); s.add_argument("--json")
    s = sp.add_parser("trim"); s.add_argument("surf"); s.add_argument("others", nargs="+"); s.add_argument("--out", required=True)
    s.add_argument("--low", type=float, default=2 / 3, help="hysteresis: keep growing the kept zone down to this support")
    a = ap.parse_args()
    if a.cmd == "select" and len(a.dirs) < 3: ap.error("select needs at least three surfaces")
    if a.cmd == "score" and len(a.dirs) < 2: ap.error("score needs at least two surfaces")
    {"seeds": cmd_seeds, "score": cmd_score, "select": cmd_select, "trim": cmd_trim}[a.cmd](a)


if __name__ == "__main__":
    main()
