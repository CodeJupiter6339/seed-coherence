"""Synthetic checks for `select` and `trim`.

Sheets are parallel planes 15 vx apart (one "wrap"), each a 60x60 tifxyz grid with 10 vx between vertices.
Every surface gets a little z noise and its own meta.json seed near the middle of the grid.
"""
import json, os, subprocess, sys
import numpy as np, tifffile

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(HERE, "..", "seedcoh.py")
N, STEP, SPACING = 60, 10.0, 15.0


def write_surface(d, sheet=0, jump_at=None, jump_sheet=1, seed_ij=(45, 30), noise=0.3, rng=None):
    """tifxyz dir: x.tif, y.tif, z.tif, meta.json. Row index i >= jump_at moves to sheet `jump_sheet`."""
    rng = rng or np.random.default_rng(0)
    os.makedirs(d, exist_ok=True)
    i, j = np.mgrid[0:N, 0:N].astype(np.float32)
    X, Y = STEP * i + 100, STEP * j + 100
    sh = np.full(X.shape, sheet, np.float32)
    if jump_at is not None:
        sh[int(jump_at):] = jump_sheet
    Z = 50 + SPACING * sh + rng.normal(0, noise, X.shape).astype(np.float32)
    for c, a in zip("xyz", (X, Y, Z)):
        tifffile.imwrite(os.path.join(d, c + ".tif"), a.astype(np.float32))
    si, sj = seed_ij
    meta = {"format": "tifxyz", "scale": [1.0, 1.0], "seed": [float(X[si, sj]), float(Y[si, sj]), float(Z[si, sj])]}
    json.dump(meta, open(os.path.join(d, "meta.json"), "w"))
    return str(d)


def run(*args):
    r = subprocess.run([sys.executable, CLI, *map(str, args)], check=True, capture_output=True, text=True)
    return r.stdout


def select(dirs, tmp_path):
    out = tmp_path / "sel.json"
    run("select", *dirs, "--json", out)
    return json.load(open(out))


def trim(surf, others, tmp_path, name="m"):
    out = tmp_path / name
    run("trim", surf, *others, "--out", out)
    return json.load(open(str(out) + ".json")), np.load(str(out) + ".npy")


def make(tmp_path, sheets, **kw):
    rng = np.random.default_rng(1)
    return [write_surface(tmp_path / n, sheet=s, rng=rng, **kw) for n, s in zip("ABCD", sheets)]


def test_four_identical_surfaces(tmp_path):
    dirs = make(tmp_path, [0, 0, 0, 0])
    r = select(dirs, tmp_path)
    assert r["choice"] in dirs
    assert r["max_support"] == 1.0
    assert all(s["support"] == 1.0 and s["median_le4"] == 1.0 for s in r["surfaces"].values())


def test_one_offset_surface_is_never_elected(tmp_path):
    dirs = make(tmp_path, [0, 0, 0, 1])
    bad = dirs[3]
    r = select(dirs, tmp_path)
    assert r["choice"] in dirs[:3]
    assert r["surfaces"][bad]["median_le4"] == 0.0 and r["surfaces"][bad]["support"] == 0.0
    assert r["surfaces"][bad]["score"] < min(r["surfaces"][g]["score"] for g in dirs[:3])
    info, keep = trim(bad, dirs[:3], tmp_path)
    assert info["kept_vertices"] == 0 and info["kept_share"] == 0.0 and not keep.any()
    # a good surface against the same three others (two good, one offset) keeps everything: support 2/3 is enough
    info, keep = trim(dirs[0], dirs[1:], tmp_path, "good")
    assert info["kept_share"] == 1.0


def test_two_against_two_abstains(tmp_path):
    dirs = make(tmp_path, [0, 0, 1, 1])
    r = select(dirs, tmp_path)
    assert r["choice"] is None
    assert r["max_support"] < 2 / 3
    assert "none" in run("select", *dirs).split()[0]


def test_surface_that_jumps_sheet_is_trimmed_to_seed_side(tmp_path):
    good = [write_surface(tmp_path / n, sheet=0, rng=np.random.default_rng(k)) for k, n in enumerate("ABC")]
    # seed at row 15, jump at row 30: rows 0-29 stay on sheet 0, rows 30-59 move one wrap over
    e = write_surface(tmp_path / "E", jump_at=N // 2, seed_ij=(15, 30), rng=np.random.default_rng(9))
    info, keep = trim(e, good, tmp_path)
    assert info["valid_vertices"] == N * N
    assert info["kept_share"] == 0.5
    assert keep[: N // 2].all() and not keep[N // 2:].any()
    # seed on the far side of the jump: the kept half is the one holding the seed, which is off the others' sheet, so nothing survives
    f = write_surface(tmp_path / "F", jump_at=N // 2, seed_ij=(45, 30), rng=np.random.default_rng(10))
    info, keep = trim(f, good, tmp_path, "mf")
    assert info["kept_share"] == 0.0
