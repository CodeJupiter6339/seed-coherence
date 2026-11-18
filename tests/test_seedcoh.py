"""Synthetic checks: a surface scored against itself agrees everywhere; a copy shifted by one sheet spacing lands in the 11-19 vx bin."""
import json, os, subprocess, sys, tempfile
import numpy as np, tifffile

HERE = os.path.dirname(os.path.abspath(__file__))


def write_mesh(d, shift=0.0, seed=(100, 100, 50)):
    os.makedirs(d)
    i, j = np.mgrid[0:60, 0:60].astype(np.float32)
    X, Y = 10 * i + 100, 10 * j + 100
    Z = np.full_like(X, 50.0) + shift
    for c, a in zip("xyz", (X, Y, Z)): tifffile.imwrite(os.path.join(d, c + ".tif"), a)
    json.dump({"seed": list(seed)}, open(os.path.join(d, "meta.json"), "w"))


def score(*dirs):
    out = os.path.join(os.path.dirname(dirs[0]), "s.json")
    subprocess.run([sys.executable, os.path.join(HERE, "..", "seedcoh.py"), "score", *dirs, "--json", out], check=True, capture_output=True)
    return json.load(open(out))


def test_self_and_one_wrap():
    t = tempfile.mkdtemp()
    a, b, c = (os.path.join(t, n) for n in "abc")
    write_mesh(a); write_mesh(b); write_mesh(c, shift=15.0)
    same = score(a, b)
    assert same["mean"]["le2"] == 1.0, same["mean"]
    assert same["pairs_starting_on_same_sheet"] == "2/2"
    off = score(a, c)
    assert off["mean"]["one_wrap_11_19"] == 1.0 and off["mean"]["le4"] == 0.0, off["mean"]
    assert off["pairs_starting_on_same_sheet"] == "0/2"


if __name__ == "__main__":
    test_self_and_one_wrap(); print("ok")
