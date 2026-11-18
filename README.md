# seed-coherence

A label-free check for sheet switches in GrowPatch surfaces.

Put four seeds on the same sheet, about 60 voxels apart. Grow each one with `vc_grow_seg_from_seed`, then measure, for every vertex of one surface, its distance to each other surface. If two surfaces followed the same sheet they overlap. If one of them switched wraps, its vertices land about one sheet spacing away (11 to 19 voxels on PHerc0490A at 8.64 µm). The scorer needs no ground truth or labels, only the meshes.

I wrote this while trying to auto-correct sheet jumps on PHerc0490A. I needed a way to tell whether a correction kept a surface on the right sheet, and "fewer bad edges" was not enough. A filter can remove bad edges and still leave the surface on the wrong wrap.

## What it found on PHerc0490A

Three regions, four seeds each (12 directed pairs per region), official m7 surface prediction and normal grids, 100 and then 140 generations. Seeds are 3 grid cells from the seed vertex of an existing mesh (`seedcoh.py seeds`).

There are two different failure modes, and the scorer separates them.

**1. The seed lands on the neighbouring sheet.** A pair "starts on the same sheet" when the two surfaces agree (≤ 4 vx) for more than 90 % of the vertices within 100 vx of their seeds.

| region (source seed) | pairs starting on the same sheet |
|---|---|
| z50 (4200 2084 5848) | 6 / 12 |
| z65 (4432 6884 7604) | 0 / 12 |
| z80 (4416 7064 9360) | 9 / 12 |

Re-seeding on a point of an existing mesh is known to sometimes land on the next sheet. Here it happened for half the pairs in one region and for all of them in another.

**2. Surfaces that start together drift apart as they grow.** For the pairs that started on the same sheet, this is the fraction of vertices within 4 vx of the other surface, by distance from the seeds (140 generations):

| distance from seeds (vx) | 0-100 | 100-200 | 200-400 | 400-800 | 800-1600 | > 1600 |
|---|---|---|---|---|---|---|
| z50 | 1.00 | 0.99 | 0.94 | 0.75 | 0.40 | 0.12 |
| z80 | 0.98 | 0.75 | 0.40 | 0.30 | 0.17 | 0.10 |

Over all pairs, 12 to 19 % of vertices sit 11 to 19 vx from the other surface. I checked that this band is a winding change and not just noise, on a separate z50 run with three seeds. I fitted an umbilicus on the level-5 CT and compared the distance to the axis for each vertex and its nearest neighbour on the other surface. For 72 to 74 % of the band points the radial difference was between 10 and 20 vx (median 12 vx), against a median of 1 vx for points that agree.

Full numbers, including the distance histogram and every pair, are in `results/`.

## Use

```bash
# 4 seeds around the seed of an existing mesh
python3 seedcoh.py seeds path/to/mesh_dir > seeds.txt
# grow them (see run_bench.sh), then
python3 seedcoh.py score A/g140 B/g140 C/g140 D/g140 --json score.json
```

`run_bench.sh PRED_ZARR NORMAL_GRIDS SOURCE_MESH OUT_DIR [THREADS]` does the seeding, growth (100 then 140 generations, `--rng-seed 490`) and scoring in one go. Each tifxyz directory needs `x.tif`, `y.tif`, `z.tif` and `meta.json`. When `meta.json` has `seed: [0,0,0]`, which resumed runs write, the scorer takes the oldest vertex in `generations.tif` as the seed.

Distances are point to surface: the other mesh is densified with bilinear points (8 × 8 per quad) before the nearest-neighbour query.

`python3 tests/test_seedcoh.py` runs two synthetic checks: a mesh against a copy of itself, and against a copy shifted by 15 vx. `python3 -m pytest tests` runs those and the `select`/`trim` tests below (needs numpy, scipy, tifffile, pytest).

## Choosing one surface and trimming it: `select` and `trim`

`score` only measures. These two subcommands act on the same agreement signal.

```bash
# pick one surface out of 3 or more grown from the same site
python3 seedcoh.py select A B C D --json select.json
# keep only the part of one surface that the others support
python3 seedcoh.py trim A B C D --out A_trim            # writes A_trim.npy and A_trim.json
python3 seedcoh.py trim A B C D --out A_trim --low 0.5  # also grow the kept zone down to 50 % support
```

**What `select` measures.** For each surface X, it takes the share of X's vertices within 4 vx of each other surface, and keeps the median over the others. It subtracts the mean share of vertices 11 to 19 vx away (the one-wrap band), so a surface that sits a wrap off from the others loses points. The highest score is the choice. It also computes a support value for each surface: the share of the other surfaces that agree (≥ 50 % of vertices within 4 vx) inside the 100 vx ring around X's own seed. If no surface reaches a support of 2/3, `select` prints `none` and chooses nothing. Abstaining is deliberate: when the surfaces split into two camps, there is no majority sheet to pick.

**What `trim` measures.** For each vertex of the surface, the support is the share of the other surfaces within 4 vx. Vertices with support of at least 2/3 are kept if they are connected (in the grid) to the seed, using the nearest supported vertex within 20 vx of the seed. Everything past the point where the others disagree is dropped. The output is a boolean mask over the grid (`.npy`) plus a small `.json` with the kept share. `trim` does not modify the surface files.

Both need the same kind of input as `score`: tifxyz directories with `x.tif`, `y.tif`, `z.tif` and a `meta.json` carrying the seed. The `scale` field is not read; vertices are taken as stored, with x and z above 0 counting as valid.

### Results on real data

These are single runs on two scrolls. Read them as case reports, not as a benchmark.

- **PHerc0490A, z50, without a seed anchor.** Seed D grew one wrap away from the others. `select` discarded it (median agreement 0, with 42 % of its mass at 11 to 19 vx) and chose B.
- **PHerc0490A, z50, with `seed_anchor_xyz`** (the option ChaseHendrick added in villa #2063). The four surfaces agree and support is 1.0 for all of them. See `docs/seed-anchor.md`.
- **PHerc0841, against the human segment w00.** No net gain. The error of the selected surface is 0.851, against 0.866 for a random pick and 0.820 for the oracle (the best of the four). The reason is that nearly all four default traces sit off w00, so there is little to choose from. `select` abstained (`none`) in 5 of the 8 sets.

In short: `select` removed the surface that hopped a wrap on PHerc0490A, and it did not help on PHerc0841, where nearly all candidates are wrong.

### Synthetic tests

`tests/test_select_trim.py` builds tifxyz surfaces as parallel planes 15 vx apart, standing in for neighbouring wraps. It checks that:

1. four identical surfaces (up to noise) give a choice and support 1;
2. three surfaces on one sheet and one 15 vx away: the odd one is never chosen, and `trim` keeps 0 % of it while the others keep everything;
3. two surfaces on each of two sheets: `select` abstains;
4. a surface that jumps to the next sheet halfway: `trim` keeps the half on the seed's side (and nothing if the seed is on the jumped side).

## Limits

- `select` and `trim` rely on the majority. If most of the surfaces left the sheet together, they agree with each other and `select` picks a surface on the wrong sheet. If all of them are wrong, as on PHerc0841, neither can help.
- `select` needs at least 3 surfaces, and with 4 a 2 against 2 split makes it abstain. The thresholds (4 vx, 2/3, the 11-19 vx band, 100 vx ring) are fixed in the code, and the 11-19 vx band comes from PHerc0490A.
- `trim` keeps only what is connected to the seed. If the seed itself is on the wrong sheet, it keeps the wrong part.
- The evidence is one PHerc0490A case, one PHerc0841 comparison and synthetic planes. The tests do not cover curved or overlapping sheets.
- One scroll, three regions, four seeds per region. It is a measurement, not a ranking of methods yet.
- Disagreement tells you that at least one of the two surfaces left the sheet, not which one.
- The tracer is not deterministic with several threads. Two identical runs gave different meshes, about 0.5 percentage points apart on edge-quality metrics, so small differences between methods need repeats.
- The 11-19 vx band is specific to this scroll and resolution. On another scroll, read the histogram against its own sheet spacing.

## Data

- Surface prediction: `s3://vesuvius-challenge-open-data/PHerc0490A/representations/predictions/surfaces/20250521151210-surface-20260413222639-surface-m7-L0-th0.2.zarr`
- CT: `PHerc0490A/volumes/20250521151210-8.640um-1.2m-116keV-masked.zarr`
- Normal grids: the published `...surface-m7-L0-th0.2.normal-grids` next to the prediction.

MIT licence.

## Notes

- [Checking `seed_anchor_xyz` (villa #2063) on PHerc0490A and PHerc0841](docs/seed-anchor.md): the anchor fixes starts on the neighbouring wrap but does not keep long traces on the sheet.
