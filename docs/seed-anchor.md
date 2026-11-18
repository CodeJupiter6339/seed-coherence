# Checking `seed_anchor_xyz` (villa #2063) on PHerc0490A and PHerc0841

ChaseHendrick added the opt-in `seed_anchor_xyz` parameter to `vc_grow_seg_from_seed` in villa #2063 (fixes #2059). It holds one vertex at the explicit seed in every solve, so a trace cannot slide onto another winding at the start. His PR shows it on PHerc0139 w031 with one seed. It is not described in `volume-cartographer/docs/` yet.

I repeated the test on two other scrolls with several seeds per site, using this repository's measures: agreement between surfaces grown from nearby seeds (`seedcoh.py score`), and distance to a human segment where one exists (`truth.py`).

Settings for every run: villa main `8e21d01`, official `m7-L0-th0.2` surface prediction and its normal grids, `--rng-seed 490`, 4 seeds per site taken from vertices of an existing mesh. The anchored arm adds `"seed_anchor_xyz": [X, Y, Z]` equal to `--seed`; nothing else changes.

## PHerc0490A, agreement between 4 seeds (no human segment)

1 thread per run, so runs are deterministic.

| site, generations | ≤ 4 vx, default → anchored | one wrap off (11-19 vx), default → anchored | pairs apart from the first 100 vx, default → anchored |
|---|---|---|---|
| z50, 30 | 0.370 → 0.758 | 0.220 → 0.020 | 6/12 → 0/12 |
| z65a90, 30 | 0.247 → 0.439 | 0.202 → 0.095 | 8/12 → 0/12 |
| z50, 100 | 0.232 → 0.340 | 0.205 → 0.145 | 6/12 → 0/12 |
| z65a90, 100 | 0.109 → 0.164 | 0.143 → 0.148 | 8/12 → 0/12 |
| z80a90, 100 | 0.211 → 0.233 | 0.146 → 0.122 | 0/12 → 0/12 |

"Apart" means a directed pair with less than 0.3 agreement in the 0-100 vx ring around the seeds. Without the anchor, at z50 one of the four seeds (D) grows one wrap away from the other three. With the anchor, every pair starts together. Past a few hundred voxels from the seed, the surfaces still drift apart.

## PHerc0841, distance to the human segment w00

16 seeds placed on w00 itself (4 centres × 4 seeds), 100 generations, 5 threads per run. The score is the share of a surface's vertices within 4 vx of w00, among the vertices w00 covers (within 60 vx).

| | default | anchored |
|---|---|---|
| seed vertex on w00 (≤ 4 vx) after growth | 14 / 16 | 16 / 16 |
| on w00, 0-100 vx from the seed | 0.53 | 0.74 |
| on w00, 100-200 vx | 0.35 | 0.45 |
| on w00, 400-800 vx | 0.16 | 0.16 |
| on w00, whole surface | 0.141 | 0.159 |
| one wrap off w00, whole surface | 0.200 | 0.208 |

The default 0-100 vx figure only counts the 14 surfaces whose seed stayed on w00. The two others had 0.0 there.

## What this shows

- The anchor does what #2063 says: the trace really starts where you put the seed. On both scrolls, the starts on the neighbouring wrap disappear.
- It does not keep the trace on the sheet further out. On PHerc0841, default and anchored surfaces both leave w00 within a few hundred voxels, and only about 15 % of either ends up on w00. This matches the limit written in #2063: "An anchor holds one point, not the sheet."
- So turn the anchor on whenever you seed by hand on a known winding. Do not trust automatic growth far from the seed. `seedcoh.py trim` keeps only the part of a surface that the other seeds agree on, but on PHerc0841 that part was empty or small for default runs.

## Limits

- 3 sites on PHerc0490A and 4 centres on PHerc0841, all with the same tracer version and prediction.
- PHerc0841 runs used 5 threads, and GrowPatch is not deterministic across thread counts, so single-run differences there are partly noise. The table gives means over 16 seeds.
- The w00 score counts as wrong any vertex more than 4 vx from w00, including vertices on a neighbouring wrap that may still be papyrus.

## Reproduce

Scripts in `scripts/`: `voxel_anchor_test.sh` (30 generations), `voxel_anchor100.sh` (100 generations), `voxel_anchor_w00.sh` (PHerc0841 against w00). Scoring uses `seedcoh.py score` and `truth.py` from this repository.
