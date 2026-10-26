# Saved certificate outputs

| file | certifier | box half-width (V0 and phi) | result |
|---|---|---|---|
| `mpmath_point.txt` | `certify_mpmath.py`, 1200 bits | 0 (the exact stroke) | CERTIFIED |
| `mpmath_box1e-90.txt` | `certify_mpmath.py`, 1200 bits | 1e-90 | CERTIFIED |
| `arb_point.txt` | `certify_arb.py`, 1200 bits | 0 | CERTIFIED |
| `arb_box1e-90.txt` | `certify_arb.py`, 1200 bits | 1e-90 | CERTIFIED |
| `box_search/` | `certify_mpmath.py` | 1e-110 to 1e-85 | certified down to 1e-95; 1e-85 fails (enclosures too wide) |

Both certifiers report the same 61 events (21 collisions, 19 cushion
rebounds, 21 pocketings) in the same order, with the same margins:
collision normal speed >= 0.502 m/s, cushion normal speed >= 0.106 m/s,
pocket entry speed >= 0.444 m/s.

The `box_search/` runs were made before the files were reorganized, so
their header names the witness `exact_runs/w_s7.json`. That file is
byte-identical to `proof/witness_21.json`.

A failure at a wider box means that the interval enclosures became too
wide to decide. It does not mean that a stroke misses: replays in the
80-digit engine still pot all 21 with V0 changed by +-0.05 m/s. The sensitive
direction is the angle phi.

## Speed ribbon

`strip/` holds 48 runs of `certify_strip.py`. Each run is a sub-strip
0.05 m/s wide in speed, by +-1e-90 rad in angle. Together they cover
V0 in V0* + [-2.3, +0.1]; see `strip/SUMMARY.md`.
`strip_single_dV0.01.txt` is a single run over +-0.01 m/s. A single
run over +-0.05 fails only because the interval enclosure is too loose
for one near-encounter (the cue and b17 just after event 47);
splitting the range into sub-strips resolves it.
