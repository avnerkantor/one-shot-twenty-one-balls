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

## Direction window (the far side)

`miss/` holds the certificates behind the two-sided bounds on the
window of directions, produced by `certify_miss.py` at 1200 bits.
They concern the set of strokes that realise the witness's event list
(the same 61 events, each ball's events in the same order).

| file | mode | result |
|---|---|---|
| `scan.txt` | `scan` | bisection at exactly represented directions: the skeleton is realised at phi* + 9.0454e-43 and phi* - 4.3587e-44, and a designated event provably fails at directions larger by a relative 1e-30 |
| `miss_plus_9.05e-43.txt` | `miss 9.05e-43 1e-90` | CERTIFIED MISS: for every direction in the box, ball 21 (`b20`) reaches the right cushion line before the corner pocket disk |
| `miss_minus_4.36e-44.txt` | `miss -4.36e-44 1e-90` | CERTIFIED MISS: for every direction in the box, `b20` reaches the bottom cushion line before the corner pocket disk |
| `point_plus_9.045e-43.txt` | `point 9.045e-43 0.01` | CERTIFIED: this direction, with V0 in V0* +/- 0.01, realises the event list and pots all 21 |
| `point_minus_4.358e-44.txt` | `point -4.358e-44 0.01` | CERTIFIED: likewise on the other side |

A miss holds for every speed, because in this witness the geometry of
the paths depends on the direction alone. A miss is not a proof that
fewer than 21 balls are pocketed: just outside the window `b20`
rebounds from the cushion inside the pocket jaws and may still drop.
Boxes wider than about 1e-80 in direction are undecided for the same
reason that stops `certify_mpmath.py` at 1e-85: a box loses far more
width through twenty thin cuts than the dynamics amplifies.

## The whole window

`window/` holds the certificates of the rectangle
V0 in V0* + [-2.3, +0.1] (which covers [V0* - 2.3, V_max]) times
phi in phi* + [-4.3586e-44, 9.0454e-43], produced by `certify_window.py`
at 400 bits: 48 sub-strips of width 0.05 m/s, each covered by an
adaptive bisection of the window into 25 pieces, the same
pieces for every sub-strip since the geometry does not depend on the
speed (see
`window/SUMMARY.md`). Every stroke in the rectangle realises the
witness's event list and pots all 21.

Why this works where boxes do not: the geometry is built in first-order
centred forms in phi, so the enclosure width is proportional to the
true sensitivity (about 3e40 m/rad) rather than to the box method's
dependency loss (about 1e80 m/rad). Pieces about 1e-43 wide are
decided in a few seconds; the widest piece is 1.2e-43 wide and the
narrowest, near the two edges, 4.6e-46.

Together with `miss/miss_plus_9.0455e-43.txt` and
`miss/miss_minus_4.3587e-44.txt`, this pins the window of directions
to five significant figures: its width is 9.4813e-43 rad.
