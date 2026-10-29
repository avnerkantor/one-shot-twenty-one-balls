# Saved certificate outputs

| file | certifier | box half-width (V0 and phi) | result |
|---|---|---|---|
| `mpmath_point.txt` | `certify_mpmath.py`, 1200 bits | 0 (the exact stroke) | CERTIFIED |
| `mpmath_box1e-90.txt` | `certify_mpmath.py`, 1200 bits | 1e-90 | CERTIFIED |
| `arb_point.txt` | `certify_arb.py`, 1200 bits | 0 | CERTIFIED |
| `arb_box1e-90.txt` | `certify_arb.py`, 1200 bits | 1e-90 | CERTIFIED |
| `box_search/` | `certify_mpmath.py` | 1e-110 to 1e-84 | certified at 1e-110, 1e-100 and 1e-95 (and 1e-90 above); 1e-88, 1e-86, 1e-85 and 1e-84 fail (enclosures too wide; the 1e-88, 1e-86 and 1e-84 runs keep V0 at +/- 1e-90) |

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
Wide boxes of directions are undecided for the same reason that stops
`certify_mpmath.py` at half-width 1e-88: a box loses far more width
through twenty thin cuts than the dynamics amplifies.

## The whole window

`window/` holds the certificates of the rectangle
V0 in V0* + [-2.3, +0.1] (which covers [V0* - 2.3, V_max]) times
phi in phi* + [-4.3586e-44, 9.0454e-43], produced by `certify_window.py`
at 400 bits: 48 sub-strips of width 0.05 m/s, each covered by an
adaptive bisection of the window into 23 pieces, the same
pieces for every sub-strip since the geometry does not depend on the
speed (see
`window/SUMMARY.md`). Every stroke in the rectangle realises the
witness's event list and pots all 21.

Why this works where boxes do not: the geometry is built in first-order
centred forms in phi, so the enclosure width is proportional to the
true sensitivity (about 3e40 m/rad) rather than to the box method's
dependency loss (about 1e80 m/rad). Pieces about 1e-43 wide are
decided in a few seconds; the widest piece is 2.4e-43 wide and the
narrowest, near the two edges, 4.6e-46.

Together with `miss/miss_plus_9.0455e-43.txt` and
`miss/miss_minus_4.3587e-44.txt`, this pins the window of directions
to five significant figures: its width is 9.4813e-43 rad.

## The lower edge of the ribbon

Along a fixed skeleton the squared speed of every ball is an affine
function of V0^2: v^2 = alpha V0^2 - beta - 2 a s (alpha, beta depend on
the direction only). A stroke can realise the event list only if every
ball is still moving at the end of each of its pieces, i.e. only if
V0 > V_c(phi) = sqrt(max (beta + 2 a ell)/alpha). `certify_lower.py`
propagates alpha and beta in the centred forms of `certify_window.py`
and encloses V_c; `lower/` holds its outputs.

| file | result |
|---|---|
| `lower/window.txt` | over phi* + [-4.3586e-44, 9.0454e-43] (400 bits, 23 pieces): V_c in [9.3203575, 9.3213142] m/s; the critical piece is everywhere the last piece of b18 (ball 19), from the left cushion to the corner pocket (0, 0) |
| `lower/point.txt` | over phi* +/- 1e-90 (1200 bits): V_c(phi*) = 9.32127016844677 m/s (V_c^2 = 86.8860775531756) |

The certified rectangle was then extended down to the edge:

- `window/sub_-2.35_-2.3.txt` ... `window/sub_-2.578685_-2.57868.txt`:
  ten sub-strips covering V0 in V0* + [-2.578685, -2.3], i.e.
  [9.321315, 9.6] m/s, over the whole window (see `window/SUMMARY.md`).
  Together with the 48 earlier sub-strips: every stroke with V0 in
  [9.321315, 12] and phi in the window realises the event list.
- `strip/edge_-2.5787298305_-2.578685.txt`: `certify_strip.py` at
  400 bits on V0 in [9.3212701695, 9.321315], phi* +/- 1e-90. Its lower
  end is 1.05e-9 m/s above V_c(phi*); b18 enters its pocket at about
  1e-5 m/s there.

So for every direction in the window the ribbon's speeds are contained
in (V_c(phi), 12] and contain [9.321315, 12]; at phi* the edge is
located to within 1.1e-9 m/s. The measure bound of the success set
improves to (12 - 9.321315) * 9.48126e-43 > 2.53e-42.

## The near-frozen cluster (open problem P4)

`cluster/window.txt` certifies a second clearance of Model M,
`construction/cluster_21.json`: the cue ball strikes twelve entry balls
over short legs and the other nine balls are struck by their neighbours
across gaps of 1.5 to 3 mm (see `construction/witness_cluster.py`). The
run is

    CERT_WITNESS=construction/cluster_21.json CERT_PREC=400 CERT_TM_DEGREE=6     CERT_MIN_WIDTH=1e-7 python proof/certify_window.py 2.5e-5 9.15e-4 -0.001 0

| file | result |
|---|---|
| `cluster/window.txt` | CERTIFIED: every stroke with phi in phi_c + [-2.5e-5, 9.15e-4] and V0 in [V0_c - 0.001, V0_c] realises the cluster's 61 events and pockets all 21 object balls; 64 pieces of width 1.469e-5 rad, 400 bits, Taylor models of degree 6; the slowest pocket entry is at 0.250018 m/s or more |

Why Taylor models: along a chain of near-frozen collisions the true
slope of every quantity with respect to the aim stays of order 1 to
10 rad/rad, but a slope computed in plain intervals loses the
cancellation that keeps it small and its width grows by about
1/(2R cos theta) at every collision. The first-order centred forms of
`certify_window.py` therefore decide pieces only 4e-12 rad wide on this
configuration (some 2.5e8 pieces). With `CERT_TM_DEGREE=6` the geometry
is carried in degree-6 Taylor models (`proof/taylor_models.py`): the
coefficients are computed from thin inputs and the inflation acts on
the remainder, which starts at order h^7, so pieces of 1.5e-5 rad are
decided in about 20 s.

The certified rectangle is 9.4e-4 rad wide in direction, against
9.48e-43 rad for the chain witness, but it is not the whole window: the
80-digit replay (`construction/tree_measure.py`) brackets the window at
phi_c + [-4.24e-5, 9.22e-4], and within about 1.5e-5 rad of either
edge the time-separation check of `certify_strip.py` cannot decide one
near pass. The speed range is narrow for the same reason: on four legs
the cue passes the ball it will strike next with 0.09 to 0.63 mm of
clearance, and the position enclosures grow with the speed range.
Nothing is claimed about the window's edges or about wider speeds.
