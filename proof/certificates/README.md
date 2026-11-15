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

The four files of `certify_mpmath.py` and `certify_arb.py` were
regenerated on 2026-11-14 from `proof/witness_21.json`. Apart from the
name of the input file in the header and the wall time they are
identical, line for line, to the earlier outputs, whose headers named
the same file by its earlier name `witness_21_exact.json`.

## The precision fix of 2026-11-14

Every certificate below this point was regenerated on 2026-11-14. The
earlier outputs, in release v2.0 and at the tag v2.1, came from
certifiers with a defect in their handling of precision, found in an
internal review. `certify_strip.py` set both mpmath contexts, the
interval context `iv` and the point context `mp`, to the precision of
the certificate, and then imported the 80-digit engine
`model/model_sim_mp.py` to obtain the plan of events. That import sets
`mp.dps = 80`, so from then on `mp` worked at 269 bits while `iv` stayed
at 400 or 1200. Three things followed.

- `lo()` and `hi()` converted interval endpoints through `mp` and so
  rounded them to nearest at 269 bits. Enclosures built from them (time
  windows, bounding boxes) could be too small by a relative 1e-81.
- `certify_window.py` and `certify_lower.py` formed the centre of each
  piece of directions as `phi* + (a + b)/2` in `mp`. The centre was
  rounded, by up to 3e-83 rad, and the domain on which a piece was
  certified was shifted off the piece. In the earlier outputs 7 of the 22
  junctions between the pieces of the chain's window, and 50 of the 63 of
  the cluster's, were left uncovered, by gaps of up to 3.3e-83 rad.
- `certify_miss.py` formed the box `centre +/- 1e-90` in `mp`, where
  adding 1e-90 to phi* changes nothing. The two misses were certified at
  single directions, not on boxes.

`certify_mpmath.py` and `certify_arb.py` were not affected. The
corrected certifiers

- copy interval endpoints without rounding (`lo`, `hi`);
- restore the precision of `mp` after the engine has been used, and
  stop with an error if the two precisions ever differ
  (`check_precision`);
- build every piece and every box of directions in interval arithmetic,
  rounded outward (`piece_domain` in `certify_window.py`), so that the
  domain on which a piece is certified contains the piece;
- read decimal arguments as intervals and round them outward, so that a
  certified set contains the set named by the decimal numbers;
- print lower bounds rounded down and upper bounds rounded up.

Two runs are now given in speeds (`CERT_V_ABS=1`) instead of offsets
from V0*. V0* is the binary64 number nearest to 11.9 and exceeds it by
3.6e-16, so V0* - 2.578685 is not 9.321315.

The new outputs certify the same statements, with the same numbers of
pieces and the same margins up to the last printed digit.

**The cue ball after its last planned event.** The plan of events ends,
for the cue ball, with its last collision. After it the cue ball rolls
on, may meet cushions, and comes to rest. These events are not in the
plan. `certify_strip.py` and `certify_window.py` follow the cue ball to
rest and check it against every other ball, as they always did, and
each output now says whether a cushion contact or the stop of the cue
ball can come before the last pocketing. For the chain witness the cue
ball meets no cushion after its last collision. It comes to rest after
the last pocketing for every stroke with V0 >= V0* - 0.3 and direction
within 1e-90 of phi*, and before the last pocketing for every such
stroke with V0 <= V0* - 0.4; over the whole window of directions the
same holds with V0* - 0.25 and V0* - 1.2 (see `strip/SUMMARY.md` and
`window/SUMMARY.md`). So over most of the ribbon the evolution up to the
last pocketing consists of the 61 planned events and one more, the stop
of the cue ball.

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
They concern the set of strokes that follow the witness's plan of
events (the same 61 events, each ball's events in the same order).

| file | mode | result |
|---|---|---|
| `scan.txt` | `scan` | bisection at exactly represented directions: the skeleton is realised geometrically at phi* + 9.0454119232e-43 and at phi* - 4.3586638168e-44, and a designated event provably fails at directions larger by a relative 1e-30. Geometry only: the timing conditions are not examined here |
| `miss_plus_9.0455e-43.txt` | `miss 9.0455e-43 1e-90` | CERTIFIED MISS on the box of half-width 1e-90: ball 21 (`b20`) reaches the right cushion line before the corner pocket disk. This is the miss of Theorem C(ii) |
| `miss_minus_4.3587e-44.txt` | `miss -4.3587e-44 1e-90` | CERTIFIED MISS on the box of half-width 1e-90: `b20` reaches the bottom cushion line before the corner pocket disk. The other miss of Theorem C(ii) |
| `miss_plus_9.05e-43.txt` | `miss 9.05e-43 1e-90` | CERTIFIED MISS, farther out |
| `miss_minus_4.36e-44.txt` | `miss -4.36e-44 1e-90` | CERTIFIED MISS, farther out |
| `point_plus_9.045e-43.txt` | `point 9.045e-43 0.01` | CERTIFIED: this direction, with V0 in V0* +/- 0.01, follows the plan and pots all 21 |
| `point_minus_4.358e-44.txt` | `point -4.358e-44 0.01` | CERTIFIED: likewise on the other side |

Each miss file prints the exact interval of offsets that its box
contains, and the width of the box, 2.0e-90 rad.

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
speed (see `window/SUMMARY.md`). Every stroke in the rectangle follows
the witness's plan of events and pots all 21.

Each output lists its pieces by their exact offsets from phi*. The
pieces are contiguous, the first begins at or below -4.3586e-44 and the
last ends at or above 9.0454e-43, and each piece is certified on a
domain of directions that contains it.

Why this works where boxes do not: the geometry is built in first-order
centred forms in phi, so the enclosure width is proportional to the
true sensitivity (about 3e40 m/rad) rather than to the box method's
dependency loss (about 1e80 m/rad). Pieces about 1e-43 wide are
decided in a few seconds; the widest piece is 2.4e-43 wide and the
narrowest, near the two edges, 4.6e-46.

Together with `miss/miss_plus_9.0455e-43.txt` and
`miss/miss_minus_4.3587e-44.txt`, this brackets the width of the window
of directions between 9.48126e-43 and 9.48137e-43 rad. The two
geometric edges found by `miss/scan.txt` sum to 9.4812783e-43 rad.

## The lower edge of the ribbon

Along a fixed skeleton the squared speed of every ball is an affine
function of V0^2: v^2 = alpha V0^2 - beta - 2 a s (alpha, beta depend on
the direction only). A stroke can follow the plan of events only if
every ball is still moving at the end of each of its pieces, i.e. only
if V0 > V_c(phi) = sqrt(max (beta + 2 a ell)/alpha). `certify_lower.py`
propagates alpha and beta in the centred forms of `certify_window.py`
and encloses V_c; `lower/` holds its outputs.

| file | result |
|---|---|
| `lower/window.txt` | over phi* + [-4.3586e-44, 9.0454e-43] (400 bits, 23 pieces): V_c in [9.3203575, 9.3213142] m/s; the critical piece is the last piece of b18 (ball 19), from the left cushion to the corner pocket (0, 0), and it attains the maximum for every direction of the window: its value of (beta + 2 a ell)/alpha is at least 86.869, that of every other piece at most 86.376 |
| `lower/point.txt` | over phi* +/- 1e-90 (1200 bits): V_c(phi*) in [9.3212701684467657, 9.3212701684467658] m/s (V_c^2 in [86.886077553175597, 86.886077553175598]); the same piece attains the maximum |

The certified rectangle was then extended down to the edge:

- `window/sub_-2.35_-2.3.txt` ... `window/sub_-2.57868_-2.5786.txt` and
  `window/abs_9.321315_9.3213201.txt`: ten runs covering V0 from
  9.321315 m/s up to V0* - 2.3, over the whole window (see
  `window/SUMMARY.md`). The lowest run is given in speeds, so that its
  lower end is the decimal number 9.321315. Together with the 48 earlier
  sub-strips: every stroke with V0 in [9.321315, 12] and phi in the
  window follows the plan of events.
- `strip/edge_9.3212701695_9.321315.txt`: `certify_strip.py` at
  400 bits on V0 in [9.3212701695, 9.321315] m/s, given in speeds,
  phi* +/- 1e-90. Its lower end is 1.05e-9 m/s above V_c(phi*); b18
  enters its pocket at 9.1e-6 m/s or more there.

So for every direction in the window the ribbon's speeds are contained
in (V_c(phi), 12] and contain [9.321315, 12]; at phi* the edge is
located to within 1.1e-9 m/s. The measure bound of the success set
improves to (12 - 9.321315) * 9.48126e-43 > 2.53e-42.

## The near-frozen cluster (open problem P4)

The ten files in `cluster/` certify a second clearance of Model M,
`construction/cluster_21.json`: the cue ball strikes twelve balls
directly, over legs of 1.5 mm to 0.54 m, and the other nine are struck
by object balls struck before them, four of them by a neighbour that
has rolled 1.5 mm and five by a ball that has rolled 0.09 to 0.46 m
(see `construction/witness_cluster.py`). The event list has 57 events:
21 ball-ball collisions, 15 cushion rebounds and 21 pocketings. The
rectangle of strokes is certified one speed interval of depth 0.001 m/s
at a time. The first run is

    CERT_WITNESS=construction/cluster_21.json CERT_PREC=400 CERT_TM_DEGREE=6 \
    CERT_MIN_WIDTH=1e-7 python proof/certify_window.py 2.5e-5 9.15e-4 -0.001 0

The ten files were regenerated on 2026-11-14 with the corrected
certifier (see "The precision fix of 2026-11-14" above). The slowest
pocket entries are lower bounds, printed rounded down.

| file | result |
|---|---|
| `cluster/window.txt` | CERTIFIED: every stroke with phi in phi_c + [-2.5e-5, 9.15e-4] and V0 in [V0_c - 0.001, V0_c] realises the cluster's 57 events and pockets all 21 object balls; 64 pieces of width 1.469e-5 rad, 400 bits, Taylor models of degree 6; the slowest pocket entry is at 0.250017 m/s or more |
| `cluster/sub_-0.002_-0.001.txt` | CERTIFIED: the same for V0 in [V0_c - 0.002, V0_c - 0.001]; 64 pieces of width 1.469e-5 rad; the slowest pocket entry is at 0.249891 m/s or more |
| `cluster/sub_-0.003_-0.002.txt` | CERTIFIED: the same for V0 in [V0_c - 0.003, V0_c - 0.002]; 64 pieces of width 1.469e-5 rad; the slowest pocket entry is at 0.249765 m/s or more |
| `cluster/sub_-0.004_-0.003.txt` | CERTIFIED: the same for V0 in [V0_c - 0.004, V0_c - 0.003]; 64 pieces of width 1.469e-5 rad; the slowest pocket entry is at 0.249639 m/s or more |
| `cluster/sub_-0.005_-0.004.txt` | CERTIFIED: the same for V0 in [V0_c - 0.005, V0_c - 0.004]; 64 pieces of width 1.469e-5 rad; the slowest pocket entry is at 0.249513 m/s or more |
| `cluster/sub_-0.006_-0.005.txt` | CERTIFIED: the same for V0 in [V0_c - 0.006, V0_c - 0.005]; 64 pieces of width 1.469e-5 rad; the slowest pocket entry is at 0.249387 m/s or more |
| `cluster/sub_-0.007_-0.006.txt` | CERTIFIED: the same for V0 in [V0_c - 0.007, V0_c - 0.006]; 64 pieces of width 1.469e-5 rad; the slowest pocket entry is at 0.24926 m/s or more |
| `cluster/sub_-0.008_-0.007.txt` | CERTIFIED: the same for V0 in [V0_c - 0.008, V0_c - 0.007]; 64 pieces of width 1.469e-5 rad; the slowest pocket entry is at 0.249134 m/s or more |
| `cluster/sub_-0.009_-0.008.txt` | CERTIFIED: the same for V0 in [V0_c - 0.009, V0_c - 0.008]; 64 pieces of width 1.469e-5 rad; the slowest pocket entry is at 0.249008 m/s or more |
| `cluster/sub_-0.010_-0.009.txt` | CERTIFIED: the same for V0 in [V0_c - 0.010, V0_c - 0.009]; 64 pieces of width 1.469e-5 rad; the slowest pocket entry is at 0.248881 m/s or more |

The other nine runs have the last two arguments `-0.002 -0.001` down to
`-0.010 -0.009`. The decimal ends of each speed interval are rounded
outward, so adjacent runs overlap and the union of the ten contains the
rectangle phi_c + [-2.5e-5, 9.15e-4] by [V0_c - 0.01, V0_c], with every
pocket entry at 0.248881 m/s or more. Each output lists its 64 pieces of
directions by their exact offsets from phi_c. In each of the ten the
pieces are contiguous, the first begins at or below -2.5e-5 and the
last ends at or above 9.15e-4, and each piece is certified on a domain
of directions that contains it. The ten runs were made four at a time
on one machine and their wall times, 10200 s to 13100 s each, reflect that.

The cue ball's cushion rebounds are events of the plan. After its last
planned event the cue ball meets no cushion before the last pocketing,
and it comes to rest after the last pocketing, for every stroke of the
rectangle (each output states this). So for the cluster the evolution
up to the last pocketing consists of the 57 planned events and no
other.

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
near pass. The speed intervals are shallow for the same reason: on
four legs the cue passes the ball it will strike next with 0.09 to
0.63 mm of clearance, and the position enclosures grow with the speed
range. On the interval [V0_c - 0.01, V0_c] taken whole the run gives
up (`cluster/whole_-0.010_0_fails.txt`, with `CERT_MIN_WIDTH=4e-6`; this
record of a failure is from the earlier code and was not regenerated): on
pieces from 5.9e-5 down to 3.7e-6 rad wide the same check fails to
separate the cue ball from ball b11, which it passes on the way to
b10. This is why the rectangle is certified as ten intervals. Nothing
is claimed about the window's edges or about speeds below V0_c - 0.01.

## Cross-check in Arb, and tests of the Taylor models

The certificates of the ribbon, the window, the lower edge and the
cluster come from one family of programs, all built on mpmath. Two
further files check them from outside. Neither is part of a proof.

| file | what it shows |
|---|---|
| `crosscheck_arb.txt` | output of `proof/crosscheck_arb.py`: the Arb certifier `certify_arb.py`, which shares no code with the mpmath certifiers, run at 52 single strokes. At 30 strokes of the chain's certified rectangle (its corners, both ends of the strips of Theorems B and D(iii), and a grid between them) and at 12 strokes of the cluster's rectangle (its corners and a grid) it certifies that the stroke realises the event list. At the two directions of Theorem C(ii), at three speeds each, it certifies that ball 21 (`b20`) meets a cushion before its pocket. At four strokes below the critical speed it certifies that ball 19 (`b18`) stops before its pocket. All 52 strokes behave as the theorems state |
| `test_taylor_models.txt` | output of `proof/test_taylor_models.py`: 3120 Taylor models of degrees 1, 2, 4 and 6 on pieces of half-width 1e-30 to 0.05, built from random expressions in the operations that `certify_window.py` uses, compared at nine points of each piece with a 1200-bit interval evaluation of the same expression. In all 28080 comparisons the true value lies inside the model. Every remainder contains zero, and 1/x and sqrt(x) refuse a model whose range contains zero. The same test fails, as it should, when the Lagrange remainder of 1/x and sqrt(x) is removed from the code |

The cross-check also shows the cue ball coming to rest before the last
pocketing at the slower strokes, as the certificates above report.
