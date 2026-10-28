# One Shot, Twenty-One Balls

**Can a single stroke of snooker clear all twenty-one object balls?**

Snooker folklore says no. This repository accompanies a short mathematics
paper showing that, in a fully specified idealized model of billiard
dynamics, the answer is *yes in principle* and *never in practice*. It
contains the computer-assisted proofs of the first half (a certified
witness, the ribbon of successful strokes, its window of directions and
its lower edge in speed) and the simulations behind the second.

<p align="center">
  <img src="figures/witness_clearance.gif" width="720"
       alt="A single cue stroke pocketing all 21 object balls, animated.">
</p>

## What is proved

**Theorem (computer-assisted).** For the table layout in
[`proof/witness_21.json`](proof/witness_21.json), every cue stroke with

* speed `V0` within `1e-90` m/s of `V0*`, and
* direction `phi` within `1e-90` rad of `phi*`

pockets all 21 object balls in Model M. So the set of winning strokes
has positive measure.

**Stronger: a long, thin ribbon.** Every stroke with speed anywhere in
`[V0* - 2.3, V0* + 0.1]` (about 9.6 to 12 m/s, the top of the model's
shot range) and direction within `1e-90` rad of `phi*` also pots all 21.
Almost any strength works, but the aim has to be right to 90 decimal
places. (The true window of directions is certified below: it is
`9.4813e-43` rad wide.) This follows from a geometric fact about this
witness. Every collision
is a moving ball striking a ball at rest, so every ball's path is the
same polyline whatever the speed; only the timing changes.
[`proof/certify_strip.py`](proof/certify_strip.py) proves that the
paths are as planned for the whole angle box. It also proves that no
two balls are ever within one diameter of each other at the same
instant for any speed in the range.

**The window of directions.** The strokes that realise the witness's
event list (the same 61 events, each ball's events in the same order)
have directions in `phi* + [-4.3586e-44, 9.0454e-43]` rad, a window
`9.4813e-43` rad wide to five significant figures, at every speed in
the ribbon. [`proof/certify_window.py`](proof/certify_window.py)
certifies the whole rectangle of speeds times directions using
first-order centred forms in the direction, which the plain box method
cannot do (boxes of half-width `1e-88` are already undecided), and
[`proof/certify_miss.py`](proof/certify_miss.py) certifies that just
outside the window, at `phi* + 9.0455e-43` and `phi* - 4.3587e-44`,
ball 21 reaches a cushion line inside the pocket jaws before the pocket
disk, at any speed.

**The lower edge.** Along a fixed skeleton every squared speed is an
affine function of `V0^2`, so the event list is lost exactly when the
slowest ball comes to rest at the lip of its pocket, at a critical
speed `V_c(phi)`. [`proof/certify_lower.py`](proof/certify_lower.py)
encloses it: `V_c(phi*) = 9.32127016844677` m/s, and `V_c` lies in
`[9.32036, 9.32131]` m/s across the window. The rectangle
`[9.321315, 12]` m/s times the window is certified, and along `phi*`
the ribbon is certified down to `9.3212701695` m/s, which is `1.05e-9`
m/s above `V_c(phi*)`. The success set has Lebesgue measure at least
`2.53e-42` (m/s times rad).
[`analysis/first_order.py`](analysis/first_order.py) explains the size
of the window by composing first-order transfer maps along the chain;
the dominant factor per thin cut is `c*l/(2R cos(theta))`, between 17
and 240 in this witness, and the composed maps reproduce the exact
sensitivity of `3.184e40` m/rad to four figures.

The box theorem is a rigorous replay of the stroke in interval
arithmetic, checked by two independent programs written against
different interval libraries (the ribbon, window and lower-edge
certificates are single implementations built on mpmath):

| | library | file |
|---|---|---|
| certifier 1 | mpmath (`mpmath.iv`) | [`proof/certify_mpmath.py`](proof/certify_mpmath.py) |
| certifier 2 | Arb (`python-flint`), written independently from the model spec | [`proof/certify_arb.py`](proof/certify_arb.py) |

At every one of the 61 events (21 collisions, 19 cushion rebounds, 21
pocketings), each certifier proves three things for every stroke in the
box:

* the event happens at a unique, transversal instant inside an explicit
  time enclosure;
* nothing else (no other ball pair, cushion, pocket or stop) happens
  first;
* the state after the event is enclosed correctly.

Saved outputs are in [`proof/certificates/`](proof/certificates/).

> **Exact numbers.** All coordinates and shot parameters are the
> binary64 (IEEE double) values that the JSON file decodes to, taken as
> exact rationals. For example, `V0* = 11.9000000000000003552713678800500929355621337890625`,
> not the decimal 11.9. With a box this thin, the difference matters.

## Reproduce the proof (about 2 minutes for certifier 1, under a second for certifier 2)

```bash
pip install -r requirements.txt

# certifier 1: the exact stroke, then the whole box
CERT_PREC=1200 python proof/certify_mpmath.py proof/witness_21.json 0 0
CERT_PREC=1200 python proof/certify_mpmath.py proof/witness_21.json 1e-90 1e-90

# certifier 2 (independent implementation)
python proof/certify_arb.py --r 0
python proof/certify_arb.py --r 1e-90

# the speed ribbon: one sub-strip (about 40 s); the 48 sub-strips covering
# V0* + [-2.3, 0.1] are listed in proof/certificates/strip/SUMMARY.md
CERT_PREC=400 python proof/certify_strip.py -0.05 0 1e-90

# the window of directions: one speed interval of the rectangle (about 30 s);
# the 58 intervals are listed in proof/certificates/window/SUMMARY.md
CERT_PREC=400 python proof/certify_window.py 4.3586e-44 9.0454e-43 -0.05 0
# the two certified misses just outside the window
CERT_PREC=1200 python proof/certify_miss.py miss 9.0455e-43 1e-90
CERT_PREC=1200 python proof/certify_miss.py miss -4.3587e-44 1e-90

# the lower edge: V_c over the whole window, then at phi* alone
CERT_PREC=400 python proof/certify_lower.py
CERT_PREC=1200 python proof/certify_lower.py 1e-90 1e-90

# the first-order (ray-transfer) analysis of the sensitivity
CERT_PREC=400 python analysis/first_order.py
```

Each run ends with `CERTIFIED` and exits with code 0, or prints
`NOT CERTIFIED: <reason>` and exits with code 1. (On Windows PowerShell,
set the variable first: `$env:CERT_PREC=1200`.)

## Why exact arithmetic is essential

A thin cut turns any small lateral error into a larger change of
direction. Errors grow by a factor of tens to hundreds per collision,
about 10^39 over the whole chain. A double-precision
simulation (16 digits) therefore carries no information about the true
trajectory after a few collisions. Earlier versions of this project
built the witness in double precision, and that witness turned out
**not** to be a witness of Model M. It is kept, with an explanation, in
[`archive/float_witness/`](archive/float_witness/). The present witness
was built in 80-digit arithmetic and then certified as above.

## Repository map

```
model/          the physics of Model M
  model_sim.py        fast double-precision engine (used for Monte Carlo)
  model_sim_mp.py     the same engine at 80 digits, no numerical nudges
  make_model_sim_mp.py  regenerates model_sim_mp.py from model_sim.py
proof/          the rigorous part
  witness_21.json     the certified layout and stroke
  certify_mpmath.py   certifier 1
  certify_arb.py      certifier 2
  certify_strip.py    the speed-ribbon certificate
  certify_window.py   the window of directions, in centred forms
  certify_miss.py     the certified misses just outside the window
  certify_lower.py    the lower edge of the ribbon in speed
  certificates/       saved outputs (box_search/, strip/, window/, miss/, lower/)
analysis/       first_order.py: the ray-transfer analysis of the sensitivity
construction/   how the witness was found
  witness_bv.py       build-and-verify constructor (grows the chain ball by ball)
  witness_bv_exact.py runs it in the 80-digit engine
  other_chains/       two more 21-ball chains it found (not certified)
  bv_robust_3ball.json  a short chain with a visibly fat success set
simulation/     how rare is rare: Monte Carlo experiments
  mc_break.py, mc_results.json, mc_sensitivity.json, fit_summary.txt
                      strokes from the opening break
  mc_random_layouts.py, random_layout_runs/
                      strokes from random layouts of the 22 balls
figures/        scripts and images
archive/        superseded approaches, kept for provenance
```

## The model in one paragraph

The table is a 3.569 m × 1.778 m rectangle. It holds 22 balls of radius
26.25 mm and has six pockets: a ball is pocketed when its centre comes
within 55 mm of a pocket point. Between events each ball decelerates at
a constant 0.10 m/s² along its direction of motion. Ball-ball collisions
are frictionless with normal restitution 0.95. A cushion reverses the
normal velocity and scales it by 0.90. There is no spin. A stroke is a
pair `(V0, phi)`. Between events the motion is explicit, so a
trajectory is a finite chain of explicit algebraic maps.

## Rebuild a witness from scratch

```bash
python construction/witness_bv_exact.py 7 my_chain.json   # ~6 minutes
CERT_PREC=1200 python proof/certify_mpmath.py my_chain.json 0 0
```

Seed 7 reproduces `proof/witness_21.json`. Other seeds give other chains,
and some seeds dead-end.

## Reproduce the simulations and figures

```bash
python simulation/mc_break.py 30000        # P(k) from the opening break (resumable)
python simulation/mc_break.py 2500 sens    # friction x restitution sensitivity grid
python simulation/mc_random_layouts.py 25000 1 run_1.json   # random layouts (seeds 1 to 8 were used)
python figures/make_figures.py             # witness layout + P(k) figure
python figures/make_gif.py                 # the animation above
```

**Are random layouts friendlier?** Not in the tail. Each trial places
all 22 balls at random and plays a random stroke. In 200,000 trials,
random layouts beat the break by a factor of two or three for two to
six balls. That advantage fades as the count grows. The best stroke
potted nine balls. Extrapolating to all twenty-one gives about 10^-12
to 10^-13, about the same as the break.

## Citation

See [`CITATION.cff`](CITATION.cff). Note that version 1 of the arXiv
preprint (arXiv:2607.12995) describes the earlier double-precision
witness; version 2 describes the certified witness and the theorems
listed above.
