# Witness result: a verified 21-ball single-stroke clearance (point witness)

## What was achieved

The build-and-verify constructor (`witness_bv.py`) produces a genuine,
simulator-certified configuration of 21 object balls and a single cue
stroke that pockets ALL 21 in one continuous run of Model M.

* File: `witness_21.json` (cue position, V0=11.9 m/s, phi, 21 ball
  coordinates, target pockets).
* Independently re-simulated end to end: 21/21 object balls pocketed,
  cue not pocketed, all six pockets used.
* Admissible configuration: minimum pairwise ball distance 15.4 cm,
  far above the 5.25 cm contact distance.
* Pocket entry speeds range 0.14 to 4.27 m/s.

This decisively refutes folk-claim (i) ("no single stroke can pot all
21") in Model M: such a stroke exists and is exhibited explicitly.

## The architecture that worked

Earlier open-loop constructors predicted the whole collision chain
analytically, then checked it; once cushion bounces accumulated the
analytic path drifted from the simulator and nothing reproduced. The
fix was BUILD-AND-VERIFY: grow the chain one collision at a time,
reading the carrier's ACTUAL simulator state after each committed
collision, placing each new ball so the struck ball departs straight at
a pocket, and re-verifying the whole prefix (including that every placed
ball actually pockets) before committing. The simulator is the single
source of truth, so there is zero drift by construction.

Key physics insight that made 21 reachable: the cue is the sole carrier
and must survive 21 collisions, so THICK cuts (near-grazing, ~75 deg)
are preferred -- they leave the cue 90%+ of its speed while the struck
ball still gets enough to reach a pocket. Randomized-greedy restarts
(`build_greedy`) with softmax tie-breaking then found a full 21-chain in
seconds.

## The open point: measure of the success set

Theorem 1 as drafted claims the success set of shots has POSITIVE
Lebesgue measure (robustness). This point witness, however, pots 21 only
at an essentially isolated shot value: bisection finds the N=21 interval
width in V0 to be below 1e-6 m/s, i.e. numerically a knife-edge. A
21-collision chain compounds sensitivity so steeply that any practical
perturbation tips at least one ball out.

Two honest options for the paper:

1. WEAKEN the theorem to existence of a point witness (N=21 is
   achievable), and state the positive-measure claim as holding by the
   piecewise-smoothness argument (the outcome map is locally constant off
   a measure-zero singular set, so SOME open set around a regular witness
   works) while noting the set is far below observable scale. The
   analytic argument is still valid; only the NUMERICAL demonstration of
   a fat neighborhood is out of floating-point reach at chain length 21.
   This is arguably the most honest framing and fits the paper's thesis
   ("measure zero vs unobservably small") perfectly.

2. BUILD a robust witness: place balls so a whole INTERVAL of cue shots
   pots all 21, by requiring generous angular tolerance at every contact
   (aim at pocket centre with large v_enter margin AND require the
   struck-ball departure cone to map into the pocket disk over a range of
   incoming directions). This is a harder construction -- effectively
   optimizing the configuration for robustness, not just feasibility --
   but the build-and-verify harness supports it: replace the per-step
   accept test with "all 9 corners of a small shot stencil still pot the
   prefix". Slower but yields a genuinely fat witness.

Option 1 is immediately publishable and consistent with the theory.
Option 2 would upgrade the numerical demonstration; recommended as a
follow-up run on stronger hardware (the stencil test multiplies sim cost
by ~9 per step).

## Files
* `witness_bv.py`  -- build-and-verify constructor (build_witness with
  backtracking; build_greedy with randomized restarts).
* `witness_21.json` -- the verified 21-ball point witness.
