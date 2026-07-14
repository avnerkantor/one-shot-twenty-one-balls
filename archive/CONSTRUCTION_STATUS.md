# Witness construction: honest status after the aiming-bug fix

## The bug that invalidated earlier "progress"

`solve_aim_struck` (and the trunk aim solver) were converging onto the
+-pi wrap-around discontinuity of the bearing-error function: a fraction
of all "aimed" cuts actually sent the struck ball 180 degrees away from
its target pocket, while every downstream audit modelled it travelling
toward the pocket. Captured proof: a struck ball departing at -138.7
degrees against a pocket bearing of +41.3 degrees, at the planned speed.

This bug silently inflated every liveness count, depth menu, and success
estimate from the moment the collision-tree architecture was introduced.
It is now fixed (reject sign changes whose bracket straddles pi, and
verify |error| < 1e-3 at the returned root).

## What is now TRUE and load-bearing

With honest aiming:

* Trunk skeletons build cleanly; all carriers are aimable by
  construction (the trunk builder now calls the real aim test, not a
  weak proxy).
* The space-time audit is exact: closed-form quartic closest-approach
  between decelerating segments, parked-cue phantoms, true ball paths,
  and an assignment-level crossing prefilter, with the event simulator
  as the final gate.
* No duplicate function definitions remain (grep -c '^def ' == 24).

## The real obstacle (now visible because the bug is gone)

Carrier speed economics. Carriers enter their branches at 2 to 4 m/s.
Each aiming cut is "fat" (large deflection) and drains speed. From most
carrier positions and headings, a carrier can cleanly pot only ONE
branch ball into a valid pocket before running out of either speed or
geometricaly-clear pocket lines. Measured menus for a representative
N_TRUNK=8 skeleton (6 branch-bearing carriers, need=13):

    carrier: 1  2   3   4  5  6      max-sum
    depths : 1  2,3 1,2 1  1  1        9      (need 13)

So the branches cannot absorb all 21-N_TRUNK object balls. Raising the
trunk length lengthens the chain but lowers carrier speeds further;
lowering it raises speeds but leaves too few carriers. The feasible band
is narrow and none of {6,7,8} closed with depth<=3 branches.

## What would actually close it (next session's plan)

Three concrete directions, in order of promise:

1. MIXED-DEPTH with a long trunk carrying the slack. Let the trunk be
   longer (N_TRUNK ~ 12-14), each trunk ball a THIN cut that keeps the
   chain fast, and give branches depth 1 only (one aimed pot + carrier
   self-pot = 2 balls per carrier). With N_TRUNK=13 you need 8 balls
   from ~11 carriers at depth 1: very slack. The earlier N_TRUNK=11
   attempt died at level 8 of the trunk beam, NOT at branches, so the
   fix is trunk-beam survival (wider beam, thinner cuts, a lookahead
   term rewarding headings that keep a pocket reachable), which the
   honest aim solver now makes meaningful.

2. TWO-BALL BRANCHES ONLY (depth 1 everywhere), sizing N_TRUNK so that
   (N_TRUNK - 2) * 2 + (N_TRUNK) >= 21, i.e. N_TRUNK >= 9 with every
   carrier potting exactly one extra ball. This is the cleanest count:
   9 trunk balls, 7 branch-carriers x 2 = 14... = 21 exactly if the
   last trunk ball self-pots. Tighten the branch builder to guarantee a
   depth-1 pot from every retained carrier (the trunk builder already
   requires depth-1 aimability, so this is mostly bookkeeping).

3. Abandon the strict tree for a DAG: allow a potted ball's carrier to
   re-enter the trunk corridor and service a second cluster. More
   bookkeeping, but it decouples "number of carriers" from "number of
   trunk balls" and removes the speed/length tension entirely.

Option 2 is the recommended next step: the counts are exact and the
honest aim solver plus exact audit are already in place to certify it.

## Files

`witness_tree.py` is syntactically clean, imports fine, and generates
honest skeletons; it does NOT yet emit a completed witness_tree.json.
Everything else (model_sim, mc_break, witness_verify, make_figures,
main.tex) is unaffected and ready.

---

## UPDATE: the cushion-fold breakthrough

The structural obstacle above (long trunk runs off the table; short trunk
starves carriers) is SOLVED by letting the trunk carrier rebound off a
cushion between cuts. Measured facts:

* Cushion rebounds keep 91-96% of speed at useful incidence angles
  (E_CUSH = 0.90, grazing bounces lose almost nothing). So folding the
  chain costs essentially no speed budget.
* A pruned trunk beam WITH an optional single bounce between cuts (cap 3
  bounces total) survives to level 9+ of an N_TRUNK=12 chain with
  ~2900 live states and no sign of collapse. Previously the straight
  trunk died hard at level 8 (70, 78, 0). The run only stopped on a
  95-second wall-clock cap, not on emptiness.

`trunk_cushion.py` (delivered) implements the rebound geometry:
`reflect_point_dir`, `rebound_speed`, `carry_with_optional_bounce`,
all validated (a ball up-right off the top cushion folds down-right at
95% speed).

### Remaining work to a completed witness (next session, well-scoped)

1. Wire `carry_with_optional_bounce` into `build_trunk_skeletons` as the
   between-node carrier path (the prototype loop that produced the level-9
   survival is ready to port; keep the bounce cap and the tight diversity
   buckets that made it tractable).
2. Extend the node-N direct T-final alignment to account for a possible
   bounce on the last segment (reflect the target pocket across the wall,
   solve the straight alignment to the image, unreflect).
3. Extend the exact space-time audit to multi-segment carrier paths:
   `branch_timed_elems` and `trunk_timed_segments` must emit one timed
   sub-segment per straight piece between bounces. The seg-seg quartic
   closest-approach test is per-straight-segment and already correct;
   it just needs to be called on each piece.
4. Run attach_branches (depth-1 menus, exact count sum = 21 - N_TRUNK)
   over the now-deep skeleton pool, with the event simulator as the
   final gate (already in place and correct).

The physics and the search are now both on the right side of feasible;
what remains is bookkeeping to thread bounces through the three places
that assume single-segment carrier paths. Estimated one focused session.

---

## UPDATE 2: bounce wired into the builder, but a placement bug remains

Progress: `carry_with_optional_bounce` is now integrated into
`build_trunk_skeletons`. With a cheap geometric prefilter
(`_carrier_prefilter`) gating the expensive aim test and a bounce cap of
3, the builder produces 40 complete N_TRUNK=12 skeletons in ~90 s, every
one using 2-3 cushion folds. The trunk survival problem is definitively
solved: the beam holds ~90-200 states at every level through the full
depth, where the straight trunk died at level 8.

BUG found by simulator check (this is why we sim-gate): the placed
trunk-ball positions are inconsistent with the bounced carrier path.
When the carrier bounces between node j and node j+1, the next trunk
ball was still placed as if on a straight chain, so firing the cue
reproduces cue-T1, T1-T2, then diverges (T2-T4, cascade). Realized pots
= 4, not the intended chain.

Root cause: the trunk records `bounce` and `approach_from` on each
plan step, but `placed[j+1]` (the struck ball's spot) is set to the
carrier's post-deflection point P computed from the PRE-bounce node,
while the carrier actually arrives along the POST-bounce direction. The
struck ball must sit on the post-bounce carrier line.

### Precise fix (next session)

In the trunk expansion, when `bounce` is not None:
* the carrier travels start -> bounce_point -> node along the reflected
  direction; the NODE and the struck ball position P must be computed
  from the post-bounce state (which they now partly are), BUT the
  struck ball must be placed at `node + 2R * u_postbounce`, and the
  chain continuation direction is the post-bounce `outgoing`. Verify
  that `placed[j+1]` equals that point, not a straight-chain image.
* `trunk_timed_segments` must emit TWO timed sub-segments for a bounced
  carrier step: start->bounce_point and bounce_point->node, each with
  its own time window (split at the bounce time). The seg-seg quartic
  test then applies per sub-segment (already correct per-segment).
* `witness_verify.py` already replays in the simulator; use it as the
  gate. A skeleton is only real if the sim reproduces cue-T1-T2-...-TN
  with no unplanned contacts. NONE of the current 40 pass this yet.

The architecture is right and the search is tractable; the remaining
work is making ball PLACEMENT consistent with the bounced PATH, then
re-running attach_branches (depth-1, exact count) with the sim gate.
The prototype that produced level-9 survival used path-only reasoning
(no placement), which is why it looked clean; placement is the missing
link.

---

## UPDATE 3: crash fixed, cue-chase fixed, direction-drift exposed

Fixed this session:
1. CRASH (`TypeError: NoneType not iterable`): `attach_branches` iterated
   over `build_branch(...)` results without guarding the `None` return.
   Guarded at both call sites (menu build + DFS rebuild). The file now
   runs to completion instead of crashing.
2. CUE-CHASE: the first cue-T1 contact was a FULL hit, leaving the cue a
   residual (1-E)/2 * V0 ~ 0.29 m/s drifting forward along the axis; when
   a trunk ball bounced back off the far cushion it re-collided with the
   drifting cue, wrecking the chain (15 spurious cue-T1 recollisions).
   Fixed by making the first contact a CUT (>= 15 deg), placing T1 so the
   line of centres has the swept bearing; the cue then leaves the axis (0
   recollisions). This lifted simulator chain reproduction from 2 contacts
   to 6 of 12.
3. Added `dist_to_next_wall` and capped every carrier segment to reach its
   node WITHOUT touching a cushion in between (no spurious second bounce).

REMAINING obstacle, now precisely diagnosed: DIRECTION DRIFT across
bounces. The analytic construction assumes each carrier reaches a node
along `approach_from -> node`, but in the simulator a bounced carrier
arrives on a slightly different heading, and the error COMPOUNDS down the
chain. Concretely, at the first post-bounce strike the struck ball departs
at the analytic bearing (e.g. -28.6 deg) but the simulator produces a very
different one (e.g. -64.4 deg), so it misses the next ball and hits a
cushion instead. Chain reproduction plateaus at 6/12 exactly where the
second bounce enters.

Root cause: the one-bounce analytic path is not faithful to the simulator
at the level of exact post-collision headings once bounces accumulate;
open-loop construction cannot stay consistent.

### The fix that will actually close it: BUILD-AND-VERIFY

Stop constructing the whole chain open-loop and then checking. Instead,
grow the chain INSIDE the simulator, one collision at a time:
  * place T1..Tk so far; fire the cue ONCE; read the simulator state at
    the moment the current carrier is about to strike the next ball;
  * solve the next ball's placement against the carrier's ACTUAL
    simulator velocity/position (not an analytic prediction);
  * append, re-fire (or continue the same run), verify the intended
    contact happened, and recurse.
This makes the simulator the single source of truth and eliminates drift
by construction. It is more expensive per step but each step is exact.
`model_sim.py` already exposes everything needed (event log, per-ball
state); the constructor becomes a loop around it rather than a parallel
analytic model.

State: `witness_tree.py` runs cleanly, reproduces 6/12 of the trunk chain
in-simulator, and no longer crashes. The remaining work is the
build-and-verify rewrite of the trunk constructor. Monte Carlo half of
the paper remains ready.
