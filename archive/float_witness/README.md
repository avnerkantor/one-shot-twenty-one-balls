# The earlier double-precision witness (superseded, not valid in Model M)

`witness_21_float.json` is the 21-ball layout described in arXiv:2607.12995
v1. It was built and checked in the double-precision engine
`model/model_sim.py`, where `witness_verify.py` still reports
"VERIFIED". **It is not a witness of Model M.**

Two things the double-precision engine does are not part of the model:

1. **Rounding.** Each thin cut multiplies a small lateral error by roughly
   10 to 70. After a handful of collisions, 16-digit arithmetic no longer
   follows the true trajectory.
2. **A separation nudge.** After every collision and cushion contact the
   engine pushes the balls apart by 1 nm (`EPS_SEP`), so that they do not
   re-detect the contact. That perturbation is amplified in the same way.

Replayed in the exact dynamics, the layout behaves differently. Both
`proof/certify_mpmath.py` (interval arithmetic) and `model/model_sim_mp.py`
(80 digits, no nudge) show the same thing: the trajectory departs from
the design at the 8th event, object ball b5 strikes b12 at the 12th
event, and only 4 balls are pocketed. Setting `EPS_SEP = 0` in the
double-precision engine already reproduces the divergence.

The construction method itself was sound. Rerun in the 80-digit engine
(`construction/witness_bv_exact.py`), it produced the certified witness in
`proof/witness_21.json`.

`model_sim_v0.py` is the independent reference engine that was used to
cross-check `model_sim.py` in double precision.
