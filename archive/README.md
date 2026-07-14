# Archive: the open-loop constructors (superseded)

These files are kept for provenance only. They implement the **earlier,
abandoned** approach to building the witness: predict the entire collision
chain analytically (a "trunk" of carrier balls with cushion "folds" and
"branch" clusters), then check it in the engine at the end.

That approach does **not** produce a verified 21-ball clearance. Once
cushion rebounds accumulate, the analytic path drifts from what the
event-driven engine actually computes, and no configuration reproduces its
intended chain end to end. Running `witness_tree.py` will report
`skeleton N failed branch attachment` for every candidate — that is the
expected (honest) failure, not a crash.

The working method is in the top-level `witness_bv.py` (build-and-verify),
which grows the chain one collision at a time *inside* the engine and so
never drifts. See the root `README.md`.

Contents:
- `witness_tree.py`    — open-loop trunk/branch constructor (does not solve it)
- `trunk_cushion.py`   — cushion-rebound geometry helper for the above
- `witness_build.py`   — earlier aiming/collision helpers
- `gen_skeletons.py`   — checkpointed skeleton generator for the above
- `attach_pool.py`     — branch-attachment experiments
- `CONSTRUCTION_STATUS.md`, `WITNESS_RESULT.md` — development notes tracing
  the bugs found (a sign wrap-around in the aim solver, cue-chase after a
  full first hit, direction drift across bounces) and the decision to
  switch to build-and-verify.
