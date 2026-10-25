"""
witness_bv_exact.py -- the build-and-verify constructor of witness_bv.py,
run against the exact Model M dynamics at 80 significant digits.

Why: a 21-collision thin-cut chain amplifies a lateral error by a factor
of order 10-70 per collision, so a double-precision trajectory carries no
information about the true one after a handful of collisions, and the
engine's 1 nm post-collision separation nudge (EPS_SEP in model_sim.py)
is itself a perturbation of that size. A witness built in double
precision is a witness for the engine's rounding, not for Model M.

Here every carrier state the constructor reads and every prefix it
verifies comes from model_sim_mp.py (no nudge, 80 digits). Placement
heuristics (free-flight look-ahead, aiming) still run in floats: they only
propose positions, which then become exact binary coordinates of the
configuration, and the 80-digit verification decides.

The output is then certified by proof/certify_mpmath.py and proof/certify_arb.py in interval arithmetic.

Run:  python construction/witness_bv_exact.py SEED [OUT.json]
(seed 7 produced proof/witness_21.json in about 6 minutes)
"""
import json
import sys
import time

import os as _os
import sys as _sys

_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(_ROOT, "model"))
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))

import model_sim_mp as MP
import witness_bv as WB


def _sim(cue0, v0, phi, plan):
    balls = [MP.Ball(cue0[0], cue0[1], "cue")]
    for i, (p, _) in enumerate(plan):
        balls.append(MP.Ball(p[0], p[1], f"b{i}"))
    MP.shoot(balls[0], v0, phi)
    return balls, MP.Simulator(balls, record_events=True)


def _carrier_after(cue0, v0, phi, plan, k):
    """As witness_bv._carrier_after, in 80 digits; returns floats."""
    balls, sim = _sim(cue0, v0, phi, plan)
    cue = balls[0]
    seen = 0
    for _ in range(4000):
        if not sim.step():
            break
        ev = sim.events[-1]
        if ev[1] == "bb":
            if seen == k:
                if cue.active and cue._s > MP.V_MIN:
                    return ("cue", float(cue.x), float(cue.y),
                            float(cue.vx), float(cue.vy), float(sim.t))
                return None
            seen += 1
    return None


def _verify_prefix(cue0, v0, phi, plan, npots_expected):
    """As witness_bv._verify_prefix, in 80 digits, stopping as soon as the
    outcome is decided (all placed balls potted, or the prefix broken)."""
    balls, sim = _sim(cue0, v0, phi, plan)
    n = len(plan)
    intended = [tuple(sorted(("cue", f"b{i}"))) for i in range(n)]
    nbb = 0
    potted = set()
    for _ in range(5000):
        if not sim.step():
            break
        ev = sim.events[-1]
        if ev[1] == "bb":
            if nbb < n and tuple(sorted(ev[2])) != intended[nbb]:
                return False
            nbb += 1
        elif ev[1] == "pocket":
            if ev[2] == "cue":
                return False
            potted.add(ev[2])
            if len(potted) == n and nbb >= n:
                return True
        if not any(b.active and b._s >= MP.V_MIN for b in balls):
            break
    return False


WB._carrier_after = _carrier_after
WB._verify_prefix = _verify_prefix


if __name__ == "__main__":
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    out = sys.argv[2] if len(sys.argv) > 2 else f"chain_seed{seed}.json"
    t0 = time.time()
    w = WB.build_greedy(seed=seed, v0=11.9, n_target=21)
    n = w["n"] if w else 0
    print(f"seed {seed}: chain of {n} in {time.time() - t0:.0f}s", flush=True)
    if w:
        w["cue"] = list(w["cue"])
        w["balls"] = [list(p) for p in w["balls"]]
        w["engine"] = "model_sim_mp (80 digits, EPS_SEP = 0)"
        json.dump(w, open(out, "w"))
