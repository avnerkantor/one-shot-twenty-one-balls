"""
tree_measure.py -- replay a tree-shaped clearance (construction/witness_tree.py
output) in the 80-digit engine and bracket its window of directions there.

The double-precision search engine nudges balls apart by 1e-9 m after each
contact; the 80-digit engine (model_sim_mp.py) does not, so the replay here
decides whether the configuration is a clearance of Model M itself (up to
rounding at 80 digits), and the bisection gives the window of directions
that a certificate would then have to confirm.

Run:  python construction/tree_measure.py TREE.json [ITERS]
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "model"))

import model_sim_mp as MP                  # noqa: E402
import witness_tree as T                   # noqa: E402


def main():
    w = json.load(open(sys.argv[1]))
    iters = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    plan = [{"pos": tuple(p), "parent": par} for p, par in zip(w["balls"], w["parents"])]
    cue0, v0, phi = tuple(w["cue"]), w["v0"], w["phi"]
    t0 = time.time()
    ok = T.check(cue0, v0, phi, plan, engine=MP)
    print(f"80-digit replay at the centre: {'realises the tree' if ok else 'FAILS: ' + T.REASON[0]} "
          f"({time.time() - t0:.1f} s)")
    if not ok:
        sys.exit(1)
    balls, sim = T.simulate(cue0, v0, phi, plan, engine=MP, max_time=60.0)
    print(f"events {len(sim.events)}, last pocketing at t = {max(float(p[0]) for p in sim.pocketed):.3f} s, "
          f"min pocket entry speed {min(float(p[3]) for p in sim.pocketed):.3f} m/s")
    depth = {f"b{i}": d for i, d in enumerate(w.get("depths", []))}
    if depth:
        print("depths:", sorted(depth.values()))
    # window of directions by bisection in the 80-digit engine
    out = []
    for sgn in (+1, -1):
        a, b = 1e-12, 0.05
        if not T.check(cue0, v0, phi + sgn * a, plan, engine=MP):
            out.append(0.0)
            continue
        for _ in range(iters):
            m = math.sqrt(a * b)
            if T.check(cue0, v0, phi + sgn * m, plan, engine=MP):
                a = m
            else:
                b = m
            if b / a < 1.05:
                break
        out.append(a)
    plus, minus = out
    print(f"window (80-digit engine): phi + [-{minus:.3e}, +{plus:.3e}] rad, width {plus + minus:.3e} rad "
          f"({time.time() - t0:.0f} s)")
    # speed range at the centre direction
    lo_v = None
    for dv in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0):
        if v0 - dv > 0 and T.check(cue0, v0 - dv, phi, plan, engine=MP):
            lo_v = dv
        else:
            break
    print(f"speed: works down to at least V0 - {lo_v} m/s" if lo_v else "speed: fails already at V0 - 0.5 m/s")
    w["window_mp"] = [minus, plus]
    json.dump(w, open(sys.argv[1], "w"), indent=1)


if __name__ == "__main__":
    main()
