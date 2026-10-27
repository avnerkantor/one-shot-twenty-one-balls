"""
mc_random_layouts.py -- P(k) for random (undesigned) layouts.

Companion to mc_break.py. Each trial places the cue ball and 21 object
balls uniformly at random on the table (rejecting overlaps and centres
inside a pocket disk), draws a stroke uniformly from
(0, V_MAX] x [0, 2 pi), and counts the object balls pocketed. Comparing
the resulting P(k) with the break's asks whether undesigned layouts are
any friendlier than the opening rack.

Run:  python simulation/mc_random_layouts.py N SEED OUT.json
"""
import json
import math
import os as _os
import random
import sys as _sys

_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(_ROOT, "model"))

import model_sim as M

V_MAX = 12.0


def random_layout(rng, n=22):
    pts = []
    while len(pts) < n:
        x = rng.uniform(M.R, M.L - M.R)
        y = rng.uniform(M.R, M.W - M.R)
        if any(math.hypot(x - px, y - py) <= M.RHO for px, py in M.POCKETS):
            continue
        if any(math.hypot(x - qx, y - qy) <= 2 * M.R for qx, qy in pts):
            continue
        pts.append((x, y))
    balls = [M.Ball(pts[0][0], pts[0][1], "cue")]
    balls += [M.Ball(x, y, f"b{i}") for i, (x, y) in enumerate(pts[1:])]
    return balls


def main():
    n = int(_sys.argv[1])
    seed = int(_sys.argv[2])
    out = _sys.argv[3]
    rng = random.Random(seed)
    counts = {}
    for _ in range(n):
        balls = random_layout(rng)
        M.shoot(balls[0], rng.uniform(1e-6, V_MAX), rng.uniform(0, 2 * math.pi))
        sim = M.Simulator(balls)
        sim.run()
        k = sum(1 for _, nm, _, _ in sim.pocketed if nm != "cue")
        counts[k] = counts.get(k, 0) + 1
    json.dump({"n": n, "seed": seed, "counts": counts}, open(out, "w"))


if __name__ == "__main__":
    main()
