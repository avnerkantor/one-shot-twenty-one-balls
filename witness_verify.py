"""
witness_verify.py -- independently verify a witness JSON file.

Loads a witness (default witness_21.json), replays the single cue stroke
in the event-driven engine of model_sim.py, and checks that:

  * the cue strikes every placed object ball, in the recorded order;
  * all object balls are pocketed and the cue is not;
  * the configuration is admissible (no two balls overlap);

then reports the margins that enter the regularity hypothesis of the
existence theorem (minimum pairwise ball separation, pocket-entry
speeds, minimum event-time gap).

Run:  python3 witness_verify.py [witness_21.json]
Exit code 0 on a fully verified clearance, 1 otherwise.
"""
import json
import math
import sys

import model_sim as M

R = M.R


def verify(path="witness_21.json"):
    w = json.load(open(path))
    cue0, v0, phi = w["cue"], w["v0"], w["phi"]
    ball_xy = w["balls"]
    pockets = w["pockets"]
    n = len(ball_xy)

    # admissibility: no two balls (or ball and cue) overlap
    pts = [tuple(cue0)] + [tuple(p) for p in ball_xy]
    min_sep = min(math.dist(pts[i], pts[j])
                  for i in range(len(pts)) for j in range(i + 1, len(pts)))
    admissible = min_sep > 2 * R

    # replay
    balls = [M.Ball(cue0[0], cue0[1], "cue")]
    for i, p in enumerate(ball_xy):
        balls.append(M.Ball(p[0], p[1], f"b{i}"))
    M.shoot(balls[0], v0, phi)
    sim = M.Simulator(balls, record_events=True)
    sim.run(max_events=8000, max_time=800)

    # collision order: cue must strike b0, b1, ... as a prefix
    bb = [tuple(sorted(e[2])) for e in sim.events if e[1] == "bb"]
    intended = [tuple(sorted(("cue", f"b{i}"))) for i in range(n)]
    order_ok = (len(bb) >= n
                and all(a == b for a, b in zip(bb[:n], intended)))

    potted = {nm for _, nm, _, _ in sim.pocketed}
    all_potted = all(f"b{i}" in potted for i in range(n))
    cue_potted = "cue" in potted

    entry_speeds = sorted(s for _, nm, _, s in sim.pocketed if nm != "cue")

    # minimum gap between consecutive event times
    times = sorted(e[0] for e in sim.events)
    gaps = [b - a for a, b in zip(times, times[1:]) if b - a > 1e-12]
    min_gap = min(gaps) if gaps else float("nan")

    n_potted = sum(1 for nm in potted if nm != "cue")

    ok = (admissible and order_ok and all_potted and not cue_potted
          and n_potted == n)

    print(f"witness file        : {path}")
    print(f"object balls        : {n}")
    print(f"cue shot            : V0={v0} m/s, phi={phi:.6f} rad, "
          f"cue at ({cue0[0]:.4f}, {cue0[1]:.4f})")
    print(f"admissible (no overlap): {admissible}  "
          f"(min separation {min_sep * 1000:.1f} mm)")
    print(f"collision order OK  : {order_ok}")
    print(f"object balls potted : {n_potted}/{n}")
    print(f"cue potted          : {cue_potted}")
    if entry_speeds:
        print(f"pocket entry speed  : min {min(entry_speeds):.3f}, "
              f"max {max(entry_speeds):.3f} m/s")
    print(f"min event-time gap  : {min_gap * 1000:.2f} ms")
    print()
    print("VERIFIED: full 21-ball single-stroke clearance"
          if ok else "NOT verified")
    return ok


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "witness_21.json"
    sys.exit(0 if verify(src) else 1)
