"""
certify_miss.py -- the far side of the direction window.

Theorems A and B certify that every stroke with direction within 1e-90
of phi* (and speed in the ribbon) follows the witness's event list and
pockets all 21 object balls. This script certifies the opposite, just
outside the window: for every direction phi in a given box, NO stroke
(V0, phi), whatever its speed, follows the witness's event list.

Why this is a finite check. In the witness every collision is a moving
ball striking a ball at rest, and such collisions, as well as cushion
rebounds, map directions to directions independently of speed. Hence
if a stroke (V0, phi) follows the event list, each ball's path is the
polyline ("skeleton") determined by phi alone: the cue leaves along
phi, meets the first designated ball, both leave along directions fixed
by the contact geometry, and so on. The skeleton is built here in
interval arithmetic over the phi box. If, for some piece, the designated
end event provably cannot be the ball's next event (the ray misses the
designated pocket disk or ball, or it reaches a cushion line or another
pocket disk first), then no stroke in the box follows the event list.

The designated event list (which ball is struck when, which cushion,
which pocket) is the one certified for the box of Theorem B; it is
re-derived here from the geometric construction over phi* +/- 1e-90.

Modes
  python certify_miss.py scan                bracket both edges of the
                                             window by bisection
  python certify_miss.py miss  DELTA HW      certify a miss on the box
                                             phi* + DELTA +/- HW
  python certify_miss.py point DELTA [DV]    certify that the exact
                                             direction phi* + DELTA, with
                                             V0 in V0* +/- DV, follows the
                                             event list and pots all 21
                                             (the full strip check at a
                                             single direction)

DELTA and HW are decimal strings; the exact binary value used is printed.
"""
import json
import os
import sys
import time

from mpmath import iv, mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import certify_strip as S  # noqa: E402  (shares constants, helpers, precision)

iv.prec = S.PREC
mp.prec = S.PREC
R, RHO = S.R, S.RHO
Fail, pos, lo, hi = S.Fail, S.pos, S.lo, S.hi


class Miss(Exception):
    """Raised with a proof that the designated event cannot happen."""


def exact(x):
    """Exact description of an mpf: mantissa * 2^exponent."""
    x = mp.mpf(x)
    return f"{x.man}*2^{x.exp}" if x.man else "0"


# ----------------------------------------------------------- targets
def targets(w, plan, rphi="1e-90"):
    """Attach to each planned wall/pocket event the cushion or pocket met,
    by the geometric construction over the certified box phi* +/- rphi."""
    PHI = iv.mpf(w["phi"]) + iv.mpf([-hi(iv.mpf(rphi)), hi(iv.mpf(rphi))])
    out = []
    try:
        geometry(w, PHI, plan, record=out)
    except Miss as exc:
        raise Fail(f"certified box reported a miss: {exc}")
    return out


# ----------------------------------------------------------- geometry
def ray_misses(P, u, Q, rad):
    """Return a reason string if the ray P + s u, s > 0, provably never
    reaches distance rad from Q; otherwise None."""
    d = (P[0] - Q[0], P[1] - Q[1])
    b = d[0] * u[0] + d[1] * u[1]
    c = d[0] ** 2 + d[1] ** 2 - rad * rad
    if not pos(c):
        raise Fail("start inside a contact disk")
    if pos(b):
        return (f"moving away from the target (start distance >= "
                f"{mp.nstr(lo(iv.sqrt(c + rad * rad)), 8)} m, radius "
                f"{mp.nstr(hi(rad), 6)} m)")
    disc = b * b - c
    if (disc < 0) is True:
        dist = iv.sqrt(d[0] ** 2 + d[1] ** 2 - b * b)
        return (f"the line of motion passes the target at distance >= "
                f"{mp.nstr(lo(dist), 12)} m > radius {mp.nstr(hi(rad), 6)} m "
                f"(margin >= {mp.nstr(lo(dist - rad), 6)} m)")
    return None


def decide(ball, ev_i, label, s_target, competitors):
    """s_target: interval arc length at which the designated event occurs
    (None if it provably never occurs). competitors: (name, d) with d the
    arc length at which a cushion line or pocket disk is reached (None if
    never). Returns normally if the designated event is provably first,
    raises Miss if it provably is not, raises Fail if undecided."""
    if s_target is None:
        raise Miss(f"event {ev_i}: {ball} never reaches {label}")
    for name, d in competitors:
        if d is None:
            continue
        if pos(iv.mpf(lo(s_target)) - d):
            raise Miss(f"event {ev_i}: {ball} reaches {name} after arc "
                       f"length <= {mp.nstr(hi(d), 10)} m, before {label} "
                       f"(arc length >= {mp.nstr(lo(s_target), 10)} m)")
        if not pos(d - s_target):
            raise Fail(f"event {ev_i}: {ball}: {name} vs {label} undecided")


def geometry(w, PHI, plan, record=None):
    """Build the skeleton over the direction box PHI following the planned
    events. If `record` is a list, the plan carries no targets yet: the
    nearest cushion/pocket is taken and appended to `record`. Otherwise
    each wall/pocket event of `plan` names its target. Returns normally
    if every designated event is provably the ball's next event (geometry
    only; timing between different balls is not examined here). Raises
    Miss with a proof that some designated event cannot happen, or Fail
    if the box is too wide to decide."""
    n = len(w["balls"])
    names = ["cue"] + [f"b{i}" for i in range(n)]
    P = {nm: (iv.mpf(p[0]), iv.mpf(p[1]))
         for nm, p in zip(names, [w["cue"]] + w["balls"])}
    U = {nm: None for nm in names}
    U["cue"] = (iv.cos(PHI), iv.sin(PHI))
    walls = {wl[0]: wl for wl in S.WALLS}
    for ev_i, ev in enumerate(plan, start=1):
        kind = ev[0]
        if kind == "bb":
            a_, b_ = ev[1], ev[2]
            mover, rest = (a_, b_) if U[a_] is not None else (b_, a_)
            if U[mover] is None or U[rest] is not None:
                raise Fail(f"event {ev_i}: not a moving-onto-resting collision")
            p, u, q = P[mover], U[mover], P[rest]
            reason = ray_misses(p, u, q, 2 * R)
            if reason is not None:
                raise Miss(f"event {ev_i}: {mover} misses {rest}: {reason}")
            s = S.first_contact(p, u, q, 2 * R)
            comps = [(f"cushion {wl[0]}", S.wall_distance(p, u, wl))
                     for wl in S.WALLS]
            comps += [(f"pocket {S.POCKET_NAMES[k]}",
                       S.first_contact(p, u, Q, RHO))
                      for k, Q in enumerate(S.POCKETS)]
            decide(mover, ev_i, f"ball {rest}", s, comps)
            C = (p[0] + s * u[0], p[1] + s * u[1])
            nx, ny = (q[0] - C[0]) / (2 * R), (q[1] - C[1]) / (2 * R)
            un = u[0] * nx + u[1] * ny
            if not pos(un):
                raise Fail(f"event {ev_i}: approach undecided")
            k = (1 + S.E) / 2 * un
            wx, wy = u[0] - k * nx, u[1] - k * ny
            f = iv.sqrt(wx * wx + wy * wy)
            if not pos(f):
                raise Fail(f"event {ev_i}: mover may stop at contact")
            P[mover], U[mover] = C, (wx / f, wy / f)
            U[rest] = (nx, ny)
            if record is not None:
                record.append(ev)
        elif kind == "wall":
            ball = ev[1]
            p, u = P[ball], U[ball]
            if u is None:
                raise Fail(f"event {ev_i}: {ball} is at rest")
            dists = {wl[0]: S.wall_distance(p, u, wl) for wl in S.WALLS}
            if record is not None:
                cands = [(hi(d), nm) for nm, d in dists.items() if d is not None]
                if not cands:
                    raise Fail(f"event {ev_i}: no cushion ahead")
                target = min(cands)[1]
                record.append(("wall", ball, target))
            else:
                target = ev[2]
            s = dists[target]
            comps = [(f"cushion {nm}", d) for nm, d in dists.items()
                     if nm != target]
            comps += [(f"pocket {S.POCKET_NAMES[k]}",
                       S.first_contact(p, u, Q, RHO))
                      for k, Q in enumerate(S.POCKETS)]
            decide(ball, ev_i, f"cushion {target}", s, comps)
            wl = walls[target]
            C = [p[0] + s * u[0], p[1] + s * u[1]]
            C[wl[1]] = wl[2]
            ux, uy = u
            if wl[1] == 0:
                ux = -S.EC * ux
            else:
                uy = -S.EC * uy
            f = iv.sqrt(ux * ux + uy * uy)
            P[ball], U[ball] = tuple(C), (ux / f, uy / f)
        else:                                                # pocket
            ball = ev[1]
            p, u = P[ball], U[ball]
            if u is None:
                raise Fail(f"event {ev_i}: {ball} is at rest")
            dists = [S.first_contact(p, u, Q, RHO) for Q in S.POCKETS]
            if record is not None:
                cands = [(hi(d), k) for k, d in enumerate(dists) if d is not None]
                if not cands:
                    raise Fail(f"event {ev_i}: no pocket ahead")
                target = min(cands)[1]
                record.append(("pocket", ball, target))
            else:
                target = ev[2]
            reason = ray_misses(p, u, S.POCKETS[target], RHO)
            if reason is not None:
                raise Miss(f"event {ev_i}: {ball} misses pocket "
                           f"{S.POCKET_NAMES[target]}: {reason}")
            s = dists[target]
            comps = [(f"cushion {wl[0]}", S.wall_distance(p, u, wl))
                     for wl in S.WALLS]
            comps += [(f"pocket {S.POCKET_NAMES[k]}", d)
                      for k, d in enumerate(dists) if k != target]
            decide(ball, ev_i, f"pocket {S.POCKET_NAMES[target]}", s, comps)
            P[ball], U[ball] = None, None
            if ball == "cue":
                raise Fail("cue ball pocketed")


# ----------------------------------------------------------- drivers
def status(w, plan_t, delta):
    """'OK', 'MISS' or 'UNDECIDED' for the exact direction phi* + delta."""
    PHI = iv.mpf(mp.mpf(w["phi"]) + delta)
    try:
        geometry(w, PHI, plan_t)
        return "OK", ""
    except Miss as exc:
        return "MISS", str(exc)
    except Fail as exc:
        return "UNDECIDED", str(exc)


def scan(w, plan_t):
    """Bracket, on each side of phi*, the last direction at which the
    skeleton is geometrically realised and the first at which a
    designated event provably fails."""
    out = {}
    for sign, name in ((1, "+"), (-1, "-")):
        inner = mp.mpf("1e-90")
        st, _ = status(w, plan_t, sign * inner)
        if st != "OK":
            raise Fail(f"side {name}: inner point 1e-90 is {st}")
        outer = mp.mpf("1e-40")
        while True:
            st, why = status(w, plan_t, sign * outer)
            if st == "MISS":
                outer_why = why
                break
            if st == "UNDECIDED":
                raise Fail(f"side {name}: {outer} undecided: {why}")
            outer *= 10
            if outer > 1:
                raise Fail(f"side {name}: no miss found")
        for _ in range(400):
            mid = (inner + outer) / 2
            st, why = status(w, plan_t, sign * mid)
            if st == "OK":
                inner = mid
            elif st == "MISS":
                outer, outer_why = mid, why
            else:
                print(f"  side {name}: undecided at {mp.nstr(mid, 20)}: {why}")
                break
            if outer - inner < outer * mp.mpf("1e-30"):
                break
        print(f"side {name}: skeleton realised at      phi* {name} {mp.nstr(inner, 25)}")
        print(f"side {name}: designated event fails at phi* {name} {mp.nstr(outer, 25)}")
        print(f"          exact: {exact(inner)} and {exact(outer)}")
        print(f"          reason: {outer_why}")
        out[name] = (inner, outer, outer_why)
    return out


def certify_miss(w, plan_t, delta, hw):
    d, h = mp.mpf(delta), mp.mpf(hw)
    centre = mp.mpf(w["phi"]) + d
    PHI = iv.mpf([centre - h, centre + h])
    print(f"direction box: phi* + {mp.nstr(d, 25)} +/- {mp.nstr(h, 6)}")
    print(f"  exact centre offset {exact(d)}, half-width {exact(h)}")
    print(f"  decimal enclosure of the offset: "
          f"[{mp.nstr(lo(iv.mpf(d) - h), 25)}, {mp.nstr(hi(iv.mpf(d) + h), 25)}]")
    print(f"  precision {S.PREC} bits")
    try:
        geometry(w, PHI, plan_t)
    except Miss as exc:
        print(f"  proof: {exc}")
        print("\nCERTIFIED MISS: no stroke with direction in this box, at any "
              "speed, follows the witness's event list.")
        return 0
    except Fail as exc:
        print(f"NOT CERTIFIED (undecided): {exc}")
        return 1
    print("NOT CERTIFIED: the skeleton is realised on this box")
    return 1


def certify_point(w, plan, plan_t, delta, dV):
    d, dv = mp.mpf(delta), mp.mpf(dV)
    centre = mp.mpf(w["phi"]) + d
    PHI = iv.mpf(centre)
    V = iv.mpf(w["v0"]) + iv.mpf([-dv, dv])
    print(f"direction: phi* + {mp.nstr(d, 25)} exactly ({exact(d)}), "
          f"V0 in V0* +/- {dV}, {S.PREC} bits")
    try:
        geometry(w, PHI, plan_t)          # designated targets are met
        pieces, same_start, designed, potted, t_end = S.build(w, V, PHI, plan)
        n_pairs, n_time = S.check_pairs(pieces, same_start, designed, t_end)
    except (Fail, Miss) as exc:
        print(f"NOT CERTIFIED: {exc}")
        return 1
    same = all(("pocket", b, S.POCKET_NAMES.index(pk)) in plan_t
               for b, pk, _, _ in potted)
    print(f"  pockets as designated       : {same}")
    print(f"  piece pairs examined        : {n_pairs} ({n_time} needed timing)")
    print(f"  min pocket entry speed      : "
          f"{mp.nstr(min(lo(p[3]) for p in potted), 6)} m/s")
    print("\nCERTIFIED: every stroke with this direction and a speed in the "
          f"range follows the event list and pockets all {len(potted)} balls.")
    return 0


def main():
    w = json.load(open(os.path.join(S.HERE, "witness_21.json")))
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        return 2
    S.check_admissible(w)
    plan = S.plan_events(w)               # proposal only
    plan_t = targets(w, plan)             # targets fixed on the certified box
    t0 = time.time()
    if args[0] == "scan":
        scan(w, plan_t)
        rc = 0
    elif args[0] == "miss":
        rc = certify_miss(w, plan_t, args[1], args[2])
    elif args[0] == "point":
        rc = certify_point(w, plan, plan_t, args[1],
                           args[2] if len(args) > 2 else "0.01")
    else:
        print(__doc__)
        return 2
    print(f"wall time: {time.time() - t0:.1f} s")
    return rc


if __name__ == "__main__":
    sys.exit(main())
