"""
certify_strip.py -- certify a strip of strokes that is wide in speed.

Claim: every stroke with V0 in [V0* - dV, V0* + dV] and phi in
[phi* - dphi, phi* + dphi] pockets all 21 object balls of
proof/witness_21.json in Model M.

Why a different method from certify_mpmath.py / certify_arb.py:
in this witness every collision is a moving ball striking a ball at
rest. Between events balls move on straight lines, and the outgoing
directions of a moving-onto-resting collision and of a cushion
rebound depend only on the incoming direction and the geometry, not
on the speed. Hence, as long as the same events happen, every ball's
path is the same polyline for every V0; the speed only changes WHEN
each ball passes each point. (Speeds are increasing functions of V0:
collisions and cushions scale speed linearly, free flight maps v to
sqrt(v^2 - 2 a s).)

The certificate therefore has two parts.

1. Skeleton (geometry). In interval arithmetic over the phi box, build
   every ball's polyline from the planned sequence of events and prove
   that along each straight piece no other cushion or pocket is met
   first, that a cue piece ending in a collision really meets the
   struck ball at the first root of the contact equation, and that every
   ball still moves at the end of every piece (checked over the whole
   speed range). Along the way, enclose for the whole speed range the time at
   which each piece starts and ends.

2. Pairs (timing). For every pair of pieces belonging to two different
   balls, prove they are never within 2R of each other at the same
   instant, for every stroke in the strip. Pieces that are contiguous at
   a common event (the collision itself) are handled exactly: the pair
   that just collided is shown to separate (the gap is t * g(t) with
   g > 0), and a ball at rest is never paired with a piece that only
   starts when it is struck.

If both parts succeed, induction over the (V0-dependent) order of events
shows the evolution of every stroke in the strip follows the skeleton,
so all 21 object balls are pocketed. Proposals (which events happen, in
which order) are taken from an 80-digit replay at the centre stroke;
they are only a plan and the proof checks everything.

Run:  python proof/certify_strip.py [dV] [dphi]      e.g. 0.01 1e-90
      python proof/certify_strip.py a b dphi         V0 in V0* + [a, b]
"""
import json
import os
import sys
import time

from mpmath import iv, mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "model"))

PREC = int(os.environ.get("CERT_PREC", "1200"))
iv.prec = PREC
mp.prec = PREC

L = iv.mpf("3.569")
W = iv.mpf("1.778")
R = iv.mpf("0.02625")
A = iv.mpf("0.10")
E = iv.mpf("0.95")
EC = iv.mpf("0.90")
RHO = iv.mpf("0.055")
Z = iv.mpf(0)
POCKETS = ((Z, Z), (L, Z), (Z, W), (L, W), (L / 2, Z), (L / 2, W))
POCKET_NAMES = ("C-bl", "C-br", "C-tl", "C-tr", "M-b", "M-t")
WALLS = (("left", 0, R), ("right", 0, L - R), ("bottom", 1, R), ("top", 1, W - R))


class Fail(Exception):
    pass


def pos(x):
    return (x > 0) is True


def lo(x):
    return mp.mpf(x.a)


def hi(x):
    return mp.mpf(x.b)


def hull(a, b):
    return iv.mpf([min(lo(a), lo(b)), max(hi(a), hi(b))])


# ---------------------------------------------------------------- plan
def plan_events(w):
    """Event list of the centre stroke from the 80-digit engine
    (heuristic input only)."""
    import model_sim_mp as M
    balls = [M.Ball(*w["cue"], "cue")] + [M.Ball(p[0], p[1], f"b{i}")
                                         for i, p in enumerate(w["balls"])]
    M.shoot(balls[0], w["v0"], w["phi"])
    s = M.Simulator(balls, record_events=True)
    n_obj = len(w["balls"])
    out, potted = [], 0
    while potted < n_obj and s.step():
        t, kind, who, extra = s.events[-1]
        if kind == "bb":
            out.append(("bb",) + tuple(who))
        elif kind == "cushion":
            out.append(("wall", who))
        elif kind == "pocket":
            out.append(("pocket", who))
            if who != "cue":
                potted += 1
    return out


# ---------------------------------------------------------------- pieces
class Piece:
    """A ball either at rest at P during [0, t1], or moving from P along
    the unit vector u for arc length ell, starting at time t0 with speed
    v0 and ending at time t1 (t0, t1, v0 are intervals over the strip)."""
    __slots__ = ("ball", "kind", "P", "u", "ell", "t0", "t1", "v0", "v1",
                 "e0", "e1")

    def __init__(self, ball, kind, P, t0, e0):
        self.ball, self.kind, self.P, self.t0, self.e0 = ball, kind, P, t0, e0
        self.u = self.ell = self.t1 = self.v0 = self.v1 = self.e1 = None


def travel(v, s):
    """Speed after free flight over arc length s (must stay positive)."""
    arg = v * v - 2 * A * s
    if not pos(arg):
        raise Fail("a ball may stop before its next event")
    return iv.sqrt(arg)


def first_contact(P, u, Q, rad):
    """Smallest s > 0 with |P + s u - Q| = rad, proved to be a simple
    (transversal) first root; None if provably no root ahead."""
    d = (P[0] - Q[0], P[1] - Q[1])
    b = d[0] * u[0] + d[1] * u[1]
    c = d[0] ** 2 + d[1] ** 2 - rad * rad
    if not pos(c):
        raise Fail("start inside a contact disk")
    if pos(b):                      # moving away: no root ahead
        return None
    disc = b * b - c
    if (disc < 0) is True:
        return None
    if not pos(disc):
        raise Fail("tangential contact cannot be excluded")
    if not pos(-b):
        raise Fail("approach direction undecided")
    return -b - iv.sqrt(disc)       # smaller root, > 0 since c > 0, b < 0


def wall_distance(P, u, wall):
    name, axis, val = wall
    p, du = P[axis], u[axis]
    if name in ("left", "bottom"):
        if pos(-du):
            return (p - val) / (-du)
    else:
        if pos(du):
            return (val - p) / du
    if pos(du) or pos(-du):
        return None                 # moving away from this wall
    raise Fail("direction parallel to a wall undecided")


def end_piece(pc, s, event):
    pc.ell = s
    pc.v1 = travel(pc.v0, s)
    # (v0 - v1)/a rewritten without cancellation: v0 and v1 both grow
    # with V0, and their interval difference would lose that
    pc.t1 = pc.t0 + 2 * s / (pc.v0 + pc.v1)
    pc.e1 = event
    return (pc.P[0] + s * pc.u[0], pc.P[1] + s * pc.u[1])


def check_single(pc, s, ignore_wall=None, ignore_pocket=None):
    """No cushion or pocket is met on the piece before arc length s."""
    for wall in WALLS:
        if wall[0] == ignore_wall:
            continue
        d = wall_distance(pc.P, pc.u, wall)
        if d is not None and not pos(d - s):
            raise Fail(f"{pc.ball}: cushion {wall[0]} may come first")
    for k, Q in enumerate(POCKETS):
        if k == ignore_pocket:
            continue
        d = first_contact(pc.P, pc.u, Q, RHO)
        if d is not None and not pos(d - s):
            raise Fail(f"{pc.ball}: pocket {POCKET_NAMES[k]} may come first")


def build(w, V, PHI, plan):
    n = len(w["balls"])
    names = ["cue"] + [f"b{i}" for i in range(n)]
    pieces = {nm: [] for nm in names}
    cur = {}
    for nm, p in zip(names, [w["cue"]] + w["balls"]):
        pc = Piece(nm, "rest", (iv.mpf(p[0]), iv.mpf(p[1])), Z, 0)
        pieces[nm].append(pc)
        cur[nm] = pc
    cue = cur["cue"]
    cue.kind, cue.u, cue.v0 = "move", (iv.cos(PHI), iv.sin(PHI)), V
    same_start, designed = [], []
    potted = []
    for ev_i, ev in enumerate(plan, start=1):
        if ev[0] == "bb":
            a_, b_ = ev[1], ev[2]
            mover, rest = (a_, b_) if cur[a_].kind == "move" else (b_, a_)
            pm, pr = cur[mover], cur[rest]
            if pm.kind != "move" or pr.kind != "rest":
                raise Fail(f"event {ev_i}: not a moving-onto-resting collision")
            s = first_contact(pm.P, pm.u, pr.P, 2 * R)
            if s is None:
                raise Fail(f"event {ev_i}: planned contact does not occur")
            check_single(pm, s)
            C = end_piece(pm, s, ev_i)
            pr.t1, pr.e1 = pm.t1, ev_i
            designed.append((pm, pr))
            nx, ny = (pr.P[0] - C[0]) / (2 * R), (pr.P[1] - C[1]) / (2 * R)
            un = pm.u[0] * nx + pm.u[1] * ny
            if not pos(un):
                raise Fail(f"event {ev_i}: not approaching")
            k = (1 + E) / 2 * un
            wx, wy = pm.u[0] - k * nx, pm.u[1] - k * ny
            f = iv.sqrt(wx * wx + wy * wy)
            if not pos(f):
                raise Fail(f"event {ev_i}: mover may stop")
            nm_ = Piece(mover, "move", C, pm.t1, ev_i)
            nm_.u, nm_.v0 = (wx / f, wy / f), f * pm.v1
            nr = Piece(rest, "move", pr.P, pm.t1, ev_i)
            nr.u, nr.v0 = (nx, ny), k * pm.v1
            pieces[mover].append(nm_)
            pieces[rest].append(nr)
            cur[mover], cur[rest] = nm_, nr
            same_start.append((nm_, nr))
        elif ev[0] == "wall":
            pc = cur[ev[1]]
            best = None
            for wall in WALLS:
                d = wall_distance(pc.P, pc.u, wall)
                if d is not None and (best is None or hi(d) < hi(best[1])):
                    best = (wall, d)
            if best is None:
                raise Fail(f"event {ev_i}: no cushion ahead")
            wall, s = best
            check_single(pc, s, ignore_wall=wall[0])
            C = list(end_piece(pc, s, ev_i))
            C[wall[1]] = wall[2]                     # exactly on the line
            ux, uy = pc.u
            if wall[1] == 0:
                ux = -EC * ux
            else:
                uy = -EC * uy
            f = iv.sqrt(ux * ux + uy * uy)
            np_ = Piece(pc.ball, "move", tuple(C), pc.t1, ev_i)
            np_.u, np_.v0 = (ux / f, uy / f), f * pc.v1
            pieces[pc.ball].append(np_)
            cur[pc.ball] = np_
        else:                                        # pocket
            pc = cur[ev[1]]
            best = None
            for k, Q in enumerate(POCKETS):
                d = first_contact(pc.P, pc.u, Q, RHO)
                if d is not None and (best is None or hi(d) < hi(best[1])):
                    best = (k, d)
            if best is None:
                raise Fail(f"event {ev_i}: no pocket ahead")
            k, s = best
            check_single(pc, s, ignore_pocket=k)
            end_piece(pc, s, ev_i)
            cur[pc.ball] = None
            if pc.ball == "cue":
                raise Fail("cue ball pocketed")
            potted.append((pc.ball, POCKET_NAMES[k], pc.t1, pc.v1))
    if len(potted) != n:
        raise Fail("plan does not pocket every object ball")
    t_end = max(hi(p[2]) for p in potted)
    # extend the cue (the only ball left) until every stroke has finished
    pc = cur["cue"]
    while True:
        stop_s = pc.v0 * pc.v0 / (2 * A)          # arc length until it stops
        best = None
        for wall in WALLS:
            d = wall_distance(pc.P, pc.u, wall)
            if d is not None and (best is None or hi(d) < hi(best[1])):
                best = (wall, d)
        wall, s = best
        if pos(s - stop_s):                       # comes to rest first
            check_single(pc, stop_s)
            pc.ell, pc.v1, pc.e1 = stop_s, Z, "x"
            pc.t1 = pc.t0 + pc.v0 / A
            if lo(pc.t1) <= t_end:
                P = (pc.P[0] + stop_s * pc.u[0], pc.P[1] + stop_s * pc.u[1])
                pieces["cue"].append(Piece("cue", "rest", P, pc.t1, "x"))
            break
        if not pos(stop_s - s):
            raise Fail("cue may stop exactly at a cushion")
        check_single(pc, s, ignore_wall=wall[0])
        C = list(end_piece(pc, s, "x"))
        if lo(pc.t1) > t_end:
            break
        C[wall[1]] = wall[2]
        ux, uy = pc.u
        if wall[1] == 0:
            ux = -EC * ux
        else:
            uy = -EC * uy
        f = iv.sqrt(ux * ux + uy * uy)
        np_ = Piece("cue", "move", tuple(C), pc.t1, "x")
        np_.u, np_.v0 = (ux / f, uy / f), f * pc.v1
        pieces["cue"].append(np_)
        pc = np_
    return pieces, same_start, designed, potted, t_end


# ---------------------------------------------------------------- pairs
def where(pc, T):
    """Enclosure of the ball's position at times T (interval), assuming
    T lies in the piece's time window."""
    if pc.kind == "rest":
        return pc.P
    tau = T - pc.t0
    tlo, thi = max(lo(tau), 0), max(hi(tau), 0)
    v0lo, v0hi = lo(pc.v0), hi(pc.v0)
    # arc length s(v0, tau) = v0 tau - a tau^2/2 up to the stop at
    # tau = v0/a, constant v0^2/(2a) after; nondecreasing in tau and v0.
    # Evaluate the quadratic only on its increasing branch.
    tlo = min(tlo, lo(iv.mpf(v0lo) / A))       # lower bound: stays valid
    slo = iv.mpf(v0lo) * tlo - A * iv.mpf(tlo) ** 2 / 2
    if thi >= lo(iv.mpf(v0hi) / A):            # may have reached its stop
        shi = iv.mpf(v0hi) ** 2 / (2 * A)
    else:
        shi = iv.mpf(v0hi) * thi - A * iv.mpf(thi) ** 2 / 2
    a_, b_ = max(lo(slo), 0), min(hi(shi), hi(pc.ell))
    if a_ > b_:
        return None     # no stroke has this ball on this piece at times T
    s = iv.mpf([a_, b_])
    return (pc.P[0] + s * pc.u[0], pc.P[1] + s * pc.u[1])


def window(pc, t_end):
    a = 0 if pc.e0 == 0 and pc.kind == "rest" else lo(pc.t0)
    b = hi(pc.t1) if pc.t1 is not None else t_end
    return a, min(b, t_end)


def apart(pa, pb, a, b, depth=0, budget=None):
    """Prove the two pieces are > 2R apart at every time in [a, b]."""
    stack = [(a, b, 0)]
    while stack:
        a_, b_, d = stack.pop()
        budget[0] -= 1
        if budget[0] < 0:
            return False
        T = iv.mpf([a_, b_])
        P, Q = where(pa, T), where(pb, T)
        if P is None or Q is None:
            continue
        dist2 = (P[0] - Q[0]) ** 2 + (P[1] - Q[1]) ** 2
        if pos(dist2 - 4 * R * R):
            continue
        if d > 48:
            return False
        m = (a_ + b_) / 2
        stack += [(a_, m, d + 1), (m, b_, d + 1)]
    return True


def separating(pa, pb):
    """The pair that has just collided separates: both pieces start at the
    same instant with centres exactly 2R apart. In local time tau the gap
    is f = |r0 + w1 tau + w2 tau^2|^2 - 4R^2 with f(0) = 0; prove f/tau > 0
    up to the end of the shorter piece."""
    r0 = (pa.P[0] - pb.P[0], pa.P[1] - pb.P[1])
    w1 = (pa.v0 * pa.u[0] - pb.v0 * pb.u[0], pa.v0 * pa.u[1] - pb.v0 * pb.u[1])
    w2 = (-A / 2 * (pa.u[0] - pb.u[0]), -A / 2 * (pa.u[1] - pb.u[1]))
    dot = lambda x, y: x[0] * y[0] + x[1] * y[1]
    g = (2 * dot(r0, w1), dot(w1, w1) + 2 * dot(r0, w2), 2 * dot(w1, w2),
         dot(w2, w2))
    tmax = min(hi(pa.t1 - pa.t0), hi(pb.t1 - pb.t0))
    stack = [(mp.mpf(0), tmax, 0)]
    while stack:
        a_, b_, d = stack.pop()
        T = iv.mpf([a_, b_])
        val = g[0] + T * (g[1] + T * (g[2] + T * g[3]))
        if pos(val):
            continue
        if d > 60:
            return False
        m = (a_ + b_) / 2
        stack += [(a_, m, d + 1), (m, b_, d + 1)]
    return True


def bbox(pc):
    if pc.kind == "rest":
        return (lo(pc.P[0]), hi(pc.P[0]), lo(pc.P[1]), hi(pc.P[1]))
    ex, ey = pc.P[0] + pc.ell * pc.u[0], pc.P[1] + pc.ell * pc.u[1]
    return (min(lo(pc.P[0]), lo(ex)), max(hi(pc.P[0]), hi(ex)),
            min(lo(pc.P[1]), lo(ey)), max(hi(pc.P[1]), hi(ey)))


def check_pairs(pieces, same_start, designed, t_end):
    skip = {(id(a), id(b)) for a, b in designed + same_start}
    skip |= {(id(b), id(a)) for a, b in designed + same_start}
    for a, b in same_start:
        if not separating(a, b):
            raise Fail(f"{a.ball}/{b.ball} may not separate after event {a.e0}")
    allp = [p for ps in pieces.values() for p in ps]
    n_pairs = n_time = 0
    for i in range(len(allp)):
        pa = allp[i]
        for pb in allp[i + 1:]:
            if pa.ball == pb.ball or (id(pa), id(pb)) in skip:
                continue
            if (pa.kind == "rest" and pb.kind == "rest"
                    and pa.e0 == 0 and pb.e0 == 0):
                continue               # both untouched: initial admissibility
            # one piece ends at the event where the other begins: the two
            # never coexist (same stroke, same instant)
            if pa.e1 is not None and pa.e1 == pb.e0:
                continue
            if pb.e1 is not None and pb.e1 == pa.e0:
                continue
            n_pairs += 1
            ba, bb = bbox(pa), bbox(pb)
            # bounding boxes more than 2R apart along x or y (in intervals)
            if any(pos(iv.mpf(p) - iv.mpf(q) - 2 * R) for p, q in
                   ((ba[0], bb[1]), (bb[0], ba[1]),
                    (ba[2], bb[3]), (bb[2], ba[3]))):
                continue
            wa, wb = window(pa, t_end), window(pb, t_end)
            a_, b_ = max(wa[0], wb[0]), min(wa[1], wb[1])
            if a_ > b_:
                continue
            n_time += 1
            if not apart(pa, pb, a_, b_, budget=[200000]):
                raise Fail(f"cannot separate {pa.ball} (from event {pa.e0}) "
                           f"and {pb.ball} (from event {pb.e0}) in time")
    return n_pairs, n_time


def check_admissible(w):
    pts = [w["cue"]] + w["balls"]
    P = [(iv.mpf(p[0]), iv.mpf(p[1])) for p in pts]
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            if not pos((P[i][0] - P[j][0]) ** 2 + (P[i][1] - P[j][1]) ** 2
                       - 4 * R * R):
                raise Fail("overlapping balls")
        for wall in WALLS:
            d = P[i][wall[1]] - wall[2]
            if not pos(d if wall[0] in ("left", "bottom") else -d):
                raise Fail("ball outside the table")
        for Q in POCKETS:
            if not pos((P[i][0] - Q[0]) ** 2 + (P[i][1] - Q[1]) ** 2 - RHO * RHO):
                raise Fail("ball inside a pocket disk")


def main():
    path = os.path.join(HERE, "witness_21.json")
    args = sys.argv[1:]
    if len(args) == 3:          # sub-strip: V0* + [a, b], phi* +/- dphi
        a, b, dphi = args
        offs = (lo(iv.mpf(a)), hi(iv.mpf(b)))
        label = f"V0 in V0* + [{a}, {b}]"
    else:                       # symmetric: V0* +/- dV, phi* +/- dphi
        dV = args[0] if args else "0.01"
        dphi = args[1] if len(args) > 1 else "1e-90"
        offs = (-hi(iv.mpf(dV)), hi(iv.mpf(dV)))
        label = f"V0 +/- {dV}"
    w = json.load(open(path))
    rphi = hi(iv.mpf(dphi))
    V = iv.mpf(w["v0"]) + iv.mpf([offs[0], offs[1]])
    PHI = iv.mpf(w["phi"]) + iv.mpf([-rphi, rphi])
    print(f"certifying strip: {label}, phi +/- {dphi}, {PREC} bits")
    t0 = time.time()
    try:
        check_admissible(w)
        plan = plan_events(w)
        pieces, same_start, designed, potted, t_end = build(w, V, PHI, plan)
        n_pairs, n_time = check_pairs(pieces, same_start, designed, t_end)
    except Fail as exc:
        print(f"NOT CERTIFIED: {exc}")
        return 1
    n_pieces = sum(len(p) for p in pieces.values())
    print(f"planned events              : {len(plan)}")
    print(f"pieces (rest + straight)    : {n_pieces}")
    print(f"collisions (moving onto rest): {len(designed)}")
    print(f"piece pairs examined        : {n_pairs} "
          f"({n_time} needed a timing argument)")
    print(f"min pocket entry speed      : "
          f"{mp.nstr(min(lo(p[3]) for p in potted), 6)} m/s (slowest stroke)")
    print(f"all pocketed by time        : {mp.nstr(t_end, 8)} s")
    print(f"wall time                   : {time.time() - t0:.1f} s")
    print("\nCERTIFIED: every stroke in the strip pockets all "
          f"{len(potted)} object balls.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
