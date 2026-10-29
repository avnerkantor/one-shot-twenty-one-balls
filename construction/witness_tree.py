"""
witness_tree.py -- a tree-shaped clearance with a wide window of directions.

The certified witness (proof/witness_21.json) is a chain: the cue ball
strikes all twenty-one object balls in turn, with thin cuts, and its
window of directions is 9.5e-43 rad wide because twenty thin cuts act in
series. Open problem (P4) of the arXiv paper asks for a configuration in
which no ball depends on more than a few earlier collisions.

This constructor grows a TREE instead of a chain. Every moving ball is a
potential carrier: it strikes a few resting balls with moderate cuts over
short legs and then drops into a pocket itself. A struck ball leaves along
the line of centres, so, as in witness_bv.py, each child is placed by the
ghost-ball rule run in reverse, on the carrier's actual simulated path.
The depth of a ball is the number of collisions between the stroke and
its own strike; the first-order factor per collision is about
c*l/(2R cos theta), so short legs l and moderate cuts theta keep it small.

The simulator (model_sim.py, double precision) is the single source of
truth during the search; every committed step is re-verified from the
stroke on a stencil of nearby strokes (speed +/- dv, direction +/- dphi),
so the finished configuration has a window of at least about 2*dphi by
construction. The result is then replayed in the 80-digit engine and its
window of directions is bracketed by bisection in both engines.

Run:  python construction/witness_tree.py SEED [OUT.json] [N_TARGET] [DPHI]
"""
import collections
import json
import math
import os
import random
import sys
import time

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "model"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import model_sim as M                      # noqa: E402
import witness_bv as WB                    # noqa: E402

R, A, E, L, W, RHO = M.R, M.A, M.E_BALL, M.L, M.W, M.RHO
POCKETS = M.POCKETS
C = (1 + E) / 2


# ------------------------------------------------------------- simulation
def simulate(cue0, v0, phi, plan, engine=M, max_events=3000, max_time=25.0):
    """Run the stroke; return (balls, sim). plan: list of dicts with 'pos'.
    The clearance itself is over within a few seconds; the cue ball's long
    roll afterwards is cut off at max_time."""
    balls = [engine.Ball(cue0[0], cue0[1], "cue")]
    for i, p in enumerate(plan):
        balls.append(engine.Ball(p["pos"][0], p["pos"][1], f"b{i}"))
    engine.shoot(balls[0], v0, phi)
    sim = engine.Simulator(balls, record_events=True)
    sim.run(max_events=max_events, max_time=max_time)
    return balls, sim


def intended_pairs(plan):
    """For each carrier, its children in order: {carrier: [child names]}."""
    kids = {}
    for i, p in enumerate(plan):
        kids.setdefault(p["parent"], []).append(f"b{i}")
    return kids


def _moved(sim):
    """Names of balls that have had an event (hence moved) so far."""
    out = set()
    for (t, kind, who, extra) in sim.events:
        if kind == "bb":
            out.update(who)
        else:
            out.add(who)
    return out


STATS = {"ranked": 0, "ray_ok": 0, "checked": 0, "passed": 0, "calls": 0,
         "centre_fail": 0, "stencil_fail": 0, "pocket_fail": 0}


REASON = [""]


def check(cue0, v0, phi, plan, engine=M, want_cue_safe=True, pending=frozenset()):
    """True if the stroke realises the tree: every collision is a planned
    (carrier, child) pair, in the planned order for each carrier, each
    child's first event is its strike, no other ball-ball contact occurs,
    the cue ball is not pocketed, and every object ball is pocketed. A
    carrier in `pending` still awaits its next child: it is removed from
    the simulation right after its last planned collision, since in the
    finished tree it will not roll free beyond that point."""
    balls, sim = simulate(cue0, v0, phi, plan, engine, max_events=0)
    by = {bl.name: bl for bl in balls}
    kids = intended_pairs(plan)
    parent = {f"b{i}": p["parent"] for i, p in enumerate(plan)}
    seen = {k: 0 for k in kids}
    touched = set()
    required = [bl for bl in balls if bl.name != "cue" and bl.name not in pending]
    n_ev = 0
    for _ in range(3000):
        if all(not bl.active for bl in required):
            return True                              # every required ball is pocketed
        if sim.t > 25.0 or not sim.step():
            REASON[0] = f'timeout/stall at t={float(sim.t):.2f}, active {[bl.name for bl in required if bl.active]}'
            return False
        for (t, kind, who, extra) in sim.events[n_ev:]:
            if kind == "bb":
                x, y = who
                if x in parent and parent[x] == y and x not in touched:
                    child, carrier = x, y
                elif y in parent and parent[y] == x and y not in touched:
                    child, carrier = y, x
                else:
                    REASON[0] = f'unplanned collision {who} at t={float(t):.3f}'
                    return False
                k = seen[carrier]
                if k >= len(kids[carrier]) or kids[carrier][k] != child:
                    REASON[0] = f'out of order {carrier}->{child} at t={float(t):.3f}'
                    return False
                seen[carrier] = k + 1
                touched.add(child)
                touched.add(carrier)
                for nm in (carrier, child):          # retire a pending ball after its last planned collision
                    if nm != "cue" and nm in pending and seen.get(nm, 0) >= len(kids.get(nm, [])):
                        bl = by[nm]
                        bl.active = False
                        bl.set_vel(0.0, 0.0)
            elif kind == "cushion":
                touched.add(who)
            elif kind == "pocket":
                if want_cue_safe and who == "cue":
                    REASON[0] = f'cue pocketed at t={float(t):.3f}'
                    return False
                touched.add(who)
        n_ev = len(sim.events)
        for bl in required:                           # a required ball at rest will never drop
            if bl.active and bl.name in touched and bl._s < engine.V_MIN:
                REASON[0] = f'{bl.name} came to rest at t={float(sim.t):.2f}'
                return False
    return False


def check_stencil(cue0, v0, phi, plan, dv, dphi, engine=M, pending=frozenset()):
    if not check(cue0, v0, phi, plan, engine, pending=pending):
        STATS["centre_fail"] += 1
        return False
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if (i or j) and not check(cue0, v0 + i * dv, phi + j * dphi, plan, engine, pending=pending):
                STATS["stencil_fail"] += 1
                return False
    return True


def state_after(cue0, v0, phi, plan, carrier, k):
    """State of `carrier` just after its k-th planned collision (k = 0: the
    moment it starts moving, i.e. the stroke for the cue or its own strike
    for an object ball). Returns (x, y, vx, vy) or None."""
    balls, sim = simulate(cue0, v0, phi, plan)
    by = {b.name: b for b in balls}
    if carrier == "cue" and k == 0:
        return (cue0[0], cue0[1], v0 * math.cos(phi), v0 * math.sin(phi), 0.0)
    # replay step by step to catch the state right after the event
    balls, sim = simulate(cue0, v0, phi, plan, max_events=0)
    by = {b.name: b for b in balls}
    kids = intended_pairs(plan).get(carrier, [])
    parent = {f"b{i}": p["parent"] for i, p in enumerate(plan)}
    n_as_carrier = 0
    for _ in range(3000):
        if not sim.step():
            return None
        t, kind, who, extra = sim.events[-1]
        if kind != "bb" or carrier not in who:
            continue
        other = who[0] if who[1] == carrier else who[1]
        if k == 0 and carrier != "cue" and parent.get(carrier) == other:
            b = by[carrier]
            return (b.x, b.y, b.vx, b.vy, sim.t)
        if other in kids:
            n_as_carrier += 1
            if n_as_carrier == k:
                b = by[carrier]
                return (b.x, b.y, b.vx, b.vy, sim.t)
    return None


# ------------------------------------------------------------- candidates
def ray_to_pocket(x, y, vx, vy, margin=0.9, horizon=4.0):
    """Does a lone ball from this state reach a pocket disk, with cushion
    rebounds under Model M, before it stops and within `horizon` metres of
    path, passing within margin*RHO of the pocket point? Returns the pocket
    index or None. Exact polyline, segment by segment (no sampling)."""
    v = math.hypot(vx, vy)
    if v < 1e-9:
        return None
    s_left = min(horizon, v * v / (2 * A))
    ux, uy = vx / v, vy / v
    for _ in range(16):
        cands = []
        if ux > 1e-12:
            cands.append(((L - R - x) / ux, 0))
        elif ux < -1e-12:
            cands.append(((R - x) / ux, 0))
        if uy > 1e-12:
            cands.append(((W - R - y) / uy, 1))
        elif uy < -1e-12:
            cands.append(((R - y) / uy, 1))
        d, axis = min(cands)
        d = max(d, 0.0)
        seg = min(d, s_left)
        ex, ey = x + seg * ux, y + seg * uy
        for gi, (gx, gy) in enumerate(POCKETS):
            t = (gx - x) * ux + (gy - y) * uy
            t = max(0.0, min(seg, t))
            if math.hypot(x + t * ux - gx, y + t * uy - gy) < margin * RHO:
                return gi
        if seg >= s_left - 1e-12:
            return None
        s_left -= seg
        x, y = ex, ey
        if axis == 0:
            ux = -M.E_CUSH * ux
        else:
            uy = -M.E_CUSH * uy
        f = math.hypot(ux, uy)
        ux, uy = ux / f, uy / f
        s_left *= f * f                      # the rebound costs kinetic energy
    return None


def path_segments(cue0, v0, phi, plan, exclude, engine=M):
    """Straight segments travelled by every ball other than `exclude` in the
    baseline simulation (positions at consecutive events)."""
    balls, sim = simulate(cue0, v0, phi, plan, engine, max_events=0)
    last = {bl.name: (0.0, bl.x, bl.y) for bl in balls}
    segs = []                                    # (t_start, (ax, ay), (bx, by))
    for _ in range(3000):
        if sim.t > 25.0 or not sim.step():
            break
        t, kind, who, extra = sim.events[-1]
        names = list(who) if kind == "bb" else [who]
        for nm in names:
            bl = next(b for b in balls if b.name == nm)
            t0, ax, ay = last[nm]
            if nm != exclude and (ax, ay) != (bl.x, bl.y):
                segs.append((t0, (ax, ay), (bl.x, bl.y)))
            last[nm] = (sim.t, bl.x, bl.y)
    for bl in balls:
        if bl.active and bl.name != exclude:
            t0, ax, ay = last[bl.name]
            segs.append((t0, (ax, ay), (bl.x, bl.y)))
    return segs


def dist_to_segments(x, y, segs, t_max):
    """Distance from (x, y) to the segments that start before t_max."""
    best = float("inf")
    for t0, (ax, ay), (bx, by) in segs:
        if t0 > t_max:
            continue
        dx, dy = bx - ax, by - ay
        d2 = dx * dx + dy * dy
        t = 0.0 if d2 == 0 else max(0.0, min(1.0, ((x - ax) * dx + (y - ay) * dy) / d2))
        d = math.hypot(ax + t * dx - x, ay + t * dy - y)
        if d < best:
            best = d
    return best


KILL = collections.Counter()
OFFSETS = (0.0, 0.2, -0.2, 0.4, -0.4, 0.6, -0.6, 0.75, -0.75)


def wrap(a):
    return (a + math.pi) % (2 * math.pi) - math.pi


def samples(cx, cy, cvx, cvy, min_leg, max_leg):
    """Points of the carrier's path where a child may be struck: before the
    carrier's own pocket, between min_leg and max_leg of path, clear of the
    cushions. Each is (t, leg, x, y, vx, vy, speed)."""
    out = []
    speed = math.hypot(cvx, cvy)
    for (ts, sx, sy, svx, svy) in WB.freeflight(cx, cy, cvx, cvy, t_end=speed / A, dt=0.002):
        if any(math.hypot(sx - gx, sy - gy) < RHO for (gx, gy) in POCKETS):
            break
        leg = math.hypot(sx - cx, sy - cy)
        if leg > max_leg:
            break
        sp = math.hypot(svx, svy)
        if sp < 0.6:
            break
        if leg < min_leg:
            continue
        if not (R + 0.08 < sx < L - R - 0.08 and R + 0.08 < sy < W - R - 0.08):
            continue
        out.append((ts, leg, sx, sy, svx, svy, sp))
    return out


def first_wall(x, y, ux, uy):
    """Distance along u to the first cushion line."""
    d = float("inf")
    if ux > 1e-12:
        d = min(d, (L - R - x) / ux)
    elif ux < -1e-12:
        d = min(d, (R - x) / ux)
    if uy > 1e-12:
        d = min(d, (W - R - y) / uy)
    elif uy < -1e-12:
        d = min(d, (R - y) / uy)
    return d


def enters_pocket(x, y, ux, uy, gi):
    """Does the ray from (x, y) along u enter the disk of pocket gi before
    its centre crosses a cushion line? In Model M a straight approach is
    accepted only within about 16.5 degrees of a corner's diagonal, or at
    more than 28.5 degrees to the cushion for a middle pocket; shallower
    balls meet the cushion line first (and may still drop after the
    rebound, which the simulation decides)."""
    gx, gy = POCKETS[gi]
    dx, dy = gx - x, gy - y
    b = dx * ux + dy * uy
    if b <= 0:
        return False
    h2 = dx * dx + dy * dy - b * b
    if h2 >= RHO * RHO:
        return False
    return b - math.sqrt(RHO * RHO - h2) < first_wall(x, y, ux, uy)


def child_geometry(sx, sy, svx, svy, sp, u):
    bx, by = sx + 2 * R * u[0], sy + 2 * R * u[1]
    cos_th = (svx * u[0] + svy * u[1]) / sp
    svx2, svy2, cvx2, cvy2 = WB.collide(svx, svy, u[0], u[1])
    return bx, by, cos_th, (svx2, svy2), (cvx2, cvy2)


def place_ok(bx, by, occupied, segs, t_max):
    if not (R + 0.02 < bx < L - R - 0.02 and R + 0.02 < by < W - R - 0.02):
        KILL["off_table"] += 1
        return False
    if any(math.hypot(bx - ox, by - oy) < 2 * R + 0.03 for (ox, oy) in occupied):
        KILL["overlap"] += 1
        return False
    if any(math.hypot(bx - px, by - py) < RHO + 2 * R for (px, py) in POCKETS):
        KILL["near_pocket"] += 1
        return False
    if dist_to_segments(bx, by, segs, t_max) < 2 * R + 0.015:
        KILL["on_path"] += 1
        return False
    return True


def interp(s0, s1, lam):
    return tuple(a + lam * (b - a) for a, b in zip(s0, s1))


def candidates(cue0, v0, phi, plan, carrier, k, mode, prefer_carrier, dv, dphi, rng,
               topm=2, max_checks=20, min_leg=0.09, max_leg=1.2, min_enter=0.30):
    """Children for `carrier` after its k-th collision. Every child is aimed
    at a point of a pocket disk, so it is a valid leaf; a child is marked
    carrier_ok if it leaves fast along a long straight run to its pocket,
    so that children of its own can be placed on the way.
    mode 'more': the carrier continues (and must keep speed); mode 'last':
    the carrier then drops into a pocket (the cue is exempt). For 'last'
    the contact point is solved so that the deflected carrier heads
    straight into some pocket, then verified in the simulator.
    Returns up to topm stencil-verified candidates, best first, as
    (score, pos, pocket_idx, info)."""
    st = state_after(cue0, v0, phi, plan, carrier, k)
    if st is None:
        return []
    cx, cy, cvx, cvy, t_state = st
    smp = samples(cx, cy, cvx, cvy, min_leg, max_leg)
    if not smp:
        return []
    occupied = [cue0] + [p["pos"] for p in plan]
    segs = path_segments(cue0, v0, phi, plan, exclude=carrier)
    need_pocket = (mode == "last" and carrier != "cue")
    is_cue = carrier == "cue"
    min_cos = 0.2 if is_cue else 0.25
    min_keep = 1.5 if is_cue else 0.9

    def evaluate(s, u, gi, extra):
        ts, leg, sx, sy, svx, svy, sp = s
        bx, by, cos_th, (svx2, svy2), (cvx2, cvy2) = child_geometry(sx, sy, svx, svy, sp, u)
        if cos_th < min_cos:
            KILL["thin_cut"] += 1
            return None
        if not place_ok(bx, by, occupied, segs, t_state + ts + 2.5):
            return None
        dep, ckeep = math.hypot(svx2, svy2), math.hypot(cvx2, cvy2)
        if not enters_pocket(bx, by, u[0], u[1], gi):
            KILL["jaw"] += 1
            return None
        gx, gy = POCKETS[gi]
        plen = math.hypot(gx - bx, gy - by)
        v_enter = math.sqrt(max(dep * dep - 2 * A * (plen - RHO), 0.0))
        if v_enter < min_enter:
            KILL["slow_entry"] += 1
            return None
        if mode == "more" and ckeep < max(min_keep, (0.8 if is_cue else 0.5) * sp):
            KILL["carrier_slow"] += 1          # the carrier must keep most of its speed
            return None
        if need_pocket:
            g2 = ray_to_pocket(sx, sy, cvx2, cvy2)
            if g2 is None or ckeep < 0.4:
                KILL["ray_miss"] += 1
                return None
            extra = dict(extra, carrier_pocket=g2)
        carrier_ok = plen >= 1.0 and dep >= 1.5
        amp = plen / (2 * R * cos_th)
        score = -math.log(amp) - 0.6 * leg + 0.05 * v_enter
        if mode == "more":
            score += 0.1 * ckeep
        if prefer_carrier:
            score += 1.5 if carrier_ok else 0.0
            if carrier_ok:
                score += 0.4 * min(plen, 2.0)
        info = {"leg": leg, "cos": cos_th, "plen": plen, "v_enter": v_enter, "ckeep": ckeep,
                "dep": dep, "carrier_ok": carrier_ok}
        info.update(extra)
        return (score, (bx, by), gi, info)

    def out_angle(s, u):
        ts, leg, sx, sy, svx, svy, sp = s
        _, _, _, _, (cvx2, cvy2) = child_geometry(sx, sy, svx, svy, sp, u)
        return math.atan2(cvy2, cvx2)

    geo = []
    for gi, (gx, gy) in enumerate(POCKETS):
        for off in OFFSETS:
            prev = None
            for s in smp:
                ts, leg, sx, sy, svx, svy, sp = s
                d0 = math.hypot(gx - sx, gy - sy)
                if d0 < 1e-6:
                    continue
                tx, ty = -(gy - sy) / d0, (gx - sx) / d0
                target = (gx + off * RHO * tx, gy + off * RHO * ty)
                u = WB.aim_ball_at_pocket(sx, sy, svx, svy, target[0], target[1])
                if u is None:
                    KILL["aim_none"] += 1
                    prev = None
                    continue
                if not need_pocket:
                    c = evaluate(s, u, gi, {"off": off})
                    if c:
                        geo.append(c)
                    continue
                c = evaluate(s, u, gi, {"off": off, "target": "grid"})   # accepts cushion routes too
                if c:
                    geo.append(c)
                ang = out_angle(s, u)
                diffs = [wrap(ang - math.atan2(py - sy, px - sx)) for (px, py) in POCKETS]
                if prev is not None:
                    for g2 in range(len(POCKETS)):
                        a0, a1 = prev[1][g2], diffs[g2]
                        if a0 * a1 < 0 and abs(a0) < 0.5 and abs(a1) < 0.5:
                            sm = interp(prev[0], s, a0 / (a0 - a1))
                            um = WB.aim_ball_at_pocket(sm[2], sm[3], sm[4], sm[5], target[0], target[1])
                            if um is None:
                                continue
                            c = evaluate(sm, um, gi, {"off": off, "target": g2})
                            if c:
                                geo.append(c)
                prev = (s, diffs)
    STATS["ranked"] += len(geo)
    geo.sort(key=lambda z: -(z[0] + rng.uniform(0, 0.3)))
    out, seen, checks = [], set(), 0
    pending = frozenset([carrier]) if mode == "more" else frozenset()
    for cand in geo:
        cell = (round(cand[1][0] / 0.03), round(cand[1][1] / 0.03), cand[2])
        if cell in seen:
            continue
        seen.add(cell)
        if checks >= max_checks:
            break
        checks += 1
        STATS["checked"] += 1
        trial = plan + [{"pos": cand[1], "parent": carrier, "pocket": cand[2]}]
        if not check_stencil(cue0, v0, phi, trial, dv, dphi, pending=pending):
            continue
        if need_pocket and not carrier_pocketed(cue0, v0, phi, trial, carrier, pending):
            STATS["pocket_fail"] += 1
            continue
        out.append(cand)
        STATS["passed"] += 1
        if len(out) >= topm:
            break
    STATS["calls"] += 1
    if os.environ.get("TREE_DEBUG"):
        print(f"    candidates({carrier}, k={k}, {mode}, prefer_carrier={prefer_carrier}): {len(geo)} geometric, "
              f"{checks} checked, {len(out)} passed; kills {dict(KILL)}", flush=True)
        KILL.clear()
    return out


def carrier_pocketed(cue0, v0, phi, plan, carrier, pending=frozenset()):
    """Is `carrier` pocketed in the simulation (pending balls retired as in check)?"""
    balls, sim = simulate(cue0, v0, phi, plan, max_events=0)
    by = {bl.name: bl for bl in balls}
    kids = intended_pairs(plan)
    seen = {k: 0 for k in kids}
    n_ev = 0
    for _ in range(3000):
        if not by[carrier].active:
            return any(name == carrier for (_, name, _, _) in sim.pocketed)
        if sim.t > 25.0 or not sim.step():
            return False
        for (t, kind, who, extra) in sim.events[n_ev:]:
            if kind == "bb":
                for nm in who:
                    if nm in kids:
                        seen[nm] += 1
                for nm in who:
                    if nm != "cue" and nm in pending and seen.get(nm, 0) >= len(kids.get(nm, [])):
                        by[nm].active = False
                        by[nm].set_vel(0.0, 0.0)
        n_ev = len(sim.events)
    return False


# ------------------------------------------------------------- growth
class Budget(Exception):
    pass


def grow(seed, n_target=21, dv=0.02, dphi=1e-6, shape=(0, 3, 2, 1), v0=11.9,
         verbose=True, budget=800, sub_budget=20):
    """Depth-first growth. The cue strikes children until n_target balls are
    placed. Every child is aimed at a pocket, so it is a valid leaf; a child
    marked carrier_ok is given up to shape[level] children of its own, the
    last of which must deflect it into a pocket, within a budget of
    sub_budget candidate calls per subtree; if that fails it stays a leaf.
    Returns the plan dict, or None (the best partial tree is saved when
    TREE_PARTIAL names a file)."""
    rng = random.Random(seed)
    cue0 = (0.35 + rng.uniform(0, 0.3), W / 2 + rng.uniform(-0.4, 0.4))
    phi = rng.uniform(-0.25, 0.25)
    plan = []
    best = []
    calls = [0]
    t0 = time.time()

    def say(msg):
        if verbose:
            print(f"  [{time.time() - t0:5.0f} s, {len(plan):2d} balls] {msg}", flush=True)

    def cands(carrier, k, mode, prefer_carrier, counter):
        if calls[0] >= budget:
            raise Budget()
        if counter is not None and counter[0] <= 0:
            return []
        calls[0] += 1
        if counter is not None:
            counter[0] -= 1
        return candidates(cue0, v0, phi, plan, carrier, k, mode, prefer_carrier, dv, dphi, rng)

    def remember():
        nonlocal best
        if len(plan) > len(best):
            best = [dict(p) for p in plan]

    def fill_sub(carrier, level, k, counter):
        """Give an object carrier children from slot k on; return True if
        the carrier ends pocketed (a 'last' child placed), False otherwise,
        in which case the plan is restored to its state at slot 0 by the
        caller."""
        remaining = n_target - len(plan)
        want = shape[min(level, len(shape) - 1)]
        if remaining == 0 or k >= want:
            return False
        modes = ["last"] if (k == want - 1 or remaining <= level + 1) else ["more", "last"]
        for mode in modes:
            want_carrier = (level + 1 < len(shape) and shape[level + 1] > 0 and remaining > 3)
            for score, pos, gi, info in cands(carrier, k, mode, want_carrier, counter):
                idx = len(plan)
                child = f"b{idx}"
                plan.append({"pos": pos, "parent": carrier, "pocket": gi, "level": level + 1, **info})
                say(f"try {child} <- {carrier} ({mode}, depth {level + 1}): leg {info['leg']:.2f} m, "
                    f"cut {math.degrees(math.acos(min(1, info['cos']))):.0f} deg, to pocket {info['plen']:.2f} m"
                    + (", carrier_ok" if info["carrier_ok"] else ""))
                if info["carrier_ok"] and remaining > 2:
                    give_children(child, level + 1)
                if mode == "last":
                    return True
                if fill_sub(carrier, level, k + 1, counter):
                    return True
                say(f"undo {child}")
                del plan[idx:]
        return False

    def give_children(carrier, level):
        """Try to make `carrier` strike children of its own; keep them only
        if the carrier is finally pocketed, else leave it a leaf."""
        n0 = len(plan)
        counter = [sub_budget]
        if not fill_sub(carrier, level, 0, counter):
            del plan[n0:]

    def fill_cue():
        k = 0
        while len(plan) < n_target:
            remaining = n_target - len(plan)
            mode = "last" if remaining == 1 else "more"
            want_carrier = remaining > 3
            got = cands("cue", k, mode, want_carrier, None)
            if not got:
                return False
            score, pos, gi, info = got[0]
            idx = len(plan)
            child = f"b{idx}"
            plan.append({"pos": pos, "parent": "cue", "pocket": gi, "level": 1, **info})
            say(f"cue child {child} ({mode}): leg {info['leg']:.2f} m, cut {math.degrees(math.acos(min(1, info['cos']))):.0f} deg, "
                f"to pocket {info['plen']:.2f} m, cue keeps {info['ckeep']:.1f} m/s" + (", carrier_ok" if info["carrier_ok"] else ""))
            if info["carrier_ok"] and n_target - len(plan) > 1:
                give_children(child, 1)
            remember()
            k += 1
        return True

    try:
        ok = fill_cue()
    except Budget:
        ok = False
        say("budget exhausted")
    if not ok or len(plan) < n_target:
        say(f"dead end: {len(plan)} balls now, best {len(best)}, {calls[0]} candidate calls")
        if os.environ.get("TREE_PARTIAL"):
            json.dump(_as_json(cue0, v0, phi, best, dv, dphi, seed, partial=True),
                      open(os.environ["TREE_PARTIAL"], "w"), indent=1)
        return None
    if not check_stencil(cue0, v0, phi, plan, dv, dphi):
        say("final stencil check failed")
        return None
    return _as_json(cue0, v0, phi, plan, dv, dphi, seed)


def _as_json(cue0, v0, phi, plan, dv, dphi, seed, partial=False):
    return {"cue": list(cue0), "v0": v0, "phi": phi,
            "balls": [list(p["pos"]) for p in plan],
            "pockets": [p["pocket"] for p in plan],
            "parents": [p["parent"] for p in plan],
            "depths": [p["level"] for p in plan],
            "n": len(plan), "dv": dv, "dphi": dphi, "seed": seed, "partial": partial,
            "engine": "model_sim (double), stencil-verified"}


# ------------------------------------------------------------- window
def window(w, engine=M, lo=1e-9, hi=0.05, iters=40):
    """Bracket the window of directions [phi - a, phi + b] on which the
    tree is realised (check passes), by bisection in each direction."""
    plan = [{"pos": tuple(p), "parent": par} for p, par in zip(w["balls"], w["parents"])]
    cue0, v0, phi = tuple(w["cue"]), w["v0"], w["phi"]
    out = []
    for sgn in (+1, -1):
        a, b = lo, hi
        if not check(cue0, v0, phi + sgn * a, plan, engine):
            out.append(0.0)
            continue
        for _ in range(iters):
            m = math.sqrt(a * b)
            if check(cue0, v0, phi + sgn * m, plan, engine):
                a = m
            else:
                b = m
            if b / a < 1.01:
                break
        out.append(a)
    return out


if __name__ == "__main__":
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    out = sys.argv[2] if len(sys.argv) > 2 else None
    n_target = int(sys.argv[3]) if len(sys.argv) > 3 else 21
    dphi = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-6
    t0 = time.time()
    w = grow(seed, n_target=n_target, dphi=dphi)
    if w is None:
        print(f"seed {seed}: dead end after {time.time() - t0:.0f} s")
        sys.exit(1)
    plus, minus = window(w)
    w["window_double"] = [minus, plus]
    print(f"seed {seed}: {w['n']} balls in {time.time() - t0:.0f} s; "
          f"window (double engine) phi + [-{minus:.3e}, +{plus:.3e}] rad; depths {w['depths']}")
    if out:
        json.dump(w, open(out, "w"), indent=1)
        print("written", out)
