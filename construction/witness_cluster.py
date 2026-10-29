"""
witness_cluster.py -- a near-frozen cluster clearance: the physical-world design.

In the first-order analysis a struck ball leaves with the direction error
delta / (2R cos theta), where delta is the striker's sideways offset at
contact. For the cue ball that offset is its aiming error times the
distance travelled, and twenty thin cuts in series give the certified
witness its window of 1e-42 rad. Trick-shot players know frozen and
near-frozen combinations as reliable: a ball that rests in its exact
place and travels only a millimetre or two before striking its neighbour
leaves along a line fixed by the resting geometry, so the aiming error
enters only through the entry balls and is not amplified down the tree.
(The design note that errors are damped by 2R cos(theta) / gap at every
collision is too optimistic: the striker's offset at contact grows with
the centre distance 2R + gap, not with the gap. What the measurement
shows is that the slopes of all quantities with respect to the aim stay
of order 1 to 10 rad/rad, against 1e42 for the chain.)

This constructor grows a branching cluster with gaps of 1.5 to 3 mm. The
cue strikes a few entry balls over short legs; every other ball first
strikes one or two neighbours across a gap along its line of centres and
then rolls, deflected, into a pocket. The 21-ball result cluster_21.json
has a window of directions 9.6e-4 rad wide in the 80-digit engine. Its
weakness is geometric: on four legs the cue passes the ball it will
strike next with 0.09 to 0.63 mm of clearance (planned pairs are not
checked by min_unplanned_clearance), so the configuration is robust in
aim but not in the placement of the balls.

The search, verification and bookkeeping are those of witness_tree.py;
only the candidate generator differs: children sit at fixed gaps, the
continuous freedom is the aim point inside the child's pocket disk, and
for a carrier's last child that offset is solved so that the deflected
carrier heads into a pocket.

Run:  python construction/witness_cluster.py SEED [OUT.json] [N_TARGET] [DPHI]
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "model"))

import model_sim as M                      # noqa: E402
import witness_bv as WB                    # noqa: E402
import witness_tree as T                   # noqa: E402

R, A, L, W, RHO, POCKETS = T.R, T.A, T.L, T.W, T.RHO, T.POCKETS
GAPS = (0.0015, 0.002, 0.0025, 0.003)
CUE_LEGS = (0.0015, 0.002, 0.003, 0.005, 0.008, 0.012, 0.02, 0.03)
ANGLES = [math.radians(a) for a in range(-76, 77, 2) if abs(a) >= 8]   # departure angles of a child
CUE_ANGLES = [math.radians(a) for a in range(-78, 79, 1) if abs(a) >= 30]  # the cue keeps speed: cuts of 30 deg or more
OVERLAP = 0.0005                                   # minimum gap between any two resting balls
PATH_MARGIN = 0.003
LAST_LEG = 0.6  
CUE_MAX = float(os.environ.get("CUE_MAX", "0.25"))      # longest cue leg between entry balls                                   # how far along its path a carrier's last child may sit


def trace_to_pocket(x, y, vx, vy, margin=0.9, horizon=6.0):
    """Like witness_tree.ray_to_pocket, but returns (pocket index, path
    length) or (None, None): the exact polyline of a lone ball under Model M
    with cushion rebounds, until it stops or has travelled `horizon`."""
    v = math.hypot(vx, vy)
    if v < 1e-9:
        return None, None
    s_left = min(horizon, v * v / (2 * A))
    ux, uy = vx / v, vy / v
    travelled = 0.0
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
        for gi, (gx, gy) in enumerate(POCKETS):
            t = (gx - x) * ux + (gy - y) * uy
            t = max(0.0, min(seg, t))
            if math.hypot(x + t * ux - gx, y + t * uy - gy) < margin * RHO:
                return gi, travelled + t
        if seg >= s_left - 1e-12:
            return None, None
        s_left -= seg
        travelled += seg
        x, y = x + seg * ux, y + seg * uy
        if axis == 0:
            ux = -M.E_CUSH * ux
        else:
            uy = -M.E_CUSH * uy
        f = math.hypot(ux, uy)
        ux, uy = ux / f, uy / f
        s_left *= f * f
    return None, None


def crosses_halfline(ax, ay, bx, by, px, py, ux, uy):
    """Does the segment a-b cross the half-line from p along u?"""
    dx, dy = bx - ax, by - ay
    den = dx * uy - dy * ux
    if abs(den) < 1e-12:
        return False
    # a + t (b - a) = p + s u
    t = ((px - ax) * uy - (py - ay) * ux) / den
    sgm = ((px - ax) * dy - (py - ay) * dx) / den
    return 0.0 <= t <= 1.0 and sgm >= 0.0


CLEAR = 0.002                                      # least clearance allowed between balls not meant to meet


def min_unplanned_clearance(cue0, v0, phi, plan, pending=frozenset(), engine=M):
    """Smallest distance beyond 2R between two balls that are not planned to
    collide, at any time of the clearance (sampled every 2 mm of motion
    inside each inter-event interval, where the motion is explicit).
    Pending carriers are retired after their last planned collision, as in
    witness_tree.check."""
    balls, sim = T.simulate(cue0, v0, phi, plan, engine, max_events=0)
    by = {b.name: b for b in balls}
    kids = T.intended_pairs(plan)
    seen = {k: 0 for k in kids}
    planned = {tuple(sorted((p["parent"], f"b{i}"))) for i, p in enumerate(plan)}
    required = [b for b in balls if b.name != "cue" and b.name not in pending]
    best = float("inf")
    graze = {}                         # (ball, pocket) -> min |distance - RHO| while moving
    n_ev = 0
    for _ in range(3000):
        if all(not b.active for b in required):
            break
        ev = sim._next_event()
        if ev is None:
            break
        dt = ev[0]
        act = [b for b in balls if b.active]
        vmax = max([b._s for b in act] + [0.0])
        n = max(1, min(60, int(dt * vmax / 0.002)))
        for j in range(n + 1):
            tau = dt * j / n
            pos = []
            for b in act:
                sp = b._s
                if sp < engine.V_MIN:
                    pos.append((b.x, b.y))
                else:
                    tt = min(tau, sp / A)
                    ux, uy = b.vx / sp, b.vy / sp
                    pos.append((b.x + b.vx * tt - 0.5 * A * tt * tt * ux,
                                b.y + b.vy * tt - 0.5 * A * tt * tt * uy))
            for i in range(len(act)):
                if act[i]._s >= engine.V_MIN:
                    for gi, (gx, gy) in enumerate(POCKETS):
                        g = abs(math.hypot(pos[i][0] - gx, pos[i][1] - gy) - RHO)
                        key = (act[i].name, gi)
                        if g < graze.get(key, float("inf")):
                            graze[key] = g
                for k in range(i + 1, len(act)):
                    if tuple(sorted((act[i].name, act[k].name))) in planned:
                        continue
                    d = math.hypot(pos[i][0] - pos[k][0], pos[i][1] - pos[k][1]) - 2 * R
                    if d < best:
                        best = d
        if sim.t > 25.0 or not sim.step():
            break
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
    entered = {(name, M.POCKET_NAMES.index(pk)) for (_, name, pk, _) in sim.pocketed}
    g_best = min([g for key, g in graze.items() if key not in entered] + [float("inf")])
    return min(best, g_best)


def straight_samples(st, legs):
    """Points at the given arc lengths along the carrier's straight path
    from its state st = (x, y, vx, vy, t); each is (t, leg, x, y, vx, vy, speed)."""
    x, y, vx, vy, t0 = st
    v = math.hypot(vx, vy)
    if v < 1e-9:
        return []
    ux, uy = vx / v, vy / v
    out = []
    for s in legs:
        sp2 = v * v - 2 * A * s
        if sp2 <= 0.25:
            break
        sp = math.sqrt(sp2)
        sx, sy = x + s * ux, y + s * uy
        if not (R + 0.08 < sx < L - R - 0.08 and R + 0.08 < sy < W - R - 0.08):
            continue
        out.append(((v - sp) / A, s, sx, sy, sp * ux, sp * uy, sp))
    return out


def place_ok(bx, by, occupied, segs, t_max):
    if not (R + 0.02 < bx < L - R - 0.02 and R + 0.02 < by < W - R - 0.02):
        T.KILL["off_table"] += 1
        return False
    if any(math.hypot(bx - ox, by - oy) < 2 * R + OVERLAP for (ox, oy) in occupied):
        T.KILL["overlap"] += 1
        return False
    if any(math.hypot(bx - px, by - py) < RHO + 2 * R for (px, py) in POCKETS):
        T.KILL["near_pocket"] += 1
        return False
    if T.dist_to_segments(bx, by, segs, t_max) < 2 * R + PATH_MARGIN:
        T.KILL["on_path"] += 1
        return False
    return True


def candidates(cue0, v0, phi, plan, carrier, k, mode, prefer_carrier, dv, dphi, rng,
               topm=2, max_checks=40, min_enter=0.25):
    """Same contract as witness_tree.candidates, for the near-frozen cluster."""
    st = T.state_after(cue0, v0, phi, plan, carrier, k)
    if st is None:
        return []
    cx, cy, cvx, cvy, t_state = st
    is_cue = carrier == "cue"
    if is_cue:
        max_checks = 60
    need_pocket = (mode == "last" and not is_cue)
    if is_cue and k == 0:
        smp = T.samples(cx, cy, cvx, cvy, 0.25, 0.9)[::3]      # the free leg to the cluster
    elif is_cue:
        smp = straight_samples(st, CUE_LEGS) + T.samples(cx, cy, cvx, cvy, 0.04, CUE_MAX)
    elif need_pocket:
        # the last child may sit farther along the path: its position is the freedom
        # that lets the deflected carrier be aimed into a pocket
        smp = straight_samples(st, GAPS) + T.samples(cx, cy, cvx, cvy, 0.004, LAST_LEG)[::2]
    else:
        smp = straight_samples(st, GAPS)
    if not smp:
        return []
    occupied = [cue0] + [p["pos"] for p in plan]
    segs = T.path_segments(cue0, v0, phi, plan, exclude=carrier)
    cue_line = None
    if not is_cue:
        kc = sum(1 for p in plan if p["parent"] == "cue")
        cst = T.state_after(cue0, v0, phi, plan, "cue", kc)
        if cst is not None:
            cv = math.hypot(cst[2], cst[3])
            if cv > 1e-9:
                cue_line = (cst[0], cst[1], cst[2] / cv, cst[3] / cv)
    min_cos = 0.2 if is_cue else 0.25
    min_keep = 1.5 if is_cue else 0.8

    def evaluate(s, u, extra):
        """Full geometric test of a child struck at sample s with line of
        centres u; its pocket, direct or after cushion rebounds, is traced."""
        ts, leg, sx, sy, svx, svy, sp = s
        bx, by, cos_th, (svx2, svy2), (cvx2, cvy2) = T.child_geometry(sx, sy, svx, svy, sp, u)
        if cos_th < min_cos:
            T.KILL["thin_cut"] += 1
            return None
        if not place_ok(bx, by, occupied, segs, t_state + ts + 2.5):
            return None
        dep, ckeep = math.hypot(svx2, svy2), math.hypot(cvx2, cvy2)
        gi, plen = trace_to_pocket(bx, by, svx2, svy2)
        if gi is None:
            T.KILL["no_pocket"] += 1
            return None
        if cue_line is not None:
            # keep the sub-cluster and its exits on its own side of the cue's path
            ux2, uy2 = svx2 / dep, svy2 / dep
            run = min(plen, T.first_wall(bx, by, ux2, uy2))
            if crosses_halfline(bx, by, bx + run * ux2, by + run * uy2, *cue_line) or                     crosses_halfline(sx, sy, bx, by, *cue_line):
                T.KILL["cue_cross"] += 1
                return None
        if is_cue and mode == "more" and T.ray_to_pocket(sx, sy, cvx2, cvy2, horizon=1.0) is not None:
            T.KILL["cue_drops"] += 1
            return None
        v_enter = math.sqrt(max(dep * dep - 2 * A * plen, 0.0))
        if v_enter < min_enter:
            T.KILL["slow_entry"] += 1
            return None
        if mode == "more" and ckeep < max(min_keep, (0.7 if is_cue else 0.5) * sp):
            T.KILL["carrier_slow"] += 1
            return None
        if need_pocket:
            g2 = T.ray_to_pocket(sx, sy, cvx2, cvy2)
            if g2 is None or ckeep < 0.4:
                T.KILL["ray_miss"] += 1
                return None
            extra = dict(extra, carrier_pocket=g2)
        carrier_ok = dep >= 1.2 and plen >= 0.3
        sens = plen * leg / (2 * R * cos_th)
        score = -math.log(max(sens, 1e-6)) + 0.05 * v_enter + ((0.3 if is_cue else 0.1) * ckeep if mode == "more" else 0.0)
        if prefer_carrier and carrier_ok and (not is_cue or sp > 6.0):
            score += 1.0 + 0.3 * min(dep, 5.0)
        info = {"leg": leg, "cos": cos_th, "plen": plen, "v_enter": v_enter, "ckeep": ckeep,
                "dep": dep, "carrier_ok": carrier_ok, "sens": sens}
        info.update(extra)
        return (score, (bx, by), gi, info)

    def out_angle(s, u):
        ts, leg, sx, sy, svx, svy, sp = s
        _, _, _, _, (cvx2, cvy2) = T.child_geometry(sx, sy, svx, svy, sp, u)
        return math.atan2(cvy2, cvx2)

    def dirs(s):
        """Lines of centres at sample s: every departure angle of the child."""
        base = math.atan2(s[5], s[4])
        angles = CUE_ANGLES if is_cue else ANGLES
        return [(a, (math.cos(base + a), math.sin(base + a))) for a in angles]

    geo = []
    for s in smp:
        prev = None
        for a, u in dirs(s):
            c = evaluate(s, u, {"angle": a})
            if c:
                geo.append(c)
            if not need_pocket:
                continue
            ang = out_angle(s, u)
            diffs = [T.wrap(ang - math.atan2(py - s[3], px - s[2])) for (px, py) in POCKETS]
            if prev is not None:
                for g2 in range(len(POCKETS)):
                    a0, a1 = prev[1][g2], diffs[g2]
                    if a0 * a1 < 0 and abs(a0) < 0.3 and abs(a1) < 0.3:
                        am = prev[0] + (a - prev[0]) * a0 / (a0 - a1)
                        base = math.atan2(s[5], s[4])
                        um = (math.cos(base + am), math.sin(base + am))
                        c = evaluate(s, um, {"angle": am, "target": g2})
                        if c:
                            geo.append(c)
            prev = (a, diffs)
    if need_pocket and len(smp) > 1:
        # solve the position along the path for each departure angle
        for a in ANGLES[::2]:
            prev = None
            for s in smp:
                base = math.atan2(s[5], s[4])
                u = (math.cos(base + a), math.sin(base + a))
                ang = out_angle(s, u)
                diffs = [T.wrap(ang - math.atan2(py - s[3], px - s[2])) for (px, py) in POCKETS]
                if prev is not None:
                    for g2 in range(len(POCKETS)):
                        a0, a1 = prev[1][g2], diffs[g2]
                        if a0 * a1 < 0 and abs(a0) < 0.5 and abs(a1) < 0.5:
                            sm = T.interp(prev[0], s, a0 / (a0 - a1))
                            basem = math.atan2(sm[5], sm[4])
                            um = (math.cos(basem + a), math.sin(basem + a))
                            c = evaluate(sm, um, {"angle": a, "target": g2, "solved": "s"})
                            if c:
                                geo.append(c)
                prev = (s, diffs)
    T.STATS["ranked"] += len(geo)
    geo.sort(key=lambda z: -(z[0] + rng.uniform(0, 0.3)))
    out, seen, checks = [], set(), 0
    pending = frozenset([carrier]) if mode == "more" else frozenset()
    for cand in geo:
        cell = (round(cand[1][0] / 0.002), round(cand[1][1] / 0.002), cand[2])
        if cell in seen:
            continue
        seen.add(cell)
        if checks >= max_checks:
            break
        checks += 1
        T.STATS["checked"] += 1
        trial = plan + [{"pos": cand[1], "parent": carrier, "pocket": cand[2]}]
        if not T.check_stencil(cue0, v0, phi, trial, dv, dphi, pending=pending):
            continue
        if need_pocket and not T.carrier_pocketed(cue0, v0, phi, trial, carrier, pending):
            T.STATS["pocket_fail"] += 1
            continue
        if min_unplanned_clearance(cue0, v0, phi, trial, pending) < CLEAR:
            T.STATS["clearance_fail"] = T.STATS.get("clearance_fail", 0) + 1
            continue
        out.append(cand)
        T.STATS["passed"] += 1
        if len(out) >= topm:
            break
    T.STATS["calls"] += 1
    if os.environ.get("TREE_DEBUG"):
        print(f"    cluster candidates({carrier}, k={k}, {mode}): {len(geo)} geometric, {checks} checked, "
              f"{len(out)} passed; kills {dict(T.KILL)}", flush=True)
        T.KILL.clear()
    return out


T.candidates = candidates            # witness_tree.grow resolves the generator at call time


if __name__ == "__main__":
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    out = sys.argv[2] if len(sys.argv) > 2 else None
    n_target = int(sys.argv[3]) if len(sys.argv) > 3 else 21
    dphi = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-4
    t0 = time.time()
    w = T.grow(seed, n_target=n_target, dv=0.05, dphi=dphi, shape=(0, 3, 2, 2, 1), v0=12.0, budget=1500, sub_budget=30)
    if w is None:
        print(f"seed {seed}: dead end after {time.time() - t0:.0f} s")
        sys.exit(1)
    w["engine"] = "witness_cluster: model_sim (double), stencil-verified"
    plus, minus = T.window(w)
    w["window_double"] = [minus, plus]
    print(f"seed {seed}: {w['n']} balls in {time.time() - t0:.0f} s; "
          f"window (double engine) phi + [-{minus:.3e}, +{plus:.3e}] rad; depths {w['depths']}")
    if out:
        json.dump(w, open(out, "w"), indent=1)
        print("written", out)
