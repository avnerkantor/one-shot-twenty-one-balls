"""
witness_bv.py -- Build-and-verify witness constructor.

The earlier open-loop constructors (witness_tree.py) predicted the whole
collision chain analytically and then checked it in the simulator; once
cushion bounces accumulated, the analytic path drifted from what the
simulator actually does, and no configuration reproduced its intended
chain.

This constructor makes the simulator the single source of truth. It
grows the chain one collision at a time:

  1. A moving "carrier" travels under the real Model M dynamics.
  2. We look ahead along the carrier's ACTUAL trajectory (as the
     simulator would integrate it, including cushion bounces) to a
     candidate strike point.
  3. We place a fresh object ball there, oriented so that after the
     (real, restitution-e) collision the struck ball departs straight
     at a chosen pocket, and the carrier is deflected to continue the
     chain.
  4. We re-simulate from scratch and confirm the intended collision
     sequence still holds exactly, then recurse.

Because every step is validated against the simulator before we commit,
there is no drift: what we build is what pots.

The public entry point is build_witness(); it returns a dict with the
cue shot, the placed object-ball coordinates, and the verified pocketing
log, or None if the greedy growth gets stuck (callers can retry with a
different seed / cue angle).
"""

import math
import random

import model_sim as M

R, A, E = M.R, M.A, M.E_BALL
L, W, RHO = M.L, M.W, M.RHO
POCKETS = M.POCKETS
PNAMES = M.POCKET_NAMES
V_MIN = M.V_MIN


# ----------------------------------------------------------------------
# trajectory look-ahead: integrate ONE ball forward exactly as the
# simulator does (constant deceleration + specular cushion rebounds),
# sampling its centre path so we can choose a strike point.
# ----------------------------------------------------------------------
def freeflight(x, y, vx, vy, t_end, dt=0.0015):
    """Yield (t, x, y, vx, vy) samples of a single ball from (x,y,v) under
    Model M cushion + friction dynamics, until it stops or t_end."""
    s = math.hypot(vx, vy)
    t = 0.0
    while t < t_end and s >= V_MIN:
        # advance dt with sub-stepping across cushion hits
        step = dt
        while step > 1e-9 and s >= V_MIN:
            ux, uy = vx / s, vy / s
            # time to stop
            t_stop = s / A
            tau = min(step, t_stop)
            # detect cushion crossing within tau
            nx = x + vx * tau - 0.5 * A * tau * tau * ux
            ny = y + vy * tau - 0.5 * A * tau * tau * uy
            hit = None
            if nx < R:
                hit = ("x", R)
            elif nx > L - R:
                hit = ("x", L - R)
            if ny < R:
                hit = ("y", R) if hit is None else hit
            elif ny > W - R:
                hit = ("y", W - R) if hit is None else hit
            if hit is None:
                x, y = nx, ny
                ns = s - A * tau
                if ns < V_MIN:
                    s = 0.0
                    vx = vy = 0.0
                else:
                    vx, vy = ux * ns, uy * ns
                    s = ns
                t += tau
                step -= tau
                yield (t, x, y, vx, vy)
            else:
                # bisect to the wall, reflect, continue
                axis, wall = hit
                lo, hi = 0.0, tau
                for _ in range(40):
                    mid = 0.5 * (lo + hi)
                    px = x + vx * mid - 0.5 * A * mid * mid * ux
                    py = y + vy * mid - 0.5 * A * mid * mid * uy
                    cur = px if axis == "x" else py
                    if (cur - wall) * ((x if axis == "x" else y) - wall) > 0:
                        lo = mid
                    else:
                        hi = mid
                th = 0.5 * (lo + hi)
                x = x + vx * th - 0.5 * A * th * th * ux
                y = y + vy * th - 0.5 * A * th * th * uy
                ns = s - A * th
                if ns < V_MIN:
                    s = 0.0
                    vx = vy = 0.0
                    yield (t + th, x, y, vx, vy)
                    return
                vx, vy = ux * ns, uy * ns
                if axis == "x":
                    vx = -vx * M.E_CUSH
                else:
                    vy = -vy * M.E_CUSH
                s = math.hypot(vx, vy)
                t += th
                step -= th
                yield (t, x, y, vx, vy)
    return


def rot(vx, vy, a):
    c, s = math.cos(a), math.sin(a)
    return (c * vx - s * vy, s * vx + c * vy)


def collide(cvx, cvy, ux, uy):
    """Real Model M elastic-with-restitution collision of a carrier moving
    (cvx,cvy) into a stationary ball, line of centres (ux,uy) (unit).
    Returns (struck_vx, struck_vy, carrier_vx, carrier_vy)."""
    vn = cvx * ux + cvy * uy            # carrier speed along line of centres
    # struck gets (1+e)/2 * vn along u ; carrier loses that normal part,
    # keeping tangential + (1-e)/2 * vn along u
    j_struck = 0.5 * (1 + E) * vn
    j_keep = 0.5 * (1 - E) * vn
    svx, svy = j_struck * ux, j_struck * uy
    cvx2 = cvx - vn * ux + j_keep * ux
    cvy2 = cvy - vn * uy + j_keep * uy
    return svx, svy, cvx2, cvy2


def aim_ball_at_pocket(cx, cy, cvx, cvy, gx, gy):
    """Given a carrier at (cx,cy) with velocity (cvx,cvy) about to strike a
    ball, choose the line of centres u so the struck ball departs straight
    at pocket (gx,gy). The struck ball sits at strike_pt + R*u on the far
    side; it departs along u, so u must point from the ball to the pocket.
    We solve: struck ball centre B = C + 2R*u (C = carrier centre at
    contact); departure dir = u; require u parallel to (G - B).
    That gives u = (G - C)/|G - C| to first order; refine by fixed point."""
    ux, uy = gx - cx, gy - cy
    n = math.hypot(ux, uy)
    if n < 1e-9:
        return None
    ux, uy = ux / n, uy / n
    for _ in range(6):
        bx, by = cx + 2 * R * ux, cy + 2 * R * uy
        dx, dy = gx - bx, gy - by
        dn = math.hypot(dx, dy)
        if dn < 1e-9:
            break
        nux, nuy = dx / dn, dy / dn
        if abs(nux - ux) + abs(nuy - uy) < 1e-12:
            ux, uy = nux, nuy
            break
        ux, uy = nux, nuy
    # carrier must actually be approaching along +u (dot with velocity > 0)
    if cvx * ux + cvy * uy <= 0.05 * math.hypot(cvx, cvy):
        return None
    return (ux, uy)


def _candidates(cue0, v0, phi, plan, step, n_target, topm=6):
    """Return up to topm verified next-strike candidates, best first. Each
    is a (score, (bx,by), pocket_idx) tuple whose placement has been
    confirmed in the simulator to (a) occur in order and (b) pocket every
    ball so far including the new one."""
    if step == 0:
        cx, cy = cue0
        cvx, cvy = v0 * math.cos(phi), v0 * math.sin(phi)
    else:
        info = _carrier_after(cue0, v0, phi, plan, step - 1)
        if info is None:
            return []
        _, cx, cy, cvx, cvy, _ = info
    speed = math.hypot(cvx, cvy)
    if speed < 0.8:
        return []
    occupied = [p for p, _ in plan]
    ranked = []
    seen_cells = set()
    for (ts, sx, sy, svx, svy) in freeflight(cx, cy, cvx, cvy,
                                             t_end=speed / A):
        if math.hypot(sx - cx, sy - cy) < 0.10:
            continue
        sp = math.hypot(svx, svy)
        if sp < 0.7:
            break
        for gi, (gx, gy) in enumerate(POCKETS):
            u = aim_ball_at_pocket(sx, sy, svx, svy, gx, gy)
            if u is None:
                continue
            bx, by = sx + 2 * R * u[0], sy + 2 * R * u[1]
            if not (R + 0.02 < bx < L - R - 0.02
                    and R + 0.02 < by < W - R - 0.02):
                continue
            svx2, svy2, cvx2, cvy2 = collide(svx, svy, u[0], u[1])
            dep = math.hypot(svx2, svy2)
            plen = math.hypot(gx - bx, gy - by)
            v_enter = math.sqrt(max(dep * dep - 2 * A * (plen - RHO), 0.0))
            if v_enter < 0.12:
                continue
            ckeep = math.hypot(cvx2, cvy2)
            if step < n_target - 1 and ckeep < 1.0:
                continue
            if any(math.hypot(bx - ox, by - oy) < 2 * R + 0.006
                   for (ox, oy) in occupied):
                continue
            margin = 8.0 if step < n_target - 1 else 0.0
            score = ckeep * margin + v_enter
            ranked.append((score, (bx, by), gi))
    # dedupe spatially and keep the best few, then verify each
    ranked.sort(key=lambda z: -z[0])
    out = []
    for cand in ranked:
        cell = (round(cand[1][0] / 0.05), round(cand[1][1] / 0.05), cand[2])
        if cell in seen_cells:
            continue
        seen_cells.add(cell)
        if _verify_prefix(cue0, v0, phi, plan + [(cand[1], cand[2])],
                          step + 1):
            out.append(cand)
        if len(out) >= topm:
            break
    return out


def build_witness(seed=0, v0=11.0, n_target=21, verbose=False,
                  topm=6, node_budget=4000):
    """Grow a verified chain with depth-first backtracking. At each step we
    take the best few simulator-verified candidates and recurse; on a dead
    end we retreat and try the next. Returns the best chain found (a dict),
    or None."""
    rng = random.Random(seed)
    cue0 = (0.5, W / 2 + rng.uniform(-0.3, 0.3))
    phi = rng.uniform(-0.05, 0.05)

    best = {"n": 0, "plan": []}
    nodes = [0]

    def dfs(plan, step):
        if nodes[0] >= node_budget:
            return
        nodes[0] += 1
        if len(plan) > best["n"]:
            best["n"] = len(plan)
            best["plan"] = list(plan)
            if verbose:
                print(f"  depth {len(plan)} (nodes {nodes[0]})", flush=True)
        if step >= n_target:
            return
        cands = _candidates(cue0, v0, phi, plan, step, n_target, topm)
        for cand in cands:
            dfs(plan + [(cand[1], cand[2])], step + 1)
            if best["n"] >= n_target:
                return

    dfs([], 0)
    if not best["plan"]:
        return None
    return {"cue": cue0, "v0": v0, "phi": phi,
            "balls": [p for p, _ in best["plan"]],
            "pockets": [g for _, g in best["plan"]],
            "n": best["n"]}


def _carrier_after(cue0, v0, phi, plan, k):
    """Simulate the plan, stop immediately after the k-th ball-ball event
    (cue-b0 is k=0, then the carrier strikes b1, b2, ...), and return the
    carrier's live state (name, x, y, vx, vy, t).

    Design: a single carrier travels the whole chain. At collision k it
    strikes a fresh ball b_k, which departs straight at a pocket; the
    carrier is deflected and continues to b_{k+1}. So the carrier is the
    SAME ball throughout: the deflected striker. For k=0 the carrier is the
    cue; after striking b0 the cue is deflected and carries on. Thus the
    carrier's name is always "cue"."""
    balls = [M.Ball(cue0[0], cue0[1], "cue")]
    for i, (p, _) in enumerate(plan):
        balls.append(M.Ball(p[0], p[1], f"b{i}"))
    bmap = {b.name: b for b in balls}
    M.shoot(balls[0], v0, phi)
    sim = M.Simulator(balls, record_events=True)
    seen = 0
    for _ in range(4000):
        if not sim.step():
            break
        ev = sim.events[-1]
        if ev[1] == "bb":
            if seen == k:
                b = bmap["cue"]
                if b.active and b._s > V_MIN:
                    return ("cue", b.x, b.y, b.vx, b.vy, sim.t)
                return None
            seen += 1
    return None


def _verify_prefix(cue0, v0, phi, plan, npots_expected):
    """Simulate the plan and confirm two things: (a) the cue strikes the
    placed balls in order (cue-b0, cue-b1, ...) as a prefix of the realised
    ball-ball events, and (b) every placed ball is actually pocketed. The
    second condition is what makes each committed step a guaranteed pot,
    not merely a planned one."""
    balls = [M.Ball(cue0[0], cue0[1], "cue")]
    for i, (p, _) in enumerate(plan):
        balls.append(M.Ball(p[0], p[1], f"b{i}"))
    M.shoot(balls[0], v0, phi)
    sim = M.Simulator(balls, record_events=True)
    sim.run(max_events=5000, max_time=500)
    intended = [tuple(sorted(("cue", f"b{i}"))) for i in range(len(plan))]
    got = [tuple(sorted(e[2])) for e in sim.events if e[1] == "bb"]
    if len(got) < len(intended):
        return False
    for a, b in zip(got[:len(intended)], intended):
        if a != b:
            return False
    # every placed ball must pocket
    potted = {nm for _, nm, _, _ in sim.pocketed}
    for i in range(len(plan)):
        if f"b{i}" not in potted:
            return False
    return True


if __name__ == "__main__":
    import time
    t0 = time.time()
    found = None
    for seed in range(40):
        w = build_witness(seed=seed, v0=11.0, n_target=6, verbose=False)
        if w:
            found = w
            print(f"seed {seed}: reached {len(w['balls'])} pots "
                  f"({time.time()-t0:.0f}s)")
            break
    if not found:
        print(f"no chain of 6 in 40 seeds ({time.time()-t0:.0f}s)")


def build_greedy(seed=0, v0=11.9, n_target=21, temp=0.4):
    """Fast greedy grower with randomized tie-breaking: among the top few
    verified candidates at each step, pick one weighted by score. Cheap
    enough to restart hundreds of times; the randomness explores the tree
    without the cost of full backtracking."""
    rng = random.Random(seed)
    cue0 = (0.5, W / 2 + rng.uniform(-0.3, 0.3))
    phi = rng.uniform(-0.05, 0.05)
    plan = []
    for step in range(n_target):
        cands = _candidates(cue0, v0, phi, plan, step, n_target, topm=4)
        if not cands:
            break
        # softmax-ish weighted choice among verified candidates
        scores = [c[0] for c in cands]
        m = max(scores)
        weights = [math.exp((s - m) / max(temp, 1e-3)) for s in scores]
        pick = rng.choices(cands, weights=weights, k=1)[0]
        plan.append((pick[1], pick[2]))
    if not plan:
        return None
    return {"cue": cue0, "v0": v0, "phi": phi,
            "balls": [p for p, _ in plan],
            "pockets": [g for _, g in plan], "n": len(plan)}


def _verify_prefix_robust(cue0, v0, phi, plan, dv=0.01, dphi=0.0015):
    """Like _verify_prefix but requires the prefix to pot every placed ball
    for a 3x3 stencil of cue shots around (v0, phi). Committing only robust
    steps builds a witness whose success set has real (positive-measure)
    width, not a knife-edge."""
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if not _verify_prefix(cue0, v0 + i * dv, phi + j * dphi,
                                  plan, len(plan)):
                return False
    return True


def _candidates_robust(cue0, v0, phi, plan, step, n_target, topm=4,
                       dv=0.01, dphi=0.0015):
    """Candidate generator using the robust verifier."""
    if step == 0:
        cx, cy = cue0
        cvx, cvy = v0 * math.cos(phi), v0 * math.sin(phi)
    else:
        info = _carrier_after(cue0, v0, phi, plan, step - 1)
        if info is None:
            return []
        _, cx, cy, cvx, cvy, _ = info
    speed = math.hypot(cvx, cvy)
    if speed < 0.8:
        return []
    occupied = [p for p, _ in plan]
    ranked = []
    for (ts, sx, sy, svx, svy) in freeflight(cx, cy, cvx, cvy,
                                             t_end=speed / A):
        if math.hypot(sx - cx, sy - cy) < 0.10:
            continue
        sp = math.hypot(svx, svy)
        if sp < 0.7:
            break
        for gi, (gx, gy) in enumerate(POCKETS):
            u = aim_ball_at_pocket(sx, sy, svx, svy, gx, gy)
            if u is None:
                continue
            bx, by = sx + 2 * R * u[0], sy + 2 * R * u[1]
            if not (R + 0.02 < bx < L - R - 0.02
                    and R + 0.02 < by < W - R - 0.02):
                continue
            svx2, svy2, cvx2, cvy2 = collide(svx, svy, u[0], u[1])
            dep = math.hypot(svx2, svy2)
            plen = math.hypot(gx - bx, gy - by)
            v_enter = math.sqrt(max(dep * dep - 2 * A * (plen - RHO), 0.0))
            # demand generous pocket-entry margin for robustness
            if v_enter < 0.35:
                continue
            ckeep = math.hypot(cvx2, cvy2)
            if step < n_target - 1 and ckeep < 1.2:
                continue
            if any(math.hypot(bx - ox, by - oy) < 2 * R + 0.02
                   for (ox, oy) in occupied):
                continue
            margin = 8.0 if step < n_target - 1 else 0.0
            score = ckeep * margin + v_enter
            ranked.append((score, (bx, by), gi))
    ranked.sort(key=lambda z: -z[0])
    out = []
    seen = set()
    for cand in ranked:
        cell = (round(cand[1][0] / 0.05), round(cand[1][1] / 0.05), cand[2])
        if cell in seen:
            continue
        seen.add(cell)
        if _verify_prefix_robust(cue0, v0, phi, plan + [(cand[1], cand[2])],
                                 dv, dphi):
            out.append(cand)
        if len(out) >= topm:
            break
    return out


def build_greedy_robust(seed=0, v0=11.9, n_target=21, temp=0.4,
                        dv=0.01, dphi=0.0015):
    """Randomized-greedy grower using the robust (stencil) verifier."""
    rng = random.Random(seed)
    cue0 = (0.5, W / 2 + rng.uniform(-0.3, 0.3))
    phi = rng.uniform(-0.05, 0.05)
    plan = []
    for step in range(n_target):
        cands = _candidates_robust(cue0, v0, phi, plan, step, n_target,
                                   topm=4, dv=dv, dphi=dphi)
        if not cands:
            break
        scores = [c[0] for c in cands]
        m = max(scores)
        weights = [math.exp((s - m) / max(temp, 1e-3)) for s in scores]
        pick = rng.choices(cands, weights=weights, k=1)[0]
        plan.append((pick[1], pick[2]))
    if not plan:
        return None
    return {"cue": cue0, "v0": v0, "phi": phi,
            "balls": [p for p, _ in plan],
            "pockets": [g for _, g in plan], "n": len(plan),
            "dv": dv, "dphi": dphi}
