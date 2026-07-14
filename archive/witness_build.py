"""
witness_build.py -- Construct the witness configuration c* and shot x*
of Theorem 1 (Lemma: witness construction).

Design: a sequential plant chain of all 21 object balls.
  * The cue ball strikes B1 full; B1 departs down-table.
  * At node k (k = 1..20), ball Bk cuts into ball B(k+1): B(k+1) departs
    along the line of centres (the next chain segment), Bk deflects into
    an assigned pocket. Aiming is solved exactly for restitution E, so
    Bk's post-collision ray points at the pocket centre.
  * B21 departs along the last line of centres, aimed directly at a
    pocket (alignment solved by tuning the last segment length).

Beam search over (segment length, pocket, side) with feasibility
filters: cut-angle window, table margins, static path clearance, speed
budget. Verification is done separately (witness_verify.py).
"""

import math
import model_sim as M

R, A, E, RHO = M.R, M.A, M.E_BALL, M.RHO
L, W = M.L, M.W
POCKETS, PNAMES = M.POCKETS, M.POCKET_NAMES

V0 = 11.5
CUE = (0.60, W / 2)
B1 = (1.00, W / 2)
THETA_MIN, THETA_MAX = 0.26, 1.00
MARGIN_WALL = 0.085
MARGIN_POCKET = RHO + 0.030
CLEAR = 2 * R + 0.010
V_CHAIN_MIN = 0.55
V_ENTER_MIN = 0.10
N_CHAIN = 21


def rot(vx, vy, a):
    c, s = math.cos(a), math.sin(a)
    return (c * vx - s * vy, s * vx + c * vy)


def ang(vx, vy):
    return math.atan2(vy, vx)


def norm_ang(a):
    while a > math.pi:
        a -= 2 * math.pi
    while a < -math.pi:
        a += 2 * math.pi
    return a


def friction_speed(v, dist):
    s2 = v * v - 2 * A * dist
    return math.sqrt(s2) if s2 > 0 else 0.0


def travel_time(v, dist):
    vf = friction_speed(v, dist)
    return (v - vf) / A


def outgoing(dx, dy, theta, sigma):
    """Striker outgoing direction/speed factor and struck-ball data."""
    ux, uy = rot(dx, dy, -sigma * theta)
    wx, wy = rot(ux, uy, sigma * math.pi / 2)
    st, ct = math.sin(theta), math.cos(theta)
    kx = st * wx + 0.5 * (1 - E) * ct * ux
    ky = st * wy + 0.5 * (1 - E) * ct * uy
    kn = math.hypot(kx, ky)
    return kx / kn, ky / kn, kn, ux, uy, 0.5 * (1 + E) * ct


def solve_theta(dx, dy, sigma, target_bearing):
    """theta such that the striker's outgoing ray has the given bearing.

    With restitution e < 1 the deflection angle is non-monotone in theta
    (a small residual forward component dominates near theta = 0), so
    there may be two roots; we take the larger one, which is the robust
    moderate-cut branch.
    """
    beta = norm_ang(target_bearing - ang(dx, dy))
    if sigma * beta <= 0 or abs(beta) >= math.pi / 2:
        return None

    def err(th):
        ox, oy = outgoing(dx, dy, th, sigma)[:2]
        return norm_ang(ang(ox, oy) - ang(dx, dy)) - beta

    n = 240
    ths = [1e-4 + (math.pi / 2 - 2e-4) * i / (n - 1) for i in range(n)]
    vals = [err(t) for t in ths]
    root = None
    for i in range(n - 1):
        if vals[i] == 0.0:
            root = ths[i]
        elif vals[i] * vals[i + 1] < 0:
            lo, hi, flo = ths[i], ths[i + 1], vals[i]
            for _ in range(70):
                mid = 0.5 * (lo + hi)
                if err(mid) * flo > 0:
                    lo = mid
                else:
                    hi = mid
            root = 0.5 * (lo + hi)   # keep scanning: last root = largest
    return root


def seg_clear(p0, p1, obstacles, exclude=()):
    x0, y0 = p0
    x1, y1 = p1
    dx, dy = x1 - x0, y1 - y0
    seg2 = dx * dx + dy * dy
    dmin = float("inf")
    for idx, (ox, oy) in enumerate(obstacles):
        if idx in exclude:
            continue
        t = ((ox - x0) * dx + (oy - y0) * dy) / seg2 if seg2 > 0 else 0.0
        t = max(0.0, min(1.0, t))
        d = math.hypot(ox - (x0 + t * dx), oy - (y0 + t * dy))
        if d < dmin:
            dmin = d
    return dmin


def in_bounds(p):
    x, y = p
    if not (MARGIN_WALL < x < L - MARGIN_WALL and
            MARGIN_WALL < y < W - MARGIN_WALL):
        return False
    return all(math.hypot(x - c[0], y - c[1]) >= MARGIN_POCKET
               for c in POCKETS)


def _expand(st, k, S_GRID):
    out = []
    px, py = st["placed"][-1]
    dx, dy = st["d"]
    for s in S_GRID:
        v_arr = friction_speed(st["v"], s)
        if v_arr < V_CHAIN_MIN:
            continue
        Lx, Ly = px + s * dx, py + s * dy
        if not in_bounds((Lx, Ly)):
            continue
        for g_idx in range(6):
            gx, gy = POCKETS[g_idx]
            pot_len = math.hypot(gx - Lx, gy - Ly)
            if pot_len > 3.30:
                continue
            pot_bearing = ang(gx - Lx, gy - Ly)
            for sigma in (+1, -1):
                th = solve_theta(dx, dy, sigma, pot_bearing)
                if th is None or not (THETA_MIN < th < THETA_MAX):
                    continue
                ox, oy, ofac, ux, uy, sfac = outgoing(dx, dy, th, sigma)
                v_pot = v_arr * ofac
                v_enter = friction_speed(v_pot, pot_len - RHO)
                if v_enter < V_ENTER_MIN:
                    continue
                v_next = v_arr * sfac
                if v_next < V_CHAIN_MIN:
                    continue
                Px, Py = Lx + 2 * R * ux, Ly + 2 * R * uy
                if not in_bounds((Px, Py)):
                    continue
                last = len(st["placed"]) - 1
                if seg_clear((px, py), (Lx, Ly), st["placed"],
                             exclude=(last,)) < CLEAR:
                    continue
                if seg_clear((Lx, Ly), (gx, gy), st["placed"],
                             exclude=(last,)) < CLEAR:
                    continue
                if min(math.hypot(Px - qx, Py - qy)
                       for qx, qy in st["placed"][:-1]) < CLEAR:
                    continue
                if any(seg_clear(pl["node"], POCKETS[pl["pocket"]],
                                 [(Px, Py)]) < CLEAR
                       for pl in st["plan"]):
                    continue
                head = ang(ux, uy)
                t_node = st["t"] + travel_time(st["v"], s)
                plan = st["plan"] + [{
                    "k": k, "node": (Lx, Ly), "s": s, "theta": th,
                    "sigma": sigma, "pocket": g_idx, "v_arr": v_arr,
                    "v_pot": v_pot, "v_enter": v_enter,
                    "v_next": v_next, "t_node": t_node,
                    "pot_len": pot_len,
                }]
                score = (v_next + 3.0 * min(v_enter, 0.5))
                out.append({
                    "placed": st["placed"] + [(Px, Py)],
                    "d": (ux, uy), "v": v_next, "t": t_node,
                    "plan": plan, "score": score,
                })
    return out


def solve_final(st, g20, sigma, g21):
    """Tune last segment length so B21's ray points at pocket g21."""
    px, py = st["placed"][-1]
    dx, dy = st["d"]
    g20x, g20y = POCKETS[g20]
    g21x, g21y = POCKETS[g21]

    def geometry(s):
        v_arr = friction_speed(st["v"], s)
        if v_arr < V_CHAIN_MIN:
            return None
        Lx, Ly = px + s * dx, py + s * dy
        if not in_bounds((Lx, Ly)):
            return None
        th = solve_theta(dx, dy, sigma, ang(g20x - Lx, g20y - Ly))
        if th is None or not (THETA_MIN < th < THETA_MAX):
            return None
        ox, oy, ofac, ux, uy, sfac = outgoing(dx, dy, th, sigma)
        Px, Py = Lx + 2 * R * ux, Ly + 2 * R * uy
        err = norm_ang(ang(g21x - Px, g21y - Py) - ang(ux, uy))
        return err, th, (Lx, Ly), (Px, Py), (ux, uy), v_arr, ofac, sfac

    prev = None
    for i in range(72):
        s = 0.14 + i * 0.01
        cur = geometry(s)
        if cur is None:
            prev = None
            continue
        if prev is not None and prev[0][0] * cur[0] < 0:
            lo, hi, flo = prev[1], s, prev[0][0]
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                gm = geometry(mid)
                if gm is None:
                    break
                if gm[0] * flo > 0:
                    lo = mid
                else:
                    hi = mid
            gm = geometry(0.5 * (lo + hi))
            if gm is not None:
                res = _check_final(st, gm, 0.5 * (lo + hi),
                                   g20, sigma, g21, (px, py))
                if res is not None:
                    return res
        prev = (cur, s)
    return None


def _check_final(st, gm, s_fin, g20, sigma, g21, prev_pos):
    err, th, node, P21, u, v_arr, ofac, sfac = gm
    g20x, g20y = POCKETS[g20]
    g21x, g21y = POCKETS[g21]
    v_pot = v_arr * ofac
    pot_len = math.hypot(g20x - node[0], g20y - node[1])
    v_enter20 = friction_speed(v_pot, pot_len - RHO)
    v21 = v_arr * sfac
    run_len = math.hypot(g21x - P21[0], g21y - P21[1])
    v_enter21 = friction_speed(v21, run_len - RHO)
    if not (v_enter20 > V_ENTER_MIN and v_enter21 > V_ENTER_MIN
            and pot_len < 3.30):
        return None
    placed = st["placed"] + [P21]
    exc20 = (len(placed) - 2,)
    c1 = seg_clear(prev_pos, node, placed[:-1], exclude=exc20)
    c2 = seg_clear(node, (g20x, g20y), placed[:-1], exclude=exc20)
    c3 = seg_clear(P21, (g21x, g21y), placed[:-1])
    c4 = min(math.hypot(P21[0] - qx, P21[1] - qy)
             for qx, qy in placed[:-1])
    if min(c1, c2, c3, c4) < CLEAR:
        return None
    if any(seg_clear(pl["node"], POCKETS[pl["pocket"]], [P21]) < CLEAR
           for pl in st["plan"]):
        return None
    t_node = st["t"] + travel_time(st["v"], s_fin)
    plan = st["plan"] + [{
        "k": N_CHAIN - 1, "node": node, "s": s_fin, "theta": th,
        "sigma": sigma, "pocket": g20, "v_arr": v_arr, "v_pot": v_pot,
        "v_enter": v_enter20, "v_next": v21, "t_node": t_node,
        "pot_len": pot_len,
    }, {
        "k": N_CHAIN, "final_run": True, "pocket": g21,
        "run_len": run_len, "v_enter": v_enter21,
        "t_node": t_node + travel_time(v21, run_len),
    }]
    return {"placed": placed, "plan": plan, "cue": CUE, "V0": V0,
            "score": v_enter21 + v_enter20 + 0.2 * min(c1, c2, c3, c4)}


def build_chain(beam_width=160):
    d0 = (1.0, 0.0)
    init = {"placed": [CUE, B1], "d": d0,
            "v": 0.5 * (1 + E) * V0,
            "t": travel_time(V0, math.dist(CUE, B1) - 2 * R),
            "plan": [], "score": 0.0}
    S_GRID = (0.16, 0.20, 0.25, 0.30, 0.36, 0.44, 0.55, 0.70, 0.85, 1.05)
    beam = [init]
    for k in range(1, N_CHAIN - 1):
        nxt = []
        for st in beam:
            nxt.extend(_expand(st, k, S_GRID))
        if not nxt:
            return None, k
        nxt.sort(key=lambda z: -z["score"])
        # diversity: cap states per spatial/heading bucket
        buckets = {}
        beam = []
        for st in nxt:
            px, py = st["placed"][-1]
            import math as _m
            key = (round(px / 0.18), round(py / 0.18),
                   round(_m.atan2(st["d"][1], st["d"][0]) / 0.5236))
            c = buckets.get(key, 0)
            if c < 4:
                buckets[key] = c + 1
                beam.append(st)
            if len(beam) >= beam_width:
                break
    comps = []
    for st in beam:
        for g20 in range(6):
            for sigma in (+1, -1):
                for g21 in range(6):
                    sol = solve_final(st, g20, sigma, g21)
                    if sol is not None:
                        comps.append(sol)
    if not comps:
        return None, N_CHAIN - 1
    comps.sort(key=lambda z: -z["score"])
    return comps[0], None


if __name__ == "__main__":
    import time, json
    t0 = time.time()
    sol, fail = build_chain()
    if sol is None:
        print(f"FAILED at node {fail}")
    else:
        print(f"built in {time.time()-t0:.1f}s  score={sol['score']:.3f}")
        for pl in sol["plan"]:
            if pl.get("final_run"):
                print(f"  B21 straight run {pl['run_len']:.3f} m -> "
                      f"{PNAMES[pl['pocket']]}  enter {pl['v_enter']:.2f}")
            else:
                print(f"  node {pl['k']:2d}  s={pl['s']:.3f}  "
                      f"theta={math.degrees(pl['theta']):5.1f}  "
                      f"-> {PNAMES[pl['pocket']]:4s}  "
                      f"v_arr={pl['v_arr']:5.2f}  "
                      f"v_enter={pl['v_enter']:5.2f}  "
                      f"v_next={pl['v_next']:5.2f}")
        with open("witness_layout.json", "w") as f:
            json.dump({"cue": sol["cue"], "V0": sol["V0"],
                       "placed": sol["placed"],
                       "plan": [{kk: vv for kk, vv in pl.items()}
                                for pl in sol["plan"]]}, f, indent=1)
        print("saved witness_layout.json")
