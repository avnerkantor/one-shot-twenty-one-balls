"""
witness_tree.py -- Construct the witness (c*, x*) as a collision tree.

Architecture:
  cue --full--> T1;  T1 cuts T2 (T1 becomes carrier of branch 1), T2 cuts
  T3, ..., T5 cuts T6 (carrier 5), T6 runs straight into a pocket.
  Each carrier services a local branch: it cuts A (A departs along the
  line of centres, aimed exactly at a pocket), deflects to B, cuts B
  (aimed), deflects to C, cuts C (aimed), and finally runs into a pocket
  itself (alignment solved by tuning the last hop length).

Counts: 6 trunk + 5 x 3 branch = 21 object balls, 21 pots.

Trunk cuts are thin (the struck ball keeps 0.975 cos(theta) of speed);
branch cuts are fat (the carrier keeps sin(theta)). All aiming is exact
for restitution E; friction enters through the speed budget only.
"""

import math
import json
import model_sim as M
import trunk_cushion as _tc
from witness_build import (rot, ang, norm_ang, friction_speed, travel_time,
                           outgoing, seg_clear, in_bounds)

R, A, E, RHO = M.R, M.A, M.E_BALL, M.RHO
L, W = M.L, M.W
POCKETS, PNAMES = M.POCKETS, M.POCKET_NAMES

V0 = 11.5
CUE = (0.40, W / 2)
T1 = (0.62, W / 2)
N_TRUNK = 12        # trunk balls; carriers = N_TRUNK - 1; branch balls = 21 - N_TRUNK
CLEAR = 2 * R + 0.007
V_ENTER_MIN = 0.10
POCKET_PATH_MARGIN = RHO + 0.008


def solve_aim_struck(node_from, dir_in, sigma, g_idx, th_lo=0.08, th_hi=1.45):
    """Solve theta so the STRUCK ball's ray (from its own centre, along
    the line of centres u) points at pocket g. Node position is the
    striker's centre at contact; struck ball sits at node + 2R u."""
    gx, gy = POCKETS[g_idx]
    dx, dy = dir_in
    nx, ny = node_from

    def err(th):
        ux, uy = rot(dx, dy, -sigma * th)
        px, py = nx + 2 * R * ux, ny + 2 * R * uy
        return norm_ang(ang(gx - px, gy - py) - ang(ux, uy))

    n = 36
    ths = [th_lo + (th_hi - th_lo) * i / (n - 1) for i in range(n)]
    vals = [err(t) for t in ths]
    for i in range(n - 1):
        if vals[i] * vals[i + 1] < 0 and \
                abs(vals[i]) + abs(vals[i + 1]) < math.pi:
            # the second condition rejects sign flips across the +-pi
            # wrap, which would converge onto an anti-aligned direction
            lo, hi, flo = ths[i], ths[i + 1], vals[i]
            for _ in range(70):
                mid = 0.5 * (lo + hi)
                if err(mid) * flo > 0:
                    lo = mid
                else:
                    hi = mid
            root = 0.5 * (lo + hi)
            if abs(err(root)) < 1e-3:
                return root
    return None



def seg_clear_x(p0, p1, pts):
    """Path clearance ignoring contact partners at the segment start."""
    import math as _m
    keep = [q for q in pts if _m.hypot(q[0] - p0[0], q[1] - p0[1]) > 2.1 * R]
    if not keep:
        return 9.9
    return seg_clear(p0, p1, keep)

def pocket_clear_path(p0, p1, skip=()):
    """Segment must stay away from all pocket capture disks."""
    for k, c in enumerate(POCKETS):
        if k in skip:
            continue
        if seg_clear(p0, p1, [c]) < POCKET_PATH_MARGIN:
            return False
    return True


def _mk_pot(node, u, v_struck, g_idx):
    px, py = node[0] + 2 * R * u[0], node[1] + 2 * R * u[1]
    gx, gy = POCKETS[g_idx]
    plen = math.hypot(gx - px, gy - py)
    v_enter = friction_speed(v_struck, plen - RHO)
    return (px, py), plen, v_enter


def _branch_expand_last(st, static_pts):
    """Last branch level: tune hop length s continuously so that after
    cutting C (aim pinned to C's pocket) the carrier's outgoing ray
    passes exactly through some pocket centre."""
    px, py = st["pos"]
    dx, dy = st["d"]
    outs = []
    for gC in range(6):
        for sigma in (+1, -1):
            for gF in range(6):
                def signed_miss(s):
                    v_arr = friction_speed(st["v"], s)
                    if v_arr < 0.45:
                        return None
                    node = (px + s * dx, py + s * dy)
                    if not in_bounds(node):
                        return None
                    th = solve_aim_struck(node, (dx, dy), sigma, gC)
                    if th is None or not (0.12 < th < 1.40):
                        return None
                    ox, oy = outgoing(dx, dy, th, sigma)[:2]
                    gx, gy = POCKETS[gF]
                    rx, ry = gx - node[0], gy - node[1]
                    if rx * ox + ry * oy <= 0.05:
                        return None
                    return ox * ry - oy * rx
                prev = None
                i = 0
                while i < 60:
                    s = 0.10 + i * 0.022
                    m = signed_miss(s)
                    if m is not None and prev is not None                             and prev[1] * m < 0:
                        lo, hi, flo = prev[0], s, prev[1]
                        for _ in range(60):
                            mid = 0.5 * (lo + hi)
                            fm = signed_miss(mid)
                            if fm is None:
                                break
                            if fm * flo > 0:
                                lo = mid
                            else:
                                hi = mid
                        s_star = 0.5 * (lo + hi)
                        cand = _finish_branch(st, static_pts, s_star,
                                              sigma, gC, gF)
                        if cand is not None:
                            outs.append(cand)
                    prev = (s, m) if m is not None else None
                    i += 1
    return outs


def _finish_branch(st, static_pts, s, sigma, gC, gF):
    px, py = st["pos"]
    dx, dy = st["d"]
    v_arr = friction_speed(st["v"], s)
    node = (px + s * dx, py + s * dy)
    th = solve_aim_struck(node, (dx, dy), sigma, gC)
    if th is None:
        return None
    ox, oy, ofac, ux, uy, sfac = outgoing(dx, dy, th, sigma)
    v_struck = v_arr * sfac
    P, plen, v_enter = _mk_pot(node, (ux, uy), v_struck, gC)
    if v_enter < V_ENTER_MIN or plen > 2.6 or not in_bounds(P):
        return None
    v_car = v_arr * ofac
    gx, gy = POCKETS[gF]
    fwd = (gx - node[0]) * ox + (gy - node[1]) * oy
    run = fwd - RHO
    v_fin = friction_speed(v_car, run)
    if v_fin < V_ENTER_MIN:
        return None
    allpts = static_pts + st["balls"]
    if seg_clear_x((px, py), node, allpts) < CLEAR:
        return None
    if seg_clear_x(P, POCKETS[gC], allpts) < CLEAR:
        return None
    if seg_clear_x(node, (gx, gy), allpts + [P]) < CLEAR:
        return None
    if min((math.hypot(P[0] - q[0], P[1] - q[1]) for q in allpts),
           default=9) < CLEAR:
        return None
    if not pocket_clear_path((px, py), node):
        return None
    if not pocket_clear_path(P, POCKETS[gC], skip=(gC,)):
        return None
    if not pocket_clear_path(node, (gx, gy), skip=(gF,)):
        return None
    t_node = st["t"] + travel_time(st["v"], s)
    stt = dict(st)
    stt["balls"] = st["balls"] + [P]
    stt["plan"] = st["plan"] + [{
        "node": node, "s": s, "theta": th, "sigma": sigma, "pocket": gC,
        "v_arr": v_arr, "v_struck": v_struck, "v_enter": v_enter,
        "t": t_node, "pot_len": plen, "ball": P}]
    stt["final"] = {"pocket": gF, "run": run, "v_enter": v_fin,
                    "miss": 0.0, "t": t_node + travel_time(v_car, run)}
    stt["score"] = st["score"] + v_fin + v_enter
    return stt


def build_branch(entry_pos, entry_dir, v_in, t_in, static_pts, depth=3):
    """Mini beam search: carrier pots `depth` balls then itself."""
    S_GRID = (0.12, 0.16, 0.21, 0.27, 0.34, 0.42, 0.52, 0.65, 0.85)
    init = {"pos": entry_pos, "d": entry_dir, "v": v_in, "t": t_in,
            "balls": [], "plan": [], "score": 0.0}
    beam = [init]
    for lvl in range(depth - 1):
        nxt = []
        for st in beam:
            px, py = st["pos"]
            dx, dy = st["d"]
            for s in S_GRID:
                v_arr = friction_speed(st["v"], s)
                if v_arr < 0.45:
                    continue
                node = (px + s * dx, py + s * dy)
                if not in_bounds(node):
                    continue
                if not pocket_clear_path((px, py), node):
                    continue
                for g in range(6):
                    for sigma in (+1, -1):
                        th = solve_aim_struck(node, (dx, dy), sigma, g)
                        if th is None or not (0.12 < th < 1.40):
                            continue
                        ox, oy, ofac, ux, uy, sfac = \
                            outgoing(dx, dy, th, sigma)
                        v_struck = v_arr * sfac
                        P, plen, v_enter = _mk_pot(node, (ux, uy),
                                                   v_struck, g)
                        if v_enter < V_ENTER_MIN or plen > 2.6:
                            continue
                        if not in_bounds(P):
                            continue
                        v_car = v_arr * ofac
                        if v_car < 0.45:
                            continue
                        allpts = static_pts + st["balls"]
                        if seg_clear_x((px, py), node, allpts) < CLEAR:
                            continue
                        if seg_clear_x(P, POCKETS[g], allpts) < CLEAR:
                            continue
                        if min((math.hypot(P[0] - q[0], P[1] - q[1])
                                for q in allpts), default=9) < CLEAR:
                            continue
                        if not pocket_clear_path(P, POCKETS[g], skip=(g,)):
                            continue
                        t_node = st["t"] + travel_time(st["v"], s)
                        nxt.append({
                            "pos": node, "d": (ox, oy), "v": v_car,
                            "t": t_node, "balls": st["balls"] + [P],
                            "plan": st["plan"] + [{
                                "node": node, "s": s, "theta": th,
                                "sigma": sigma, "pocket": g,
                                "v_arr": v_arr, "v_struck": v_struck,
                                "v_enter": v_enter, "t": t_node,
                                "pot_len": plen, "ball": P}],
                            "score": v_car + 2.0 * min(v_enter, 0.6),
                        })
        if not nxt:
            return None
        nxt.sort(key=lambda z: -z["score"])
        seen = {}
        beam = []
        for cand in nxt:
            key = (round(cand["pos"][0] / 0.15),
                   round(cand["pos"][1] / 0.15),
                   round(math.atan2(cand["d"][1], cand["d"][0]) / 0.52),
                   cand["plan"][-1]["pocket"])
            c = seen.get(key, 0)
            if c < 2:
                seen[key] = c + 1
                beam.append(cand)
            if len(beam) >= 40:
                break
    # last level: continuous alignment of the carrier's final ray
    finals = []
    for st in beam:
        finals.extend(_branch_expand_last(st, static_pts))
    finals.sort(key=lambda z: -z["score"])
    out, seen = [], {}
    for f in finals:
        b0 = f["balls"][0]
        key = (f["final"]["pocket"], round(b0[0] / 0.2),
               round(b0[1] / 0.2))
        if seen.get(key, 0) < 2:
            seen[key] = seen.get(key, 0) + 1
            f = dict(f)
            f["entry"] = entry_pos
            out.append(f)
        if len(out) >= 6:
            break
    return out


def carrier_dirs_exact(dx, dy, theta, sigma):
    return outgoing(dx, dy, theta, sigma)




def solve_selfpot_theta(prev, d, sigma, g_idx, s):
    """theta such that the deflected carrier's ray from the trunk node
    (at distance s along d from prev) passes through pocket g."""
    node = (prev[0] + s * d[0], prev[1] + s * d[1])
    gx, gy = POCKETS[g_idx]

    def err(th):
        ox, oy = outgoing(d[0], d[1], th, sigma)[:2]
        rx, ry = gx - node[0], gy - node[1]
        if rx * ox + ry * oy <= 0.05:
            return None
        return ox * ry - oy * rx

    prev_v = None
    for i in range(70):
        th = 0.14 + i * (1.20 / 69)
        m = err(th)
        if m is not None and prev_v is not None and prev_v[1] * m < 0:
            lo, hi, flo = prev_v[0], th, prev_v[1]
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                fm = err(mid)
                if fm is None:
                    break
                if fm * flo > 0:
                    lo = mid
                else:
                    hi = mid
            return 0.5 * (lo + hi), node
        prev_v = (th, m) if m is not None else None
    return None

def carrier_has_target(entry, d):
    """A carrier is viable only if, within a few hop lengths, some
    pocket is reachable by a first cut (theta window not degenerate)."""
    for hop in (0.18, 0.30, 0.45, 0.62, 0.85, 1.10, 1.40):
        node = (entry[0] + hop * d[0], entry[1] + hop * d[1])
        if not in_bounds(node):
            continue
        for g in range(6):
            for sigma in (+1, -1):
                th = solve_aim_struck(node, d, sigma, g)
                if th is not None and 0.12 < th < 1.40:
                    return True
    return False


_CDP_CACHE = {}


def _carrier_prefilter(node, cd, v_car):
    """Cheap necessary condition: some pocket lies ahead within a plausible
    cone and reach, so a depth-1 branch is at least geometrically possible."""
    reach = v_car * v_car / (2 * A)      # max roll distance
    for gx, gy in POCKETS:
        rx, ry = gx - node[0], gy - node[1]
        dist = math.hypot(rx, ry)
        if dist > reach + 0.4 or dist < 0.10:
            continue
        # within a wide forward cone of the carrier heading?
        if rx * cd[0] + ry * cd[1] > -0.30 * dist:
            return True
    return False


def carrier_depth1_possible(entry, d, v):
    """Coarse existence test for a depth-1 branch: some (cut ball pocket,
    side, final pocket) shows a sign change of the aligned-self-pot miss
    over hop length. Cached by rounded state."""
    ck = (round(entry[0], 2), round(entry[1], 2),
          round(math.atan2(d[1], d[0]), 2), round(v, 1))
    hit = _CDP_CACHE.get(ck)
    if hit is not None:
        return hit
    res = _carrier_depth1_core(entry, d, v)
    _CDP_CACHE[ck] = res
    return res


def _carrier_depth1_core(entry, d, v):
    signs = {}
    for i in range(9):
        s = 0.12 + i * 0.15
        v_arr = friction_speed(v, s)
        if v_arr < 0.45:
            break
        node = (entry[0] + s * d[0], entry[1] + s * d[1])
        if not in_bounds(node):
            continue
        for gC in range(6):
            for sigma in (+1, -1):
                th = solve_aim_struck(node, d, sigma, gC)
                if th is None or not (0.12 < th < 1.40):
                    continue
                ox, oy = outgoing(d[0], d[1], th, sigma)[:2]
                for gF in range(6):
                    gx, gy = POCKETS[gF]
                    rx, ry = gx - node[0], gy - node[1]
                    if rx * ox + ry * oy <= 0.05:
                        continue
                    m = 1 if ox * ry - oy * rx > 0 else -1
                    key = (gC, sigma, gF)
                    if key in signs and signs[key] != m:
                        return True
                    signs[key] = m
    return False

def forward_room(p, d):
    """Distance from p along heading d until the trunk ray exits the
    in-bounds region (table minus wall/pocket margins)."""
    x, y = p
    best = 0.0
    for i in range(1, 40):
        t = i * 0.06
        q = (x + d[0] * t, y + d[1] * t)
        if in_bounds(q):
            best = t
        else:
            break
    return best


def build_trunk_skeletons(v0=V0, want=120, filter_live=True):
    """Trunk beam through four cut nodes, then a direct solve of the
    fifth node: for each level-4 state, side, cut angle, and pocket,
    the segment length aligning T6's departure ray with the pocket
    centre is closed form (the miss is linear in s)."""
    # The cue strikes T1 as a CUT (not full): T1 departs along the line of
    # centres, the cue deflects to the opposite side and leaves the chain
    # axis, so it never chases a bounced trunk ball back into the pack. We
    # sweep the initial cut angle and side. For each, T1 is PLACED so that
    # the line of centres at contact has the required bearing (the cue
    # travels +x from CUE and contacts T1 after a fixed run).
    TH_GRID = (0.20, 0.30, 0.42, 0.56, 0.72)
    S_GRID = (0.18, 0.26, 0.36, 0.48, 0.62)
    cue_run = T1[0] - CUE[0]          # cue centre travel before contact
    beam = []
    for th0 in (0.30, 0.42, 0.56, 0.72, 0.90):
        for sg0 in (+1, -1):
            ux, uy = rot(1.0, 0.0, -sg0 * th0)
            ox, oy, ofac = outgoing(1.0, 0.0, th0, sg0)[:3]
            Xc = CUE[0] + cue_run                 # cue centre at contact
            T1p = (Xc + 2 * R * ux, CUE[1] + 2 * R * uy)
            if not in_bounds(T1p):
                continue
            v_t1 = v0 * 0.5 * (1 + E) * math.cos(th0)
            if v_t1 < 3.0:
                continue
            t1 = travel_time(v0, cue_run)
            beam.append({
                "pos": T1p, "d": (ux, uy), "v": v_t1, "t": t1,
                "placed": [CUE, T1p], "trunk_plan": [], "score": 0.0,
                "nbounce": 0,
                "cue_out": (ox, oy), "cue_v": v0 * ofac})
    for k in range(1, N_TRUNK - 1):
        nxt = []
        for st in beam:
            nb = st.get("nbounce", 0)
            # start states: straight from current carrier, or one bounce
            starts = [(st["pos"], st["d"], st["v"], None)]
            if nb < 3:
                for o in _tc.carry_with_optional_bounce(
                        st["pos"], st["d"], st["v"], 0):
                    starts.append((o["pos"], o["d"], o["v"],
                                   (st["pos"], o["bounce"])))
            for (px, py), (dx, dy), vst, bounce in starts:
                # cap the segment so the carrier reaches the node WITHOUT
                # touching a cushion in between; after a bounce this also
                # prevents a spurious second bounce. Keeps the analytic path
                # exact and simulator-consistent.
                wall_cap = _tc.dist_to_next_wall((px, py), (dx, dy)) \
                    - 2 * R - 0.02
                for s in S_GRID:
                    if s > wall_cap:
                        continue
                    v_arr = friction_speed(vst, s)
                    if v_arr < 1.2:
                        continue
                    node = (px + s * dx, py + s * dy)
                    if not in_bounds(node):
                        continue
                    if not pocket_clear_path((px, py), node):
                        continue
                    for th in TH_GRID:
                        for sigma in (+1, -1):
                            ox, oy, ofac, ux, uy, sfac = \
                                outgoing(dx, dy, th, sigma)
                            v_next = v_arr * sfac
                            v_car = v_arr * ofac
                            if v_next < 1.0 or v_car < 1.2:
                                continue
                            P = (node[0] + 2 * R * ux,
                                 node[1] + 2 * R * uy)
                            if not in_bounds(P):
                                continue
                            allpts = list(st["placed"])
                            # clearance from the (possibly post-bounce)
                            # approach point to the node
                            if seg_clear_x((px, py), node,
                                           allpts) < CLEAR:
                                continue
                            if min(math.hypot(P[0] - q[0], P[1] - q[1])
                                   for q in allpts) < CLEAR:
                                continue
                            if not _carrier_prefilter(
                                    node, (ox, oy), v_car):
                                continue
                            if not carrier_depth1_possible(
                                    node, (ox, oy), v_car):
                                continue
                            t_node = st["t"] + travel_time(vst, s)
                            nxt.append({
                                "pos": P, "d": (ux, uy), "v": v_next,
                                "t": t_node,
                                "placed": st["placed"] + [P],
                                "nbounce": nb + (1 if bounce else 0),
                                "trunk_plan": st["trunk_plan"] + [{
                                    "node": node, "s": s, "theta": th,
                                    "sigma": sigma, "v_arr": v_arr,
                                    "v_next": v_next, "v_car": v_car,
                                    "t": t_node, "selfpot": None,
                                    "bounce": bounce,
                                    "approach_from": (px, py),
                                    "carrier_dir": (ox, oy)}],
                                "score": v_next
                                + 0.5 * min(v_car, 3.0)
                                - 0.6 * abs(node[1] - W / 2) / W,
                            })
        if not nxt:
            return []
        nxt.sort(key=lambda z: -z["score"])
        seen = {}
        beam = []
        for st in nxt:
            key = (round(st["pos"][0] / 0.24), round(st["pos"][1] / 0.24),
                   round(math.atan2(st["d"][1], st["d"][0]) / 0.5),
                   st.get("nbounce", 0))
            c = seen.get(key, 0)
            if c < 3:
                seen[key] = c + 1
                beam.append(st)
            if len(beam) >= 90:
                break

    # final node + straight-in alignment, from the best level-(N-2) states
    outs = []
    TH_FINE = [0.22 + i * (1.08 / 17) for i in range(18)]
    beam.sort(key=lambda z: -z["score"])
    for st in beam[:80]:
        prev = st["pos"]
        d5 = st["d"]
        v_pre = st["v"]
        for sig in (+1, -1):
            for th in TH_FINE:
                ox, oy, ofac, ux, uy, sfac = \
                    outgoing(d5[0], d5[1], th, sig)
                cross_ud = ux * d5[1] - uy * d5[0]
                if abs(cross_ud) < 1e-9:
                    continue
                for g in range(6):
                    gx, gy = POCKETS[g]
                    rx, ry = gx - prev[0], gy - prev[1]
                    s_star = (ux * ry - uy * rx) / cross_ud
                    if not (0.20 <= s_star <= 1.30):
                        continue
                    v_arr = friction_speed(v_pre, s_star)
                    if v_arr < 1.2:
                        continue
                    node = (prev[0] + s_star * d5[0],
                            prev[1] + s_star * d5[1])
                    if not in_bounds(node):
                        continue
                    P6 = (node[0] + 2 * R * ux, node[1] + 2 * R * uy)
                    if not in_bounds(P6):
                        continue
                    fwd = (gx - P6[0]) * ux + (gy - P6[1]) * uy
                    if fwd <= 0.2:
                        continue
                    v_next = v_arr * sfac
                    v_car = v_arr * ofac
                    if v_next < 1.0 or v_car < 1.5:
                        continue
                    run = fwd - RHO
                    v_enter = friction_speed(v_next, run)
                    if v_enter < V_ENTER_MIN:
                        continue
                    statics = st["placed"]
                    if seg_clear_x(prev, node, statics) < CLEAR:
                        continue
                    if seg_clear_x(P6, (gx, gy), statics) < CLEAR:
                        continue
                    if not pocket_clear_path(P6, (gx, gy), skip=(g,)):
                        continue
                    if not pocket_clear_path(prev, node):
                        continue
                    if min(math.hypot(P6[0] - q[0], P6[1] - q[1])
                           for q in statics) < CLEAR:
                        continue
                    # T_N is a pure straight-in pot; no branch hangs off it
                    t_node = st["t"] + travel_time(v_pre, s_star)
                    s2 = dict(st)
                    s2["trunk_plan"] = st["trunk_plan"] + [{
                        "node": node, "s": s_star, "theta": th,
                        "sigma": sig, "v_arr": v_arr, "v_next": v_next,
                        "v_car": v_car, "t": t_node, "selfpot": None,
                        "carrier_dir": (ox, oy), "no_branch": True}]
                    s2["placed"] = st["placed"] + [P6]
                    s2["pos"] = P6
                    s2["t6_final"] = {"pocket": g, "run": run,
                                      "v_enter": v_enter, "miss": 0.0}
                    s2["score2"] = st["score"] + v_enter
                    outs.append(s2)
    outs.sort(key=lambda z: -z["score2"])
    seen = {}
    final = []
    for st in outs:
        n5 = st["trunk_plan"][-1]["node"]
        key = (st["t6_final"]["pocket"],
               tuple(tp["sigma"] for tp in st["trunk_plan"]),
               round(n5[0] / 0.15), round(n5[1] / 0.15))
        c = seen.get(key, 0)
        if c < 2:
            seen[key] = c + 1
            if (not filter_live) or all(
                    carrier_depth1_possible(tp["node"], tp["carrier_dir"],
                                            tp["v_car"])
                    for tp in st["trunk_plan"]):
                final.append(st)
        if len(final) >= want:
            break
    return final



def seg_near_frac(a, b, p):
    """(distance, fraction along a->b) of the closest point to p."""
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    if L2 <= 0:
        return math.hypot(p[0] - ax, p[1] - ay), 0.0
    f = ((p[0] - ax) * dx + (p[1] - ay) * dy) / L2
    f = max(0.0, min(1.0, f))
    return math.hypot(p[0] - (ax + f * dx), p[1] - (ay + f * dy)), f


def branch_timed_elems(br):
    """(segments with time windows, balls with departure times)."""
    segs = []
    t_prev = br["entry_t"]
    prev = br["entry"]
    v = br["entry_v"]
    for pl in br["plan"]:
        segs.append((prev, pl["node"], t_prev, pl["t"]))
        prev, t_prev = pl["node"], pl["t"]
        segs.append((pl["ball"], POCKETS[pl["pocket"]], pl["t"],
                     pl["t"] + travel_time(pl["v_struck"],
                                           max(pl["pot_len"] - RHO, 0.01))))
    segs.append((prev, POCKETS[br["final"]["pocket"]], t_prev,
                 br["final"]["t"]))
    balls = [(pl["ball"], pl["t"]) for pl in br["plan"]]
    return segs, balls


T_MARGIN = 0.05     # seconds of required temporal separation



SEG_SEG_CLEAR = 2 * R + 0.002


def _seg_kin(seg):
    """Exact kinematics of a timed segment under deceleration A:
    returns (a, unit_dir, v0, t0, t1)."""
    a, c, t0, t1 = seg
    dx, dy = c[0] - a[0], c[1] - a[1]
    d = math.hypot(dx, dy)
    T = max(t1 - t0, 1e-9)
    v0 = d / T + 0.5 * A * T
    if d > 0:
        u = (dx / d, dy / d)
    else:
        u = (0.0, 0.0)
    return a, u, v0, t0, t1, d


def _seg_pos(kin, t):
    a, u, v0, t0, t1, d = kin
    tau = min(max(t, t0), t1) - t0
    s = min(v0 * tau - 0.5 * A * tau * tau, d)
    return (a[0] + u[0] * s, a[1] + u[1] * s)


def seg_seg_conflict(s1, s2):
    """Exact space-time closest approach of two timed segments under
    constant deceleration: relative motion is quadratic in t, so the
    squared distance is a quartic minimized in closed form."""
    lo = max(s1[2], s2[2])
    hi = min(s1[3], s2[3])
    if lo >= hi:
        return False
    from model_sim import _cubic_real_roots

    def coeffs(seg):
        a, u, v0, t0, t1, d = _seg_kin(seg)
        c0 = (a[0] + u[0] * (-v0 * t0 - 0.5 * A * t0 * t0),
              a[1] + u[1] * (-v0 * t0 - 0.5 * A * t0 * t0))
        c1 = (u[0] * (v0 + A * t0), u[1] * (v0 + A * t0))
        c2 = (-0.5 * A * u[0], -0.5 * A * u[1])
        return c0, c1, c2

    a0, a1, a2 = coeffs(s1)
    b0, b1, b2 = coeffs(s2)
    P = (a0[0] - b0[0], a0[1] - b0[1])
    V = (a1[0] - b1[0], a1[1] - b1[1])
    Q = (a2[0] - b2[0], a2[1] - b2[1])
    # f(t) = |P + V t + Q t^2|^2 ; f'(t)/2 = cubic
    c3 = 2 * (Q[0] * Q[0] + Q[1] * Q[1])
    c2_ = 3 * (V[0] * Q[0] + V[1] * Q[1])
    c1_ = (V[0] * V[0] + V[1] * V[1]) + 2 * (P[0] * Q[0] + P[1] * Q[1])
    c0_ = P[0] * V[0] + P[1] * V[1]
    cand = [lo, hi] + [t for t in _cubic_real_roots(c3, c2_, c1_, c0_)
                       if lo < t < hi]
    thr2 = SEG_SEG_CLEAR * SEG_SEG_CLEAR
    for t in cand:
        dx = P[0] + V[0] * t + Q[0] * t * t
        dy = P[1] + V[1] * t + Q[1] * t * t
        if dx * dx + dy * dy < thr2:
            return True
    return False


def branch_conflicts(br, other_balls, other_segs):
    """Static conflict test with timing: a segment only conflicts with a
    ball that is still on its spot when the segment is traversed."""
    segs, balls = branch_timed_elems(br)
    def passage_time(seg, frac):
        a, u, v0, t0, t1, dd = _seg_kin(seg)
        dist = frac * dd
        disc = v0 * v0 - 2 * A * dist
        tau = (v0 - math.sqrt(max(disc, 0.0))) / A
        return t0 + tau

    for (bp, bt) in balls:
        for (qp, qt) in other_balls:
            if math.hypot(bp[0] - qp[0], bp[1] - qp[1]) < CLEAR:
                return True                    # coexist from t = 0
        for seg in other_segs:
            a, c, t0, t1 = seg
            d, f = seg_near_frac(a, c, bp)
            if d < CLEAR and math.hypot(bp[0] - a[0], bp[1] - a[1]) \
                    > 2.1 * R:
                if passage_time(seg, f) < bt + T_MARGIN:
                    return True
    for seg in segs:
        a, c, t0, t1 = seg
        for (qp, qt) in other_balls:
            d, f = seg_near_frac(a, c, qp)
            if d < CLEAR and math.hypot(qp[0] - a[0], qp[1] - a[1]) \
                    > 2.1 * R:
                if passage_time(seg, f) < qt + T_MARGIN:
                    return True
    return False


def sim_check(skel, got):
    """Replay a candidate assembly in the event simulator; accept only a
    clean 21-ball clearance with exactly the planned contacts."""
    import model_sim as MS
    trunk_pos = skel["placed"][1:]
    balls = [MS.Ball(CUE[0], CUE[1], "cue")]
    for i, p in enumerate(trunk_pos, 1):
        balls.append(MS.Ball(p[0], p[1], f"T{i}"))
    planned = {("T1", "cue")}
    for i in range(1, len(trunk_pos)):
        planned.add(tuple(sorted((f"T{i}", f"T{i+1}"))))
    for i, br in got:
        for ci, p in enumerate(br["balls"]):
            nm = f"B{i+1}{'ABCDE'[ci]}"
            balls.append(MS.Ball(p[0], p[1], nm))
            planned.add(tuple(sorted((f"T{i+1}", nm))))
    MS.shoot(balls[0], V0, 0.0)
    sim = MS.Simulator(balls, record_events=True)
    sim.run(max_events=3000, max_time=400.0)
    n_obj = sum(1 for _, nm, _, _ in sim.pocketed if nm != "cue")
    realised = {tuple(sorted(ev[2])) for ev in sim.events if ev[1] == "bb"}
    return n_obj == 21 and realised == planned



def cue_phantoms():
    """The cue ball creeps forward after the full hit and parks on the
    trunk line; model its path as permanent phantom obstacles."""
    v_res = 0.5 * (1 - E) * V0
    x0 = T1[0] - 2 * R                      # cue centre at contact
    travel = v_res * v_res / (2 * A)
    pts = []
    k = 0
    while True:
        x = x0 + min(k * 0.08, travel)
        pts.append((x, T1[1]))
        if k * 0.08 >= travel:
            break
        k += 1
    return pts


def trunk_timed_segments(skel):
    """True approach segments of the trunk balls (from each ball's
    initial position to its collision node) and the last ball's
    straight run, with traversal time windows."""
    segs = []
    t_prev = travel_time(V0, math.dist(CUE, T1) - 2 * R)
    for j, tp in enumerate(skel["trunk_plan"]):
        start = skel["placed"][j + 1]        # T_{j+1} initial position
        segs.append((start, tp["node"], t_prev, tp["t"]))
        t_prev = tp["t"]
    tf = skel["t6_final"]
    gp = POCKETS[tf["pocket"]]
    last_t = skel["trunk_plan"][-1]["t"]
    segs.append((skel["placed"][-1], gp, last_t,
                 last_t + travel_time(skel["trunk_plan"][-1]["v_next"],
                                      max(tf["run"], 0.01))))
    return segs


def assignment_crossings_ok(got):
    """Reject assignments where two branches' moving pieces would meet:
    space-time closest approach between segments of different branches."""
    seg_sets = []
    for _, v in got:
        seg_sets.append(branch_timed_elems(v)[0])
    for i in range(len(seg_sets)):
        for j in range(i + 1, len(seg_sets)):
            for s1 in seg_sets[i]:
                for s2 in seg_sets[j]:
                    if seg_seg_conflict(s1, s2):
                        return False
    return True

def attach_branches(skel, max_sim_tries=400):
    trunk = skel["trunk_plan"]
    phantoms = cue_phantoms()
    base_static = list(skel["placed"]) + phantoms
    need = 21 - N_TRUNK
    carrier_idx = [i for i, tp in enumerate(trunk)
                   if not tp.get("no_branch")]
    variants = {}
    menus = []
    for i in carrier_idx:
        tp = trunk[i]
        menu = []
        for d in (1, 2, 3):
            vs = build_branch(tp["node"], tp["carrier_dir"], tp["v_car"],
                              tp["t"], base_static, depth=d)
            if not vs:
                continue
            for v in vs:
                v["entry_t"] = tp["t"]
                v["entry_v"] = tp["v_car"]
            menu.append(d)
            variants[(i, d)] = vs
        if not menu:
            return None
        menus.append(menu)
    if sum(max(m) for m in menus) < need or \
            sum(min(m) for m in menus) > need:
        return None

    order = sorted(range(len(carrier_idx)), key=lambda m: len(menus[m]))
    rebuilds = [0]
    MAX_REBUILDS = 30

    def dfs(pos, got, balls_acc, segs_acc, remaining):
        if pos == len(order):
            if remaining == 0:
                yield got
            return
        m = order[pos]
        i = carrier_idx[m]
        rest = order[pos + 1:]
        lo = sum(min(menus[j]) for j in rest)
        hi = sum(max(menus[j]) for j in rest)
        tp = trunk[i]
        for d in sorted(menus[m], key=lambda x: -x):
            if not (lo <= remaining - d <= hi):
                continue
            cands = [v for v in variants[(i, d)]
                     if not branch_conflicts(v, balls_acc, segs_acc)]
            if not cands and rebuilds[0] < MAX_REBUILDS:
                rebuilds[0] += 1
                fresh = build_branch(tp["node"], tp["carrier_dir"],
                                     tp["v_car"], tp["t"],
                                     base_static +
                                     [b for b, _ in balls_acc],
                                     depth=d)
                if fresh:
                    for v in fresh:
                        v["entry_t"] = tp["t"]
                        v["entry_v"] = tp["v_car"]
                    cands = [v for v in fresh
                             if not branch_conflicts(v, balls_acc,
                                                     segs_acc)]
            for v in cands[:4]:
                segs, balls = branch_timed_elems(v)
                yield from dfs(pos + 1, got + [(i, v)],
                               balls_acc + balls, segs_acc + segs,
                               remaining - d)
        return

    seed_balls = [(p, 1e9) for p in phantoms]
    seed_segs = trunk_timed_segments(skel)
    tried = 0
    generated = 0
    for got in dfs(0, [], seed_balls, seed_segs, need):
        generated += 1
        if generated > 20000:
            return None
        if not assignment_crossings_ok(got):
            continue
        tried += 1
        cand = sorted(got, key=lambda z: z[0])
        if sim_check(skel, cand):
            branches = []
            for i, br in cand:
                tp = trunk[i]
                branches.append({"balls": br["balls"], "plan": br["plan"],
                                 "final": br["final"],
                                 "entry_v": tp["v_car"],
                                 "entry_t": tp["t"]})
            out = dict(skel)
            out["branches"] = branches
            return out
        if tried >= max_sim_tries:
            return None
    return None


def build_all(v0=V0):
    skels = build_trunk_skeletons(v0)
    print(f"skeletons: {len(skels)}")
    for i, sk in enumerate(skels):
        sol = attach_branches(sk)
        if sol is not None:
            print(f"skeleton {i} succeeded")
            return sol, None
        print(f"skeleton {i} failed branch attachment")
    return None, "no skeleton admitted branches"


if __name__ == "__main__":
    import time
    t0 = time.time()
    sol, fail = build_all()
    if sol is None:
        print(f"FAILED at trunk node {fail}")
    else:
        print(f"BUILT in {time.time()-t0:.0f}s")
        nballs = 2 + len(sol["placed"]) - 2 + \
            sum(len(br["balls"]) for br in sol["branches"])
        print("trunk nodes:", len(sol["trunk_plan"]),
              " branches:", len(sol["branches"]),
              " total object balls:",
              len(sol["placed"]) - 1 +
              sum(len(br["balls"]) for br in sol["branches"]))
        for i, tp in enumerate(sol["trunk_plan"], 1):
            print(f"  trunk {i}: s={tp['s']:.2f} th={math.degrees(tp['theta']):.0f} "
                  f"v_next={tp['v_next']:.2f} v_car={tp['v_car']:.2f}")
        for i, br in enumerate(sol["branches"], 1):
            pots = [PNAMES[p["pocket"]] for p in br["plan"]]
            fin = br["final"]
            print(f"  branch {i}: pots {pots} + self->{PNAMES[fin['pocket']]}"
                  f" (enter {fin['v_enter']:.2f}, miss {fin['miss']*1000:.0f}mm)")
        tf = sol["t6_final"]
        print(f"  T6 -> {PNAMES[tf['pocket']]} enter {tf['v_enter']:.2f} "
              f"miss {tf['miss']*1000:.0f}mm")
        with open("witness_tree.json", "w") as f:
            json.dump(sol, f, indent=1, default=list)
        print("saved witness_tree.json")
