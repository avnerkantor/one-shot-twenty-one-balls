"""
trunk_cushion.py -- Long trunk that FOLDS via cushion rebounds.

Key idea (resolves the "trunk runs off the table" collapse): between two
cut nodes, the moving trunk ball may rebound off one cushion. Rebounds
keep 91-96% of speed at grazing incidence, so the chain can be arbitrarily
long while staying on the table. Aiming and self-pot geometry are unchanged
from the straight-trunk case; only the carrier PATH between nodes changes.

We reuse model_sim for cushion physics ground truth, but compute rebound
geometry analytically here for the search.
"""
import math
import model_sim as M

L, W, R, A, E_CUSH = M.L, M.W, M.R, M.A, M.E_CUSH
LO, HI_X, HI_Y = R, L - R, W - R   # centre bounds


def reflect_point_dir(p, d, wall):
    """Given a ball at p heading d (unit), find where it first meets the
    given wall centre-line and the post-rebound direction. wall in
    {'x_lo','x_hi','y_lo','y_hi'}. Returns (hit_point, new_dir, dist) or
    None if it doesn't reach that wall going forward."""
    x, y = p
    dx, dy = d
    if wall == 'x_lo':
        if dx >= -1e-9: return None
        t = (LO - x) / dx
        ny = y + dy * t
        if not (LO <= ny <= HI_Y): return None
        return (LO, ny), (-dx, dy), t
    if wall == 'x_hi':
        if dx <= 1e-9: return None
        t = (HI_X - x) / dx
        ny = y + dy * t
        if not (LO <= ny <= HI_Y): return None
        return (HI_X, ny), (-dx, dy), t
    if wall == 'y_lo':
        if dy >= -1e-9: return None
        t = (LO - y) / dy
        nx = x + dx * t
        if not (LO <= nx <= HI_X): return None
        return (nx, LO), (dx, -dy), t
    if wall == 'y_hi':
        if dy <= 1e-9: return None
        t = (HI_Y - y) / dy
        nx = x + dx * t
        if not (LO <= nx <= HI_X): return None
        return (nx, HI_Y), (dx, -dy), t
    return None


def rebound_speed(v_in, d, wall):
    """Speed after a cushion rebound with normal restitution E_CUSH."""
    dx, dy = d
    if wall in ('x_lo', 'x_hi'):
        vn, vt = abs(dx) * v_in * E_CUSH, abs(dy) * v_in
    else:
        vn, vt = abs(dy) * v_in * E_CUSH, abs(dx) * v_in
    return math.hypot(vn, vt)


def carry_with_optional_bounce(p, d, v, target_reach):
    """Yield reachable (endpoint, dir, v_at_end, path_segments) options:
    (a) straight to a point target_reach away, if in bounds;
    (b) one rebound off each reachable wall, then continue.
    Each option is a dict with keys pos,d,v,segs (list of (a,b) pts)."""
    opts = []
    # straight option handled by caller's node grid; here we do bounces.
    for wall in ('x_lo', 'x_hi', 'y_lo', 'y_hi'):
        r = reflect_point_dir(p, d, wall)
        if r is None:
            continue
        hit, nd, dist = r
        if dist < 0.05:
            continue
        v_hit = math.sqrt(max(v * v - 2 * A * dist, 0.0))
        if v_hit < 1.2:
            continue
        v_out = rebound_speed(v_hit, d, wall)
        opts.append({"bounce": hit, "pos": hit, "d": nd, "v": v_out,
                     "pre_dist": dist})
    return opts


if __name__ == "__main__":
    # sanity: a ball heading up-right off the top cushion folds down-right
    p = (2.0, 1.5); d = (0.5, 0.5); dn = math.hypot(*d); d = (d[0]/dn, d[1]/dn)
    for o in carry_with_optional_bounce(p, d, 5.0, 0.5):
        print(f"bounce at ({o['bounce'][0]:.2f},{o['bounce'][1]:.2f}) "
              f"new_dir=({o['d'][0]:.2f},{o['d'][1]:.2f}) v={o['v']:.2f}")


def dist_to_next_wall(p, d):
    """Forward distance from p (unit dir d) until the ball centre reaches
    any cushion centre-line. Used to cap a post-bounce segment so it does
    not incur a SECOND bounce before the next collision."""
    x, y = p
    dx, dy = d
    best = float("inf")
    if dx > 1e-9:
        best = min(best, (HI_X - x) / dx)
    elif dx < -1e-9:
        best = min(best, (LO - x) / dx)
    if dy > 1e-9:
        best = min(best, (HI_Y - y) / dy)
    elif dy < -1e-9:
        best = min(best, (LO - y) / dy)
    return best
