"""
model_sim.py -- Event-driven simulator for Model M (the dynamics of Section 2).

Model M:
  * Table: rectangle [0,L] x [0,W], six pocket points on the boundary
    (four corners, two side-cushion midpoints).
  * Balls: identical disks of radius R. Between events each moving ball
    decelerates at constant rate A antiparallel to its velocity
    (straight-line motion, quadratic in t).
  * Ball-ball collision: instantaneous, frictionless, equal masses,
    normal restitution E_BALL. Tangential components unchanged.
  * Cushion: normal component reversed and scaled by E_CUSH.
  * Pocket: a ball is captured (removed) when its centre enters the open
    disk of radius RHO about a pocket point.
  * Singular events (simultaneous contacts, grazing) are measure zero;
    the integrator resolves near-ties in time order with tolerances.

All quantities SI.
"""

import numpy as np

# ---------------- model parameters (the "model instance") ----------------
L = 3.569          # playing length (m)
W = 1.778          # playing width  (m)
R = 0.026250       # ball radius (m)  (52.5 mm diameter)
A = 0.10           # deceleration (m/s^2)
E_BALL = 0.95      # ball-ball restitution
E_CUSH = 0.90      # ball-cushion restitution
RHO = 0.055        # pocket capture radius (m)
V_MIN = 1e-4       # speed below which a ball is stopped

POCKETS = np.array([
    [0.0, 0.0], [L, 0.0], [0.0, W], [L, W],   # corners
    [L / 2.0, 0.0], [L / 2.0, W],             # middles
])
POCKET_NAMES = ["C-bl", "C-br", "C-tl", "C-tr", "M-b", "M-t"]

TOL = 1e-10        # time tolerance
EPS_SEP = 1e-9     # post-collision separation nudge


def _real_roots_quartic(c4, c3, c2, c1, c0):
    """Real nonnegative roots of c4 t^4 + ... + c0 = 0."""
    coeffs = np.array([c4, c3, c2, c1, c0], dtype=float)
    nz = np.nonzero(np.abs(coeffs) > 1e-300)[0]
    if len(nz) == 0:
        return []
    coeffs = coeffs[nz[0]:]
    if len(coeffs) == 1:
        return []
    r = np.roots(coeffs)
    out = [x.real for x in r if abs(x.imag) < 1e-8 and x.real > TOL]
    return out


class Ball:
    __slots__ = ("pos", "vel", "active", "name")

    def __init__(self, x, y, name):
        self.pos = np.array([x, y], float)
        self.vel = np.zeros(2)
        self.active = True
        self.name = name

    @property
    def speed(self):
        return float(np.hypot(self.vel[0], self.vel[1]))


class Simulator:
    def __init__(self, balls, record_events=False):
        self.balls = balls
        self.t = 0.0
        self.record = record_events
        self.events = []          # (t, kind, names, extra)
        self.pocketed = []        # (t, ball name, pocket name, entry speed)

    # ---- kinematics on [0, horizon]: quadratic motion coefficients ----
    def _coeffs(self, b):
        """pos(t) = p + v t - 0.5 A t^2 uhat  (valid until stop)."""
        s = b.speed
        if s < V_MIN:
            return b.pos.copy(), np.zeros(2), np.zeros(2), 0.0
        uhat = b.vel / s
        return b.pos.copy(), b.vel.copy(), -0.5 * A * uhat, s / A

    def _pos_at(self, b, dt):
        p, v, q, tstop = self._coeffs(b)
        tau = min(dt, tstop) if tstop > 0 else 0.0
        return p + v * tau + q * tau * tau

    def _vel_at(self, b, dt):
        s = b.speed
        if s < V_MIN:
            return np.zeros(2)
        tau = min(dt, s / A)
        return b.vel * (1.0 - A * tau / s)

    # ---------------- event search ----------------
    def _next_event(self):
        """Return (dt, kind, payload) of earliest event, horizon-capped."""
        moving = [b for b in self.balls if b.active and b.speed >= V_MIN]
        if not moving:
            return None
        # horizon: earliest stop time
        horizon = min(b.speed / A for b in moving)
        best = (horizon, "stop", None)

        active = [b for b in self.balls if b.active]
        # ball-ball
        for i, bi in enumerate(active):
            for bj in active[i + 1:]:
                if bi.speed < V_MIN and bj.speed < V_MIN:
                    continue
                dt = self._pair_time(bi, bj, horizon)
                if dt is not None and dt < best[0]:
                    best = (dt, "bb", (bi, bj))
        # cushion + pocket for moving balls
        for b in moving:
            dt = self._cushion_time(b, horizon)
            if dt is not None and dt < best[0]:
                best = (dt, "cushion", (b, dt))
            dtp = self._pocket_time(b, horizon)
            if dtp is not None and dtp[0] < best[0]:
                best = (dtp[0], "pocket", (b, dtp[1]))
        return best

    def _pair_time(self, bi, bj, horizon):
        pi, vi, qi, ti = self._coeffs(bi)
        pj, vj, qj, tj = self._coeffs(bj)
        # valid joint quadratic window
        tmax = horizon
        for tstop, b in ((ti, bi), (tj, bj)):
            if b.speed >= V_MIN:
                tmax = min(tmax, tstop)
        if tmax <= TOL:
            return None
        dp = pi - pj
        # quick prune
        dist = np.hypot(dp[0], dp[1]) - 2 * R
        vmax = bi.speed + bj.speed
        if dist > vmax * tmax + 1e-12:
            return None
        dv = vi - vj
        dq = qi - qj
        # |dp + dv t + dq t^2|^2 = (2R)^2
        c4 = dq @ dq
        c3 = 2 * dv @ dq
        c2 = dv @ dv + 2 * dp @ dq
        c1 = 2 * dp @ dv
        c0 = dp @ dp - (2 * R) ** 2
        roots = _real_roots_quartic(c4, c3, c2, c1, c0)
        roots = [t for t in roots if t <= tmax + TOL]
        if not roots:
            return None
        t = min(roots)
        # approaching check at contact
        rel_p = dp + dv * t + dq * t * t
        rel_v = dv + 2 * dq * t
        if rel_p @ rel_v >= -1e-12:
            return None
        return t

    def _cushion_time(self, b, horizon):
        p, v, q, tstop = self._coeffs(b)
        tmax = min(horizon, tstop)
        best = None
        for axis, lo, hi in ((0, R, L - R), (1, R, W - R)):
            for wall in (lo, hi):
                # p + v t + q t^2 = wall
                roots = np.roots([q[axis], v[axis], p[axis] - wall]) \
                    if abs(q[axis]) > 1e-300 else \
                    ([-(p[axis] - wall) / v[axis]] if abs(v[axis]) > 1e-300 else [])
                for r_ in np.atleast_1d(roots):
                    r_ = complex(r_)
                    if abs(r_.imag) > 1e-9:
                        continue
                    t = r_.real
                    if TOL < t <= tmax + TOL:
                        # moving toward the wall?
                        vel_ax = v[axis] + 2 * q[axis] * t
                        if (wall == hi and vel_ax > 0) or (wall == lo and vel_ax < 0):
                            if best is None or t < best:
                                best = t
        return best

    def _pocket_time(self, b, horizon):
        p, v, q, tstop = self._coeffs(b)
        tmax = min(horizon, tstop)
        best = None
        for k, c in enumerate(POCKETS):
            dp = p - c
            dist = np.hypot(dp[0], dp[1]) - RHO
            if dist > b.speed * tmax + 1e-12:
                continue
            c4 = q @ q
            c3 = 2 * v @ q
            c2 = v @ v + 2 * dp @ q
            c1 = 2 * dp @ v
            c0 = dp @ dp - RHO ** 2
            roots = _real_roots_quartic(c4, c3, c2, c1, c0)
            roots = [t for t in roots if t <= tmax + TOL]
            if roots:
                t = min(roots)
                if best is None or t < best[0]:
                    best = (t, k)
        return best

    # ---------------- event resolution ----------------
    def _advance(self, dt):
        for b in self.balls:
            if b.active and b.speed >= V_MIN:
                newp = self._pos_at(b, dt)
                newv = self._vel_at(b, dt)
                b.pos, b.vel = newp, newv
                if b.speed < V_MIN:
                    b.vel[:] = 0.0
        self.t += dt

    def step(self):
        ev = self._next_event()
        if ev is None:
            return False
        dt, kind, payload = ev
        self._advance(dt)
        if kind == "stop":
            for b in self.balls:
                if b.active and b.speed < 2 * V_MIN:
                    b.vel[:] = 0.0
            if self.record:
                self.events.append((self.t, "stop", None, None))
        elif kind == "bb":
            bi, bj = payload
            n = bj.pos - bi.pos
            n /= np.hypot(n[0], n[1])
            vi_n = bi.vel @ n
            vj_n = bj.vel @ n
            # equal masses, restitution E_BALL along normal
            vi_n2 = 0.5 * (1 - E_BALL) * vi_n + 0.5 * (1 + E_BALL) * vj_n
            vj_n2 = 0.5 * (1 + E_BALL) * vi_n + 0.5 * (1 - E_BALL) * vj_n
            bi.vel += (vi_n2 - vi_n) * n
            bj.vel += (vj_n2 - vj_n) * n
            # nudge apart to avoid re-detection
            bi.pos -= EPS_SEP * n
            bj.pos += EPS_SEP * n
            if self.record:
                self.events.append((self.t, "bb", (bi.name, bj.name), None))
        elif kind == "cushion":
            b, _ = payload
            for axis, lo, hi in ((0, R, L - R), (1, R, W - R)):
                if b.pos[axis] <= lo + 1e-7 and b.vel[axis] < 0:
                    b.vel[axis] *= -E_CUSH
                    b.pos[axis] = lo + EPS_SEP
                elif b.pos[axis] >= hi - 1e-7 and b.vel[axis] > 0:
                    b.vel[axis] *= -E_CUSH
                    b.pos[axis] = hi - EPS_SEP
            if self.record:
                self.events.append((self.t, "cushion", b.name, None))
        elif kind == "pocket":
            b, k = payload
            b.active = False
            spd = b.speed
            b.vel[:] = 0.0
            self.pocketed.append((self.t, b.name, POCKET_NAMES[k], spd))
            if self.record:
                self.events.append((self.t, "pocket", b.name, POCKET_NAMES[k]))
        return True

    def run(self, max_events=5000, max_time=200.0):
        n = 0
        while n < max_events and self.t < max_time:
            moving = any(b.active and b.speed >= V_MIN for b in self.balls)
            if not moving:
                break
            if not self.step():
                break
            n += 1
        return n


# ---------------- configurations ----------------
def opening_rack(eps=5e-4, cue_pos=(0.65, W / 2.0)):
    """Standard opening configuration, reds separated by eps (Section 4)."""
    balls = []
    cue = Ball(cue_pos[0], cue_pos[1], "cue")
    balls.append(cue)
    # colours on their spots
    spots = {
        "yellow": (0.737, W / 2 - 0.292),
        "green": (0.737, W / 2 + 0.292),
        "brown": (0.737, W / 2),
        "blue": (L / 2, W / 2),
        "pink": (2.667, W / 2),
        "black": (L - 0.324, W / 2),
    }
    for name, (x, y) in spots.items():
        balls.append(Ball(x, y, name))
    # reds: triangle, apex just up-table of pink
    pitch = 2 * R + eps
    apex_x = 2.667 + 2 * R + eps
    k = 0
    for row in range(5):
        x = apex_x + row * pitch * np.sqrt(3) / 2
        for j in range(row + 1):
            y = W / 2 + (j - row / 2.0) * pitch
            k += 1
            balls.append(Ball(x, y, f"red{k}"))
    return balls, cue


def shoot(cue, v0, phi):
    cue.vel = np.array([v0 * np.cos(phi), v0 * np.sin(phi)])


if __name__ == "__main__":
    import time
    balls, cue = opening_rack()
    shoot(cue, 8.0, 0.0)
    s = Simulator(balls)
    t0 = time.time()
    n = s.run()
    print(f"events={n}  wall={time.time()-t0:.2f}s  sim_t={s.t:.1f}s  "
          f"pocketed={[(b, p) for _, b, p, _ in s.pocketed]}")
