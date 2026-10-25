"""
model_sim.py -- Fast event-driven simulator for Model M (Section 2 dynamics).

Model M:
  * Table: rectangle [0,L] x [0,W]; six pocket points on the boundary
    (four corners, two long-cushion midpoints).
  * Balls: identical disks of radius R. Between events each moving ball
    decelerates at constant rate A antiparallel to its velocity, so its
    centre follows p(t) = p0 + v0 t - (A/2) t^2 vhat until it stops.
  * Ball-ball collision: instantaneous, frictionless, equal masses,
    normal restitution E_BALL; tangential components unchanged.
  * Cushion: normal velocity component reversed, scaled by E_CUSH.
  * Pocket: a ball is captured (removed) when its centre enters the open
    disk of radius RHO about a pocket point.

Root finding: contact conditions are quartic in t on a bounded window;
we split the window at the (closed-form) critical points of the quartic
and bisect on monotone pieces. No linear-algebra calls in the hot path.
"""

import math

# ---------------- model instance ----------------
L = 3.569
W = 1.778
R = 0.026250
A = 0.10
E_BALL = 0.95
E_CUSH = 0.90
RHO = 0.055
V_MIN = 1e-4

POCKETS = ((0.0, 0.0), (L, 0.0), (0.0, W), (L, W), (L / 2, 0.0), (L / 2, W))
POCKET_NAMES = ("C-bl", "C-br", "C-tl", "C-tr", "M-b", "M-t")

TOL = 1e-10
EPS_SEP = 1e-9
FOUR_R2 = (2 * R) ** 2
RHO2 = RHO ** 2


def _cubic_real_roots(a, b, c, d):
    """Real roots of a x^3 + b x^2 + c x + d = 0."""
    if abs(a) < 1e-300:
        if abs(b) < 1e-300:
            if abs(c) < 1e-300:
                return ()
            return (-d / c,)
        disc = c * c - 4 * b * d
        if disc < 0:
            return ()
        s = math.sqrt(disc)
        return ((-c + s) / (2 * b), (-c - s) / (2 * b))
    b, c, d = b / a, c / a, d / a
    p = c - b * b / 3.0
    q = 2 * b ** 3 / 27.0 - b * c / 3.0 + d
    off = -b / 3.0
    disc = (q / 2) ** 2 + (p / 3) ** 3
    if disc > 0:
        s = math.sqrt(disc)
        u = math.copysign(abs(-q / 2 + s) ** (1 / 3), -q / 2 + s)
        v = math.copysign(abs(-q / 2 - s) ** (1 / 3), -q / 2 - s)
        return (u + v + off,)
    if abs(p) < 1e-300:
        return (off,)
    m = 2 * math.sqrt(-p / 3.0)
    arg = 3 * q / (p * m)
    arg = max(-1.0, min(1.0, arg))
    theta = math.acos(arg) / 3.0
    return (m * math.cos(theta) + off,
            m * math.cos(theta - 2 * math.pi / 3) + off,
            m * math.cos(theta + 2 * math.pi / 3) + off)


def _first_root_quartic(c4, c3, c2, c1, c0, tmax):
    """Smallest t in (TOL, tmax] with f(t)=0, f(0)>0 assumed; None if none.

    f(t) = c4 t^4 + c3 t^3 + c2 t^2 + c1 t + c0.
    Splits [0, tmax] at critical points of f, bisects monotone pieces.
    """
    def f(t):
        return (((c4 * t + c3) * t + c2) * t + c1) * t + c0

    crits = [t for t in _cubic_real_roots(4 * c4, 3 * c3, 2 * c2, c1)
             if TOL < t < tmax]
    crits.sort()
    knots = [TOL] + crits + [tmax]
    fa = f(knots[0])
    for i in range(len(knots) - 1):
        a_, b_ = knots[i], knots[i + 1]
        fb = f(b_)
        if fa > 0.0 and fb <= 0.0:
            lo, hi = a_, b_
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if f(mid) > 0.0:
                    lo = mid
                else:
                    hi = mid
                if hi - lo < 1e-14:
                    break
            return hi
        fa = fb
    return None


class Ball:
    __slots__ = ("x", "y", "vx", "vy", "active", "name", "_s")

    def __init__(self, x, y, name):
        self.x, self.y = x, y
        self.vx = self.vy = 0.0
        self.active = True
        self.name = name
        self._s = 0.0

    def set_vel(self, vx, vy):
        self.vx, self.vy = vx, vy
        self._s = math.hypot(vx, vy)

    @property
    def speed(self):
        return self._s


class Simulator:
    def __init__(self, balls, record_events=False):
        self.balls = balls
        self.t = 0.0
        self.record = record_events
        self.events = []
        self.pocketed = []

    def _advance(self, dt):
        for b in self.balls:
            s = b._s
            if not b.active or s < V_MIN:
                continue
            tau = dt if dt < s / A else s / A
            ux, uy = b.vx / s, b.vy / s
            b.x += b.vx * tau - 0.5 * A * tau * tau * ux
            b.y += b.vy * tau - 0.5 * A * tau * tau * uy
            ns = s - A * tau
            if ns < V_MIN:
                b.set_vel(0.0, 0.0)
            else:
                b.set_vel(ux * ns, uy * ns)
        self.t += dt

    def _next_event(self):
        moving = [b for b in self.balls if b.active and b._s >= V_MIN]
        if not moving:
            return None
        horizon = min(b._s / A for b in moving)
        best_t, best_kind, best_pay = horizon, "stop", None

        active = [b for b in self.balls if b.active]
        n = len(active)
        for i in range(n):
            bi = active[i]
            si = bi._s
            mi = si >= V_MIN
            for j in range(i + 1, n):
                bj = active[j]
                sj = bj._s
                if not mi and sj < V_MIN:
                    continue
                # window on which both quadratics are valid
                tmax = best_t
                if mi and si / A < tmax:
                    tmax = si / A
                if sj >= V_MIN and sj / A < tmax:
                    tmax = sj / A
                if tmax <= TOL:
                    continue
                dpx, dpy = bi.x - bj.x, bi.y - bj.y
                gap = math.hypot(dpx, dpy) - 2 * R
                if gap > (si + sj) * tmax:
                    continue
                dvx, dvy = bi.vx - bj.vx, bi.vy - bj.vy
                qix = qiy = qjx = qjy = 0.0
                if mi:
                    qix, qiy = -0.5 * A * bi.vx / si, -0.5 * A * bi.vy / si
                if sj >= V_MIN:
                    qjx, qjy = -0.5 * A * bj.vx / sj, -0.5 * A * bj.vy / sj
                dqx, dqy = qix - qjx, qiy - qjy
                c4 = dqx * dqx + dqy * dqy
                c3 = 2 * (dvx * dqx + dvy * dqy)
                c2 = dvx * dvx + dvy * dvy + 2 * (dpx * dqx + dpy * dqy)
                c1 = 2 * (dpx * dvx + dpy * dvy)
                c0 = dpx * dpx + dpy * dpy - FOUR_R2
                if c0 <= 0.0:      # already overlapping: skip (post-nudge)
                    continue
                t = _first_root_quartic(c4, c3, c2, c1, c0, tmax)
                if t is None:
                    continue
                # approaching?
                rpx = dpx + dvx * t + dqx * t * t
                rpy = dpy + dvy * t + dqy * t * t
                rvx = dvx + 2 * dqx * t
                rvy = dvy + 2 * dqy * t
                if rpx * rvx + rpy * rvy >= -1e-12:
                    continue
                best_t, best_kind, best_pay = t, "bb", (bi, bj)

        for b in moving:
            s = b._s
            tmax = min(best_t, s / A)
            if tmax <= TOL:
                continue
            ux, uy = b.vx / s, b.vy / s
            qx, qy = -0.5 * A * ux, -0.5 * A * uy
            # cushions: axis-aligned quadratics
            for p0, v0, q0, lo, hi in ((b.x, b.vx, qx, R, L - R),
                                       (b.y, b.vy, qy, R, W - R)):
                for wall in (lo, hi):
                    cc, bb_, aa = p0 - wall, v0, q0
                    if abs(aa) < 1e-300:
                        roots = (-cc / bb_,) if abs(bb_) > 1e-300 else ()
                    else:
                        disc = bb_ * bb_ - 4 * aa * cc
                        if disc < 0:
                            roots = ()
                        else:
                            sq = math.sqrt(disc)
                            roots = ((-bb_ + sq) / (2 * aa),
                                     (-bb_ - sq) / (2 * aa))
                    for t in roots:
                        if TOL < t <= tmax:
                            vel_ax = v0 + 2 * q0 * t
                            if (wall == hi and vel_ax > 0) or \
                               (wall == lo and vel_ax < 0):
                                if t < best_t:
                                    best_t, best_kind, best_pay = \
                                        t, "cushion", b
                                tmax = t
            # pockets
            for k in range(6):
                cx, cy = POCKETS[k]
                dpx, dpy = b.x - cx, b.y - cy
                if math.hypot(dpx, dpy) - RHO > s * tmax:
                    continue
                c4 = qx * qx + qy * qy
                c3 = 2 * (b.vx * qx + b.vy * qy)
                c2 = b.vx * b.vx + b.vy * b.vy + 2 * (dpx * qx + dpy * qy)
                c1 = 2 * (dpx * b.vx + dpy * b.vy)
                c0 = dpx * dpx + dpy * dpy - RHO2
                if c0 <= 0.0:
                    return (TOL, "pocket", (b, k))
                t = _first_root_quartic(c4, c3, c2, c1, c0, min(tmax, best_t))
                if t is not None and t < best_t:
                    best_t, best_kind, best_pay = t, "pocket", (b, k)
        return best_t, best_kind, best_pay

    def step(self):
        ev = self._next_event()
        if ev is None:
            return False
        dt, kind, pay = ev
        self._advance(dt)
        if kind == "stop":
            for b in self.balls:
                if b.active and b._s < 2 * V_MIN:
                    b.set_vel(0.0, 0.0)
        elif kind == "bb":
            bi, bj = pay
            nx, ny = bj.x - bi.x, bj.y - bi.y
            d = math.hypot(nx, ny)
            nx, ny = nx / d, ny / d
            vin = bi.vx * nx + bi.vy * ny
            vjn = bj.vx * nx + bj.vy * ny
            vin2 = 0.5 * (1 - E_BALL) * vin + 0.5 * (1 + E_BALL) * vjn
            vjn2 = 0.5 * (1 + E_BALL) * vin + 0.5 * (1 - E_BALL) * vjn
            bi.set_vel(bi.vx + (vin2 - vin) * nx, bi.vy + (vin2 - vin) * ny)
            bj.set_vel(bj.vx + (vjn2 - vjn) * nx, bj.vy + (vjn2 - vjn) * ny)
            bi.x -= EPS_SEP * nx
            bi.y -= EPS_SEP * ny
            bj.x += EPS_SEP * nx
            bj.y += EPS_SEP * ny
            if self.record:
                self.events.append((self.t, "bb", (bi.name, bj.name), None))
        elif kind == "cushion":
            b = pay
            if b.x <= R + 1e-7 and b.vx < 0:
                b.set_vel(-E_CUSH * b.vx, b.vy)
                b.x = R + EPS_SEP
            elif b.x >= L - R - 1e-7 and b.vx > 0:
                b.set_vel(-E_CUSH * b.vx, b.vy)
                b.x = L - R - EPS_SEP
            if b.y <= R + 1e-7 and b.vy < 0:
                b.set_vel(b.vx, -E_CUSH * b.vy)
                b.y = R + EPS_SEP
            elif b.y >= W - R - 1e-7 and b.vy > 0:
                b.set_vel(b.vx, -E_CUSH * b.vy)
                b.y = W - R - EPS_SEP
            if self.record:
                self.events.append((self.t, "cushion", b.name, None))
        elif kind == "pocket":
            b, k = pay
            b.active = False
            spd = b._s
            b.set_vel(0.0, 0.0)
            self.pocketed.append((self.t, b.name, POCKET_NAMES[k], spd))
            if self.record:
                self.events.append((self.t, "pocket", b.name, POCKET_NAMES[k]))
        return True

    def run(self, max_events=6000, max_time=300.0):
        n = 0
        while n < max_events and self.t < max_time:
            if not any(b.active and b._s >= V_MIN for b in self.balls):
                break
            if not self.step():
                break
            n += 1
        return n


# ---------------- configurations ----------------
def opening_rack(eps=5e-4, cue_pos=(0.65, W / 2.0)):
    """Standard opening configuration with reds separated by eps."""
    balls = [Ball(cue_pos[0], cue_pos[1], "cue")]
    spots = {"yellow": (0.737, W / 2 - 0.292), "green": (0.737, W / 2 + 0.292),
             "brown": (0.737, W / 2), "blue": (L / 2, W / 2),
             "pink": (2.667, W / 2), "black": (L - 0.324, W / 2)}
    for name, (x, y) in spots.items():
        balls.append(Ball(x, y, name))
    pitch = 2 * R + eps
    apex_x = 2.667 + 2 * R + eps
    k = 0
    for row in range(5):
        x = apex_x + row * pitch * math.sqrt(3) / 2
        for j in range(row + 1):
            y = W / 2 + (j - row / 2.0) * pitch
            k += 1
            balls.append(Ball(x, y, f"red{k}"))
    return balls, balls[0]


def shoot(cue, v0, phi):
    cue.set_vel(v0 * math.cos(phi), v0 * math.sin(phi))


if __name__ == "__main__":
    import time
    balls, cue = opening_rack()
    shoot(cue, 8.0, 0.0)
    s = Simulator(balls)
    t0 = time.time()
    n = s.run()
    print(f"events={n}  wall={time.time()-t0:.3f}s  sim_t={s.t:.1f}s  "
          f"pocketed={[(b, p) for _, b, p, _ in s.pocketed]}")
