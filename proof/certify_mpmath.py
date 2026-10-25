"""
witness_certify.py -- rigorous interval-arithmetic certificate for a witness.

Replays the witness of proof/witness_21.json in Model M itself (no numerical
nudges, no stopping threshold: the exact dynamics of Definition 2), with
every quantity enclosed in an interval (mpmath.iv, outward rounding).
The initial shot may be a whole box  V0 in [V0* - dV, V0* + dV],
phi in [phi* - dphi, phi* + dphi];  every shot in the box is covered.

At each step the script
  1. proposes the next event from a high-precision midpoint computation;
  2. proves, in interval arithmetic, that
       * the proposed clock has a unique, simple (transversal) zero in an
         explicit time interval T, with strictly negative derivative there,
         and is strictly positive on [0, inf T];
       * every other clock (all ball pairs, all cushions, all pockets,
         all stopping times) is strictly positive on [0, sup T];
  3. advances the enclosure to T and applies the event map.
Proposals are only heuristics; soundness rests on step 2 alone. A failed
proof aborts with an explanation, it never certifies.

The run stops at the 21st object-ball pocketing: from then on N = 21
regardless of what the cue ball does, and N <= 21 always.

Run:  CERT_PREC=1200 python proof/certify_mpmath.py [proof/witness_21.json] [dV] [dphi]

Coordinates and shot parameters are the binary64 values the JSON decodes
to, taken as exact rationals (e.g. V0* = 11.9000000000000003552713678800500929355621337890625).
The box half-widths are rounded outward.
Exit code 0 iff the certificate is complete.
"""
import json
import os
from os import path as _os_path

_HERE = _os_path.dirname(_os_path.abspath(__file__))
import sys
import time

from mpmath import iv, mp

PREC = int(os.environ.get("CERT_PREC", "256"))
iv.prec = PREC
mp.prec = PREC

# ---------------- Model M, decimal parameters enclosed exactly ----------
L = iv.mpf("3.569")
W = iv.mpf("1.778")
R = iv.mpf("0.02625")
A = iv.mpf("0.10")
E = iv.mpf("0.95")
EC = iv.mpf("0.90")
RHO = iv.mpf("0.055")
POCKETS = ((iv.mpf(0), iv.mpf(0)), (L, iv.mpf(0)), (iv.mpf(0), W),
           (L, W), (L / 2, iv.mpf(0)), (L / 2, W))
POCKET_NAMES = ("C-bl", "C-br", "C-tl", "C-tr", "M-b", "M-t")
FOUR_R2 = 4 * R * R
RHO2 = RHO * RHO
# walls: (axis, kind) with kind 0 -> coordinate = R, 1 -> coordinate = L-R / W-R
WALLS = ((0, 0), (0, 1), (1, 0), (1, 1))
WALL_NAMES = ("left", "right", "bottom", "top")

ZERO = iv.mpf(0)


class CertificateError(Exception):
    pass


def lo(x):
    return mp.mpf(x.a)


def hi(x):
    return mp.mpf(x.b)


def mid(x):
    return mp.mpf(x.mid)


def width(x):
    return hi(x) - lo(x)


def pos(x):
    return (x > 0) is True


def neg(x):
    return (x < 0) is True


# ---------------- polynomials with interval coefficients ----------------
def horner(c, t):
    """c = (c0, c1, ..., cn); evaluate sum c_k t^k."""
    acc = c[-1]
    for ck in reversed(c[:-1]):
        acc = acc * t + ck
    return acc


def deriv(c):
    return tuple(k * c[k] for k in range(1, len(c))) or (ZERO,)


def prove_positive(c, a, b, depth=60, budget=None):
    """Prove p(t) > 0 for all t in [a, b] and all coefficient values."""
    if budget is None:
        budget = [20000]
    stack = [(mp.mpf(a), mp.mpf(b), depth)]
    dc = deriv(c)
    while stack:
        a_, b_, d = stack.pop()
        budget[0] -= 1
        if budget[0] < 0:
            return False
        T = iv.mpf([a_, b_])
        if pos(horner(c, T)):
            continue
        m = (a_ + b_) / 2
        centred = horner(c, iv.mpf(m)) + horner(dc, T) * (T - m)
        if pos(centred):
            continue
        if not pos(horner(c, iv.mpf(a_))) or d == 0:
            return False
        stack.append((a_, m, d - 1))
        stack.append((m, b_, d - 1))
    return True


# ---------------- midpoint root proposals (heuristic only) --------------
def first_down_root(c, tmax):
    """Smallest t in (0, tmax] where the polynomial crosses zero downward,
    computed in mp at the coefficient midpoints. Heuristic."""
    cm = [mid(x) for x in c]
    while len(cm) > 1 and cm[-1] == 0:
        cm.pop()
    if len(cm) < 2:
        return None
    dm = [k * cm[k] for k in range(1, len(cm))]
    f = lambda t: sum(ck * t ** k for k, ck in enumerate(cm))
    df = lambda t: sum(ck * t ** k for k, ck in enumerate(dm))
    roots = []
    try:
        for r in mp.polyroots(list(reversed(cm)), maxsteps=200, extraprec=200):
            if abs(mp.im(r)) < mp.mpf(10) ** -30 * (1 + abs(r)):
                roots.append(mp.re(r))
    except mp.NoConvergence:
        # fall back to sampling + bisection
        n = 4000
        prev_t, prev_f = mp.mpf(0), f(mp.mpf(0))
        for k in range(1, n + 1):
            t = tmax * k / n
            ft = f(t)
            if prev_f > 0 and ft <= 0:
                roots.append(mp.findroot(f, (prev_t, t), solver="bisect"))
                break
            prev_t, prev_f = t, ft
    best = None
    for r in roots:
        if r <= mp.mpf(10) ** -40 or r > tmax:
            continue
        if df(r) >= 0:
            continue
        if best is None or r < best:
            best = r
    if best is None:
        return None
    for _ in range(6):                      # polish
        d = df(best)
        if d == 0:
            break
        best = best - f(best) / d
    return best


# ---------------- the state -----------------------------------------------
class State:
    def __init__(self, cue, v0, phi, balls):
        self.names = ["cue"] + [f"b{i}" for i in range(len(balls))]
        pts = [cue] + balls
        self.px = [iv.mpf(p[0]) for p in pts]
        self.py = [iv.mpf(p[1]) for p in pts]
        self.vx = [ZERO for _ in pts]
        self.vy = [ZERO for _ in pts]
        self.vx[0] = v0 * iv.cos(phi)
        self.vy[0] = v0 * iv.sin(phi)
        self.active = [True] * len(pts)
        self.moving = [False] * len(pts)
        self.moving[0] = True
        self.fresh = set()       # clocks known to vanish exactly at t = 0
        self.t = ZERO

    def n(self):
        return len(self.px)

    def speed(self, i):
        return iv.sqrt(self.vx[i] ** 2 + self.vy[i] ** 2)

    def quad(self, i):
        """Coefficient of t^2 in the position: -(a/2) v/|v|."""
        if not self.moving[i]:
            return ZERO, ZERO
        s = self.speed(i)
        return -A / 2 * self.vx[i] / s, -A / 2 * self.vy[i] / s

    # clock polynomials ----------------------------------------------------
    def clock_bb(self, i, j):
        qix, qiy = self.quad(i)
        qjx, qjy = self.quad(j)
        dpx, dpy = self.px[i] - self.px[j], self.py[i] - self.py[j]
        dvx, dvy = self.vx[i] - self.vx[j], self.vy[i] - self.vy[j]
        dqx, dqy = qix - qjx, qiy - qjy
        return (dpx * dpx + dpy * dpy - FOUR_R2,
                2 * (dpx * dvx + dpy * dvy),
                dvx * dvx + dvy * dvy + 2 * (dpx * dqx + dpy * dqy),
                2 * (dvx * dqx + dvy * dqy),
                dqx * dqx + dqy * dqy)

    def clock_wall(self, i, w):
        axis, kind = WALLS[w]
        p = self.px[i] if axis == 0 else self.py[i]
        v = self.vx[i] if axis == 0 else self.vy[i]
        q = self.quad(i)[axis]
        top = (L if axis == 0 else W) - R
        if kind == 0:
            return (p - R, v, q)
        return (top - p, -v, -q)

    def clock_pocket(self, i, k):
        qx, qy = self.quad(i)
        cx, cy = POCKETS[k]
        dpx, dpy = self.px[i] - cx, self.py[i] - cy
        vx, vy = self.vx[i], self.vy[i]
        return (dpx * dpx + dpy * dpy - RHO2,
                2 * (dpx * vx + dpy * vy),
                vx * vx + vy * vy + 2 * (dpx * qx + dpy * qy),
                2 * (vx * qx + vy * qy),
                qx * qx + qy * qy)

    def clocks(self):
        """All clocks that can vanish: (key, coefficients, fresh?)."""
        idx = [i for i in range(self.n()) if self.active[i]]
        out = []
        for a_ in range(len(idx)):
            i = idx[a_]
            for j in idx[a_ + 1:]:
                if not (self.moving[i] or self.moving[j]):
                    continue
                key = ("bb", i, j)
                out.append((key, self.clock_bb(i, j)))
            if self.moving[i]:
                for w in range(4):
                    out.append((("wall", i, w), self.clock_wall(i, w)))
                for k in range(6):
                    out.append((("pocket", i, k), self.clock_pocket(i, k)))
        res = []
        for key, c in out:
            if key in self.fresh:
                res.append((key, c[1:], True))   # clock = t * g(t)
            else:
                res.append((key, c, False))
        return res

    def stop_times(self):
        return {i: self.speed(i) / A for i in range(self.n())
                if self.active[i] and self.moving[i]}


# ---------------- quick rigorous exclusion --------------------------------
def far_apart(st, key, tb):
    """Cheap sufficient condition that a clock stays positive on [0, tb]:
    distances shrink at most at the sum of the speeds."""
    kind = key[0]
    if kind == "bb":
        _, i, j = key
        d = iv.sqrt((st.px[i] - st.px[j]) ** 2 + (st.py[i] - st.py[j]) ** 2)
        return pos(d - 2 * R - (st.speed(i) + st.speed(j)) * tb)
    if kind == "pocket":
        _, i, k = key
        cx, cy = POCKETS[k]
        d = iv.sqrt((st.px[i] - cx) ** 2 + (st.py[i] - cy) ** 2)
        return pos(d - RHO - st.speed(i) * tb)
    return False


# ---------------- one certified step --------------------------------------
def certified_step(st):
    stops = st.stop_times()
    if not stops:
        raise CertificateError("all balls at rest before the clearance")
    horizon = min(mid(s) for s in stops.values())

    clocks = st.clocks()
    # 1. proposal
    best_t, best = horizon, ("stop", min(stops, key=lambda i: mid(stops[i])))
    for key, c, fresh in clocks:
        if key[0] in ("bb", "pocket") and far_apart(st, key, iv.mpf(best_t)):
            continue
        r = first_down_root(c, best_t)
        if r is not None and r < best_t:
            best_t, best = r, key

    # 2. enclose the event time
    if best[0] == "stop":
        T = stops[best[1]]
        ev_c = None
    else:
        ev_c = next(c for k, c, _ in clocks if k == best)
        dc = deriv(ev_c)
        delta = mp.mpf(2) ** (-PREC + 40) * (1 + best_t)
        while True:
            a_, b_ = best_t - delta, best_t + delta
            if (a_ > 0 and pos(horner(ev_c, iv.mpf(a_)))
                    and neg(horner(ev_c, iv.mpf(b_)))
                    and neg(horner(dc, iv.mpf([a_, b_])))):
                break
            delta *= 2
            if delta > mp.mpf("1e-4") * (1 + best_t):
                raise CertificateError(f"cannot isolate a simple zero for "
                                       f"{best} near t={mp.nstr(best_t, 12)}")
        T = iv.mpf([a_, b_])
        # Positive on [0, a_]: f(a_) > 0 and f' < 0 on a monotone stretch
        # [a_ - eta, a_], then a direct proof on [0, a_ - eta].
        eta = a_
        while eta > delta:
            if neg(horner(dc, iv.mpf([max(a_ - eta, 0), a_]))):
                break
            eta /= 2
        else:
            raise CertificateError(f"{best}: no monotone stretch before event")
        if a_ - eta > 0 and not prove_positive(ev_c, 0, a_ - eta):
            raise CertificateError(f"{best}: clock not positive before event")

    tb = hi(T)
    # 3. every other clock strictly positive on [0, sup T]
    for key, c, fresh in clocks:
        if key == best:
            continue
        if far_apart(st, key, iv.mpf(tb)):
            continue
        if not prove_positive(c, 0, tb):
            raise CertificateError(
                f"cannot exclude {key} before {best} at t~{mp.nstr(best_t, 12)}")
    for i, s in stops.items():
        if best == ("stop", i):
            continue
        if not pos(s - T):
            raise CertificateError(f"stop of {st.names[i]} not separated")

    # 4. advance every moving ball to T
    for i in range(st.n()):
        if not (st.active[i] and st.moving[i]):
            continue
        s = st.speed(i)
        ux, uy = st.vx[i] / s, st.vy[i] / s
        if best == ("stop", i):
            st.px[i] = st.px[i] + ux * s * s / (2 * A)
            st.py[i] = st.py[i] + uy * s * s / (2 * A)
            st.vx[i] = st.vy[i] = ZERO
            st.moving[i] = False
            continue
        st.px[i] = st.px[i] + st.vx[i] * T - A / 2 * ux * T * T
        st.py[i] = st.py[i] + st.vy[i] * T - A / 2 * uy * T * T
        ns = s - A * T
        st.vx[i], st.vy[i] = ux * ns, uy * ns
    st.t = st.t + T
    st.fresh = set()

    # 5. event map
    info = {"kind": best[0], "t": st.t}
    if best[0] == "bb":
        _, i, j = best
        nx = (st.px[j] - st.px[i]) / (2 * R)   # exact unit normal at contact
        ny = (st.py[j] - st.py[i]) / (2 * R)
        vin = st.vx[i] * nx + st.vy[i] * ny
        vjn = st.vx[j] * nx + st.vy[j] * ny
        info["approach"] = vin - vjn            # normal approach speed
        vin2 = (1 - E) / 2 * vin + (1 + E) / 2 * vjn
        vjn2 = (1 + E) / 2 * vin + (1 - E) / 2 * vjn
        st.vx[i] += (vin2 - vin) * nx
        st.vy[i] += (vin2 - vin) * ny
        st.vx[j] += (vjn2 - vjn) * nx
        st.vy[j] += (vjn2 - vjn) * ny
        for b in (i, j):
            if not pos(st.speed(b)):
                raise CertificateError(f"{st.names[b]} may be at rest after "
                                       f"collision")
            st.moving[b] = True
        st.fresh.add(best)
        info["balls"] = (st.names[i], st.names[j])
    elif best[0] == "wall":
        _, i, w = best
        axis, kind = WALLS[w]
        val = R if kind == 0 else (L if axis == 0 else W) - R
        if axis == 0:
            info["approach"] = abs(st.vx[i])
            st.px[i] = val
            st.vx[i] = -EC * st.vx[i]
        else:
            info["approach"] = abs(st.vy[i])
            st.py[i] = val
            st.vy[i] = -EC * st.vy[i]
        st.fresh.add(best)
        info["balls"] = (st.names[i], WALL_NAMES[w])
    elif best[0] == "pocket":
        _, i, k = best
        info["entry_speed"] = st.speed(i)
        st.active[i] = False
        st.moving[i] = False
        info["balls"] = (st.names[i], POCKET_NAMES[k])
    else:
        info["balls"] = (st.names[best[1]],)
    info["T"] = T
    return info


def check_admissible(st):
    n = st.n()
    for i in range(n):
        for j in range(i + 1, n):
            if not pos((st.px[i] - st.px[j]) ** 2 + (st.py[i] - st.py[j]) ** 2
                       - FOUR_R2):
                raise CertificateError(f"balls {i},{j} overlap")
        for w in range(4):
            if not pos(st.clock_wall(i, w)[0]):
                raise CertificateError(f"ball {i} not inside the table")
        for k in range(6):
            if not pos(st.clock_pocket(i, k)[0]):
                raise CertificateError(f"ball {i} inside a pocket disk")


def certify(path=_os_path.join(_HERE, "witness_21.json"), dV="0", dphi="0", verbose=True):
    w = json.load(open(path))
    # half-widths rounded outward, so the box contains the stated one
    rV, rphi = hi(iv.mpf(dV)), hi(iv.mpf(dphi))
    V0 = iv.mpf(w["v0"]) + iv.mpf([-rV, rV])
    PHI = iv.mpf(w["phi"]) + iv.mpf([-rphi, rphi])
    st = State(w["cue"], V0, PHI, w["balls"])
    check_admissible(st)
    n_obj = st.n() - 1
    need = int(os.environ.get("CERT_NEED", n_obj))  # tests only
    log = []
    t0 = time.time()
    potted = 0
    while potted < need:
        info = certified_step(st)
        log.append(info)
        if info["kind"] == "pocket":
            if info["balls"][0] == "cue":
                raise CertificateError("cue ball pocketed")
            potted += 1
        if verbose:
            extra = ""
            if "approach" in info:
                extra = f"  normal speed >= {mp.nstr(lo(info['approach']), 4)}"
            if "entry_speed" in info:
                extra = f"  entry speed >= {mp.nstr(lo(info['entry_speed']), 4)}"
            print(f"{len(log):3d}  t in [{mp.nstr(lo(info['t']), 10)}, "
                  f"+{mp.nstr(width(info['t']), 2)}]  "
                  f"{info['kind']:6s} {' / '.join(info['balls']):14s}{extra}")
        if len(log) > 2000:
            raise CertificateError("too many events")
    return st, log, time.time() - t0


def summary(st, log):
    approach = [lo(e["approach"]) for e in log if e["kind"] == "bb"]
    cush = [lo(e["approach"]) for e in log if e["kind"] == "wall"]
    entry = [lo(e["entry_speed"]) for e in log if e["kind"] == "pocket"]
    worst = max(max(width(st.px[i]), width(st.py[i]), width(st.vx[i]),
                    width(st.vy[i])) for i in range(st.n()))
    return {
        "events": len(log),
        "collisions": len(approach),
        "cushions": len(cush),
        "pocketings": len(entry),
        "min_collision_normal_speed": mp.nstr(min(approach), 6),
        "min_cushion_normal_speed": mp.nstr(min(cush), 6) if cush else None,
        "min_pocket_entry_speed": mp.nstr(min(entry), 6),
        "final_enclosure_width": mp.nstr(worst, 3),
        "final_time": mp.nstr(mid(st.t), 10),
    }


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else _os_path.join(_HERE, "witness_21.json")
    dV = sys.argv[2] if len(sys.argv) > 2 else "0"
    dphi = sys.argv[3] if len(sys.argv) > 3 else "0"
    print(f"certifying {src}: V0 +/- {dV}, phi +/- {dphi}, "
          f"{PREC}-bit interval arithmetic")
    try:
        st, log, secs = certify(src, dV, dphi)
    except CertificateError as exc:
        print(f"\nNOT CERTIFIED: {exc}")
        sys.exit(1)
    print()
    for k, v in summary(st, log).items():
        print(f"{k:28s}: {v}")
    print(f"wall time                   : {secs:.1f} s")
    print(f"\nCERTIFIED: every shot in the box pockets all "
          f"{st.n() - 1} object balls.")
    sys.exit(0)
