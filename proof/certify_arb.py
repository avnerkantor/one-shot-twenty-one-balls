#!/usr/bin/env python3
"""
witness_certify_arb.py -- independent rigorous certificate checker (Arb ball arithmetic).

Claim: for every shot (V0, phi) in [v0-r, v0+r] x [phi-r, phi+r] the initial configuration of
proof/witness_21.json is admissible and model M pockets all 21 object balls (the cue ball not
pocketed before the 21st object ball).

Method (event-driven, everything in python-flint arb balls; nothing relies on floating point):
  * The state (positions, velocities of every ball) is a vector of balls that contains the true
    state for every parameter in the box, at the (parameter dependent) time of the last event.
  * Between events every moving ball moves on a straight line p + d*g(t), g(t) = s t - a t^2/2,
    s = |v|, d = v/s, until its stop time s/a.
  * Event functions:
      stop                    : t = s/a                                   (closed form)
      wall, pocket, bb with one ball at rest : polynomial in the arc length g of the moving ball,
                                degree 1 or 2; first root by closed form, convexity/sign argument
      bb with both balls moving: quartic in t; first root by rigorous left-first bisection
                                (positivity proofs) + sign change + f'<0 + interval Newton.
  * Each step: all closed-form candidate enclosures are computed; H = min of their upper ends is a
    time by which some event provably happens. Every remaining event function is shown either
    positive on [0,H] or its first root is isolated. The earliest event E must have an enclosure
    whose upper end is < the lower end of every other candidate. Then the state is advanced to
    the enclosure of E's time and E is applied.
  * The pair (or ball/cushion) that just made contact has an exactly-zero clock; it is handled by
    factoring the root at 0 (pair) or by the strict sign of the normal velocity (cushion).
"""
import sys, json, time, argparse
from fractions import Fraction
from flint import arb, ctx

# ----------------------------------------------------------------------------- utilities
class NotCertified(Exception):
    pass

def fail(msg):
    raise NotCertified(msg)

def hull(lo, hi):
    """Ball containing [lower(lo), upper(hi)]."""
    return arb(lo.lower()).union(arb(hi.upper()))

def lowb(x):
    return x.lower()

def upb(x):
    return x.upper()

def amin(a, b):
    return a if (a.upper() <= b.upper()) else b

def fmt_lower(x, sig=10):
    """Decimal string with `sig` significant digits that is <= every point of x (rounded down)."""
    lo = arb(x.lower())
    if not lo.is_finite():
        return "-inf"
    if not lo.is_exact():
        raise RuntimeError("lower bound not exact")
    man, ex = lo.man_exp()
    q = Fraction(int(man)) * (Fraction(2) ** int(ex))
    if q == 0:
        return "0"
    neg = q < 0
    aq = -q if neg else q
    e = 0
    while aq >= 10 ** (e + 1):
        e += 1
    while aq < 10 ** e:
        e -= 1
    scale = Fraction(10) ** (sig - 1 - e)
    v = aq * scale
    if neg:
        n = -(-(-v.numerator // v.denominator))   # -ceil(|v|): rounds toward -inf
    else:
        n = v.numerator // v.denominator          # floor(|v|)
    s = str(abs(n)).rjust(sig, "0")
    expo = e
    out = s[0] + "." + s[1:] + ("e%+d" % expo if expo != 0 else "")
    return ("-" if n < 0 else "") + out

def fmt_rad(x):
    """Upper bound of the radius of x, 4 digits (display only)."""
    r = x.rad()
    if r.is_zero():
        return "0"
    return arb(r.upper()).str(4, radius=False)

def fmt_up(x):
    """Upper bound of a nonnegative quantity, 4 digits (display only)."""
    if x.is_zero():
        return "0"
    return arb(x.upper()).str(4, radius=False)

def peval(c, t):
    r = c[-1]
    for k in range(len(c) - 2, -1, -1):
        r = r * t + c[k]
    return r

def pderiv(c):
    return [c[k] * k for k in range(1, len(c))] or [arb(0)]

def prove_pos(c, dc, lo, hi):
    """True if c(t) > 0 provably for all t in [lo, hi]."""
    I = hull(lo, hi)
    if peval(c, I) > 0:
        return True
    m = (lo + hi) / 2
    mv = peval(c, m) + peval(dc, I) * (I - m)
    return mv > 0

def first_root(c, X, maxdepth=200, maxnodes=200000):
    """Polynomial c (arb coefficients) on [0, X] (X exact, X>0).
    Returns ('none',) if c > 0 on [0,X] for every realisation,
            ('root', J, dJ) if every realisation is > 0 on [0, lo J), has a unique root in J and
                             c' < 0 on J (dJ = enclosure of c' on J),
            ('fail', msg) otherwise."""
    dc = pderiv(c)
    stack = [(arb(0), X, 0)]
    nodes = 0
    while stack:
        lo, hi, dep = stack.pop()
        nodes += 1
        if nodes > maxnodes:
            return ('fail', 'too many subdivisions')
        if prove_pos(c, dc, lo, hi):
            continue
        I = hull(lo, hi)
        dI = peval(dc, I)
        if dI < 0 and peval(c, lo) > 0 and peval(c, hi) < 0:
            J = newton_refine(c, dc, lo, hi)
            return ('root', J, peval(dc, J))
        if peval(c, lo) < 0:
            return ('fail', 'polynomial negative at left end of an unexcluded interval')
        if dep >= maxdepth:
            return ('fail', 'bisection depth exceeded')
        m = (lo + hi) / 2
        stack.append((m, hi, dep + 1))
        stack.append((lo, m, dep + 1))
    return ('none',)

def newton_refine(c, dc, lo, hi, maxit=200):
    """Precondition: c(lo)>0, c(hi)<0, c'<0 on [lo,hi] (all realisations). Interval Newton
    contraction; every root in [lo,hi] stays in the returned enclosure."""
    J = hull(lo, hi)
    for _ in range(maxit):
        m = J.mid()
        N = m - peval(c, m) / peval(dc, J)
        if not N.is_finite():
            break
        K = J.intersection(N)
        if K.rad() >= J.rad() * arb("0.75"):
            J = K
            break
        J = K
    return J

# ----------------------------------------------------------------------------- model constants
def constants():
    C = {}
    C['L'] = arb("3.569"); C['W'] = arb("1.778"); C['R'] = arb("0.02625")
    C['rho'] = arb("0.055"); C['a'] = arb("0.10"); C['e'] = arb("0.95"); C['ec'] = arb("0.90")
    L, W = C['L'], C['W']
    C['pockets'] = [("C-bl", arb(0), arb(0)), ("C-br", L, arb(0)), ("C-tl", arb(0), W),
                    ("C-tr", L, W), ("M-b", L / 2, arb(0)), ("M-t", L / 2, W)]
    return C

WALLS = ["left", "right", "bottom", "top"]

def bname(i):
    return "cue" if i == 0 else "b%d" % (i - 1)

# ----------------------------------------------------------------------------- the checker
class Ball:
    __slots__ = ("x", "y", "vx", "vy", "moving", "alive")

def wall_data(C, b, w, dx, dy):
    """(distance to wall line >=0 expected, component of d toward the wall)."""
    R, L, W = C['R'], C['L'], C['W']
    if w == 0: return b.x - R, -dx
    if w == 1: return (L - R) - b.x, dx
    if w == 2: return b.y - R, -dy
    return (W - R) - b.y, dy

def g_to_t(g, s, a):
    disc = s * s - 2 * a * g
    if not disc > 0:
        fail("arc length not provably before stop")
    return 2 * g / (s + disc.sqrt())

def run(V, PHI, data, C, out, label=""):
    R, a, e, ec, rho = C['R'], C['a'], C['e'], C['ec'], C['rho']
    L, W = C['L'], C['W']
    R2x = 2 * R
    D2 = R2x * R2x
    balls = []
    pts = [data['cue']] + data['balls']
    for k, (px, py) in enumerate(pts):
        b = Ball(); b.x = arb(float(px)); b.y = arb(float(py)); b.alive = True
        if k == 0:
            b.vx = V * PHI.cos(); b.vy = V * PHI.sin(); b.moving = True
        else:
            b.vx = arb(0); b.vy = arb(0); b.moving = False
        balls.append(b)
    n = len(balls)
    # ---------------- admissibility
    for i in range(n):
        b = balls[i]
        if not (b.x > R and b.x < L - R and b.y > R and b.y < W - R):
            fail("initial centre of %s not strictly inside the table" % bname(i))
        for (pn, PX, PY) in C['pockets']:
            if not ((b.x - PX) ** 2 + (b.y - PY) ** 2 > rho * rho):
                fail("initial centre of %s not outside pocket %s" % (bname(i), pn))
        for j in range(i + 1, n):
            c = balls[j]
            if not ((b.x - c.x) ** 2 + (b.y - c.y) ** 2 > D2):
                fail("initial distance %s-%s not > 2R" % (bname(i), bname(j)))
    if not (V > 0):
        fail("V0 not provably positive")

    T = arb(0)
    contact = None
    pocketed = 0
    ev = 0
    stats = {'bb': 0, 'wall': 0, 'pocket': 0, 'stop': 0}
    minsp = {'bb': None, 'wall': None, 'pocket': None, 'pocketspeed': None}
    maxtw = arb(0)
    def upd_min(k, v):
        if minsp[k] is None or v.lower() < minsp[k].lower():
            minsp[k] = v

    while True:
        mv = [i for i in range(n) if balls[i].alive and balls[i].moving]
        if not mv:
            fail("all balls at rest with %d object balls pocketed" % pocketed)
        S, DX, DY, GM, TS = {}, {}, {}, {}, {}
        for i in mv:
            b = balls[i]
            s = (b.vx * b.vx + b.vy * b.vy).sqrt()
            if not s > 0:
                fail("speed of moving ball %s not provably positive" % bname(i))
            S[i] = s; DX[i] = b.vx / s; DY[i] = b.vy / s
            GM[i] = s * s / (2 * a); TS[i] = s / a
        cands = []        # (kind, ids, tenc)
        generic_g = []    # (kind, ids, i, coeffs-in-g)
        # ---- stops
        for i in mv:
            cands.append(('stop', (i,), TS[i]))
        # ---- cushions
        for i in mv:
            b = balls[i]
            for w in range(4):
                dist, tow = wall_data(C, b, w, DX[i], DY[i])
                if contact == ('wall', i, w):
                    if tow < 0:
                        continue
                    fail("ball %s not provably leaving cushion %s" % (bname(i), WALLS[w]))
                if tow < 0:
                    continue
                if not dist > 0:
                    fail("ball %s not provably inside cushion %s" % (bname(i), WALLS[w]))
                if tow > 0:
                    gs = dist / tow
                    if gs > GM[i]:
                        continue
                    if gs < GM[i]:
                        cands.append(('wall', (i, w), g_to_t(gs, S[i], a)))
                        continue
                generic_g.append(('wall', (i, w), i, [dist, -tow]))
        # ---- pockets and bb with one ball at rest: monic quadratic in g, g^2 + 2 b g + c
        def quad(kind, ids, i, qx, qy, rr, zero):
            c = qx * qx + qy * qy - rr * rr
            bb = DX[i] * qx + DY[i] * qy
            if zero:
                # c is exactly 0 for the true state: g (g + 2b)
                if bb > 0:
                    return
                generic_g.append((kind, ids, i, [2 * bb, arb(1)]))   # h(g) = g + 2b
                return
            if not c > 0:
                fail("%s %s: initial clearance not provably positive" % (kind, ids))
            if bb > 0:
                return                       # both roots negative
            disc = bb * bb - c
            if disc < 0:
                return                       # no real root
            if bb < 0 and disc > 0:
                gs = c / (-bb + disc.sqrt())  # smaller (positive) root
                if gs > GM[i]:
                    return
                if gs < GM[i]:
                    cands.append((kind, ids, g_to_t(gs, S[i], a)))
                    return
            generic_g.append((kind, ids, i, [c, 2 * bb, arb(1)]))
        for i in mv:
            b = balls[i]
            for k, (pn, PX, PY) in enumerate(C['pockets']):
                quad('pocket', (i, k), i, b.x - PX, b.y - PY, rho, False)
        quartic = []
        for i in range(n):
            if not balls[i].alive: continue
            for j in range(i + 1, n):
                if not balls[j].alive: continue
                mi, mj = balls[i].moving, balls[j].moving
                if not (mi or mj): continue
                zero = (contact == ('bb', i, j))
                if mi and mj:
                    quartic.append((i, j, zero))
                else:
                    k, l = (i, j) if mi else (j, i)
                    quad('bb', (i, j), k, balls[k].x - balls[l].x, balls[k].y - balls[l].y, R2x, zero)
        if not cands:
            fail("no closed-form candidate")
        H = cands[0][2].upper()
        for cd in cands[1:]:
            u = cd[2].upper()
            if u < H: H = u
        H = arb(H)
        # ---- generic single-ball event functions (in arc length g)
        for (kind, ids, i, cf) in generic_g:
            if H < TS[i].lower():
                gH = S[i] * H - a * H * H / 2
            else:
                gH = GM[i]
            G = arb(gH.upper())
            rr = first_root(cf, G)
            if rr[0] == 'none':
                continue
            if rr[0] == 'fail':
                fail("%s %s: %s" % (kind, ids, rr[1]))
            cands.append((kind, ids, g_to_t(rr[1], S[i], a)))
        # ---- both-moving pairs (quartic in t)
        for (i, j, zero) in quartic:
            bi, bj = balls[i], balls[j]
            Ax = bi.x - bj.x; Ay = bi.y - bj.y
            Bx = bi.vx - bj.vx; By = bi.vy - bj.vy
            if not zero:
                dist = (Ax * Ax + Ay * Ay).sqrt()
                if dist - (S[i] + S[j]) * H > R2x:
                    continue
            Cx = -(a / 2) * (DX[i] - DX[j]); Cy = -(a / 2) * (DY[i] - DY[j])
            c0 = Ax * Ax + Ay * Ay - D2
            c1 = 2 * (Ax * Bx + Ay * By)
            c2 = Bx * Bx + By * By + 2 * (Ax * Cx + Ay * Cy)
            c3 = 2 * (Bx * Cx + By * Cy)
            c4 = Cx * Cx + Cy * Cy
            if zero:
                if not c1 > 0:
                    fail("pair %s-%s not provably separating after contact" % (bname(i), bname(j)))
                cf = [c1, c2, c3, c4]          # f(t) = t * h(t), true c0 = 0
            else:
                if not c0 > 0:
                    fail("pair %s-%s clearance not provably positive" % (bname(i), bname(j)))
                cf = [c0, c1, c2, c3, c4]
            H1 = H; clipped = False
            for tl in (TS[i].lower(), TS[j].lower()):
                if tl < H1: H1 = arb(tl); clipped = True
            if not H1 > 0:
                fail("pair horizon degenerate")
            rr = first_root(cf, H1)
            if rr[0] == 'fail':
                fail("pair %s-%s: %s" % (bname(i), bname(j), rr[1]))
            if rr[0] == 'root':
                J = rr[1]
                if zero and not J.lower() > 0:
                    fail("re-contact root not separated from 0")
                cands.append(('bb', (i, j), J))
                continue
            if clipped:
                # window [H1, H]: positions enclosed as boxes (accounts for possible stops)
                def box(k):
                    bk = balls[k]
                    g1 = S[k] * H1 - a * H1 * H1 / 2
                    if H < TS[k].lower():
                        g2 = S[k] * H - a * H * H / 2
                    else:
                        g2 = GM[k]
                    g = hull(g1, g2)
                    return bk.x + DX[k] * g, bk.y + DY[k] * g
                xi, yi = box(i); xj, yj = box(j)
                if not ((xi - xj) ** 2 + (yi - yj) ** 2 > D2):
                    fail("pair %s-%s: window before a stop not cleared" % (bname(i), bname(j)))
        # ---- choose the earliest event and separate it from all others
        best = min(range(len(cands)), key=lambda q: float(cands[q][2].mid()))
        E = cands[best]
        tE = E[2]
        if not tE.lower() > 0:
            fail("event time not provably positive")
        for q, cd in enumerate(cands):
            if q == best: continue
            if not cd[2].lower() > tE.upper():
                fail("cannot separate event %s %s from %s %s" % (E[0], E[1], cd[0], cd[1]))
        if tE.rad() > maxtw: maxtw = arb(tE.rad())
        # ---- advance every moving ball to time tE
        kind, ids = E[0], E[1]
        for i in mv:
            b = balls[i]
            if kind == 'stop' and ids[0] == i:
                b.x = b.x + DX[i] * GM[i]; b.y = b.y + DY[i] * GM[i]
                b.vx = arb(0); b.vy = arb(0); b.moving = False
                continue
            g = S[i] * tE - a * tE * tE / 2
            sp = S[i] - a * tE
            b.x = b.x + DX[i] * g; b.y = b.y + DY[i] * g
            b.vx = DX[i] * sp; b.vy = DY[i] * sp
        T = T + tE
        ev += 1
        stats[kind] += 1
        line = "%3d  %-6s" % (ev, kind)
        if kind == 'stop':
            line += " %-12s t>=%s" % (bname(ids[0]), fmt_lower(T))
            contact = None
        elif kind == 'wall':
            i, w = ids
            b = balls[i]
            if w == 0: b.x = R; vn = -b.vx; b.vx = -ec * b.vx
            elif w == 1: b.x = L - R; vn = b.vx; b.vx = -ec * b.vx
            elif w == 2: b.y = R; vn = -b.vy; b.vy = -ec * b.vy
            else: b.y = W - R; vn = b.vy; b.vy = -ec * b.vy
            if not vn > 0:
                fail("cushion approach speed not positive")
            upd_min('wall', vn)
            line += " %-12s t>=%s  vn>=%s" % (bname(i) + " " + WALLS[w], fmt_lower(T), fmt_lower(vn))
            contact = ('wall', i, w)
        elif kind == 'pocket':
            i, k = ids
            b = balls[i]
            pn, PX, PY = C['pockets'][k]
            vr = -((b.x - PX) * b.vx + (b.y - PY) * b.vy) / rho
            sp = (b.vx * b.vx + b.vy * b.vy).sqrt()
            if not vr > 0:
                fail("pocket entry speed not positive")
            upd_min('pocket', vr); upd_min('pocketspeed', sp)
            line += " %-12s t>=%s  vr>=%s  |v|>=%s" % (bname(i) + " " + pn, fmt_lower(T),
                                                       fmt_lower(vr), fmt_lower(sp))
            b.alive = False; b.moving = False
            contact = None
            if i == 0:
                out(line)
                fail("cue ball pocketed")
            pocketed += 1
        else:  # bb
            i, j = ids
            bi, bj = balls[i], balls[j]
            nx = (bj.x - bi.x) / R2x; ny = (bj.y - bi.y) / R2x
            ui = bi.vx * nx + bi.vy * ny; uj = bj.vx * nx + bj.vy * ny
            appr = ui - uj
            if not appr > 0:
                fail("bb approach speed not positive")
            upd_min('bb', appr)
            dv = (1 + e) / 2 * (uj - ui)
            bi.vx = bi.vx + dv * nx; bi.vy = bi.vy + dv * ny
            bj.vx = bj.vx - dv * nx; bj.vy = bj.vy - dv * ny
            bi.moving = True; bj.moving = True
            line += " %-12s t>=%s  vn>=%s" % (bname(i) + " " + bname(j), fmt_lower(T), fmt_lower(appr))
            contact = ('bb', i, j)
        out(line)
        if pocketed == 21:
            break
    # final enclosure width
    wmax = arb(0)
    for b in balls:
        if b.alive:
            for z in (b.x, b.y, b.vx, b.vy):
                if z.rad() > wmax: wmax = arb(z.rad())
    return dict(ev=ev, stats=stats, minsp=minsp, T=T, wmax=wmax, maxtw=maxtw)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r", default="0", help="half-width of the (V0, phi) box (decimal string)")
    ap.add_argument("--prec", type=int, default=1200)
    ap.add_argument("--json", default=__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "witness_21.json"))
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    ctx.prec = args.prec
    lines = []
    fh = open(args.out, "w") if args.out else None
    def out(s):
        print(s, flush=True)
        if fh: fh.write(s + "\n"); fh.flush()
    t0 = time.time()
    status = 1
    try:
        with open(args.json) as f:
            data = json.load(f)
        v0 = data['v0']; phi = data['phi']
        if not (isinstance(v0, float) and isinstance(phi, float)):
            v0 = float(v0); phi = float(phi)
        rq = Fraction(args.r)
        if rq < 0:
            fail("negative r")
        r_up = (arb(rq.numerator) / arb(rq.denominator)).upper() if rq != 0 else arb(0)
        V = arb(v0, r_up) if rq != 0 else arb(v0)
        PHI = arb(phi, r_up) if rq != 0 else arb(phi)
        out("witness_certify_arb.py  model M, Arb ball arithmetic, prec=%d bits" % args.prec)
        out("data: %s   v0=%r  phi=%r  r=%s (radius used: %s)" % (args.json, v0, phi, args.r,
                                                                   fmt_rad(V) if rq != 0 else "0"))
        out("events: idx kind names  t>=elapsed-time-lower-bound  speeds (lower bounds)")
        C = constants()
        res = run(V, PHI, data, C, out)
        st, ms = res['stats'], res['minsp']
        out("")
        out("summary: %d events (bb %d, wall %d, pocket %d, stop %d)" %
            (res['ev'], st['bb'], st['wall'], st['pocket'], st['stop']))
        out("  min bb normal approach speed      >= %s" % (fmt_lower(ms['bb']) if ms['bb'] else "n/a"))
        out("  min cushion normal approach speed >= %s" % (fmt_lower(ms['wall']) if ms['wall'] else "n/a"))
        out("  min pocket radial entry speed     >= %s" % fmt_lower(ms['pocket']))
        out("  min pocket entry speed |v|        >= %s" % fmt_lower(ms['pocketspeed']))
        out("  final time enclosure  [%s, +%s]" % (fmt_lower(res['T']), fmt_rad(res['T'])))
        out("  final state enclosure max radius   %s ; max event-time radius %s" %
            (fmt_up(res['wmax']), fmt_up(res['maxtw'])))
        out("  wall time %.1f s" % (time.time() - t0))
        out("CERTIFIED")
        status = 0
    except NotCertified as ex:
        out("NOT CERTIFIED: %s" % ex)
        out("  wall time %.1f s" % (time.time() - t0))
    if fh: fh.close()
    sys.exit(status)

if __name__ == "__main__":
    main()
