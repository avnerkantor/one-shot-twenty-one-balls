"""
certify_window.py -- certify the whole window of directions at once.

Claim: every stroke with direction phi in [phi* - DMINUS, phi* + DPLUS]
and speed V0 in [V0* - DV, V0* + DV] realises the witness's event list
(the same 61 events, each ball's events in the same order) and pockets
all 21 object balls.

Why a plain box of directions cannot do this. Through twenty thin cuts
the true dynamics amplifies a change of direction by about 3e40 m/rad,
but a box of directions evaluated in interval arithmetic loses about
1e40 more (the dependency problem), which is why certify_mpmath.py and
certify_strip.py stop at boxes of width about 1e-85. The window is
about 1e-42 wide: 1e47 boxes.

The remedy is a first-order centered form in the one variable phi.
Every quantity q along the geometric construction of the skeleton is
represented over a piece [phic - h, phic + h] as

    q(phi) in c + D * (phi - phic),        |phi - phic| <= h,

with c a thin interval (the value at the centre, with rounding) and D
an interval enclosing the slope. Sums, differences and products
propagate exactly to first order; square roots and divisions use the
mean value theorem on the range; cos and sin of phi itself use their
derivatives over the piece. The box width enters only through the
second-order terms, so the enclosures stay as sharp as the dynamics
itself across pieces far wider than 1e-85.

The skeleton (geometry) is built in this arithmetic, exactly as in
certify_strip.py, including the checks that each designated event is
the first along its ray. The timing checks (pieces of different balls
never within 2R at the same instant; the pair that has just collided
separates) reuse certify_strip.py on the ranges of the centred forms.

An adaptive bisection covers the window: a piece that cannot be decided
is split, down to a minimum width. The output lists every certified
piece; the pieces are contiguous and cover the window exactly.

Run:  CERT_PREC=400 python certify_window.py [DMINUS DPLUS [DV]]
      CERT_PREC=400 python certify_window.py DMINUS DPLUS VA VB
      (V0 in V0* +/- DV, or V0 in V0* + [VA, VB]; defaults
      4.358e-44 9.045e-43 0.01)
"""
import json
import os
import sys
import time

from mpmath import iv, mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import certify_strip as S  # noqa: E402

iv.prec = S.PREC
mp.prec = S.PREC
R, RHO, A, E, EC, Z = S.R, S.RHO, S.A, S.E, S.EC, S.Z
Fail, pos, lo, hi = S.Fail, S.pos, S.lo, S.hi


# ------------------------------------------------------- centred forms
class CF:
    """q(phi) in c + D*(phi - phic) on the current piece, T = [-h, h]."""
    __slots__ = ("c", "D")
    T = iv.mpf(0)

    def __init__(self, c, D=None):
        self.c = c if isinstance(c, type(Z)) else iv.mpf(c)
        self.D = Z if D is None else D

    @staticmethod
    def lift(x):
        return x if isinstance(x, CF) else CF(x)

    def range(self):
        return self.c + self.D * CF.T

    def __neg__(self):
        return CF(-self.c, -self.D)

    def __add__(self, o):
        o = CF.lift(o)
        return CF(self.c + o.c, self.D + o.D)

    __radd__ = __add__

    def __sub__(self, o):
        o = CF.lift(o)
        return CF(self.c - o.c, self.D - o.D)

    def __rsub__(self, o):
        return CF.lift(o) - self

    def __mul__(self, o):
        o = CF.lift(o)
        # (c1 + D1 t)(c2 + D2 t) = c1 c2 + (c1 D2 + c2 D1 + D1 D2 t) t
        return CF(self.c * o.c, self.c * o.D + o.c * self.D + self.D * o.D * CF.T)

    __rmul__ = __mul__

    def inv(self):
        rng = self.range()
        if not (pos(rng) or pos(-rng)):
            raise Fail("division by a form whose range contains zero")
        # 1/(c + D t) = 1/c - D t / (c (c + D t))
        return CF(1 / self.c, -self.D / (self.c * rng))

    def __truediv__(self, o):
        return self * CF.lift(o).inv()

    def __rtruediv__(self, o):
        return CF.lift(o) * self.inv()

    def sqrt(self):
        rng = self.range()
        if not pos(rng):
            raise Fail("square root of a form whose range is not positive")
        # sqrt(c + D t) = sqrt(c) + D t / (sqrt(c + D t) + sqrt(c))
        sc = iv.sqrt(self.c)
        return CF(sc, self.D / (iv.sqrt(rng) + sc))

    @classmethod
    def set_piece(cls, h):
        cls.T = iv.mpf([-h, h])

    @classmethod
    def cos_sin(cls, phic):
        """Forms of cos(phic + t) and sin(phic + t) over the piece."""
        phic_iv = iv.mpf(phic)
        PHI = phic_iv + cls.T
        return CF(iv.cos(phic_iv), -iv.sin(PHI)), CF(iv.sin(phic_iv), iv.cos(PHI))


# With CERT_TM_DEGREE=N (N >= 2) the geometry is carried in Taylor models
# of degree N instead of first-order centred forms (taylor_models.py);
# the two classes share the interface used below. Needed when the slopes
# computed in intervals inflate along a chain of near-frozen collisions.
TM_DEGREE = int(os.environ.get("CERT_TM_DEGREE", "0"))
if TM_DEGREE:
    import taylor_models as _TM  # noqa: E402
    _TM.Fail = Fail
    _TM.TM.N = TM_DEGREE
    CF = _TM.TM  # noqa: F811


def posr(x):
    return pos(x.range())


def negr(x):
    return pos(-x.range())


def rng(x):
    return x.range() if isinstance(x, CF) else x


# ----------------------------------------------------------- geometry
class Piece:
    """As in certify_strip.Piece, with geometry in centred forms. In
    addition al, be, be1 (centred forms in phi) record the speed scaling
    along the skeleton: on this piece the squared speed at arc length s is
    al * V0^2 - (be + 2 a s), so be1 = be + 2 a ell is the value at the end
    of the piece (Lemma on speed scaling)."""
    __slots__ = ("ball", "kind", "P", "u", "ell", "t0", "t1", "v0", "v1",
                 "e0", "e1", "al", "be", "be1")

    def __init__(self, ball, kind, P, t0, e0):
        self.ball, self.kind, self.P, self.t0, self.e0 = ball, kind, P, t0, e0
        self.u = self.ell = self.t1 = self.v0 = self.v1 = self.e1 = None
        self.al = self.be = self.be1 = None


def first_contact(P, u, Q, rad):
    """Smallest s > 0 with |P + s u - Q| = rad as a centred form, proved a
    simple first root; None if provably no root ahead."""
    d = (P[0] - Q[0], P[1] - Q[1])
    b = d[0] * u[0] + d[1] * u[1]
    c = d[0] * d[0] + d[1] * d[1] - rad * rad
    if not posr(c):
        raise Fail("start inside a contact disk")
    if posr(b):
        return None
    disc = b * b - c
    if negr(disc):
        return None
    if not posr(disc):
        raise Fail("tangential contact cannot be excluded")
    if not posr(-b):
        raise Fail("approach direction undecided")
    return -b - disc.sqrt()


def wall_distance(P, u, wall):
    name, axis, val = wall
    p, du = P[axis], u[axis]
    if name in ("left", "bottom"):
        if posr(-du):
            return (p - val) / (-du)
    else:
        if posr(du):
            return (CF(val) - p) / du
    if posr(du) or negr(du):
        return None
    raise Fail("direction parallel to a wall undecided")


def end_piece(pc, s, event):
    """s: centred form (or interval) arc length of the piece."""
    sr = rng(s)
    pc.ell = s
    if pc.al is not None:
        pc.be1 = pc.be + s * (2 * A)
    pc.v1 = S.travel(pc.v0, sr)
    pc.t1 = pc.t0 + 2 * sr / (pc.v0 + pc.v1)
    pc.e1 = event
    return (pc.P[0] + s * pc.u[0], pc.P[1] + s * pc.u[1])


def check_single(pc, s, ignore_wall=None, ignore_pocket=None):
    """No cushion or pocket is met on the piece before arc length s."""
    sr = rng(s)
    for wall in S.WALLS:
        if wall[0] == ignore_wall:
            continue
        d = wall_distance(pc.P, pc.u, wall)
        if d is not None and not pos(d.range() - sr):
            raise Fail(f"{pc.ball}: cushion {wall[0]} may come first")
    for k, Q in enumerate(S.POCKETS):
        if k == ignore_pocket:
            continue
        try:
            d = first_contact(pc.P, pc.u, Q, RHO)
        except Fail as exc:
            if "tangential" in str(exc) and passes_beyond(pc.P, pc.u, Q, RHO, sr):
                continue                    # the near-tangent lies beyond the end of the piece
            raise Fail(f"{pc.ball}: pocket {S.POCKET_NAMES[k]}: {exc}")
        if d is not None and not pos(d.range() - sr):
            raise Fail(f"{pc.ball}: pocket {S.POCKET_NAMES[k]} may come first")


def passes_beyond(P, u, Q, rad, sr):
    """True if any entry of the ray into the disk of radius rad about Q
    provably happens beyond the arc length range sr (centred-form inputs)."""
    d = (P[0] - Q[0], P[1] - Q[1])
    b = (d[0] * u[0] + d[1] * u[1]).range()
    disc = b * b - (d[0] * d[0] + d[1] * d[1]).range() + rad * rad
    root = iv.sqrt(iv.mpf([0, max(S.hi(disc), 0)]))
    return pos(-b - root - sr)


def reflect(pc, wall, s, event):
    C = list(end_piece(pc, s, event))
    C[wall[1]] = CF(wall[2])                     # exactly on the line
    ux, uy = pc.u
    if wall[1] == 0:
        ux = ux * (-EC)
    else:
        uy = uy * (-EC)
    f2 = ux * ux + uy * uy
    f = f2.sqrt()
    np_ = Piece(pc.ball, "move", tuple(C), pc.t1, event)
    np_.u, np_.v0 = (ux / f, uy / f), f.range() * pc.v1
    np_.al, np_.be = f2 * pc.al, f2 * pc.be1
    return np_


def build(w, V, phic, h, plan_t):
    """Build the skeleton over the piece phic +/- h (phic an mpf) with
    speeds V (interval). plan_t carries the designated cushions and
    pockets. Returns what certify_strip.build returns, with geometry as
    centred forms."""
    CF.set_piece(h)
    n = len(w["balls"])
    names = ["cue"] + [f"b{i}" for i in range(n)]
    pieces = {nm: [] for nm in names}
    cur = {}
    for nm, p in zip(names, [w["cue"]] + w["balls"]):
        pc = Piece(nm, "rest", (CF(p[0]), CF(p[1])), Z, 0)
        pieces[nm].append(pc)
        cur[nm] = pc
    cue = cur["cue"]
    cue.kind = "move"
    cue.u = CF.cos_sin(phic)
    cue.v0 = V
    cue.al, cue.be = CF(1), CF(0)
    walls = {wl[0]: wl for wl in S.WALLS}
    same_start, designed, potted = [], [], []
    for ev_i, ev in enumerate(plan_t, start=1):
        if ev[0] == "bb":
            a_, b_ = ev[1], ev[2]
            mover, rest = (a_, b_) if cur[a_].kind == "move" else (b_, a_)
            pm, pr = cur[mover], cur[rest]
            if pm.kind != "move" or pr.kind != "rest":
                raise Fail(f"event {ev_i}: not a moving-onto-resting collision")
            s = first_contact(pm.P, pm.u, pr.P, 2 * R)
            if s is None:
                raise Fail(f"event {ev_i}: planned contact does not occur")
            check_single(pm, s)
            C = end_piece(pm, s, ev_i)
            pr.t1, pr.e1 = pm.t1, ev_i
            designed.append((pm, pr))
            nx, ny = (pr.P[0] - C[0]) / (2 * R), (pr.P[1] - C[1]) / (2 * R)
            un = pm.u[0] * nx + pm.u[1] * ny
            if not posr(un):
                raise Fail(f"event {ev_i}: not approaching")
            k = un * ((1 + E) / 2)
            wx, wy = pm.u[0] - k * nx, pm.u[1] - k * ny
            f2 = wx * wx + wy * wy
            f = f2.sqrt()
            k2 = k * k
            nm_ = Piece(mover, "move", C, pm.t1, ev_i)
            nm_.u, nm_.v0 = (wx / f, wy / f), f.range() * pm.v1
            nm_.al, nm_.be = f2 * pm.al, f2 * pm.be1
            nr = Piece(rest, "move", pr.P, pm.t1, ev_i)
            nr.u, nr.v0 = (nx, ny), k.range() * pm.v1
            nr.al, nr.be = k2 * pm.al, k2 * pm.be1
            pieces[mover].append(nm_)
            pieces[rest].append(nr)
            cur[mover], cur[rest] = nm_, nr
            same_start.append((nm_, nr))
        elif ev[0] == "wall":
            pc = cur[ev[1]]
            wall = walls[ev[2]]
            s = wall_distance(pc.P, pc.u, wall)
            if s is None:
                raise Fail(f"event {ev_i}: designated cushion not ahead")
            check_single(pc, s, ignore_wall=wall[0])
            np_ = reflect(pc, wall, s, ev_i)
            pieces[pc.ball].append(np_)
            cur[pc.ball] = np_
        else:                                                # pocket
            pc = cur[ev[1]]
            k = ev[2]
            s = first_contact(pc.P, pc.u, S.POCKETS[k], RHO)
            if s is None:
                raise Fail(f"event {ev_i}: designated pocket not ahead")
            check_single(pc, s, ignore_pocket=k)
            end_piece(pc, s, ev_i)
            cur[pc.ball] = None
            if pc.ball == "cue":
                raise Fail("cue ball pocketed")
            potted.append((pc.ball, S.POCKET_NAMES[k], pc.t1, pc.v1))
    if len(potted) != n:
        raise Fail("plan does not pocket every object ball")
    t_end = max(hi(p[2]) for p in potted)
    # extend the cue (the only ball left) until every stroke has finished
    pc = cur["cue"]
    while True:
        stop_s = pc.v0 * pc.v0 / (2 * A)
        best = None
        for wall in S.WALLS:
            d = wall_distance(pc.P, pc.u, wall)
            if d is not None and (best is None or hi(d.range()) < hi(best[1].range())):
                best = (wall, d)
        wall, s = best
        if pos(s.range() - stop_s):               # comes to rest first
            check_single(pc, stop_s)
            pc.ell, pc.v1, pc.e1 = stop_s, Z, "x"
            pc.t1 = pc.t0 + pc.v0 / A
            if lo(pc.t1) <= t_end:
                P = (pc.P[0] + pc.u[0] * stop_s, pc.P[1] + pc.u[1] * stop_s)
                pieces["cue"].append(Piece("cue", "rest", P, pc.t1, "x"))
            break
        if not pos(stop_s - s.range()):
            raise Fail("cue may stop exactly at a cushion")
        check_single(pc, s, ignore_wall=wall[0])
        np_ = reflect(pc, wall, s, "x")
        if lo(pc.t1) > t_end:
            break
        pieces["cue"].append(np_)
        pc = np_
    return pieces, same_start, designed, potted, t_end


def to_intervals(pieces, same_start, designed):
    """Range-valued copies of the pieces for certify_strip.check_pairs."""
    conv = {}
    for nm, ps in pieces.items():
        out = []
        for pc in ps:
            q = S.Piece(pc.ball, pc.kind, (rng(pc.P[0]), rng(pc.P[1])), pc.t0, pc.e0)
            q.u = None if pc.u is None else (rng(pc.u[0]), rng(pc.u[1]))
            q.ell = None if pc.ell is None else rng(pc.ell)
            q.t1, q.v0, q.v1, q.e1 = pc.t1, pc.v0, pc.v1, pc.e1
            conv[id(pc)] = q
            out.append(q)
        pieces[nm] = out
    ss = [(conv[id(a)], conv[id(b)]) for a, b in same_start]
    ds = [(conv[id(a)], conv[id(b)]) for a, b in designed]
    return pieces, ss, ds


def certify_piece(w, V, phic, h, plan_t):
    pieces, same_start, designed, potted, t_end = build(w, V, phic, h, plan_t)
    pieces, same_start, designed = to_intervals(pieces, same_start, designed)
    n_pairs, n_time = S.check_pairs(pieces, same_start, designed, t_end)
    return potted, n_pairs, n_time


# ------------------------------------------------------------ driver
def exact(x):
    x = mp.mpf(x)
    return f"{'-' if x < 0 else ''}{x.man}*2^{x.exp}" if x.man else "0"


VERBOSE = os.environ.get("CERT_VERBOSE") == "1"


def main():
    path = os.environ.get("CERT_WITNESS", os.path.join(S.HERE, "witness_21.json"))
    w = json.load(open(path))
    args = sys.argv[1:]
    dminus = mp.mpf(args[0] if len(args) > 0 else "4.358e-44")
    dplus = mp.mpf(args[1] if len(args) > 1 else "9.045e-43")
    if len(args) > 3:                        # V0 in V0* + [VA, VB]
        va, vb = mp.mpf(args[2]), mp.mpf(args[3])
        vlabel = f"V0 in V0* + [{args[2]}, {args[3]}]"
    else:                                    # V0 in V0* +/- DV
        dV = mp.mpf(args[2] if len(args) > 2 else "0.01")
        va, vb = -dV, dV
        vlabel = f"V0 in V0* +/- {mp.nstr(dV, 6)}"
    min_width = mp.mpf(os.environ.get("CERT_MIN_WIDTH", "1e-91"))   # give up below this piece width
    S.check_admissible(w)
    import certify_miss as M
    plan = S.plan_events(w)                  # proposal only
    plan_t = M.targets(w, plan)              # designated cushions and pockets
    phistar = mp.mpf(w["phi"])
    V = iv.mpf(w["v0"]) + iv.mpf([va, vb])
    print(f"window: phi* + [-{mp.nstr(dminus, 10)}, {mp.nstr(dplus, 10)}], "
          f"{vlabel}, {S.PREC} bits")
    print(f"  exact offsets: -{exact(dminus)}, +{exact(dplus)}")
    rel = os.path.relpath(os.path.abspath(path), os.path.dirname(S.HERE)).replace(os.sep, "/")
    print(f"  witness: {rel}" + (f", Taylor models of degree {TM_DEGREE}" if TM_DEGREE else ", first-order centred forms"))
    t0 = time.time()
    stack = [(-dminus, dplus)]
    done = []
    n_try = 0
    min_pot = None
    while stack:
        a, b = stack.pop()
        phic, h = phistar + (a + b) / 2, (b - a) / 2
        n_try += 1
        t1 = time.time()
        try:
            potted, n_pairs, n_time = certify_piece(w, V, phic, h, plan_t)
        except Fail as exc:
            if VERBOSE:
                print(f"  piece [{mp.nstr(a, 8)}, {mp.nstr(b, 8)}] width {mp.nstr(b - a, 4)}: "
                      f"FAIL after {time.time() - t1:.1f} s: {exc}", flush=True)
            if b - a < min_width:
                print(f"\nNOT CERTIFIED: piece [{mp.nstr(a, 12)}, {mp.nstr(b, 12)}] "
                      f"below minimum width: {exc}")
                return 1
            m = (a + b) / 2
            stack += [(m, b), (a, m)]
            continue
        v = min(lo(p[3]) for p in potted)
        min_pot = v if min_pot is None else min(min_pot, v)
        done.append((a, b, n_time))
        if VERBOSE:
            print(f"  piece [{mp.nstr(a, 8)}, {mp.nstr(b, 8)}] width {mp.nstr(b - a, 4)}: "
                  f"ok in {time.time() - t1:.1f} s ({n_time} timing pairs)", flush=True)
        if len(done) % 20 == 0:
            print(f"  {len(done)} pieces certified, {len(stack)} pending, "
                  f"{time.time() - t0:.0f} s", flush=True)
    done.sort()
    # the pieces are contiguous and cover the window exactly
    assert done[0][0] == -dminus and done[-1][1] == dplus
    assert all(done[i][1] == done[i + 1][0] for i in range(len(done) - 1))
    widths = [b - a for a, b, _ in done]
    print(f"\npieces certified            : {len(done)} (of {n_try} tried)")
    print(f"piece widths                : {mp.nstr(min(widths), 4)} to "
          f"{mp.nstr(max(widths), 4)}")
    print(f"min pocket entry speed      : {mp.nstr(min_pot, 6)} m/s")
    print(f"wall time                   : {time.time() - t0:.1f} s")
    print("\npieces (offsets from phi*, exact):")
    for a, b, nt in done:
        print(f"  [{exact(a)}, {exact(b)}]  timing pairs {nt}")
    print(f"\nCERTIFIED: every stroke with direction in phi* + "
          f"[-{mp.nstr(dminus, 10)}, {mp.nstr(dplus, 10)}] and {vlabel} "
          f"realises the event list and pockets all {len(w['balls'])} "
          f"object balls.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
