"""
first_order.py -- the first-order (ray-transfer) analysis of the witness.

Explains the size of the sensitivity of ball 21's position at its pocket
to the direction of the stroke (Section "The width of the ribbon" of the
arXiv paper). This is an analysis of the certified skeleton, not a proof
file: the skeleton is built with the centred forms of
proof/certify_window.py at phi* (half-width 1e-90), and all quantities
below are read off its centre values and slopes.

A perturbation of the cue ball's line is described, as in ray optics, by
its lateral offset delta at a reference point and its direction error eps.
  free path s:        (delta, eps) -> (delta + s eps, eps)
  collision, cut theta, at the contact point:
      db   = -delta - 2R cos(theta) eps            (impact parameter)
      delta' = -g(theta) db,  eps' = eps - psi'(theta) db / (2R cos theta)
      g = (1-c)/sqrt(sin^2 + (1-c)^2 cos^2),  c = (1+e)/2,
      psi' = -c (1 - (2-c) cos^2) / (sin^2 + (1-c)^2 cos^2)
      the struck ball's direction turns by -delta/(2R cos theta)
  cushion, incidence alpha to the normal, exit alpha' (tan a' = tan a / e_c):
      (delta, eps) -> (-delta cos a'/cos a, -eps cos^2 a'/(e_c cos^2 a))
Composed along the skeleton and followed by ball 21's path to its pocket,
this reproduces the exact derivative of ball 21's lateral position at the
pocket, which the centred forms carry as a slope. The dominant factor per
collision is A = c ell/(2R cos theta).

Run:  CERT_PREC=400 python analysis/first_order.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "proof"))
from mpmath import iv, mp  # noqa: E402
import certify_strip as S  # noqa: E402
import certify_window as Wd  # noqa: E402
import certify_miss as M  # noqa: E402

mp.prec = iv.prec = S.PREC
R, EC = mp.mpf("0.02625"), mp.mpf("0.90")
C = (1 + mp.mpf("0.95")) / 2


def mid(x):
    return mp.mpf(x.range().mid) if isinstance(x, Wd.CF) else mp.mpf(x.mid)


def slope(x):
    return mp.mpf(x.D.mid)


def psi_prime(th):
    return -C * (1 - (2 - C) * mp.cos(th) ** 2) / (mp.sin(th) ** 2 + (1 - C) ** 2 * mp.cos(th) ** 2)


def g(th):
    return (1 - C) / mp.sqrt(mp.sin(th) ** 2 + (1 - C) ** 2 * mp.cos(th) ** 2)


def main():
    w = json.load(open(os.path.join(S.HERE, "witness_21.json")))
    plan_t = M.targets(w, S.plan_events(w))
    pieces = Wd.build(w, iv.mpf(w["v0"]), mp.mpf(w["phi"]), mp.mpf("1e-90"), plan_t)[0]
    ev = {i + 1: e for i, e in enumerate(plan_t)}
    cue = [p for p in pieces["cue"] if p.kind == "move"]
    normal = {}                       # outgoing unit normal of each collision
    for nm, ps in pieces.items():
        for p in ps:
            if nm != "cue" and p.kind == "move" and isinstance(p.e0, int) and ev[p.e0][0] == "bb":
                normal[p.e0] = (mid(p.u[0]), mid(p.u[1]))

    # walk the cue's pieces: cut angles, legs, and the transfer chain
    rows, legs = [], []               # rows: (event, ball, theta); legs: path length to the next collision
    state = (mp.mpf(0), mp.mpf(1))    # d/dphi: offset 0, direction error 1
    leg = mp.mpf(0)
    for k, p in enumerate(cue):
        e = p.e1
        s = mid(p.ell)
        leg += s
        state = (state[0] + s * state[1], state[1])
        if not isinstance(e, int):
            break
        if ev[e][0] == "wall":
            ax = 0 if ev[e][2] in ("left", "right") else 1
            cin, cout = abs(mid(p.u[ax])), abs(mid(cue[k + 1].u[ax]))
            state = (-(cout / cin) * state[0], -(cout ** 2 / (EC * cin ** 2)) * state[1])
            continue
        n = normal[e]
        th = mp.acos(mid(p.u[0]) * n[0] + mid(p.u[1]) * n[1])
        rows.append((e, ev[e][2], th))
        legs.append(leg)
        leg = mp.mpf(0)
        if len(rows) == len(w["balls"]):          # the last collision: the struck ball's direction
            dth = -state[0] / (2 * R * mp.cos(th))
            break
        db = -(state[0] + 2 * R * mp.cos(th) * state[1])
        state = (-g(th) * db, state[1] - psi_prime(th) * db / (2 * R * mp.cos(th)))

    last_ball = rows[-1][1]
    b21 = [p for p in pieces[last_ball] if p.kind == "move"]
    ell21 = sum(mid(p.ell) for p in b21)
    chain = abs(dth) * ell21
    # exact derivative from the centred forms: lateral position of ball 21 at the pocket
    p = b21[-1]
    ux, uy, ell = mid(p.u[0]), mid(p.u[1]), mid(p.ell)
    dPx, dPy = slope(p.P[0]) + ell * slope(p.u[0]), slope(p.P[1]) + ell * slope(p.u[1])
    exact = abs(-dPx * uy + dPy * ux)

    print(f"{'event':>5} {'ball':>5} {'cut angle':>10} {'next leg (m)':>13} {'A = c l/(2R cos)':>17}")
    dom = []
    for i, (e, ball, th) in enumerate(rows[:-1]):
        A = C * legs[i + 1] / (2 * R * mp.cos(th))
        dom.append(A)
        print(f"{e:>5} {ball:>5} {mp.nstr(th * 180 / mp.pi, 5):>10} {mp.nstr(legs[i + 1], 5):>13} {mp.nstr(A, 5):>17}")
    e, ball, th = rows[-1]
    print(f"{e:>5} {ball:>5} {mp.nstr(th * 180 / mp.pi, 4):>10} {'(ball 21 to pocket)':>13} "
          f"{mp.nstr(ell21 / (2 * R * mp.cos(th)), 5):>17}")
    ds = sorted(dom)
    prod = mp.mpf(1)
    for A in dom:
        prod *= A
    ths = [r[2] * 180 / mp.pi for r in rows[:-1]]
    pps = [abs(psi_prime(r[2])) / C for r in rows[:-1]]
    print(f"\ncut angles of the first {len(dom)} collisions: {mp.nstr(min(ths), 4)} to {mp.nstr(max(ths), 4)} deg")
    print(f"|dpsi/dtheta| / c over them           : {mp.nstr(min(pps), 5)} to {mp.nstr(max(pps), 5)}")
    print(f"dominant factors A                   : {mp.nstr(ds[0], 4)} to {mp.nstr(ds[-1], 4)}, "
          f"median {mp.nstr((ds[len(ds) // 2 - 1] + ds[len(ds) // 2]) / 2, 4)}, product {mp.nstr(prod, 4)}")
    print(f"stroke to first contact              : {mp.nstr(legs[0], 5)} m")
    print(f"ball 21 path to pocket               : {mp.nstr(ell21, 5)} m, last factor {mp.nstr(ell21 / (2 * R * mp.cos(th)), 5)}")
    print(f"product of dominant terms            : {mp.nstr(legs[0] * prod * ell21 / (2 * R * mp.cos(th)), 4)} m/rad")
    print(f"full first-order chain               : {mp.nstr(chain, 6)} m/rad")
    print(f"exact derivative (centred forms)     : {mp.nstr(exact, 6)} m/rad")
    print(f"ratio chain / exact                  : {mp.nstr(chain / exact, 6)}")


if __name__ == "__main__":
    main()
