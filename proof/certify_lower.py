"""
certify_lower.py -- the exact lower speed limit of the ribbon.

Along a fixed skeleton the squared speed of a ball at any point of its
path is an affine function of the squared speed of the stroke:

    v^2 = alpha * V0^2 - beta,

where alpha > 0 and beta >= 0 depend on the direction phi alone. The
cue ball starts with alpha = 1, beta = 0; free flight over arc length s
adds 2 a s to beta; a collision multiplies (alpha, beta) of the struck
ball by k^2 and of the striking ball by |u - k n|^2; a cushion contact
multiplies them by |u'|^2 (notation of the speed-independence lemma).
A stroke can realise the event list only if every ball is still moving
when it reaches the designated end of each of its pieces, that is only
if V0 > V_c(phi), where

    V_c(phi)^2 = max over moving pieces of beta_end / alpha.

This program encloses V_c over a window of directions
[phi* - DMINUS, phi* + DPLUS] with the centred forms of
certify_window.py (the window is bisected until the geometry is
decided), and names the critical piece. The enclosure is rigorous: the
skeleton is the one of the window certificates, and alpha, beta are
propagated in the same arithmetic.

Run:  CERT_PREC=400  python certify_lower.py [DMINUS DPLUS]
      (default: the certified window 4.3586e-44 9.0454e-43)
      CERT_PREC=1200 python certify_lower.py 1e-90 1e-90   (about phi*)
"""
import json
import os
import sys
import time

from mpmath import iv, mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import certify_strip as S  # noqa: E402
import certify_window as Wd  # noqa: E402
import certify_miss as M  # noqa: E402

lo, hi = S.lo, S.hi


def critical(pieces):
    """(lo, hi) of max_j beta_end/alpha over the designated moving pieces,
    with the name of the piece attaining the upper bound."""
    best_lo, best_hi, who = None, None, None
    for ps in pieces.values():
        for pc in ps:
            if pc.kind != "move" or pc.al is None or not isinstance(pc.e1, int):
                continue                    # rest pieces, cue after its last event
            r = (pc.be1 / pc.al).range()
            if best_lo is None or lo(r) > best_lo:
                best_lo = lo(r)
            if best_hi is None or hi(r) > best_hi:
                best_hi, who = hi(r), (pc.ball, pc.e0, pc.e1)
    return best_lo, best_hi, who


def main():
    w = json.load(open(os.path.join(S.HERE, "witness_21.json")))
    args = sys.argv[1:]
    dminus = mp.mpf(args[0] if len(args) > 0 else "4.3586e-44")
    dplus = mp.mpf(args[1] if len(args) > 1 else "9.0454e-43")
    S.check_admissible(w)
    plan = S.plan_events(w)                  # proposal only
    plan_t = M.targets(w, plan)              # designated cushions and pockets
    phistar = mp.mpf(w["phi"])
    V = iv.mpf(w["v0"])                      # any speed at which no ball stops
    print(f"lower speed limit over phi* + [-{mp.nstr(dminus, 10)}, "
          f"{mp.nstr(dplus, 10)}], {S.PREC} bits")
    print(f"  exact offsets: -{Wd.exact(dminus)}, +{Wd.exact(dplus)}")
    t0 = time.time()
    stack, done = [(-dminus, dplus)], []
    while stack:
        a, b = stack.pop()
        phic, h = phistar + (a + b) / 2, (b - a) / 2
        try:
            pieces = Wd.build(w, V, phic, h, plan_t)[0]
        except S.Fail as exc:
            if b - a < mp.mpf("1e-91"):
                print(f"NOT DECIDED: piece [{mp.nstr(a, 12)}, {mp.nstr(b, 12)}]: {exc}")
                return 1
            m = (a + b) / 2
            stack += [(m, b), (a, m)]
            continue
        done.append((a, b) + critical(pieces))
    done.sort()
    assert done[0][0] == -dminus and done[-1][1] == dplus
    assert all(done[i][1] == done[i + 1][0] for i in range(len(done) - 1))
    vc_lo = min(d[2] for d in done)
    vc_hi = max(d[3] for d in done)
    whos = {d[4] for d in done}
    print(f"pieces of the window        : {len(done)}")
    print(f"critical piece              : {', '.join(f'{b} from event {e0} to event {e1}' for b, e0, e1 in sorted(whos))}")
    print(f"V_c^2 (max beta_end/alpha)  : [{mp.nstr(mp.mpf(vc_lo), 20)}, {mp.nstr(mp.mpf(vc_hi), 20)}]")
    r = iv.sqrt(iv.mpf([vc_lo, vc_hi]))
    print(f"V_c                         : [{mp.nstr(lo(r), 20)}, {mp.nstr(hi(r), 20)}] m/s")
    print(f"V_c - V0*                   : [{mp.nstr(lo(r) - mp.mpf(w['v0']), 15)}, "
          f"{mp.nstr(hi(r) - mp.mpf(w['v0']), 15)}] m/s")
    print(f"wall time                   : {time.time() - t0:.1f} s")
    print("\nPROVED: for every direction phi in the window, a stroke (V0, phi) with "
          "V0 <= V_c(phi) does not realise the event list: the critical ball "
          "comes to rest at or before the designated end of its piece. Here "
          f"V_c(phi) lies in the interval above.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
