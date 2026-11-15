#!/usr/bin/env python3
"""
test_taylor_models.py -- containment tests for proof/taylor_models.py.

A Taylor model of a function f on a piece phi = phic + t, |t| <= h, is sound
if f(phic + t) lies in a_0 + a_1 t + ... + a_N t^N + R for every t of the
piece. This script checks that property on random expressions built from the
operations that certify_window.py uses (cos, sin, +, -, *, /, sqrt):

  * the expression is evaluated once in Taylor models, at 400 bits;
  * at nine points t of the piece, including both ends, the same expression
    is evaluated in plain interval arithmetic at 1200 bits, which gives a
    much thinner enclosure J of the true value;
  * the test passes at t if J lies inside the model's enclosure at t, and
    fails if the two are disjoint.

It also checks, for every model that appears, that the remainder R contains
zero (the expansions of 1/x and sqrt(x) are taken about the midpoint of the
constant coefficient, which then lies in the range of the model), and that
1/x and sqrt(x) refuse a model whose range contains zero.

These are tests, not proofs: they exercise the code on finitely many inputs.

usage: python proof/test_taylor_models.py
"""
import os, random, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mpmath import iv, mp
from taylor_models import TM, Fail

LOW, HIGH = 400, 1200


def const(rng):
    return iv.mpf(rng.randint(-64, 64)) / 32          # exact dyadic number


def build(rng, steps):
    """A random straight-line program over two inputs (cos, sin)."""
    prog = []
    for n in range(2, 2 + steps):
        prog.append((rng.choice(("add", "sub", "mul", "lin", "div", "sqrt", "norm")),
                     rng.randrange(n), rng.randrange(n), rng.randint(-64, 64), rng.randint(1, 64)))
    return prog


def run(prog, c, s, sqrt, check=None):
    """Evaluate the program; sqrt is the square root of the number type."""
    v = [c, s]
    for op, i, j, k1, k2 in prog:
        a, b = v[i], v[j]
        p, q = iv.mpf(k1) / 32, iv.mpf(k2) / 32          # q > 0
        if op == "add":
            r = a + b
        elif op == "sub":
            r = a - b * p
        elif op == "mul":
            r = a * b
        elif op == "lin":
            r = a * p + q
        elif op == "div":
            r = a / (b * b + q)                           # denominator >= q > 0
        elif op == "sqrt":
            r = sqrt(a * a + q)                           # argument >= q > 0
        else:
            r = a / sqrt(b * b + 1)                       # |r| <= |a|
        if check is not None:
            check(r)
        v.append(r)
    return v


def main():
    t0 = time.time()
    rng = random.Random(20261115)
    stats = {"points": 0, "inside": 0, "disjoint": 0, "partial": 0, "models": 0, "rem_without_zero": 0,
             "skipped": 0}
    for N in (1, 2, 4, 6):
        for h_s in ("1e-30", "1e-12", "1e-6", "1e-3", "0.05"):
            for trial in range(6):
                mp.prec = iv.prec = LOW
                TM.N = N
                h = mp.mpf(h_s)
                TM.set_piece(h)
                phic = mp.mpf(rng.randint(-3000, 3000)) / 1024
                prog = build(rng, 24)

                def check(r):
                    stats["models"] += 1
                    if not (r.R.a <= 0 <= r.R.b):
                        stats["rem_without_zero"] += 1

                c, s = TM.cos_sin(phic)
                check(c); check(s)
                try:
                    models = run(prog, c, s, lambda x: x.sqrt(), check)
                except Fail:
                    stats["skipped"] += 1                 # enclosure too wide to take 1/x or sqrt
                    continue
                pts = [h * k / 4 for k in range(-4, 5)]
                encl = []
                for t in pts:
                    T = iv.mpf(t)
                    row = []
                    for m in models:
                        acc = m.a[N]
                        for k in range(N - 1, -1, -1):
                            acc = m.a[k] + acc * T
                        row.append(acc + m.R)
                    encl.append(row)
                mp.prec = iv.prec = HIGH
                for t, row in zip(pts, encl):
                    P = iv.mpf(phic) + iv.mpf(t)
                    exact = run(prog, iv.cos(P), iv.sin(P), iv.sqrt)
                    for E, J in zip(row, exact):
                        stats["points"] += 1
                        if E.a <= J.a and J.b <= E.b:
                            stats["inside"] += 1
                        elif J.b < E.a or E.b < J.a:
                            stats["disjoint"] += 1
                        else:
                            stats["partial"] += 1
    # 1/x and sqrt(x) must refuse a range that contains zero
    mp.prec = iv.prec = LOW
    TM.N = 4
    TM.set_piece(mp.mpf("0.5"))
    c, s = TM.cos_sin(mp.mpf(0))
    refused = 0
    for f in (lambda: s.inv(), lambda: s.sqrt(), lambda: (s * s - 1).sqrt(), lambda: 1 / (c - 1)):
        try:
            f()
        except Fail:
            refused += 1
    print("test_taylor_models.py  degrees 1, 2, 4, 6; half-widths 1e-30 to 0.05; %d and %d bits" % (LOW, HIGH))
    print("models built                 : %d (remainder without zero: %d)" % (stats["models"], stats["rem_without_zero"]))
    print("programs skipped             : %d (an enclosure was too wide for 1/x or sqrt)" % stats["skipped"])
    print("point checks                 : %d" % stats["points"])
    print("  true value inside the model: %d" % stats["inside"])
    print("  disjoint from the model    : %d" % stats["disjoint"])
    print("  overlapping, not inside    : %d" % stats["partial"])
    print("zero ranges refused          : %d of 4" % refused)
    print("wall time                    : %.1f s" % (time.time() - t0))
    ok = (stats["disjoint"] == 0 and stats["partial"] == 0 and stats["rem_without_zero"] == 0
          and refused == 4 and stats["inside"] > 0)
    print("ALL CONTAINMENT TESTS PASSED" if ok else "CONTAINMENT TESTS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
