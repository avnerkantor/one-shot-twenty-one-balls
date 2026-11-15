#!/usr/bin/env python3
"""
crosscheck_arb.py -- cross-check of Theorems B to E with the Arb certifier.

Theorems B to E are certified by certify_strip.py, certify_window.py,
certify_miss.py and certify_lower.py, which share code and all use mpmath.
This script checks their statements at single strokes with certify_arb.py,
which shares no code with them and uses Arb ball arithmetic.

For each stroke below, certify_arb.run() certifies the evolution event by
event. The certified events are then compared, ball by ball, with the event
list of the configuration (the events that certify_arb.py itself certifies
at the stroke stored in the JSON file):

  realizes   every object ball undergoes exactly its events of the list, in
             order, and the cue ball undergoes first its events of the list;
             after them the cue ball only meets cushions or stops
  deviates   some ball provably undergoes an event that is not its next
             event of the list (the run is certified up to that event)

Each stroke is a single point of the stroke plane (the speed is an exact
decimal, enclosed in a ball of Arb), so this is a check of the theorems at
finitely many strokes, not a second proof of them.

usage: python proof/crosscheck_arb.py [--prec 1200]
"""
import argparse, json, os, re, sys, time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import certify_arb as CA
from flint import arb, ctx

EVENT = re.compile(r"\s*(\d+)\s+(bb|wall|pocket|stop)\s+(\S+)(?:\s+(\S+))?\s+t>=")


def q(s):
    """Exact rational from a decimal string, enclosed in a ball."""
    f = Fraction(s)
    return arb(f.numerator) / arb(f.denominator)


def events(data, V, PHI, C):
    """Events certified by certify_arb at one stroke, and how the run ended."""
    lines = []
    try:
        CA.run(V, PHI, data, C, lines.append)
        end = "all 21 pocketed"
    except CA.NotCertified as ex:
        end = "stopped: %s" % ex
    ev = []
    for ln in lines:
        m = EVENT.match(ln)
        if m:
            ev.append((m.group(2), m.group(3), m.group(4)))
    return ev, end


def by_ball(ev):
    seq = {}
    for kind, a, b in ev:
        if kind == "bb":
            seq.setdefault(a, []).append(("bb", b))
            seq.setdefault(b, []).append(("bb", a))
        elif kind == "stop":
            seq.setdefault(a, []).append(("stop",))
        else:
            seq.setdefault(a, []).append((kind, b))
    return seq


def compare(ev, end, ref):
    """('realizes', note) or ('deviates', note) or ('undecided', note)."""
    seq, want = by_ball(ev), by_ball(ref)
    # first event, in the certified order, that is not the ball's next event of the list
    pos = {}
    for kind, a, b in ev:
        parts = [(a, ("bb", b)), (b, ("bb", a))] if kind == "bb" else \
                [(a, ("stop",))] if kind == "stop" else [(a, (kind, b))]
        for name, item in parts:
            k = pos.get(name, 0)
            rest = want.get(name, [])
            if k < len(rest):
                if item != rest[k]:
                    return "deviates", "%s: %s in place of %s" % (name, " ".join(item), " ".join(rest[k]))
                pos[name] = k + 1
            elif name != "cue" or item[0] not in ("wall", "stop"):
                return "deviates", "%s: extra event %s" % (name, " ".join(item))
    if end != "all 21 pocketed":
        return "undecided", end
    for name, rest in want.items():
        if pos.get(name, 0) != len(rest):
            return "undecided", "%s has events left" % name
    extra = len(seq.get("cue", [])) - len(want["cue"])
    return "realizes", ("cue ball: %d later event(s) before the last pocketing" % extra) if extra else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prec", type=int, default=1200)
    args = ap.parse_args()
    ctx.prec = args.prec
    C = CA.constants()
    t0 = time.time()
    bad = 0
    print("crosscheck_arb.py  model M, Arb ball arithmetic, prec=%d bits, single strokes" % args.prec)

    def block(title, path, counts, strokes):
        nonlocal bad
        data = json.load(open(os.path.join(ROOT, path)))
        v0, phi = arb(float(data["v0"])), arb(float(data["phi"]))
        ref, end = events(data, v0, phi, C)
        n = (len(ref), sum(e[0] == "bb" for e in ref), sum(e[0] == "wall" for e in ref),
             sum(e[0] == "pocket" for e in ref))
        print("")
        print("%s (%s)" % (title, path))
        print("  event list at the stored stroke: %d events (bb %d, wall %d, pocket %d); %s" % (n + (end,)))
        if n != counts or end != "all 21 pocketed":
            print("  UNEXPECTED event list")
            bad += 1
            return
        for speed, dphi, expect, why in strokes:
            V = v0 + q(speed[1:]) if speed.startswith("@") else q(speed)
            got, note = compare(*events(data, V, phi + q(dphi), C), ref)
            ok = got == expect
            bad += not ok
            sp = ("V0* %s" % speed[1:]) if speed.startswith("@") else speed
            print("  V0 = %-14s phi = phi* %+.5e : %-9s %s%s" % (
                sp, float(Fraction(dphi)), got, ("(" + note + ") ") if note else "",
                why if ok else "UNEXPECTED, expected " + expect))

    lo, hi = "-4.3586e-44", "9.0454e-43"
    chain = []
    for v in ("9.321315", "9.4", "@-2.3", "10", "10.5", "11", "11.5", "12"):
        for d in (lo, "0", hi):
            chain.append((v, d, "realizes", "Theorems C(i), D(ii)"))
    for d in ("-1e-90", "0", "1e-90"):
        chain.append(("@-2.3", d, "realizes", "Theorem B"))
        chain.append(("9.3212701695", d, "realizes", "Theorem D(iii)"))
    for v in ("9.321315", "10.5", "12"):
        chain.append((v, "9.0455e-43", "deviates", "Theorem C(ii)"))
        chain.append((v, "-4.3587e-44", "deviates", "Theorem C(ii)"))
    for d in (lo, "0", hi):
        chain.append(("9.32035", d, "deviates", "Theorem D(i), below every V_c"))
    chain.append(("9.32127016", "0", "deviates", "Theorem D(i), below V_c(phi*)"))
    block("chain witness", "proof/witness_21.json", (61, 21, 19, 21), chain)

    cluster = []
    for v in ("@-0.01", "@-0.005", "@0"):
        for d in ("-2.5e-5", "0", "4.45e-4", "9.15e-4"):
            cluster.append((v, d, "realizes", "Theorem E"))
    block("cluster", "construction/cluster_21.json", (57, 21, 15, 21), cluster)

    print("")
    print("wall time %.1f s" % (time.time() - t0))
    print("CROSS-CHECK PASSED: every stroke behaves as the theorems state." if not bad
          else "CROSS-CHECK FAILED at %d stroke(s)." % bad)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
