"""
make_figures.py -- regenerate the figures used in the papers.

  pk_decay.pdf/png       : estimated P(k) with Wilson intervals + log-linear
                           fit + extrapolation to k=21  (from mc_results.json)
  witness_layout.pdf/png : the verified 21-ball clearance, cue path in red
                           (from proof/witness_21.json, replayed in the
                           80-digit exact engine model_sim_mp)

  cluster_layout.pdf/png : the near-frozen cluster of Theorem E, replayed
                           the same way (from construction/cluster_21.json),
                           with an enlarged view of the first contacts
  amplification.pdf/png  : the twenty factors A = c l/(2R cos theta) of the
                           first-order analysis (from the saved output of
                           analysis/first_order.py)

Run:  python figures/make_figures.py   (writes into figures/)
Dependencies: matplotlib, mpmath.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "model"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

import model_sim_mp as M


def replay(w):
    """Replay a witness in the 80-digit engine up to the last object-ball
    pocketing. Returns every ball's position at every event (straight
    lines in between, so these vertices draw the exact routes) and the
    number of object balls pocketed."""
    balls = [M.Ball(w["cue"][0], w["cue"][1], "cue")]
    for i, p in enumerate(w["balls"]):
        balls.append(M.Ball(p[0], p[1], f"b{i}"))
    M.shoot(balls[0], w["v0"], w["phi"])
    sim = M.Simulator(balls, record_events=True)
    paths = {b.name: [(float(b.x), float(b.y))] for b in balls}
    n_obj = len(balls) - 1
    potted = 0
    for _ in range(9000):
        if not sim.step():
            break
        for b in balls:
            xy = (float(b.x), float(b.y))
            if xy != paths[b.name][-1] and (b.active or
                                            sim.events[-1][2] == b.name):
                paths[b.name].append(xy)
        potted = sum(1 for _, nm, _, _ in sim.pocketed if nm != "cue")
        if potted == n_obj:
            break
    return paths, potted


def draw_clearance(ax, w, paths, labels=True, lw_scale=1.0):
    """Table, pockets, balls at their initial positions and all routes."""
    L, W, R, RHO = (float(v) for v in (M.L, M.W, M.R, M.RHO))
    pockets = [(float(x), float(y)) for x, y in M.POCKETS]
    ax.add_patch(Rectangle((0, 0), L, W, fill=False, lw=1.5))
    for cx, cy in pockets:
        ax.add_patch(Circle((cx, cy), RHO, fc="0.85", ec="0.5",
                            lw=0.8, zorder=0))
    ax.add_patch(Circle(w["cue"], R, fc="white", ec="black", lw=1.3,
                        zorder=3))
    for i, p in enumerate(w["balls"]):
        ax.add_patch(Circle(p, R, fc="0.4", ec="black", zorder=3))
        if labels == "inside":
            ax.text(p[0], p[1], str(i + 1), ha="center", va="center_baseline",
                    color="white", fontsize=5, fontweight="bold", zorder=4)
        elif labels:
            ax.annotate(str(i + 1), p, textcoords="offset points",
                        xytext=(3, 4), fontsize=6)
    for name, pts in paths.items():
        if name == "cue":
            continue
        ax.plot([q[0] for q in pts], [q[1] for q in pts],
                lw=0.6 * lw_scale, ls=":", color="0.45", zorder=1)
    pts = paths["cue"]
    ax.plot([q[0] for q in pts], [q[1] for q in pts],
            lw=0.8 * lw_scale, color="C3", alpha=0.7, zorder=2)
    ax.set_aspect("equal")
    ax.axis("off")


def fig_cluster(src=os.path.join(ROOT, "construction", "cluster_21.json"),
                out=os.path.join(HERE, "cluster_layout.pdf"),
                zoom=None):
    """The cluster of Theorem E: the whole table, and an enlarged view of
    the region `zoom` = (x0, x1, y0, y1), by default the bounding box of
    the object balls."""
    if not os.path.exists(src):
        print(f"skip {out}: {src} not found")
        return
    w = json.load(open(src))
    paths, potted = replay(w)
    assert potted == len(w["balls"]), f"replay pockets {potted} balls"
    L, W, R = (float(v) for v in (M.L, M.W, M.R))
    if zoom is None:
        xs = [p[0] for p in w["balls"]] + [w["cue"][0]]
        ys = [p[1] for p in w["balls"]] + [w["cue"][1]]
        m = 3 * R
        zoom = (min(xs) - m, max(xs) + m, min(ys) - m, max(ys) + m)
    x0, x1, y0, y1 = zoom
    fig, (ax, az) = plt.subplots(
        1, 2, figsize=(9, 3.3),
        gridspec_kw={"width_ratios": [(L + 0.2) / (W + 0.2),
                                      (x1 - x0) / (y1 - y0)]})
    draw_clearance(ax, w, paths, labels=False)
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                           lw=0.6, ls="--", ec="0.3", zorder=4))
    ax.set_xlim(-0.1, L + 0.1)
    ax.set_ylim(-0.1, W + 0.1)
    draw_clearance(az, w, paths, labels="inside", lw_scale=1.2)
    az.annotate("cue", w["cue"], textcoords="offset points",
                xytext=(-4, 6), fontsize=7)
    az.set_xlim(x0, x1)
    az.set_ylim(y0, y1)
    az.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                           lw=0.6, ls="--", ec="0.3", zorder=4,
                           clip_on=False))
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
    print(f"wrote {out}")


def read_amplification(src):
    """The cut angles and factors A of the first twenty collisions, from
    the table printed by analysis/first_order.py."""
    rows = []
    for line in open(src):
        f = line.split()
        if len(f) == 5 and f[0].isdigit() and f[1].startswith("b"):
            rows.append((int(f[1][1:]) + 1, float(f[2]), float(f[4])))
    return rows


def fig_amplification(
        src=os.path.join(ROOT, "analysis", "first_order_output.txt"),
        out=os.path.join(HERE, "amplification.pdf")):
    if not os.path.exists(src):
        print(f"skip {out}: {src} not found "
              "(CERT_PREC=400 python analysis/first_order.py > it)")
        return
    rows = read_amplification(src)
    assert len(rows) == 20, len(rows)
    ks = [r[0] for r in rows]
    As = [r[2] for r in rows]
    med = sorted(As)[len(As) // 2 - 1: len(As) // 2 + 1]
    med = math.sqrt(med[0] * med[1])
    fig, ax = plt.subplots(figsize=(5.6, 2.7))
    ax.bar(ks, [a - 1 for a in As], bottom=1, width=0.62, color="0.4",
           zorder=2)
    ax.set_yscale("log")
    ax.set_ylim(1, 400)
    ax.set_xlim(0.3, 20.7)
    ax.set_xticks(range(1, 21))
    ax.set_yticks([1, 10, 100])
    ax.set_yticklabels(["1", "10", "100"])
    ax.tick_params(axis="both", labelsize=8, length=2.5, color="0.5")
    ax.yaxis.grid(True, which="major", lw=0.5, color="0.85", zorder=0)
    ax.axhline(med, lw=0.8, ls="--", color="0.15", zorder=3)
    ax.annotate(f"median {med:.0f}", (20.6, med), ha="right", va="bottom",
                textcoords="offset points", xytext=(0, 2), fontsize=7.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("0.5")
    ax.set_xlabel("collision", fontsize=9)
    ax.set_ylabel("$A(\\theta, \\ell)$", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
    print(f"wrote {out}")


def fig_witness(src=os.path.join(ROOT, "proof", "witness_21.json"),
                out=os.path.join(HERE, "witness_layout.pdf")):
    if not os.path.exists(src):
        print(f"skip {out}: {src} not found")
        return
    w = json.load(open(src))
    fig, ax = plt.subplots(figsize=(9, 4.7))
    L, W, R, RHO = (float(v) for v in (M.L, M.W, M.R, M.RHO))
    pockets = [(float(x), float(y)) for x, y in M.POCKETS]
    ax.add_patch(Rectangle((0, 0), L, W, fill=False, lw=1.5))
    for (cx, cy), nm in zip(pockets, M.POCKET_NAMES):
        ax.add_patch(Circle((cx, cy), RHO, fc="0.85", ec="0.5",
                            lw=0.8, zorder=0))

    ax.add_patch(Circle(w["cue"], R, fc="white", ec="black", lw=1.3,
                        zorder=3))
    ax.annotate("cue", w["cue"], textcoords="offset points",
                xytext=(-4, 6), fontsize=7)

    for i, p in enumerate(w["balls"]):
        ax.add_patch(Circle(p, R, fc="0.4", ec="black", zorder=3))
        ax.annotate(str(i + 1), p, textcoords="offset points",
                    xytext=(3, 4), fontsize=6)

    # replay up to the last object-ball pocketing (the certified part),
    # recording every ball's position at every event: straight lines in
    # between, so these vertices draw the exact routes
    balls = [M.Ball(w["cue"][0], w["cue"][1], "cue")]
    for i, p in enumerate(w["balls"]):
        balls.append(M.Ball(p[0], p[1], f"b{i}"))
    M.shoot(balls[0], w["v0"], w["phi"])
    sim = M.Simulator(balls, record_events=True)
    paths = {b.name: [(float(b.x), float(b.y))] for b in balls}
    n_obj = len(balls) - 1
    for _ in range(9000):
        if not sim.step():
            break
        for b in balls:
            xy = (float(b.x), float(b.y))
            if xy != paths[b.name][-1] and (b.active or
                                            sim.events[-1][2] == b.name):
                paths[b.name].append(xy)
        if sum(1 for _, nm, _, _ in sim.pocketed if nm != "cue") == n_obj:
            break
    for b in balls[1:]:
        pts = paths[b.name]
        ax.plot([q[0] for q in pts], [q[1] for q in pts], lw=0.6, ls=":",
                color="0.45", zorder=1)
    pts = paths["cue"]
    ax.plot([q[0] for q in pts], [q[1] for q in pts],
            lw=0.8, color="C3", alpha=0.7, zorder=2)

    ax.set_xlim(-0.1, L + 0.1)
    ax.set_ylim(-0.1, W + 0.1)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    den = 1 + z * z / n
    ctr = (p + z * z / (2 * n)) / den
    hw = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(ctr - hw, 1e-12), min(ctr + hw, 1.0))


def fig_pk(src=os.path.join(ROOT, "simulation", "mc_results.json"),
           out=os.path.join(HERE, "pk_decay.pdf")):
    if not os.path.exists(src):
        print(f"skip {out}: {src} not found")
        return
    d = json.load(open(src))
    n = d["n"]
    counts = {int(k): v for k, v in d["counts"].items()}
    ks = sorted(counts)
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for k in ks:
        c = counts[k]
        p = c / n
        lo, hi = wilson(c, n)
        ax.errorbar([k], [max(p, 1e-12)],
                    yerr=[[max(p - lo, 0)], [max(hi - p, 0)]],
                    fmt="o", ms=4, color="0.15", capsize=2, lw=0.9)
    fit_pts = [(k, counts[k] / n) for k in ks if k >= 1
               and counts[k] >= 5]
    if len(fit_pts) >= 3:
        xs = [k for k, _ in fit_pts]
        ys = [math.log10(p) for _, p in fit_pts]
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / \
            sum((x - mx) ** 2 for x in xs)
        a = my - b * mx
        grid = [xs[0], 21]
        ax.plot(grid, [10 ** (a + b * x) for x in grid], ls="--",
                lw=1.0, color="0.4",
                label=f"$\\log_{{10}}P = {a:.2f} {b:+.2f}k$")
        ax.plot([21], [10 ** (a + b * 21)], marker="s", mfc="white",
                mec="0.2", ms=6)
        ax.legend(frameon=False, fontsize=8)
    ax.set_yscale("log")
    ax.set_xlabel("$k$ (object balls pocketed)")
    ax.set_ylabel("$\\hat P(k)$")
    ax.set_xticks(range(0, 22, 3))
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    fig_witness()
    fig_pk()
    fig_cluster()
    fig_amplification()
