"""
make_gif.py -- animate the verified 21-ball single-stroke clearance.

Replays witness_21.json in the event-driven engine, sampling every ball's
position at a fixed frame rate, and renders an animated GIF: the cue
traces its path in red, each ball vanishes as it is pocketed, and a
counter tracks the running total.

Run:  python3 make_gif.py [witness_21.json] [out.gif]
Deps: matplotlib, imageio, model_sim.py
"""
import json
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import imageio.v2 as imageio

import model_sim as M

R, L, W, RHO, A = M.R, M.L, M.W, M.RHO, M.A
V_MIN = M.V_MIN

FPS = 30
DT = 1.0 / FPS
TAIL_SECONDS = 1.2          # how long the red trail lingers behind the cue
HOLD_START_FRAMES = 25      # freeze on the full setup before the break
HOLD_END_FRAMES = 40        # freeze on the final cleared table


def sample_run(w):
    """Replay the witness, advancing event by event, and capture a frame
    snapshot at each fixed frame time by interpolating along the current
    decelerating segment. Returns a list of (snapshot_dict, n_potted)."""
    balls = [M.Ball(w["cue"][0], w["cue"][1], "cue")]
    for i, p in enumerate(w["balls"]):
        balls.append(M.Ball(p[0], p[1], f"b{i}"))
    M.shoot(balls[0], w["v0"], w["phi"])
    sim = M.Simulator(balls, record_events=False)

    frames = []
    t_next = 0.0
    guard = 0
    while guard < 200000:
        guard += 1
        # peek at the next event time without applying it
        ev = sim._next_event()
        if ev is None:
            break
        dt_ev = ev[0]
        t_event = sim.t + dt_ev
        # emit every frame whose timestamp falls in [sim.t, t_event]
        while t_next <= t_event + 1e-12:
            tau = t_next - sim.t
            if tau >= -1e-9:
                frames.append(_snapshot(balls, max(tau, 0.0)))
            t_next += DT
        # now actually advance one event (recomputes the same event and
        # applies it); positions become exact at t_event
        if not sim.step():
            break
        # stop once everything is at rest or pocketed
        if not any(b.active and b._s >= V_MIN for b in balls):
            break
    frames.append(_snapshot(balls, 0.0))
    return frames


def _snapshot(balls, tau):
    """Positions of all balls tau seconds into the current straight-line
    decelerating segment, without mutating the simulator."""
    out = {}
    n_potted = 0
    for b in balls:
        if not b.active:
            out[b.name] = None
            if b.name != "cue":
                n_potted += 1
            continue
        s = b._s
        if s < V_MIN or tau <= 0:
            out[b.name] = (b.x, b.y)
            continue
        ux, uy = b.vx / s, b.vy / s
        t_stop = s / A
        t = min(tau, t_stop)
        x = b.x + b.vx * t - 0.5 * A * t * t * ux
        y = b.y + b.vy * t - 0.5 * A * t * t * uy
        out[b.name] = (x, y)
    return out, n_potted


STRIKE_FRAMES = 18          # frames of the stick thrusting into the cue
STICK_LEN = 1.15            # drawn length of the cue stick (m)
STICK_PULLBACK = 0.28       # how far behind the ball the tip starts (m)


def _fig_to_rgb(fig):
    import numpy as np
    fig.canvas.draw()
    buf = fig.canvas.buffer_rgba()
    img = np.asarray(buf)[:, :, :3].copy()
    plt.close(fig)
    return img


def _new_axes():
    fig, ax = plt.subplots(figsize=(9, 4.7), dpi=80)
    ax.add_patch(Rectangle((0, 0), L, W, fill=False, lw=1.5))
    for (cx, cy), nm in zip(M.POCKETS, M.POCKET_NAMES):
        ax.add_patch(Circle((cx, cy), RHO, fc="0.85", ec="0.5",
                            lw=0.8, zorder=0))
    ax.set_title("One stroke, twenty-one balls", fontsize=11)
    ax.set_xlim(-0.1, L + 0.1)
    ax.set_ylim(-0.1, W + 0.1)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax


def _draw_setup(ax, w, snap):
    """Draw the object balls still on the table, numbered, with their
    dotted routes to their target pockets."""
    n = len(w["balls"])
    pockets = w["pockets"]
    for i in range(n):
        p = snap.get(f"b{i}")
        if p is None:
            continue
        ax.add_patch(Circle(p, R, fc="0.4", ec="black", zorder=3))
        ax.annotate(str(i + 1), p, textcoords="offset points",
                    xytext=(3, 4), fontsize=6)
        gx, gy = M.POCKETS[pockets[i]]
        ax.plot([p[0], gx], [p[1], gy], lw=0.4, ls=":",
                color="0.7", zorder=1)


def _draw_cue_stick(ax, cue_xy, phi, tip_gap, alpha=1.0):
    """Draw a cue stick pointing at the cue ball along direction phi, with
    its tip a distance tip_gap behind the ball surface (tip_gap=0 means
    touching). The stick recedes in the -phi direction."""
    dx, dy = math.cos(phi), math.sin(phi)
    tip_x = cue_xy[0] - (R + tip_gap) * dx
    tip_y = cue_xy[1] - (R + tip_gap) * dy
    butt_x = tip_x - STICK_LEN * dx
    butt_y = tip_y - STICK_LEN * dy
    # shaft
    ax.plot([butt_x, tip_x], [butt_y, tip_y], lw=3.2, color="#b5822e",
            solid_capstyle="round", zorder=5, alpha=alpha)
    # ferrule + tip
    fx = tip_x + 0.05 * dx
    fy = tip_y + 0.05 * dy
    ax.plot([tip_x, fx], [tip_y, fy], lw=3.2, color="0.25",
            solid_capstyle="round", zorder=6, alpha=alpha)


def render(w, frames, out):
    tail_len = int(TAIL_SECONDS * FPS)
    phi = w["phi"]
    cue0 = tuple(w["cue"])

    # snapshot of the initial setup (all balls at rest, cue at start)
    setup_snap = frames[0][0]

    images = []

    # ---- opening hold on the full setup ----
    for _ in range(HOLD_START_FRAMES):
        fig, ax = _new_axes()
        _draw_setup(ax, w, setup_snap)
        ax.add_patch(Circle(cue0, R, fc="white", ec="black", lw=1.3,
                            zorder=4))
        ax.text(0.02, 1.06, "potted: 0/21", transform=ax.transAxes,
                fontsize=11, va="top", fontweight="bold")
        fig.tight_layout()
        images.append(_fig_to_rgb(fig))

    # ---- strike: stick appears behind the ball and thrusts in ----
    for s in range(STRIKE_FRAMES):
        frac = s / (STRIKE_FRAMES - 1)          # 0 -> 1
        tip_gap = STICK_PULLBACK * (1 - frac)   # recedes to contact
        alpha = min(1.0, 0.3 + 1.4 * frac)      # fades in as it thrusts
        fig, ax = _new_axes()
        _draw_setup(ax, w, setup_snap)
        _draw_cue_stick(ax, cue0, phi, tip_gap, alpha=alpha)
        ax.add_patch(Circle(cue0, R, fc="white", ec="black", lw=1.3,
                            zorder=4))
        ax.text(0.02, 1.06, "potted: 0/21", transform=ax.transAxes,
                fontsize=11, va="top", fontweight="bold")
        fig.tight_layout()
        images.append(_fig_to_rgb(fig))

    # ---- the break and clearance (stick is gone the instant it launches) ----
    cue_pts = [snap.get("cue") for snap, _ in frames]
    for fi, (snap, n_potted) in enumerate(frames):
        fig, ax = _new_axes()
        _draw_setup(ax, w, snap)

        lo = max(0, fi - tail_len)
        seg = [q for q in cue_pts[lo:fi + 1] if q is not None]
        if len(seg) >= 2:
            ax.plot([q[0] for q in seg], [q[1] for q in seg],
                    lw=1.1, color="C3", alpha=0.85, zorder=2)

        cue = snap.get("cue")
        if cue is not None:
            ax.add_patch(Circle(cue, R, fc="white", ec="black",
                                lw=1.3, zorder=4))

        ax.text(0.02, 1.06, f"potted: {n_potted}/21",
                transform=ax.transAxes, fontsize=11, va="top",
                fontweight="bold")
        fig.tight_layout()
        images.append(_fig_to_rgb(fig))

    # ---- hold on the cleared table ----
    for _ in range(HOLD_END_FRAMES):
        images.append(images[-1])

    imageio.mimsave(out, images, fps=FPS, loop=0)
    print(f"wrote {out}  ({len(images)} frames, {len(images)/FPS:.1f}s)")


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "witness_21.json"
    out = sys.argv[2] if len(sys.argv) > 2 else "witness_clearance.gif"
    w = json.load(open(src))
    frames = sample_run(w)
    print(f"sampled {len(frames)} frames from the run")
    render(w, frames, out)
