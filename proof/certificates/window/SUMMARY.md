# Direction window x speed ribbon: phi in phi* + [-4.3586e-44, 9.0454e-43], V0 in V0* + [-2.578685, +0.1]

48 sub-strips of width 0.05 m/s, each certified over the whole window by proof/certify_window.py (400 bits, centred forms in phi).

| V0 offset | result | pieces |
|---|---|---|
| [-2.3, -2.25] | CERTIFIED | 23 (of 45 tried) |
| [-2.25, -2.2] | CERTIFIED | 23 (of 45 tried) |
| [-2.2, -2.15] | CERTIFIED | 23 (of 45 tried) |
| [-2.15, -2.1] | CERTIFIED | 23 (of 45 tried) |
| [-2.1, -2.05] | CERTIFIED | 23 (of 45 tried) |
| [-2.05, -2] | CERTIFIED | 23 (of 45 tried) |
| [-2, -1.95] | CERTIFIED | 23 (of 45 tried) |
| [-1.95, -1.9] | CERTIFIED | 23 (of 45 tried) |
| [-1.9, -1.85] | CERTIFIED | 23 (of 45 tried) |
| [-1.85, -1.8] | CERTIFIED | 23 (of 45 tried) |
| [-1.8, -1.75] | CERTIFIED | 23 (of 45 tried) |
| [-1.75, -1.7] | CERTIFIED | 23 (of 45 tried) |
| [-1.7, -1.65] | CERTIFIED | 23 (of 45 tried) |
| [-1.65, -1.6] | CERTIFIED | 23 (of 45 tried) |
| [-1.6, -1.55] | CERTIFIED | 23 (of 45 tried) |
| [-1.55, -1.5] | CERTIFIED | 23 (of 45 tried) |
| [-1.5, -1.45] | CERTIFIED | 23 (of 45 tried) |
| [-1.45, -1.4] | CERTIFIED | 23 (of 45 tried) |
| [-1.4, -1.35] | CERTIFIED | 23 (of 45 tried) |
| [-1.35, -1.3] | CERTIFIED | 23 (of 45 tried) |
| [-1.3, -1.25] | CERTIFIED | 23 (of 45 tried) |
| [-1.25, -1.2] | CERTIFIED | 23 (of 45 tried) |
| [-1.2, -1.15] | CERTIFIED | 23 (of 45 tried) |
| [-1.15, -1.1] | CERTIFIED | 23 (of 45 tried) |
| [-1.1, -1.05] | CERTIFIED | 23 (of 45 tried) |
| [-1.05, -1] | CERTIFIED | 23 (of 45 tried) |
| [-1, -0.95] | CERTIFIED | 23 (of 45 tried) |
| [-0.95, -0.9] | CERTIFIED | 23 (of 45 tried) |
| [-0.9, -0.85] | CERTIFIED | 23 (of 45 tried) |
| [-0.85, -0.8] | CERTIFIED | 23 (of 45 tried) |
| [-0.8, -0.75] | CERTIFIED | 23 (of 45 tried) |
| [-0.75, -0.7] | CERTIFIED | 23 (of 45 tried) |
| [-0.7, -0.65] | CERTIFIED | 23 (of 45 tried) |
| [-0.65, -0.6] | CERTIFIED | 23 (of 45 tried) |
| [-0.6, -0.55] | CERTIFIED | 23 (of 45 tried) |
| [-0.55, -0.5] | CERTIFIED | 23 (of 45 tried) |
| [-0.5, -0.45] | CERTIFIED | 23 (of 45 tried) |
| [-0.45, -0.4] | CERTIFIED | 23 (of 45 tried) |
| [-0.4, -0.35] | CERTIFIED | 23 (of 45 tried) |
| [-0.35, -0.3] | CERTIFIED | 23 (of 45 tried) |
| [-0.3, -0.25] | CERTIFIED | 23 (of 45 tried) |
| [-0.25, -0.2] | CERTIFIED | 23 (of 45 tried) |
| [-0.2, -0.15] | CERTIFIED | 23 (of 45 tried) |
| [-0.15, -0.1] | CERTIFIED | 23 (of 45 tried) |
| [-0.1, -0.05] | CERTIFIED | 23 (of 45 tried) |
| [-0.05, 0] | CERTIFIED | 23 (of 45 tried) |
| [0, 0.05] | CERTIFIED | 23 (of 45 tried) |
| [0.05, 0.1] | CERTIFIED | 23 (of 45 tried) |

`window_dV0.01.txt` is a single run over V0* +/- 0.01.

## Down to the lower edge

Ten further sub-strips, run with the same window after the lower speed
limit was enclosed by `certify_lower.py` (V_c in [9.3203575, 9.3213142]
m/s over the window, see `../lower/window.txt`). Their lower end,
V0* - 2.578685 = 9.321315 m/s, exceeds every V_c on the window by at
least 8e-7 m/s. Last column: lower bound of the pocket entry speed of
the slowest ball (b18, the paper's ball 19).

| V0 offset | result | pieces | min pocket entry speed (m/s) |
|---|---|---|---|
| [-2.35, -2.3] | CERTIFIED | 23 (of 45 tried) | 0.135431 |
| [-2.4, -2.35] | CERTIFIED | 23 (of 45 tried) | 0.119554 |
| [-2.45, -2.4] | CERTIFIED | 23 (of 45 tried) | 0.101323 |
| [-2.5, -2.45] | CERTIFIED | 23 (of 45 tried) | 0.0791247 |
| [-2.55, -2.5] | CERTIFIED | 23 (of 45 tried) | 0.0477108 |
| [-2.575, -2.55] | CERTIFIED | 23 (of 45 tried) | 0.0170906 |
| [-2.578, -2.575] | CERTIFIED | 23 (of 45 tried) | 0.00737133 |
| [-2.5786, -2.578] | CERTIFIED | 23 (of 45 tried) | 0.00260683 |
| [-2.57868, -2.5786] | CERTIFIED | 23 (of 45 tried) | 0.000676002 |
| [-2.578685, -2.57868] | CERTIFIED | 23 (of 45 tried) | 0.000246617 |
