"""3D manifold panel per regime (equilibrium branches, separatrix surface, limit-cycle manifold, unstable-cycle manifold) with a movable current plane."""
import numpy as np, pickle, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in ("/usr/share/fonts/opentype/inter/Inter-Regular.otf", "/usr/share/fonts/opentype/inter/Inter-SemiBold.otf"): font_manager.fontManager.addfont(f)
plt.rcParams["font.family"] = "Inter"
INK, PAPER, MUTED, GOLD, TEAL, VIOLET, MAGENTA, RED, BLUE = "#0B0F19", "#F4F1EA", "#A9B1C1", "#E3B23C", "#35B8A6", "#9C86E6", "#E6329B", "#E8453C", "#5B8DEF"
Vmin, Vmax, nmin = -90, 30, -0.05

def resample_closed(v, n, k=80):
    """resample a closed curve by angle about its centroid (like the MATLAB code) so rings line up across I."""
    cx, cy = v.mean(), n.mean(); sx = v.max() - v.min() + 1e-9; sy = n.max() - n.min() + 1e-9
    ang = np.arctan2((n - cy) / sy, (v - cx) / sx); o = np.argsort(ang); ang, v, n = ang[o], v[o], n[o]
    ang, u = np.unique(ang, return_index=True); v, n = v[u], n[u]
    q = np.linspace(-np.pi, np.pi, k)
    return np.interp(q, ang, v, period=2 * np.pi), np.interp(q, ang, n, period=2 * np.pi)

def resample_open(v, n, k=60):
    L = np.concatenate([[0], np.cumsum(np.hypot(np.diff(v), np.diff(n)))]); q = np.linspace(0, L[-1], k)
    return np.interp(q, L, v), np.interp(q, L, n)

def surfaces(d):
    """-> list of (X, Y, Z, color, alpha) surfaces + line branches."""
    Il = d["Ilist"]; out = []
    for key, col, al, k in (("lc", MAGENTA, 0.55, 80), ("ulc", RED, 0.45, 80)):
        rows = [(i, resample_closed(d[key][i][0], d[key][i][1], k)) for i in sorted(d[key])]
        # split into contiguous runs of I index
        runs, cur = [], []
        for i, r in rows:
            if cur and i != cur[-1][0] + 1: runs.append(cur); cur = []
            cur.append((i, r))
        if cur: runs.append(cur)
        for run in runs:
            if len(run) < 2: continue
            X = np.array([r[1][0] for r in run]); Y = np.array([r[1][1] for r in run]); Z = np.array([[Il[r[0]]] * k for r in run])
            out.append((X, Y, Z, col, al))
    for side in (0, 1):
        rows = [(i, resample_open(*d["sep"][i][side])) for i in range(len(Il)) if len(d["sep"][i]) > side]
        runs, cur = [], []
        for i, r in rows:
            if cur and i != cur[-1][0] + 1: runs.append(cur); cur = []
            cur.append((i, r))
        if cur: runs.append(cur)
        for run in runs:
            if len(run) < 2: continue
            X = np.array([r[1][0] for r in run]); Y = np.array([r[1][1] for r in run]); Z = np.array([[Il[r[0]]] * 60 for r in run])
            out.append((X, Y, Z, VIOLET, 0.35))
    return out

def branches(d):
    Il = d["Ilist"]; B = {"stable": [], "unstable": [], "saddle": []}
    for i, eqs in enumerate(d["eq"]):
        for (V, n, kind, ev) in eqs: B[kind].append((V, n, Il[i]))
    return B

def render(d, Icur, out_png, size=(680, 570), elev=22, azim=-55, surf=None, br=None):
    Il = d["Ilist"]; nmax = d["nmax"]
    fig = plt.figure(figsize=(size[0] / 100, size[1] / 100), dpi=100, facecolor=INK)
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.96], projection="3d", facecolor=INK)
    ax.set_box_aspect((1.25, 1.0, 1.0))
    surf = surf or surfaces(d); br = br or branches(d)
    for X, Y, Z, col, al in surf:
        ax.plot_surface(X, Y, Z, color=col, alpha=al, linewidth=0, antialiased=False, shade=True)
    for kind, col, ls, lw in (("stable", BLUE, "-", 2.6), ("unstable", RED, ":", 2.6), ("saddle", VIOLET, "--", 2.2)):
        pts = br[kind]
        if pts:
            P = np.array(pts); ax.plot(P[:, 0], P[:, 1], P[:, 2], ls, color=col, lw=lw)
    # bifurcation plane (gold, faint) and the current plane (teal)
    for I, col, al in ((d["Ibif"], GOLD, 0.10), (Icur, TEAL, 0.28)):
        X, Y = np.meshgrid([Vmin, Vmax], [nmin, nmax]); ax.plot_surface(X, Y, np.full_like(X, I), color=col, alpha=al, linewidth=0, shade=False)
        ax.plot([Vmin, Vmax, Vmax, Vmin, Vmin], [nmin, nmin, nmax, nmax, nmin], [I] * 5, color=col, lw=1.2 if col == TEAL else 0.8, ls="-" if col == TEAL else ":")
    ax.set_xlim(Vmin, Vmax); ax.set_ylim(nmin, nmax); ax.set_zlim(Il[0], Il[-1])
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel("V (mV)", color=MUTED, fontsize=9, labelpad=-2); ax.set_ylabel("n", color=MUTED, fontsize=9, labelpad=-2); ax.set_zlabel("I", color=MUTED, fontsize=9, labelpad=-4)
    ax.tick_params(colors=MUTED, labelsize=7, pad=-2)
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane): pane.set_facecolor(INK); pane.set_edgecolor("#2A3347")
    ax.grid(True); ax.xaxis._axinfo["grid"]["color"] = "#1E2638"; ax.yaxis._axinfo["grid"]["color"] = "#1E2638"; ax.zaxis._axinfo["grid"]["color"] = "#1E2638"
    fig.savefig(out_png, dpi=100, facecolor=INK); plt.close(fig)

if __name__ == "__main__":
    name = sys.argv[1]; d = pickle.load(open(f"data/{name}.pkl", "rb"))
    levels = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    os.makedirs(f"p3d/{name}", exist_ok=True)
    surf = surfaces(d); br = branches(d)
    Is = np.linspace(d["Ilist"][0], d["Ilist"][-1], levels) if levels > 1 else [float(sys.argv[3])]
    for k, I in enumerate(Is):
        render(d, I, f"p3d/{name}/L{k:03d}.png", surf=surf, br=br)
        if k % 10 == 0: print(name, k, flush=True)
    print("done", name)
