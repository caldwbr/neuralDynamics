"""Precompute, per regime and per current level: equilibria, stable limit cycle, unstable cycle (subH), separatrix. -> data/<regime>.pkl"""
import numpy as np, pickle, os, sys, time
from model import REGIMES, equilibria, rk4, separatrix, ninf

RANGES = {  # I grid for the 3D manifold and the animation, and the bifurcation current
    "SNIC": dict(I=(0.0, 8.0), Ibif=4.6, nmax=0.75),
    "SN":   dict(I=(0.0, 8.0), Ibif=4.6, nmax=0.75),
    "subH": dict(I=(35.0, 65.0), Ibif=48.9, nmax=1.0),
    "supH": dict(I=(0.0, 22.0), Ibif=14.6, nmax=0.75),
}
NI = 121
Vmin, Vmax, nmin = -90, 30, -0.05

def cycles_vectorized(p, Ilist, nmax, dt=0.01, settle=600.0, sample=250.0, backward=False, min_std=0.8):
    """stable LC (forward) or unstable LC (backward) at every I in Ilist at once. Returns dict I-index -> (V, n, period)."""
    seeds_per_I = []
    for I in Ilist:
        eq = equilibria(p, I, Vmin, Vmax)
        if not backward:
            s = [(10, 0.5), (-10, 0.6), (0, 0.4), (-20, 0.7), (10, 0.85), (-50, 0.05), (-70, 0.10)]
            for (Ve, ne, kind, ev) in eq:
                if kind == "unstable": s += [(Ve + 0.2, ne + 0.003), (Ve + 0.8, ne + 0.015), (Ve + 2.0, ne + 0.04)]
            s = s[:10] + [(10, 0.5)] * (10 - len(s[:10]))
        else:
            st = [e for e in eq if e[2] == "stable"]
            if st:
                Vs, ns = st[0][0], st[0][1]
                s = [(Vs + a, ns + b) for a, b in [(0.3, 0.003), (0.8, 0.008), (1.5, 0.02), (3.0, 0.045), (5.0, 0.08)]]
            else:
                s = [(np.nan, np.nan)] * 5
        seeds_per_I.append(s)
    nS = len(seeds_per_I[0])
    V = np.array([[s[0] for s in row] for row in seeds_per_I], float); n = np.array([[s[1] for s in row] for row in seeds_per_I], float)
    Iarr = np.array(Ilist, float)[:, None]
    t0 = time.time()
    V, n, _ = rk4(p, Iarr, V, n, dt, int(settle / dt), backward=backward)
    # clamp runaways so they don't blow up the whole array
    V = np.clip(np.nan_to_num(V, nan=0.0), -200, 200); n = np.clip(np.nan_to_num(n, nan=0.0), -1, 2)
    V, n, rec = rk4(p, Iarr, V, n, dt, int(sample / dt), record_every=4, backward=backward)
    RV = np.array([r[0] for r in rec]); Rn = np.array([r[1] for r in rec])  # (T, nI, nS)
    print(f"   integrated {len(Ilist)} levels x {nS} seeds in {time.time()-t0:.0f}s", flush=True)
    out = {}
    for i in range(len(Ilist)):
        best = None
        for s in range(nS):
            v = RV[:, i, s]; nn = Rn[:, i, s]
            if not np.all(np.isfinite(v)): continue
            if not (np.all(v >= Vmin - 2) and np.all(v <= Vmax + 2) and np.all(nn >= nmin - 0.05) and np.all(nn <= nmax + 0.05)): continue
            if np.std(v) <= min_std: continue
            amp = v.max() - v.min()
            if best is None or amp > best[0]:
                m = v.mean(); up = np.nonzero((v[:-1] < m) & (v[1:] >= m))[0]
                per = np.mean(np.diff(up)) * dt * 4 if len(up) > 2 else np.nan
                best = (amp, v, nn, per)
        if best is not None:
            amp, v, nn, per = best
            if np.isfinite(per):
                npts = int(round(per / (dt * 4))); v = v[-npts:]; nn = nn[-npts:]
            out[i] = (v.astype(np.float32), nn.astype(np.float32), float(per))
    return out

def run(name):
    p = REGIMES[name]; R = RANGES[name]; Ilist = np.linspace(R["I"][0], R["I"][1], NI)
    print(name, flush=True)
    data = dict(name=name, p=p, Ilist=Ilist, Ibif=R["Ibif"], nmax=R["nmax"])
    data["eq"] = [equilibria(p, I, Vmin, Vmax) for I in Ilist]
    data["lc"] = cycles_vectorized(p, Ilist, R["nmax"])
    print(f"   stable cycles at {len(data['lc'])} of {NI} levels", flush=True)
    data["ulc"] = cycles_vectorized(p, Ilist, R["nmax"], settle=400.0, sample=100.0, backward=True, min_std=1.0) if name == "subH" else {}
    if name == "subH": print(f"   unstable cycles at {len(data['ulc'])} levels", flush=True)
    t0 = time.time(); data["sep"] = [separatrix(p, I, Vmin, Vmax, nmin, R["nmax"]) for I in Ilist]
    print(f"   separatrices in {time.time()-t0:.0f}s", flush=True)
    os.makedirs("data", exist_ok=True)
    pickle.dump(data, open(f"data/{name}.pkl", "wb"))

if __name__ == "__main__":
    for name in (sys.argv[1:] or ["SNIC", "SN", "subH", "supH"]): run(name)
