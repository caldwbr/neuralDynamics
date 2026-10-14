"""Izhikevich's I_Na,p + I_K model (persistent sodium + potassium), the four excitability regimes.
State (V, n). dV/dt = (I - gL(V-EL) - gNa m_inf(V)(V-ENa) - gK n (V-EK))/C ; dn/dt = (n_inf(V) - n)/tau
"""
import numpy as np
from scipy.optimize import brentq

REGIMES = {
    "SNIC": dict(C=1, E_L=-80, E_Na=60, E_K=-90, g_L=8, g_Na=20, g_K=10, V_m_half=-20, k_m=15, V_n_half=-25, k_n=5, tau=1.0),
    "SN":   dict(C=1, E_L=-80, E_Na=60, E_K=-90, g_L=8, g_Na=20, g_K=10, V_m_half=-20, k_m=15, V_n_half=-25, k_n=5, tau=0.152),
    "subH": dict(C=1, E_L=-78, E_Na=60, E_K=-90, g_L=1, g_Na=4, g_K=4, V_m_half=-30, k_m=7, V_n_half=-45, k_n=5, tau=1.0),
    "supH": dict(C=1, E_L=-78, E_Na=60, E_K=-90, g_L=8, g_Na=20, g_K=10, V_m_half=-20, k_m=15, V_n_half=-45, k_n=5, tau=1.0),
}

def minf(p, V): return 1.0 / (1.0 + np.exp((p["V_m_half"] - V) / p["k_m"]))
def ninf(p, V): return 1.0 / (1.0 + np.exp((p["V_n_half"] - V) / p["k_n"]))
def dV(p, I, V, n): return (I - p["g_L"] * (V - p["E_L"]) - p["g_Na"] * minf(p, V) * (V - p["E_Na"]) - p["g_K"] * n * (V - p["E_K"])) / p["C"]
def dn(p, V, n): return (ninf(p, V) - n) / p["tau"]
def vnull(p, I, V):
    """n on the V-nullcline"""
    return (I - p["g_L"] * (V - p["E_L"]) - p["g_Na"] * minf(p, V) * (V - p["E_Na"])) / (p["g_K"] * (V - p["E_K"]))
def g_eq(p, I, V): return I - p["g_L"] * (V - p["E_L"]) - p["g_Na"] * minf(p, V) * (V - p["E_Na"]) - p["g_K"] * ninf(p, V) * (V - p["E_K"])

def jacobian(p, I, V, n, h=1e-4):
    return np.array([[(dV(p, I, V + h, n) - dV(p, I, V - h, n)) / (2 * h), (dV(p, I, V, n + h) - dV(p, I, V, n - h)) / (2 * h)],
                     [(dn(p, V + h, n) - dn(p, V - h, n)) / (2 * h), (dn(p, V, n + h) - dn(p, V, n - h)) / (2 * h)]])

def equilibria(p, I, Vmin=-90, Vmax=30):
    """list of (V, n, kind, eigenvalues) ; kind in stable/unstable/saddle (+ 'focus' flag via complex eigs)"""
    Vs = np.linspace(Vmin, Vmax, 3000); g = g_eq(p, I, Vs)
    out = []
    for k in np.nonzero(np.diff(np.sign(g)) != 0)[0]:
        try: Ve = brentq(lambda v: g_eq(p, I, v), Vs[k], Vs[k + 1])
        except ValueError: continue
        ne = ninf(p, Ve); ev = np.linalg.eigvals(jacobian(p, I, Ve, ne)); re = ev.real
        kind = "saddle" if re[0] * re[1] < 0 else ("stable" if np.all(re < 0) else "unstable")
        out.append((Ve, ne, kind, ev))
    return out

def rk4(p, I, V, n, dt, steps, record_every=0, backward=False):
    """vectorized RK4 over arrays V, n (any shape); I scalar or array broadcastable. Returns final (V, n) and optional recording."""
    s = -1.0 if backward else 1.0
    rec = []
    for k in range(steps):
        k1v = s * dV(p, I, V, n); k1n = s * dn(p, V, n)
        k2v = s * dV(p, I, V + 0.5 * dt * k1v, n + 0.5 * dt * k1n); k2n = s * dn(p, V + 0.5 * dt * k1v, n + 0.5 * dt * k1n)
        k3v = s * dV(p, I, V + 0.5 * dt * k2v, n + 0.5 * dt * k2n); k3n = s * dn(p, V + 0.5 * dt * k2v, n + 0.5 * dt * k2n)
        k4v = s * dV(p, I, V + dt * k3v, n + dt * k3n); k4n = s * dn(p, V + dt * k3v, n + dt * k3n)
        V = V + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v); n = n + dt / 6 * (k1n + 2 * k2n + 2 * k3n + k4n)
        if record_every and (k % record_every == 0): rec.append((V.copy(), n.copy()))
    return V, n, rec

def limit_cycle(p, I, Vmin=-90, Vmax=30, nmin=-0.05, nmax=1.0, dt=0.005, settle=600.0, sample=200.0, seeds=None, min_std=0.8):
    """Find the largest stable limit cycle at current I by long integration from several seeds. Returns (Vs, ns, period_ms) or None."""
    base = [(10, 0.5), (-10, 0.6), (0, 0.4), (-20, 0.7), (10, 0.85), (-50, 0.05), (-70, 0.10)]
    for (Ve, ne, kind, ev) in equilibria(p, I, Vmin, Vmax):
        if kind == "unstable": base += [(Ve + 0.2, ne + 0.003), (Ve + 0.8, ne + 0.015), (Ve + 2.0, ne + 0.04)]
    if seeds: base = seeds + base
    V0 = np.array([b[0] for b in base], float); n0 = np.array([b[1] for b in base], float)
    V, n, _ = rk4(p, I, V0, n0, dt, int(settle / dt))
    V, n, rec = rk4(p, I, V, n, dt, int(sample / dt), record_every=4)
    RV = np.array([r[0] for r in rec]); Rn = np.array([r[1] for r in rec])   # (T, seeds)
    best = None
    for s in range(RV.shape[1]):
        v = RV[:, s]; nn = Rn[:, s]
        if not (np.all(v >= Vmin - 2) and np.all(v <= Vmax + 2) and np.all(nn >= nmin - 0.05) and np.all(nn <= nmax + 0.05)): continue
        if np.std(v) <= min_std: continue
        amp = v.max() - v.min()
        if best is None or amp > best[0]:
            # period from upward crossings of the mean
            m = v.mean(); up = np.nonzero((v[:-1] < m) & (v[1:] >= m))[0]
            per = np.mean(np.diff(up)) * dt * 4 if len(up) > 2 else np.nan
            best = (amp, v, nn, per)
    if best is None: return None
    amp, v, nn, per = best
    # cut one period (last) for a clean closed curve
    if np.isfinite(per):
        npts = int(round(per / (dt * 4))); v = v[-npts:]; nn = nn[-npts:]
    return v, nn, per

def unstable_cycle(p, I, Vmin=-90, Vmax=30, nmin=-0.05, nmax=1.0, dt=0.005, settle=400.0, sample=80.0):
    """Unstable (saddle) limit cycle around a stable focus, found by backward-time integration from near the focus."""
    eq = equilibria(p, I, Vmin, Vmax); st = [e for e in eq if e[2] == "stable"]
    if not st: return None
    Vs, ns = st[0][0], st[0][1]
    d = [(0.3, 0.003), (0.8, 0.008), (1.5, 0.02), (3.0, 0.045), (5.0, 0.08)]
    V0 = np.array([Vs + a for a, b in d]); n0 = np.array([ns + b for a, b in d])
    V, n, _ = rk4(p, I, V0, n0, dt, int(settle / dt), backward=True)
    V, n, rec = rk4(p, I, V, n, dt, int(sample / dt), record_every=4, backward=True)
    RV = np.array([r[0] for r in rec]); Rn = np.array([r[1] for r in rec])
    best = None
    for s in range(RV.shape[1]):
        v = RV[:, s]; nn = Rn[:, s]
        if not (np.all(v >= Vmin - 2) and np.all(v <= Vmax + 2) and np.all(nn >= nmin - 0.05) and np.all(nn <= nmax + 0.05)): continue
        if np.std(v) <= 1.0: continue
        amp = v.max() - v.min()
        if best is None or amp > best[0]: best = (amp, v, nn)
    if best is None: return None
    amp, v, nn = best
    m = v.mean(); up = np.nonzero((v[:-1] < m) & (v[1:] >= m))[0]
    if len(up) > 2:
        npts = int(round(np.mean(np.diff(up)))); v = v[-npts:]; nn = nn[-npts:]
    return v, nn

def separatrix(p, I, Vmin=-90, Vmax=30, nmin=-0.05, nmax=1.0, dt=0.002, T=25.0):
    """Stable manifold of the saddle (threshold curve): backward-time integration from the saddle along its stable eigenvector. Returns list of two (V, n) arrays."""
    eq = equilibria(p, I, Vmin, Vmax); sd = [e for e in eq if e[2] == "saddle"]
    if not sd: return []
    Vs, ns = sd[0][0], sd[0][1]
    J = jacobian(p, I, Vs, ns); w, v = np.linalg.eig(J); i = np.argmin(w.real); vs = v[:, i].real; vs = vs / np.linalg.norm(vs)
    out = []
    for sgn in (+1, -1):
        V = np.array([Vs + sgn * 5e-3 * vs[0]]); n = np.array([ns + sgn * 5e-3 * vs[1]])
        _, _, rec = rk4(p, I, V, n, dt, int(T / dt), record_every=5, backward=True)
        v = np.array([r[0][0] for r in rec]); nn = np.array([r[1][0] for r in rec])
        oob = np.nonzero((v < Vmin) | (v > Vmax) | (nn < nmin) | (nn > nmax))[0]
        if len(oob): v = v[: oob[0]]; nn = nn[: oob[0]]
        if len(v) > 5: out.append((v, nn))
    return out

if __name__ == "__main__":
    for name, Ilist in (("SNIC", np.arange(2, 7.01, 0.25)), ("SN", np.arange(2, 7.01, 0.25)), ("subH", np.arange(40, 55.1, 1.0)), ("supH", np.arange(0, 40.1, 2.0))):
        p = REGIMES[name]; print(f"\n{name}")
        for I in Ilist:
            eq = equilibria(p, I)
            desc = ", ".join(f"{k}@{V:.1f}{'(focus)' if abs(ev.imag[0])>0 else ''} re={ev.real[0]:+.3f}" for V, n, k, ev in eq)
            print(f"  I={I:5.2f}: {len(eq)} eq | {desc}")
