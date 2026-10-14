"""Tuning for the rate-code -> burst-code video: SNIC period curve to high current, and a burster search
with a time-stretched SN fast subsystem (C and tau_n scaled together = pure time rescaling of the fast dynamics)."""
import numpy as np, itertools, sys
from model import REGIMES, minf, ninf, limit_cycle

if sys.argv[1] == "snic":
    p = REGIMES["SNIC"]
    for I in [4.52, 4.54, 4.56, 4.58, 4.6, 4.65, 4.7, 4.8, 5.0, 5.5, 6, 7, 8, 10, 12, 14, 17, 20, 25, 30, 40, 50]:
        r = limit_cycle(p, I, settle=400.0, sample=150.0)
        print(f"I={I:6.2f}  period={r[2]:7.2f}  amp={r[0].max()-r[0].min():5.1f}" if r else f"I={I:6.2f}  no cycle", flush=True)

if sys.argv[1] == "burst":
    s = float(sys.argv[2])   # time stretch of the fast subsystem
    p = dict(REGIMES["SN"]); p["C"] *= s; p["tau"] *= s
    def hinf(V, Vh, k): return 1 / (1 + np.exp((V - Vh) / k))
    def sim(Iapp, gh, Vh, k, Eh, tau, T=1500.0, dt=0.005):
        V, n, h = -70.0, 0.0, 0.2; VV = []; HH = []
        for i in range(int(T / dt)):
            Ih = gh * h * (Eh - V)
            dVdt = (Iapp + Ih - p["g_L"] * (V - p["E_L"]) - p["g_Na"] * minf(p, V) * (V - p["E_Na"]) - p["g_K"] * n * (V - p["E_K"])) / p["C"]
            dndt = (ninf(p, V) - n) / p["tau"]; dhdt = (hinf(V, Vh, k) - h) / tau
            V += dt * dVdt; n += dt * dndt; h += dt * dhdt
            if i % 20 == 0: VV.append(V); HH.append(h)
        return np.array(VV), np.array(HH), dt * 20
    res = []
    for Vh, k, gh, Eh, Iapp, tau in itertools.product((-55, -60, -65), (3, 5), (0.5, 0.8, 1.2), (-30, -20), (-3, -6, -9), (30 * s,)):
        V, h, ds = sim(Iapp, gh, Vh, k, Eh, tau); t = np.arange(len(V)) * ds
        V = V[len(V) // 3:]; h = h[len(h) // 3:]; t = t[len(t) // 3:]
        sp = np.nonzero((V[1:] > -25) & (V[:-1] <= -25))[0]
        if len(sp) < 4: continue
        ts = t[sp]; isi = np.diff(ts); med = np.median(isi); gaps = isi[isi > 4 * med]
        if len(gaps) < 2: continue
        nb = len(gaps) + 1; spb = len(ts) / nb; intra = float(np.median(isi[isi <= 4 * med])); per = float(np.median(gaps)) + intra * (spb - 1)
        duty = intra * (spb - 1) / per
        silentV = np.median(V[V < -50]) if np.any(V < -50) else np.nan
        res.append((Vh, k, gh, Eh, Iapp, round(tau, 1), nb, round(spb, 1), round(intra, 2), round(per, 0), round(duty, 2), round(float(h.min()), 2), round(float(h.max()), 2), round(float(silentV), 1)))
    print("Vh k gh Eh Iapp tau | bursts spk/burst intraISI period duty hmin hmax silentV")
    for r in sorted(res, key=lambda r: r[10]): print(r, flush=True)
