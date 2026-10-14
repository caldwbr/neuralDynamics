"""Vectorized burster search at tau_h = 76 (RK4), SN fast subsystem stretched x1.5."""
import numpy as np, itertools
from rb_net import burster_params, SCALE_MS
p = burster_params()
import sys
combos = list(itertools.product((-55, -60, -65, -70), (3, 5, 7), (0.5, 0.8, 1.2, 1.8), (-30, -20), (-3, -6, -9, -12)))
Vh = np.array([c[0] for c in combos], float); kh = np.array([c[1] for c in combos], float); gh = np.array([c[2] for c in combos], float)
Eh = np.array([c[3] for c in combos], float); Ia = np.array([c[4] for c in combos], float); tauh = float(sys.argv[1])
M = len(combos); V = np.full(M, -70.0); n = np.zeros(M); h = np.full(M, 0.4)
def rhs(V, n, h):
    m = 1 / (1 + np.exp((p["V_m_half"] - V) / p["k_m"])); ni = 1 / (1 + np.exp((p["V_n_half"] - V) / p["k_n"])); hi = 1 / (1 + np.exp((V - Vh) / kh))
    dV = (Ia + gh * h * (Eh - V) - p["g_L"] * (V - p["E_L"]) - p["g_Na"] * m * (V - p["E_Na"]) - p["g_K"] * n * (V - p["E_K"])) / p["C"]
    return dV, (ni - n) / p["tau"], (hi - h) / tauh
dt = 0.01; T = 2000.0; N = int(T / dt); above = V > -20; spikes = [[] for _ in range(M)]; hmin = np.full(M, 9.0); hmax = np.zeros(M); vsil = np.zeros(M); nsil = np.zeros(M)
for k in range(N):
    k1 = rhs(V, n, h); k2 = rhs(V + 0.5 * dt * k1[0], n + 0.5 * dt * k1[1], h + 0.5 * dt * k1[2]); k3 = rhs(V + 0.5 * dt * k2[0], n + 0.5 * dt * k2[1], h + 0.5 * dt * k2[2]); k4 = rhs(V + dt * k3[0], n + dt * k3[1], h + dt * k3[2])
    V = V + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); n = n + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]); h = np.clip(h + dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]), 0, 1)
    a = V > -20; new = a & ~above; above = a; t = (k + 1) * dt
    if t > 700:
        for i in np.nonzero(new)[0]: spikes[i].append(t)
        hmin = np.minimum(hmin, h); hmax = np.maximum(hmax, h); sl = V < -50; vsil += np.where(sl, V, 0); nsil += sl
res = []
for i, c in enumerate(combos):
    sp = np.array(spikes[i])
    if len(sp) < 4: res.append((c, "silent" if len(sp) == 0 else f"{len(sp)} spikes")); continue
    isi = np.diff(sp); med = np.median(isi); gaps = isi[isi > 4 * med]; intra = isi[isi <= 4 * med]
    if len(gaps) < 2: res.append((c, f"tonic {len(sp)/1300*1000/SCALE_MS:.0f} Hz")); continue
    nb = len(gaps) + 1; spb = len(sp) / nb; per = np.median(gaps) + np.median(intra) * (spb - 1)
    res.append((c, dict(spb=round(spb, 1), isi=round(float(np.median(intra)), 2), hz=round(1000 / (np.median(intra) * SCALE_MS)), per=round(float(per)), bhz=round(1000 / (per * SCALE_MS), 2), duty=round(1 - float(np.median(gaps)) / per, 2), h=(round(float(hmin[i]), 2), round(float(hmax[i]), 2)), vsil=round(float(vsil[i] / max(nsil[i], 1)), 1))))
good = [r for r in res if isinstance(r[1], dict)]
print(len(good), "bursting combos of", M)
for c, r in sorted(good, key=lambda r: abs(r[1]["hz"] - 80) + 40 * abs(r[1]["duty"] - 0.25) + 10 * abs(r[1]["bhz"] - 2)): print(c, r)
print("--- non bursting summary:", sum(1 for r in res if r[1] == "silent"), "silent;", sum(1 for r in res if isinstance(r[1], str) and r[1].startswith("tonic")), "tonic")
