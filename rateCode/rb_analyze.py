"""h-parametrized geometry of the burster's fast subsystem (SN regime + HCN conductance folded into the leak) -> data/burst_h.pkl
Same dict layout as data/<regime>.pkl so render3d.surfaces/branches work: 'Ilist' holds the h grid."""
import numpy as np, pickle, time
from model import REGIMES, equilibria, separatrix
from analyze import cycles_vectorized
from rb_net import HCN, I_BURST, STRETCH
Vmin, Vmax, nmin, nmax = -90, 30, -0.05, 0.75
base = dict(REGIMES["SN"])            # unstretched: identical orbits, faster convergence
def p_of_h(h):
    gL = base["g_L"] + HCN["g_h"] * h; EL = (base["g_L"] * base["E_L"] + HCN["g_h"] * h * HCN["E_h"]) / gL
    p = dict(base); p["g_L"] = gL; p["E_L"] = EL; return p
H = np.linspace(0.26, 0.54, 71)
t0 = time.time()
eq = [equilibria(p_of_h(h), I_BURST, Vmin, Vmax) for h in H]
# fold: last h with a stable rest state
rest = [any(e[2] == "stable" for e in E) for E in eq]
h_fold = H[max(i for i, r in enumerate(rest) if r)] if any(rest) else None
print("fold (rest disappears) near h =", h_fold, "| rest exists at", sum(rest), "of", len(H), "levels", flush=True)
pv = dict(base); gL = base["g_L"] + HCN["g_h"] * H[:, None]; pv["g_L"] = gL; pv["E_L"] = (base["g_L"] * base["E_L"] + HCN["g_h"] * H[:, None] * HCN["E_h"]) / gL
lc = cycles_vectorized(pv, [I_BURST] * len(H), nmax, settle=800.0, sample=300.0)
# cycles_vectorized indexes by position; our 'I' is constant so it computed per-level with per-level p through broadcasting
for i in lc: lc[i] = (lc[i][0], lc[i][1], lc[i][2] * STRETCH)
print("stable cycles at", len(lc), "levels; h range", (H[min(lc)], H[max(lc)]) if lc else None, flush=True)
sep = [separatrix(p_of_h(h), I_BURST, Vmin, Vmax, nmin, nmax) for h in H]
d = dict(name="burst_h", p=base, Ilist=H, Ibif=float(h_fold), nmax=nmax, eq=eq, lc=lc, ulc={}, sep=sep, h_hc=float(H[min(lc)]) if lc else None)
pickle.dump(d, open("data/burst_h.pkl", "wb"))
print("done in", round(time.time() - t0), "s; homoclinic near h =", d["h_hc"])
