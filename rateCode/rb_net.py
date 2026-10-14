"""Network of I_Na,p + I_K neurons with an optional HCN (h) gate and population-mean synaptic coupling.
Used for the single burster (scene 2) and the RSP -> DMN -> cortex chain (scene 3) of the rate-code / burst-code video.
All times in model units; SCALE_MS converts to real milliseconds for the on-screen clock.
"""
import numpy as np
from model import REGIMES

SCALE_MS = 3.6          # real ms per model time unit (chosen so the burster's intra-burst rate is ~80 Hz and the SNIC cell spans 10-60 Hz)
STRETCH = 2.1           # the burster's fast subsystem: SN regime with C and tau_n both x2.1 (pure time rescaling of the fast dynamics)
HCN = dict(Vh=-60.0, k=7.0, g_h=1.2, E_h=-30.0, tau_h=100.0)   # grid search (RK4, stretch 2.1): ~11 spikes/burst at ~78 Hz, ~2.1 bursts/s, duty 0.26
I_BURST = -12.0         # the hyperpolarized burster's bias current (silent-phase V about -59 mV)
SPIKE_TH = -30.0

def burster_params():
    p = dict(REGIMES["SN"]); p["C"] *= STRETCH; p["tau"] *= STRETCH; return p

class Net:
    """pops: list of dicts with keys name, N, p (regime dict), I (bias), hcn (dict or None), noise (mV/sqrt(unit)), tau_s.
    G[i][j]: synaptic gain from population j onto population i (current = G * S_j(t - D[i][j]) * (E_syn - V)).
    """
    def __init__(self, pops, G=None, D=None, dt=0.01, seed=0, E_syn=0.0):
        self.rng = np.random.default_rng(seed); self.dt = dt; self.t = 0.0; self.E_syn = E_syn
        self.pops = pops; self.K = len(pops); self.N = sum(p["N"] for p in pops)
        self.idx = []; a = 0
        for p in pops: self.idx.append(slice(a, a + p["N"])); a += p["N"]
        def arr(key, default=None):
            out = np.zeros(self.N)
            for p, sl in zip(pops, self.idx): out[sl] = p["p"].get(key, default) if default is not None else p["p"][key]
            return out
        self.C = arr("C"); self.tau = arr("tau"); self.gL = arr("g_L"); self.EL = arr("E_L"); self.gNa = arr("g_Na"); self.gK = arr("g_K")
        self.Vm = arr("V_m_half"); self.km = arr("k_m"); self.Vn = arr("V_n_half"); self.kn = arr("k_n"); self.ENa = arr("E_Na"); self.EK = arr("E_K")
        self.gh = np.zeros(self.N); self.Vh = np.full(self.N, -60.0); self.kh = np.ones(self.N); self.Eh = np.full(self.N, -30.0); self.tauh = np.ones(self.N)
        self.I = np.zeros(self.N); self.noise = np.zeros(self.N); self.tau_s = np.ones(self.N); self.pop_of = np.zeros(self.N, int)
        for k, (p, sl) in enumerate(zip(pops, self.idx)):
            self.I[sl] = p["I"]; self.noise[sl] = p.get("noise", 0.0); self.tau_s[sl] = p.get("tau_s", 1.5); self.pop_of[sl] = k
            if p.get("hcn"):
                h = p["hcn"]; n = p["N"]
                self.gh[sl] = h["g_h"]; self.Vh[sl] = h["Vh"]; self.kh[sl] = h["k"]; self.Eh[sl] = h["E_h"]
                self.tauh[sl] = h["tau_h"] * (1 + p.get("tau_h_jitter", 0.0) * self.rng.uniform(-1, 1, n))
        self.G = np.array(G if G is not None else np.zeros((self.K, self.K)), float)
        self.D = np.array(D if D is not None else np.zeros((self.K, self.K)), float)
        self.maxdelay = int(np.ceil(self.D.max() / dt)) + 2
        self.shist = np.zeros((self.maxdelay, self.N))   # ring buffer of every neuron's s
        self.hp = 0
        # connectivity: W[i, j] weight from neuron j onto neuron i (default: all-to-all within the population blocks of G, weight G/N_j)
        self.W = np.zeros((self.N, self.N))
        for i in range(self.K):
            for j in range(self.K):
                if self.G[i, j] != 0: self.W[self.idx[i], self.idx[j]] = self.G[i, j] / pops[j]["N"]
        self.blocks = [(i, j) for i in range(self.K) for j in range(self.K) if self.G[i, j] != 0]
        # state
        self.V = np.full(self.N, -65.0); self.n = np.zeros(self.N); self.h = np.zeros(self.N); self.s = np.zeros(self.N)
        self.Iext = np.zeros(self.N)   # per-step external pulses (set by caller)
        self.gmul = np.ones(self.N)    # per-neuron multiplier on the synaptic drive (heterogeneity)
        self.spikes = []               # (t, neuron)
        self.above = self.V > SPIKE_TH

    def minf(self, V): return 1 / (1 + np.exp((self.Vm - V) / self.km))
    def ninf(self, V): return 1 / (1 + np.exp((self.Vn - V) / self.kn))
    def hinf(self, V): return 1 / (1 + np.exp((V - self.Vh) / self.kh))
    def syn_current(self, V):
        if not self.blocks: return 0.0
        drive = np.zeros(self.N)
        for (i, j) in self.blocks:
            back = int(round(self.D[i, j] / self.dt)); sj = self.shist[(self.hp - back) % self.maxdelay, self.idx[j]]
            drive[self.idx[i]] += self.W[self.idx[i], self.idx[j]] @ sj
        return self.gmul * drive * (self.E_syn - V)
    def sparse(self, i, j, K, rng, g=None):
        """rewire block (i <- j) so each postsynaptic neuron gets K random presynaptic cells, weight g/K (g defaults to G[i, j])."""
        g = self.G[i, j] if g is None else g; Ni = self.pops[i]["N"]; Nj = self.pops[j]["N"]
        B = np.zeros((Ni, Nj))
        for r in range(Ni): B[r, rng.choice(Nj, K, replace=False)] = g / K
        self.W[self.idx[i], self.idx[j]] = B
    def rhs(self, V, n, h, Isyn):
        Ih = self.gh * h * (self.Eh - V)
        dV = (self.I + self.Iext + Ih + Isyn - self.gL * (V - self.EL) - self.gNa * self.minf(V) * (V - self.ENa) - self.gK * n * (V - self.EK)) / self.C
        dn = (self.ninf(V) - n) / self.tau
        dh = np.where(self.gh > 0, (self.hinf(V) - h) / self.tauh, 0.0)
        return dV, dn, dh
    def step(self, nsteps, record_every=0):
        rec = []
        dt = self.dt
        for k in range(nsteps):
            Isyn = self.syn_current(self.V)
            V, n, h = self.V, self.n, self.h
            k1 = self.rhs(V, n, h, Isyn)
            k2 = self.rhs(V + 0.5 * dt * k1[0], n + 0.5 * dt * k1[1], h + 0.5 * dt * k1[2], Isyn)
            k3 = self.rhs(V + 0.5 * dt * k2[0], n + 0.5 * dt * k2[1], h + 0.5 * dt * k2[2], Isyn)
            k4 = self.rhs(V + dt * k3[0], n + dt * k3[1], h + dt * k3[2], Isyn)
            self.V = V + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) + self.noise * np.sqrt(dt) * self.rng.standard_normal(self.N)
            self.n = n + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            self.h = np.clip(h + dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]), 0, 1)
            self.t += dt
            above = self.V > SPIKE_TH; new = above & ~self.above; self.above = above
            self.s *= np.exp(-dt / self.tau_s)
            if new.any():
                ids = np.nonzero(new)[0]; self.s[ids] += 1.0
                self.spikes += [(self.t, int(i)) for i in ids]
            self.hp = (self.hp + 1) % self.maxdelay; self.shist[self.hp] = self.s
            if record_every and (k % record_every == 0): rec.append((self.t, self.V.copy(), self.n.copy(), self.h.copy()))
        return rec
    def set_rest(self, Vr):
        self.V[:] = Vr; self.n[:] = self.ninf(self.V); self.h[:] = np.where(self.gh > 0, self.hinf(self.V), 0.0)
        self.above = self.V > SPIKE_TH
