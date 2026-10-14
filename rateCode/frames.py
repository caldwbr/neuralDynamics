"""Frame renderer for "Four ways a neuron starts to fire".
usage: python3 frames.py <job>   where job in intro, SNIC, SN, subH, supH, outro   (writes frames/<job>/f_%05d.png and frames/<job>.json)
"""
import numpy as np, pickle, json, sys, os, glob, math
import soundfile as sf
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import font_manager
from PIL import Image, ImageDraw, ImageFont
from model import REGIMES, equilibria, vnull, ninf, dV, dn, rk4

for f in glob.glob("/usr/share/fonts/opentype/inter/Inter-*.otf"): font_manager.fontManager.addfont(f)
plt.rcParams["font.family"] = "Inter"
FPS = 30; W, H = 1920, 1080
INK, CARD, PAPER, MUTED, GOLD, TEAL, VIOLET, MAGENTA, RED, BLUE, ORANGE, CYAN = "#0B0F19", "#151C2C", "#F4F1EA", "#A9B1C1", "#E3B23C", "#35B8A6", "#9C86E6", "#E6329B", "#E8453C", "#5B8DEF", "#F2A33A", "#4CC9F0"
Vmin, Vmax, nmin = -90, 30, -0.05
HERE = os.path.dirname(os.path.abspath(__file__))
FONT = lambda w, s: ImageFont.truetype(f"/usr/share/fonts/opentype/inter/Inter-{w}.otf", s)

CFG = {  # per-regime animation settings
    "SNIC": dict(title="Regime 1 of 4 — saddle-node on an invariant circle (SNIC)", short="SNIC", Ilow=2.0, Ihigh=5.2, Iend=8.0, Idown=None, rate=12.0, win=160.0,
                 kicks={"A": [(0.55, 12.0)]}, reset_B=False, labels=True),
    "SN":   dict(title="Regime 2 of 4 — saddle-node off the invariant circle (SN)", short="SN", Ilow=2.5, Ihigh=6.0, Iend=6.0, Idown=3.0, rate=1.6, win=16.0,
                 kicks={"A": [(0.72, 12.0)]}, reset_B=True, labels=True),
    "subH": dict(title="Regime 3 of 4 — subcritical Andronov–Hopf", short="sub-Hopf", Ilow=44.0, Ihigh=55.0, Iend=55.0, Idown=45.0, rate=5.0, win=50.0,
                 kicks={"A": [(0.22, 7.0), (0.84, 34.0)]}, reset_B=True, labels=True),
    "supH": dict(title="Regime 4 of 4 — supercritical Andronov–Hopf", short="super-Hopf", Ilow=8.0, Ihigh=17.0, Iend=22.0, Idown=None, rate=3.0, win=30.0,
                 kicks={"A": [(0.5, 10.0)]}, reset_B=False, labels=True),
}
ORDER = ["SNIC", "SN", "subH", "supH"]
def wav_len(k): a, sr = sf.read(f"{HERE}/tts/{k}.wav"); return len(a) / sr
def nearest_level(d, I):
    Il = d["Ilist"]; return int(np.argmin(np.abs(Il - I)))
def p3d_path(name, I, d):
    k = int(round((I - d["Ilist"][0]) / (d["Ilist"][-1] - d["Ilist"][0]) * 60)); k = max(0, min(60, k))
    return f"{HERE}/p3d/{name}/L{k:03d}.png"

# ---------------------------------------------------------------- simulation
class Sim:
    def __init__(self, p, I0, nmax, rate, seed=0):
        self.p = p; self.nmax = nmax; self.rate = rate; self.rng = np.random.default_rng(seed)
        self.dt = 0.01; self.t = 0.0
        eq = equilibria(p, I0, Vmin, Vmax); st = [e for e in eq if e[2] == "stable"]
        self.rest = (st[0][0], st[0][1]) if st else (-60.0, 0.1)
        self.nV, self.nn = self.rest
        self.M = 14
        self.pV = self.rng.uniform(Vmin + 5, Vmax - 5, self.M); self.pn = self.rng.uniform(nmin + 0.02, nmax - 0.02, self.M)
        self.age = self.rng.uniform(0, 6, self.M)
        self.tails = [[] for _ in range(self.M)]; self.ntail = []
        self.hist_t = []; self.hist_V = []
    def step(self, I, seconds):
        ms = self.rate * seconds; steps = max(1, int(round(ms / self.dt))); dt = ms / steps
        V = np.concatenate([[self.nV], self.pV]); n = np.concatenate([[self.nn], self.pn])
        # a little voltage noise on the neuron so it leaves an unstable focus; record sub-frame positions for the tails
        V, n, rec = rk4(self.p, I, V, n, dt, steps, record_every=max(1, steps // 10))
        V[0] += self.rng.normal(0, 0.25 * math.sqrt(max(ms, 1e-6)))
        V = np.clip(V, Vmin - 5, Vmax + 5); n = np.clip(n, nmin - 0.02, self.nmax + 0.02)
        self.nV, self.nn = float(V[0]), float(n[0]); self.pV, self.pn = V[1:], n[1:]
        self.t += ms; self.hist_t.append(self.t); self.hist_V.append(self.nV)
        self.age += seconds
        for i in range(self.M):
            self.tails[i] += [(float(r[0][i + 1]), float(r[1][i + 1])) for r in rec]; self.tails[i] = self.tails[i][-200:]
            if self.age[i] > 7.5:
                self.age[i] = 0; self.pV[i] = self.rng.uniform(Vmin + 5, Vmax - 5); self.pn[i] = self.rng.uniform(nmin + 0.02, self.nmax - 0.02); self.tails[i] = []
        self.ntail += [(float(r[0][0]), float(r[1][0])) for r in rec]; self.ntail = self.ntail[-400:]
    def kick(self, dv): self.nV += dv; self.ntail = []
    def reset(self, I):
        eq = equilibria(self.p, I, Vmin, Vmax); st = [e for e in eq if e[2] == "stable"]
        if st: self.nV, self.nn = st[0][0], st[0][1]
        self.ntail = []

# ---------------------------------------------------------------- drawing
def draw_plane(ax, p, d, I, sim, nmax, labels, show_field=True, intro_labels=False):
    ax.cla(); ax.set_facecolor(INK)
    ax.set_xlim(Vmin, Vmax); ax.set_ylim(nmin, nmax)
    Vf = np.linspace(Vmin, Vmax, 700)
    if show_field:
        gx, gy = np.meshgrid(np.linspace(Vmin + 3, Vmax - 3, 26), np.linspace(nmin + 0.02, nmax - 0.02, 16))
        u = dV(p, I, gx, gy); v = dn(p, gx, gy)
        # normalise in screen-ish units (V spans 120 mV, n spans ~0.8)
        us = u / 120.0; vs = v / (nmax - nmin); m = np.hypot(us, vs) + 1e-9
        L = 0.022; ex = gx + us / m * L * 120; ey = gy + vs / m * L * (nmax - nmin)
        segs = np.stack([np.stack([gx.ravel(), gy.ravel()], 1), np.stack([ex.ravel(), ey.ravel()], 1)], 1)
        ax.add_collection(LineCollection(segs, colors=[(0.66, 0.69, 0.76, 0.22)], linewidths=1.2))
        ax.scatter(ex.ravel(), ey.ravel(), s=4, color=(0.66, 0.69, 0.76, 0.35), linewidths=0)
    nv = vnull(p, I, Vf); mask = (nv > nmin - 0.2) & (nv < nmax + 0.2)
    ax.plot(Vf[mask], nv[mask], color=TEAL, lw=2.2, label="V-nullcline  (dV/dt = 0)")
    ax.plot(Vf, ninf(p, Vf), color=ORANGE, lw=2.2, label="n-nullcline  (dn/dt = 0)")
    k = nearest_level(d, I)
    for (v, nn) in d["sep"][k]: ax.plot(v, nn, color=VIOLET, lw=1.6, ls="--", alpha=0.9)
    if k in d["ulc"]: v, nn, _ = d["ulc"][k]; ax.plot(np.append(v, v[0]), np.append(nn, nn[0]), color=RED, lw=2.0, ls="--", alpha=0.95)
    if k in d["lc"]: v, nn, _ = d["lc"][k]; ax.plot(np.append(v, v[0]), np.append(nn, nn[0]), color=MAGENTA, lw=2.8, alpha=0.95)
    # particles with fading tails
    for tail in sim.tails:
        if len(tail) > 2:
            T = np.array(tail); segs = np.stack([T[:-1], T[1:]], 1); a = np.linspace(0.05, 0.6, len(segs))
            ax.add_collection(LineCollection(segs, colors=[(0.75, 0.78, 0.85, x) for x in a], linewidths=1.6))
    ax.scatter(sim.pV, sim.pn, s=22, color=(0.85, 0.87, 0.92, 0.9), linewidths=0, zorder=5)
    if len(sim.ntail) > 2:
        T = np.array(sim.ntail); segs = np.stack([T[:-1], T[1:]], 1); a = np.linspace(0.05, 0.9, len(segs))
        ax.add_collection(LineCollection(segs, colors=[(0.89, 0.70, 0.24, x) for x in a], linewidths=3.0))
    ax.scatter([sim.nV], [sim.nn], s=260, color=GOLD, alpha=0.25, linewidths=0, zorder=6)
    ax.scatter([sim.nV], [sim.nn], s=70, color=GOLD, edgecolors=PAPER, linewidths=1.2, zorder=7)
    # equilibria
    for (Ve, ne, kind, ev) in equilibria(p, I, Vmin, Vmax):
        if kind == "stable": ax.scatter([Ve], [ne], s=130, color=BLUE, edgecolors=PAPER, linewidths=1.4, zorder=8)
        elif kind == "saddle": ax.scatter([Ve], [ne], s=140, marker="x", color=VIOLET, linewidths=2.6, zorder=8)
        else: ax.scatter([Ve], [ne], s=130, facecolors="none", edgecolors=RED, linewidths=2.2, zorder=8)
        if labels:
            lab = {"stable": "rest (stable)", "saddle": "saddle — threshold", "unstable": "unstable focus"}[kind]
            if kind == "stable" and abs(ev.imag[0]) > 0: lab = "rest (stable focus)"
            dx, dy, ha, va = {"stable": (-3, -0.045, "right", "top"), "saddle": (4, 0.06, "left", "bottom"), "unstable": (4, 0.03, "left", "bottom")}[kind]
            ax.annotate(lab, (Ve, ne), (Ve + dx, ne + dy), color={"stable": BLUE, "saddle": VIOLET, "unstable": RED}[kind], fontsize=12, ha=ha, va=va)
    if labels and k in d["lc"]:
        v, nn, _ = d["lc"][k]; i = int(np.argmax(nn)); ax.annotate("limit cycle — spiking", (v[i], nn[i]), (v[i], nn[i] + 0.02), color=MAGENTA, fontsize=12, ha="center", va="bottom")
    if labels and k in d["ulc"]:
        v, nn, _ = d["ulc"][k]; i = int(np.argmin(nn)); ax.annotate("unstable cycle — threshold", (v[i], nn[i]), (v[i], nn[i] - 0.02), color=RED, fontsize=12, ha="center", va="top")
    if labels and d["sep"][k]:
        v, nn = d["sep"][k][0]; j = len(v) // 2; ax.annotate("threshold (separatrix)", (v[j], nn[j]), (v[j] + 2, nn[j] + 0.03), color=VIOLET, fontsize=12, alpha=0.9)
    ax.set_xlabel("V  (membrane voltage, mV)", color=MUTED, fontsize=14); ax.set_ylabel("n  (K⁺ gate)", color=MUTED, fontsize=14)
    ax.tick_params(colors=MUTED, labelsize=12)
    for s in ax.spines.values(): s.set_color("#2A3347")
    if intro_labels:
        leg = ax.legend(loc="upper left", fontsize=13, facecolor=CARD, edgecolor="#2A3347", labelcolor=PAPER); leg.get_frame().set_alpha(0.9)

def draw_trace(ax, sim, win):
    ax.cla(); ax.set_facecolor(CARD)
    if sim.hist_t:
        t = np.array(sim.hist_t); v = np.array(sim.hist_V); t0 = t[-1]
        m = t > t0 - win; ax.plot(t[m] - t0, v[m], color=GOLD, lw=2.0)
        ax.set_xlim(-win, 0)
    ax.set_ylim(-95, 35); ax.set_ylabel("V (mV)", color=MUTED, fontsize=12); ax.set_xlabel("time (ms, model units)", color=MUTED, fontsize=11)
    ax.tick_params(colors=MUTED, labelsize=10)
    for s in ax.spines.values(): s.set_color("#2A3347")
    ax.text(0.01, 0.95, "the gold cell's voltage", transform=ax.transAxes, color=MUTED, fontsize=11, va="top")

class Canvas:
    def __init__(self):
        self.fig = plt.figure(figsize=(W / 100, H / 100), dpi=100, facecolor=INK)
        self.ax = self.fig.add_axes([70 / W, 1 - 1000 / H, 1070 / W, 870 / H])
        self.ax2 = self.fig.add_axes([1270 / W, 1 - 380 / H, 610 / W, 240 / H])
    def render(self, p3d, title, I, Irange, Ibif, sub=None):
        self.fig.canvas.draw()
        buf = np.asarray(self.fig.canvas.buffer_rgba())[..., :3]
        im = Image.fromarray(buf.copy()); d = ImageDraw.Draw(im)
        if p3d:
            im.paste(Image.open(p3d).convert("RGB"), (1200, 440))
            d.text((1215, 445), "the stack of phase planes along I", font=FONT("Regular", 15), fill=MUTED)
        d.text((60, 28), title, font=FONT("SemiBold", 34), fill=PAPER)
        if sub: d.text((60, 72), sub, font=FONT("Regular", 20), fill=MUTED)
        # current readout + slider
        x0, x1, y = 1270, 1800, 62
        d.text((x0, 18), f"I = {I:5.2f} pA", font=FONT("SemiBold", 30), fill=GOLD)
        d.rounded_rectangle([x0, y, x1, y + 12], radius=6, fill="#2A3347")
        u = (I - Irange[0]) / (Irange[1] - Irange[0]); u = max(0, min(1, u))
        d.rounded_rectangle([x0, y, x0 + int(u * (x1 - x0)), y + 12], radius=6, fill=GOLD)
        xb = x0 + int((Ibif - Irange[0]) / (Irange[1] - Irange[0]) * (x1 - x0))
        d.line([xb, y - 8, xb, y + 20], fill=PAPER, width=2); d.text((xb, y + 24), f"bifurcation {Ibif:g}", font=FONT("Regular", 14), fill=PAPER, anchor="ma")
        d.text((x0, y + 24), f"{Irange[0]:g}", font=FONT("Regular", 14), fill=MUTED); d.text((x1, y + 24), f"{Irange[1]:g}", font=FONT("Regular", 14), fill=MUTED, anchor="ra")
        return im

def montage(labels_on=True, cur=None):
    """2x2 of the four 3D panels at their bifurcation level"""
    im = Image.new("RGB", (W, H), INK); d = ImageDraw.Draw(im)
    names = {"SNIC": "1 · saddle-node on an invariant circle", "SN": "2 · saddle-node off the circle", "subH": "3 · subcritical Hopf", "supH": "4 · supercritical Hopf"}
    for k, name in enumerate(ORDER):
        dd = pickle.load(open(f"{HERE}/data/{name}.pkl", "rb"))
        pan = Image.open(p3d_path(name, dd["Ilist"][-1], dd)).convert("RGB").resize((800, 670))
        x = 80 + (k % 2) * 960; y = 140 + (k // 2) * 470
        pan = pan.resize((760, 420)); im.paste(pan, (x + 40, y))
        d.text((x + 60, y + 400), names[name], font=FONT("SemiBold", 26), fill=PAPER)
    return im

# ---------------------------------------------------------------- jobs
def job_regime(name):
    cfg = CFG[name]; p = REGIMES[name]; d = pickle.load(open(f"{HERE}/data/{name}.pkl", "rb")); nmax = d["nmax"]
    dA, dB, dC = wav_len(f"{name}_A") + 1.0, wav_len(f"{name}_B") + 1.0, wav_len(f"{name}_C") + 2.0
    T = dA + dB + dC; nF = int(round(T * FPS))
    Irange = (d["Ilist"][0], d["Ilist"][-1])
    def I_of(t):
        if t < dA: return cfg["Ilow"]
        if t < dA + dB: u = (t - dA) / dB; return cfg["Ilow"] + (cfg["Ihigh"] - cfg["Ilow"]) * u
        u = (t - dA - dB) / dC
        if cfg["Idown"] is not None:
            if u < 0.45: return cfg["Ihigh"]
            if u < 0.8: return cfg["Ihigh"] + (cfg["Idown"] - cfg["Ihigh"]) * (u - 0.45) / 0.35
            return cfg["Idown"]
        return cfg["Ihigh"] + (cfg["Iend"] - cfg["Ihigh"]) * min(1, u / 0.85)
    sim = Sim(p, cfg["Ilow"], nmax, cfg["rate"], seed=ORDER.index(name) + 1)
    events = [(dA * f, "kick", dv) for f, dv in cfg["kicks"].get("A", [])]
    if cfg["reset_B"]: events.append((dA + 0.3, "reset", 0))
    if name in ("subH", "supH"): events.append((dA + dB + 0.2, "kick", 4.0))   # a nudge off the now-unstable focus
    events.sort()
    cv = Canvas(); out = f"{HERE}/frames/{name}"; os.makedirs(out, exist_ok=True)
    phase_sub = {0: "low current: the resting state and its neighbours", 1: "the current rises", 2: "past the bifurcation"}
    ei = 0
    for f in range(nF):
        t = f / FPS; I = I_of(t)
        while ei < len(events) and events[ei][0] <= t:
            _, kind, dv = events[ei]; ei += 1
            if kind == "kick": sim.kick(dv)
            else: sim.reset(I)
        sim.step(I, 1 / FPS)
        draw_plane(cv.ax, p, d, I, sim, nmax, cfg["labels"]); draw_trace(cv.ax2, sim, cfg["win"])
        ph = 0 if t < dA else (1 if t < dA + dB else 2)
        im = cv.render(p3d_path(name, I, d), cfg["title"], I, Irange, d["Ibif"], sub=phase_sub[ph])
        im.save(f"{out}/f_{f:05d}.png", compress_level=1)
        if f % 150 == 0: print(name, f, nF, flush=True)
    json.dump(dict(frames=nF, dur=T, audio=[(f"{name}_A", 0.0), (f"{name}_B", dA), (f"{name}_C", dA + dB)]), open(f"{HERE}/frames/{name}.json", "w"))

def job_intro():
    name = "SNIC"; p = REGIMES[name]; d = pickle.load(open(f"{HERE}/data/{name}.pkl", "rb")); nmax = d["nmax"]
    La = wav_len("intro"); T = La + 1.5; nF = int(round(T * FPS)); Irange = (d["Ilist"][0], d["Ilist"][-1])
    tR = 0.55 * La   # "raise it" — lift the current for ~5 s and bring it back
    tM = La - 9.0    # the last nine seconds: the four regimes as a montage
    def I_of(t):
        if t < tR: return 2.0
        u = (t - tR) / 7.0
        if u < 0.5: return 2.0 + 2.3 * (u / 0.5)
        if u < 1.0: return 4.3 - 2.3 * ((u - 0.5) / 0.5)
        return 2.0
    sim = Sim(p, 2.0, nmax, 12.0, seed=9)
    cv = Canvas(); out = f"{HERE}/frames/intro"; os.makedirs(out, exist_ok=True)
    mont = montage()
    for f in range(nF):
        t = f / FPS; I = I_of(t); sim.step(I, 1 / FPS)
        if t < tM:
            draw_plane(cv.ax, p, d, I, sim, nmax, labels=True, intro_labels=True); draw_trace(cv.ax2, sim, 160.0)
            im = cv.render(p3d_path(name, I, d), "Four ways a neuron starts to fire", I, Irange, d["Ibif"], sub="the phase plane of a two-dimensional neuron: voltage V and potassium gate n")
            if t > tM - 1.0:   # dissolve into the montage
                a = (t - (tM - 1.0)) / 1.0; im = Image.blend(im, mont, a)
        else:
            im = mont.copy(); dd = ImageDraw.Draw(im); dd.text((60, 28), "Four ways a neuron starts to fire", font=FONT("SemiBold", 34), fill=PAPER)
            dd.text((60, 72), "four bifurcations, four shapes in the stack of phase planes — one at a time", font=FONT("Regular", 20), fill=MUTED)
        im.save(f"{out}/f_{f:05d}.png", compress_level=1)
        if f % 150 == 0: print("intro", f, nF, flush=True)
    json.dump(dict(frames=nF, dur=T, audio=[("intro", 0.6)]), open(f"{HERE}/frames/intro.json", "w"))

def job_outro():
    La = wav_len("outro"); T = La + 3.0; nF = int(round(T * FPS)); out = f"{HERE}/frames/outro"; os.makedirs(out, exist_ok=True)
    mont = montage()
    rows = [("SNIC", "integrator · rate from zero · no bistability", GOLD), ("SN", "integrator · finite rate · bistable (hysteresis)", GOLD),
            ("subH", "resonator · full-size spikes at once · bistable (hysteresis)", GOLD), ("supH", "resonator · small oscillations that grow · no bistability", GOLD)]
    for f in range(nF):
        t = f / FPS; im = mont.copy(); dd = ImageDraw.Draw(im)
        dd.text((60, 28), "Four doorways from rest to rhythm", font=FONT("SemiBold", 34), fill=PAPER)
        dd.text((60, 72), "Izhikevich's I_Na,p + I_K model · the same two equations, four parameter sets", font=FONT("Regular", 20), fill=MUTED)
        for k, (name, text, col) in enumerate(rows):
            if t > 0.25 * La * (k + 1) / 4 + 2.0:
                x = 80 + (k % 2) * 960; y = 140 + (k // 2) * 470
                dd.text((x + 60, y + 436), text, font=FONT("Regular", 20), fill=col)
        if t > La - 6: dd.text((W - 60, 34), "Brad Caldwell  ·  theoriesofconsciousness.com", font=FONT("Regular", 22), fill=MUTED, anchor="ra")
        im.save(f"{out}/f_{f:05d}.png", compress_level=1)
    json.dump(dict(frames=nF, dur=T, audio=[("outro", 0.6)]), open(f"{HERE}/frames/outro.json", "w"))

if __name__ == "__main__":
    job = sys.argv[1]
    if job == "intro": job_intro()
    elif job == "outro": job_outro()
    else: job_regime(job)
    print("done", job)
