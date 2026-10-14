"""Frame renderer for "Rate code, burst code" (RSP/HCN1 video).
usage: python3 rb_frames.py <job>   job in intro, s1, s2, s3, outro  -> rb_frames/<job>/f_%05d.png + rb_frames/<job>.json
"""
import numpy as np, pickle, json, sys, os, math
import soundfile as sf
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import proj3d
from PIL import Image, ImageDraw, ImageFont, ImageFilter
from frames import INK, CARD, PAPER, MUTED, GOLD, TEAL, VIOLET, MAGENTA, RED, BLUE, ORANGE, CYAN, FONT, W, H, FPS
from model import REGIMES, equilibria, vnull, ninf, dV, dn, limit_cycle
from rb_net import Net, burster_params, HCN, I_BURST, SCALE_MS
import render3d

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/rb_frames"
SKIP = int(os.environ.get("SKIP", "1"))   # test mode: draw only every SKIP-th frame (simulation still runs every frame)
Vmin, Vmax, nmin, nmax = -90, 30, -0.05, 0.75
SLOW1 = 10.0; UPS1 = 1000.0 / SCALE_MS / SLOW1          # model units per video second in scenes 1-2 (27.8)
SLOW3 = 4.0;  UPS3 = 1000.0 / SCALE_MS / SLOW3          # scene 3 (69.4)
WIN1 = 1000.0 / SCALE_MS                                 # trace window: 1 s real = 278 units
GREY = "#2A3347"
def wav_len(k): a, sr = sf.read(f"{HERE}/rb_tts/{k}.wav"); return len(a) / sr
def real_s(units): return units * SCALE_MS / 1000.0
def hexrgb(h): return tuple(int(h[i:i + 2], 16) for i in (1, 3, 5))

# ---------------------------------------------------------------- phase plane
def draw_pp(ax, p, I, tail, cell, lc=None, sep=(), eqs=None, labels=(), field=True, ghost=None):
    ax.cla(); ax.set_facecolor(INK); ax.set_xlim(Vmin, Vmax); ax.set_ylim(nmin, nmax)
    Vf = np.linspace(Vmin, Vmax, 700)
    if field:
        gx, gy = np.meshgrid(np.linspace(Vmin + 3, Vmax - 3, 24), np.linspace(nmin + 0.02, nmax - 0.02, 15))
        u = dV(p, I, gx, gy); v = dn(p, gx, gy); us = u / 120.0; vs = v / (nmax - nmin); m = np.hypot(us, vs) + 1e-9
        L = 0.022; ex = gx + us / m * L * 120; ey = gy + vs / m * L * (nmax - nmin)
        segs = np.stack([np.stack([gx.ravel(), gy.ravel()], 1), np.stack([ex.ravel(), ey.ravel()], 1)], 1)
        ax.add_collection(LineCollection(segs, colors=[(0.66, 0.69, 0.76, 0.2)], linewidths=1.1))
        ax.scatter(ex.ravel(), ey.ravel(), s=3.5, color=(0.66, 0.69, 0.76, 0.3), linewidths=0)
    nv = vnull(p, I, Vf); mask = (nv > nmin - 0.2) & (nv < nmax + 0.2)
    ax.plot(Vf[mask], nv[mask], color=TEAL, lw=2.2); ax.plot(Vf, ninf(p, Vf), color=ORANGE, lw=2.2)
    for (v, nn) in sep: ax.plot(v, nn, color=VIOLET, lw=1.6, ls="--", alpha=0.9)
    if lc is not None:
        v, nn = lc[0], lc[1]; ax.plot(np.append(v, v[0]), np.append(nn, nn[0]), color=MAGENTA, lw=2.8, alpha=0.95)
    if len(tail) > 2:
        T = np.array(tail); segs = np.stack([T[:-1], T[1:]], 1); a = np.linspace(0.04, 0.9, len(segs))
        ax.add_collection(LineCollection(segs, colors=[(0.89, 0.70, 0.24, x) for x in a], linewidths=3.0))
    ax.scatter([cell[0]], [cell[1]], s=260, color=GOLD, alpha=0.25, linewidths=0, zorder=6)
    ax.scatter([cell[0]], [cell[1]], s=70, color=GOLD, edgecolors=PAPER, linewidths=1.2, zorder=7)
    eqs = equilibria(p, I, Vmin, Vmax) if eqs is None else eqs
    for (Ve, ne, kind, ev) in eqs:
        if kind == "stable": ax.scatter([Ve], [ne], s=130, color=BLUE, edgecolors=PAPER, linewidths=1.4, zorder=8)
        elif kind == "saddle": ax.scatter([Ve], [ne], s=140, marker="x", color=VIOLET, linewidths=2.6, zorder=8)
        else: ax.scatter([Ve], [ne], s=130, facecolors="none", edgecolors=RED, linewidths=2.2, zorder=8)
    if ghost is not None: ax.scatter([ghost[0]], [ghost[1]], s=170, facecolors="none", edgecolors=MUTED, linewidths=1.6, ls=":", zorder=8)
    for (x, y, text, col, ha, va) in labels: ax.annotate(text, (x, y), color=col, fontsize=12.5, ha=ha, va=va)
    ax.set_xlabel("V  (membrane voltage, mV)", color=MUTED, fontsize=14); ax.set_ylabel("n  (K⁺ gate)", color=MUTED, fontsize=14)
    ax.tick_params(colors=MUTED, labelsize=12)
    for s in ax.spines.values(): s.set_color(GREY)

def style_axes(ax, ylabel, xlabel=None):
    ax.set_facecolor(CARD); ax.set_ylabel(ylabel, color=MUTED, fontsize=12); ax.tick_params(colors=MUTED, labelsize=10)
    if xlabel: ax.set_xlabel(xlabel, color=MUTED, fontsize=11)
    for s in ax.spines.values(): s.set_color(GREY)

def draw_vtrace(ax, hist_t, hist_V, spikes, t_now, win, label, marks=(), ylim=(-95, 35)):
    """hist_t in units; shown in real seconds relative to now"""
    ax.cla(); style_axes(ax, "V (mV)")
    k = 1000.0 / SCALE_MS
    if len(hist_t):
        t = np.array(hist_t); v = np.array(hist_V); m = t > t_now - win
        ax.plot((t[m] - t_now) / k, v[m], color=GOLD, lw=1.6)
    sp = np.array([s for s in spikes if s > t_now - win]);
    if len(sp): ax.vlines((sp - t_now) / k, ylim[1] - 9, ylim[1] - 1, color=PAPER, lw=1.2)
    for (te, text, col) in marks:
        if t_now - win < te <= t_now:
            ax.axvline((te - t_now) / k, color=col, lw=1.4, ls="--", alpha=0.9); ax.text((te - t_now) / k, ylim[0] + 4, text, color=col, fontsize=10.5, ha="left", va="bottom", rotation=90)
    ax.set_xlim(-win / k, 0); ax.set_ylim(*ylim)
    ax.text(0.01, 0.04, label, transform=ax.transAxes, color=MUTED, fontsize=11, va="bottom")

def chrome(im, title, sub, footer, clock_units=None, slow=None):
    d = ImageDraw.Draw(im)
    d.text((60, 28), title, font=FONT("SemiBold", 34), fill=PAPER)
    if sub: d.text((60, 72), sub, font=FONT("Regular", 20), fill=MUTED)
    if footer: d.text((W - 60, H - 26), footer, font=FONT("Regular", 17), fill=MUTED, anchor="rd")
    if clock_units is not None:
        d.text((W - 60, 30), f"t = {real_s(clock_units):6.2f} s real", font=FONT("SemiBold", 26), fill=PAPER, anchor="ra")
        if slow: d.text((W - 60, 66), f"slowed {slow:g}×  ·  1 model unit = {SCALE_MS} ms", font=FONT("Regular", 17), fill=MUTED, anchor="ra")
    return d

def fig_to_image(fig):
    fig.canvas.draw(); return Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy())

def spike_rate_hz(spikes, t_now):
    """instantaneous rate from the last ISI (real Hz), 0 if the cell has been quiet longer than 1.6 ISI or > 0.3 s"""
    if len(spikes) < 2: return 0.0
    isi = spikes[-1] - spikes[-2]; since = t_now - spikes[-1]
    if since > max(1.6 * isi, 8.0) or since > 0.3 * WIN1: return 0.0
    return 1000.0 / (isi * SCALE_MS)

# ================================================================ scene 1: rate code
I_BASE, I_WEAK, I_STRONG = 4.63, 6.7, 22.0
def job_s1():
    p = REGIMES["SNIC"]
    dA, dB, dC = wav_len("s1_A") + 1.0, wav_len("s1_B") + 1.0, wav_len("s1_C") + 2.0
    T = dA + dB + dC; nF = int(round(T * FPS))
    # events (video s): the strong step 3.0 s into B for 3.3 s (= 0.33 s real), the weak step at B + 14.5
    tS0, tS1 = dA + 3.0, dA + 6.3; tW0, tW1 = dA + 15.0, dA + 18.3
    def I_of(t):
        if tS0 <= t < tS1: return I_STRONG
        if tW0 <= t < tW1: return I_WEAK
        return I_BASE
    def status(t):
        if tS0 <= t < tS1: return ("EVENT", RED)
        if tW0 <= t < tW1: return ("weaker event", ORANGE)
        return ("baseline", MUTED)
    lcs = {}
    for I in (I_BASE, I_WEAK, I_STRONG):
        r = limit_cycle(p, I, settle=400.0, sample=120.0); lcs[I] = (r[0], r[1], r[2]); print(f"LC at I={I}: period {r[2]:.2f} units = {r[2]*SCALE_MS:.1f} ms -> {1000/(r[2]*SCALE_MS):.0f} Hz", flush=True)
    # ghost of the saddle-node at baseline: the slowest point of the baseline cycle
    v, nn, _ = lcs[I_BASE]; sp = np.hypot(np.diff(v) / 120, np.diff(nn) / 0.8); gi = int(np.argmin(sp)); ghost = (v[gi], nn[gi])
    net = Net([dict(name="c", N=1, p=p, I=I_BASE, noise=0.25)], dt=0.01, seed=11); net.set_rest(-60.0); net.V[:] = ghost[0] + 3; net.n[:] = ghost[1]
    steps = int(round(UPS1 / FPS / net.dt))
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100, facecolor=INK)
    ax = fig.add_axes([60 / W, 1 - 1000 / H, 1000 / W, 880 / H]); axV = fig.add_axes([1130 / W, 1 - 400 / H, 740 / W, 280 / H]); axI = fig.add_axes([1130 / W, 1 - 560 / H, 740 / W, 100 / H])
    tail = []; hist_t = []; hist_V = []; hist_I = []; rate_s = 0.0
    os.makedirs(f"{OUT}/s1", exist_ok=True)
    for f in range(nF):
        t = f / FPS; I = I_of(t); net.I[:] = I
        rec = net.step(steps, record_every=10)
        for (tt, V, n, h) in rec: tail.append((float(V[0]), float(n[0]))); hist_t.append(tt); hist_V.append(float(V[0])); hist_I.append(I)
        tail = tail[-320:]
        spikes = [s for s, i in net.spikes]
        r = spike_rate_hz(spikes, net.t); rate_s += (r - rate_s) * (0.5 if r > rate_s else 0.12)
        if f % SKIP: continue
        labels = [(ghost[0] + 7, ghost[1] + 0.05, "slow passage — the ghost of the saddle-node", MUTED, "left", "bottom")] if I == I_BASE else []
        vv, nn, _ = lcs[I]; j = int(np.argmax(nn)); labels.append((vv[j], nn[j] + 0.02, "spiking loop — one lap = one spike", MAGENTA, "center", "bottom"))
        draw_pp(ax, p, I, tail, (net.V[0], net.n[0]), lc=lcs[I], labels=labels, ghost=ghost if I == I_BASE else None)
        draw_vtrace(axV, hist_t, hist_V, spikes, net.t, WIN1, "the cell's voltage — ticks above: spikes")
        axI.cla(); style_axes(axI, "input I", "time (s, real) — the last second")
        tt = np.array(hist_t); m = tt > net.t - WIN1; k = 1000.0 / SCALE_MS
        axI.plot((tt[m] - net.t) / k, np.array(hist_I)[m], color=TEAL, lw=2.0); axI.set_xlim(-WIN1 / k, 0); axI.set_ylim(0, 26); axI.set_yticks([5, 15, 25])
        im = fig_to_image(fig)
        d = chrome(im, "1 · The rate code — an event is a change of firing rate", "a class-1 (SNIC) neuron: the input current steps up for a third of a second, the rate follows", "Izhikevich I_Na,p + I_K model · Brad Caldwell", clock_units=net.t, slow=SLOW1)
        # rate meter
        x0, y0 = 1130, 620; st, col = status(t)
        d.text((x0, y0), "firing rate", font=FONT("Regular", 22), fill=MUTED)
        d.text((x0, y0 + 30), f"{rate_s:4.0f}", font=FONT("SemiBold", 120), fill=GOLD); d.text((x0 + 300, y0 + 102), "Hz", font=FONT("SemiBold", 44), fill=GOLD)
        bx0, bx1, by = x0, x0 + 740, y0 + 190
        d.rounded_rectangle([bx0, by, bx1, by + 22], radius=11, fill=GREY)
        u = max(0, min(1, rate_s / 100.0)); d.rounded_rectangle([bx0, by, bx0 + int(u * (bx1 - bx0)), by + 22], radius=11, fill=GOLD if rate_s < 45 else RED)
        for hz in range(0, 101, 20):
            xx = bx0 + int(hz / 100 * (bx1 - bx0)); d.line([xx, by + 26, xx, by + 34], fill=MUTED, width=2); d.text((xx, by + 38), f"{hz}", font=FONT("Regular", 15), fill=MUTED, anchor="ma")
        d.text((x0, by + 70), st, font=FONT("SemiBold", 30), fill=col)
        d.text((x0, by + 112), f"I = {I:4.1f}  ·  {'the loop is tight and fast' if I > 10 else ('a smaller step, a middling rate' if I > 5 else 'the trajectory crawls past the ghost')}", font=FONT("Regular", 18), fill=MUTED)
        im.save(f"{OUT}/s1/f_{f:05d}.png", compress_level=1)
        if f % 300 == 0: print("s1", f, nF, flush=True)
    json.dump(dict(frames=nF, dur=T, audio=[("s1_A", 0.0), ("s1_B", dA), ("s1_C", dA + dB)]), open(f"{OUT}/s1.json", "w"))

# ================================================================ scene 2: burst code
def p_of_h(h):
    base = burster_params(); gL = base["g_L"] + HCN["g_h"] * h; EL = (base["g_L"] * base["E_L"] + HCN["g_h"] * h * HCN["E_h"]) / gL
    p = dict(base); p["g_L"] = gL; p["E_L"] = EL; return p

class Proj3D:
    def __init__(self, d, size=(740, 460)):
        fig = plt.figure(figsize=(size[0] / 100, size[1] / 100), dpi=100)
        ax = fig.add_axes([0.02, 0.02, 0.96, 0.96], projection="3d"); ax.set_box_aspect((1.25, 1.0, 1.0))
        ax.set_xlim(Vmin, Vmax); ax.set_ylim(nmin, d["nmax"]); ax.set_zlim(d["Ilist"][0], d["Ilist"][-1]); ax.view_init(elev=22, azim=-55)
        fig.canvas.draw(); self.M = ax.get_proj(); self.ax = ax; self.size = size; self.fig = fig
    def __call__(self, x, y, z):
        xs, ys, _ = proj3d.proj_transform(np.asarray(x, float), np.asarray(y, float), np.asarray(z, float), self.M)
        P = self.ax.transData.transform(np.c_[xs, ys]); return P[:, 0], self.size[1] - P[:, 1]

def job_s2():
    d = pickle.load(open(f"{HERE}/data/burst_h.pkl", "rb")); Hs = d["Ilist"]; h_fold, h_hc = d["Ibif"], d["h_hc"]
    proj = Proj3D(d)
    keys = ["s2_A", "s2_B", "s2_C", "s2_D", "s2_E"]; durs = [wav_len(k) + 1.0 for k in keys]; durs[-1] += 2.0
    starts = np.concatenate([[0], np.cumsum(durs)[:-1]]); T = float(np.sum(durs)); nF = int(round(T * FPS))
    tB, tD = starts[1], starts[3]
    p = burster_params()
    net = Net([dict(name="b", N=1, p=p, I=I_BURST, hcn=HCN, noise=0.0)], dt=0.01, seed=2); net.set_rest(-80.0); net.h[:] = 0.0
    gh_full = net.gh.copy(); net.gh[:] = 0.0            # the HCN gate is switched on at the start of B
    steps = int(round(UPS1 / FPS / net.dt))
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100, facecolor=INK)
    ax = fig.add_axes([60 / W, 1 - 1000 / H, 900 / W, 880 / H])
    axV = fig.add_axes([1040 / W, 1 - 350 / H, 830 / W, 230 / H]); axH = fig.add_axes([1040 / W, 1 - 520 / H, 830 / W, 130 / H])
    tail = []; hist_t = []; hist_V = []; hist_h = []; traj3 = []; marks = []
    # event scheduling in D: (fraction of the silent phase or 'mid'/'early', amplitude, label)
    plan = [(0.7, 12.0, "strong shove", RED, 2.5), (0.7, 5.0, "weak shove", ORANGE, 7.5), ("mid", 12.0, "shove mid-burst", VIOLET, 13.5), (0.15, 12.0, "shove too early", CYAN, 18.0)]
    plan_i = 0; last_gap_used = None; SIL = 96.5; GAP = 15.0
    os.makedirs(f"{OUT}/s2", exist_ok=True)
    for f in range(nF):
        t = f / FPS
        if t >= tB and net.gh[0] == 0: net.gh[:] = gh_full
        # events: instantaneous voltage kicks ("shoves"), placed by burst phase and not before a narration cue
        spikes = [s for s, i in net.spikes]
        if plan_i < len(plan) and spikes and t >= tD + plan[plan_i][4]:
            rule, dv, lab, col, _ = plan[plan_i]
            last = spikes[-1]; gap = net.t - last
            if rule == "mid":
                onset = None
                for i in range(len(spikes) - 1, 0, -1):
                    if spikes[i] - spikes[i - 1] > GAP: onset = spikes[i]; break
                if onset is not None and gap < 6.0 and net.t >= onset + 18.0 and (last_gap_used is None or onset > last_gap_used):
                    net.V[:] += dv; marks.append((net.t, lab, col)); plan_i += 1; last_gap_used = onset
            elif gap > GAP and net.t >= last + rule * SIL and (last_gap_used is None or last > last_gap_used):
                net.V[:] += dv; marks.append((net.t, lab, col)); plan_i += 1; last_gap_used = last
        rec = net.step(steps, record_every=10)
        for (tt, V, n, h) in rec:
            tail.append((float(V[0]), float(n[0]))); hist_t.append(tt); hist_V.append(float(V[0])); hist_h.append(float(h[0])); traj3.append((float(V[0]), float(n[0]), float(h[0])))
        tail = tail[-200:]; traj3 = traj3[-int(150 / (steps * net.dt / 10)) :]     # about one burst cycle of 3-D trail
        spikes = [s for s, i in net.spikes]
        if f % SKIP: continue
        h = float(net.h[0]); ph = p_of_h(h); k = int(np.argmin(np.abs(Hs - h)))
        lc = d["lc"].get(k) if h >= h_hc - 0.003 else None
        eqs = equilibria(ph, I_BURST, Vmin, Vmax)
        labels = []
        for (Ve, ne, kind, ev) in eqs:
            if kind == "stable": labels.append((Ve - 3, ne - 0.045, "rest", BLUE, "right", "top"))
            if kind == "saddle": labels.append((Ve + 4, ne + 0.05, "saddle — threshold", VIOLET, "left", "bottom"))
        if lc is not None: j = int(np.argmax(lc[1])); labels.append((lc[0][j], lc[1][j] + 0.02, "spiking loop (70–80 Hz)", MAGENTA, "center", "bottom"))
        draw_pp(ax, ph, I_BURST, tail, (net.V[0], net.n[0]), lc=lc, sep=d["sep"][k], eqs=eqs, labels=labels)
        draw_vtrace(axV, hist_t, hist_V, spikes, net.t, WIN1, "voltage — ticks: spikes", marks=marks)
        axH.cla(); style_axes(axH, "h (HCN1 gate)", "time (s, real) — the last second")
        tt = np.array(hist_t); m = tt > net.t - WIN1; kk = 1000.0 / SCALE_MS
        if m.any(): axH.plot((tt[m] - net.t) / kk, np.array(hist_h)[m], color=CYAN, lw=2.0)
        axH.axhline(h_fold, color=GOLD, lw=1.0, ls="--", alpha=0.8); axH.axhline(h_hc, color=RED, lw=1.0, ls="--", alpha=0.8)
        axH.text(-WIN1 / kk + 0.01, h_fold + 0.006, "fold — burst starts", color=GOLD, fontsize=9.5, va="bottom"); axH.text(-WIN1 / kk + 0.01, h_hc - 0.006, "homoclinic — burst ends", color=RED, fontsize=9.5, va="top")
        axH.set_xlim(-WIN1 / kk, 0); axH.set_ylim(-0.02 if t < tB else 0.30, 0.52)
        Ih = HCN["g_h"] * h * (HCN["E_h"] - float(net.V[0])) if net.gh[0] > 0 else 0.0
        axH.text(0.99, 0.93, f"h = {h:5.3f}    I_h = {Ih:4.1f}    I + I_h = {I_BURST + Ih:5.1f}" + ("    (HCN1 not yet open)" if net.gh[0] == 0 else ""), transform=axH.transAxes, color=CYAN, fontsize=11, ha="right", va="top")
        im = fig_to_image(fig)
        ph_name = ["hyperpolarize it: the resting input is now below threshold", "the HCN1 gate h opens in the dark and closes in the light: a slow loop", "one closed loop in three dimensions", "events through a sedated synapse: shoves, not pushes", "a switch and a clock, not a dial"][int(np.searchsorted(starts, t, side="right") - 1)]
        dr = chrome(im, "2 · The burst code — hyperpolarization + HCN1 turn the dial into a switch", ph_name, "SN fast subsystem (bistable) + HCN1 gate h · I_h = g_h · h · (E_h − V) · Brad Caldwell", clock_units=net.t, slow=SLOW1)
        # 3-D panel with the trajectory
        pan = Image.open(f"{HERE}/p3d/burst_h/L{k:03d}.png").convert("RGBA")
        if len(traj3) > 2:
            P = np.array(traj3); px, py = proj(P[:, 0], P[:, 1], P[:, 2])
            ov = Image.new("RGBA", pan.size, (0, 0, 0, 0)); od = ImageDraw.Draw(ov); n = len(px)
            for i in range(1, n):
                a = int(255 * (0.08 + 0.92 * i / n)); od.line([(px[i - 1], py[i - 1]), (px[i], py[i])], fill=(227, 178, 60, a), width=2)
            od.ellipse([px[-1] - 6, py[-1] - 6, px[-1] + 6, py[-1] + 6], fill=(227, 178, 60, 255), outline=(244, 241, 234, 255))
            pan = Image.alpha_composite(pan, ov)
        im.paste(pan.convert("RGB"), (1080, 560))
        dr.text((1090, 566), "the stack of fast phase planes along h — the burst is one loop", font=FONT("Regular", 15), fill=MUTED)
        im.save(f"{OUT}/s2/f_{f:05d}.png", compress_level=1)
        if f % 300 == 0: print("s2", f, nF, "events", len(marks), flush=True)
    print("event marks:", [(round(m[0], 1), m[1]) for m in marks])
    json.dump(dict(frames=nF, dur=T, audio=[(k, float(s)) for k, s in zip(keys, starts)]), open(f"{OUT}/s2.json", "w"))

# ================================================================ scene 3: the spread
RING_LABELS_DMN = ["posterior cingulate", "precuneus", "medial prefrontal", "hippocampal formation", "angular gyrus", "lateral temporal"]
RING_LABELS_CTX = ["visual", "auditory", "somatosensory", "motor", "parietal", "prefrontal", "insula", "temporal"]
def job_s3():
    from rb_pop import simulate, N1, N2, N3, UPS
    keys = ["s3_A", "s3_B", "s3_C", "s3_D"]; durs = [wav_len(k) + 1.0 for k in keys]; durs[-1] += 2.5
    starts = np.concatenate([[0], np.cumsum(durs)[:-1]]); T = float(np.sum(durs)); nF = int(round(T * FPS))
    BURN = 2.0; t_couple = starts[1] + 3.0; t_sed = starts[2] + 1.5; sed_dur = 10.0
    cache = f"{HERE}/data/pop_scene.pkl"
    if os.path.exists(cache) and pickle.load(open(cache, "rb")).get("T_video") == T + BURN:
        pop = pickle.load(open(cache, "rb")); print("loaded cached population run", flush=True)
    else:
        pop = simulate(T + BURN, sed_start=t_sed + BURN, sed_dur=sed_dur, couple_at=t_couple + BURN, g_rsp=0.045, g_12=0.3, g_23=0.5)
        pickle.dump(pop, open(cache, "wb"))
    sp = pop["spikes"]; sp_t = sp[:, 0] / UPS - BURN; sp_i = sp[:, 1].astype(int)   # video seconds
    N = N1 + N2 + N3; Tm = pop["T"] / UPS - BURN
    # smoothed population activity (EEG-ish): exponential smoothing of mean s, tau = 4 units
    S = pop["S"]; dtm = np.median(np.diff(pop["T"])); a = dtm / 4.0; E = np.zeros_like(S)
    for i in range(1, len(S)): E[i] = E[i - 1] + a * (S[i] - E[i - 1])
    # per-neuron visual glow: decaying with tau 4 units, from the spike list, sampled per frame
    order = np.argsort(sp_t); sp_t, sp_i = sp_t[order], sp_i[order]
    glow = np.zeros(N); ptr = 0; decay = math.exp(-(UPS / FPS) / 5.0)
    # ring geometry
    cx, cy = 520, 585; rng = np.random.default_rng(7)
    pos = np.zeros((N, 2)); cols = [GOLD] * N1 + [TEAL] * N2 + [VIOLET] * N3
    for i in range(N1):
        r, th = (42 if i < 10 else 88), (i / 10 if i < 10 else (i - 10) / 20) * 2 * math.pi + (0 if i < 10 else 0.15); pos[i] = (cx + r * math.cos(th), cy + r * math.sin(th))
    for i in range(N2): th = i / N2 * 2 * math.pi - math.pi / 2; pos[N1 + i] = (cx + 228 * math.cos(th), cy + 228 * math.sin(th))
    for i in range(N3): th = i / N3 * 2 * math.pi - math.pi / 2 + 0.05; pos[N1 + N2 + i] = (cx + 392 * math.cos(th), cy + 392 * math.sin(th))
    spokes = [(int(rng.integers(0, N1)), N1 + int(rng.integers(0, N2))) for _ in range(70)] + [(N1 + int(rng.integers(0, N2)), N1 + N2 + int(rng.integers(0, N3))) for _ in range(110)]
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100, facecolor=INK)
    axR = fig.add_axes([1090 / W, 1 - 640 / H, 780 / W, 520 / H]); axE = fig.add_axes([1090 / W, 1 - 980 / H, 780 / W, 250 / H])
    WINV = 3.0   # video seconds shown in the raster (= 0.75 s real)
    os.makedirs(f"{OUT}/s3", exist_ok=True)
    for f in range(nF):
        t = f / FPS
        glow *= decay
        while ptr < len(sp_t) and sp_t[ptr] <= t: glow[sp_i[ptr]] += 1.0; ptr += 1
        if f % SKIP: continue
        # raster
        axR.cla(); axR.set_facecolor(CARD)
        m = (sp_t > t - WINV) & (sp_t <= t)
        if m.any(): axR.scatter((sp_t[m] - t) / SLOW3, sp_i[m], s=4, c=[cols[i] for i in sp_i[m]], linewidths=0)
        for y in (N1 - 0.5, N1 + N2 - 0.5): axR.axhline(y, color=GREY, lw=1.0)
        axR.set_xlim(-WINV / SLOW3, 0); axR.set_ylim(-0.5, N - 0.5); axR.set_yticks([N1 / 2, N1 + N2 / 2, N1 + N2 + N3 / 2]); axR.set_yticklabels(["RSP L5\n(30)", "DMN\n(40)", "cortex\n(60)"])
        axR.tick_params(colors=MUTED, labelsize=10.5); axR.set_xlabel("time (s, real) — the last ¾ second", color=MUTED, fontsize=11)
        for s in axR.spines.values(): s.set_color(GREY)
        for (te, lab, col) in ((t_couple, "coupling on", GOLD), (t_sed, "sedation deepens", CYAN), (t_sed + sed_dur, "", CYAN)):
            if t - WINV < te <= t:
                axR.axvline((te - t) / SLOW3, color=col, lw=1.2, ls="--"); axR.text((te - t) / SLOW3 + 0.005, N - 2, lab, color=col, fontsize=10, va="top")
        axR.text(0.01, 0.985, "spike raster — one row per neuron", transform=axR.transAxes, color=MUTED, fontsize=11, va="top")
        # EEG-ish traces
        axE.cla(); axE.set_facecolor(CARD)
        m2 = (Tm > t - WINV) & (Tm <= t)
        if m2.any():
            x = (Tm[m2] - t) / SLOW3
            for kpop, (col, off, lab) in enumerate(((GOLD, 2.0, "RSP"), (TEAL, 1.0, "DMN"), (VIOLET, 0.0, "cortex"))):
                axE.plot(x, off + E[m2, kpop] * 1.6, color=col, lw=1.6); axE.text(-WINV / SLOW3 + 0.01, off + 0.42, lab, color=col, fontsize=10.5, va="top")
            axE.plot(x, -1.1 + (E[m2, 0] * 0.3 + E[m2, 1] * 0.5 + E[m2, 2] * 1.0) * 1.3, color=PAPER, lw=1.8); axE.text(-WINV / SLOW3 + 0.01, -0.62, "sum (the EEG sees this)", color=PAPER, fontsize=10.5, va="top")
        axE.set_xlim(-WINV / SLOW3, 0); axE.set_ylim(-1.3, 3.0); axE.set_yticks([]); axE.tick_params(colors=MUTED, labelsize=10)
        axE.set_xlabel("time (s, real)", color=MUTED, fontsize=11)
        for s in axE.spines.values(): s.set_color(GREY)
        axE.text(0.01, 0.985, "population activity (smoothed synaptic output)", transform=axE.transAxes, color=MUTED, fontsize=11, va="top")
        im = fig_to_image(fig)
        ph = int(np.searchsorted(starts, t, side="right") - 1)
        sub = ["the retrosplenial hub: thirty layer-5 cells with HCN1, each on its own clock", "weak mutual excitation: the thirty clocks become one volley", "the next stages, hyperpolarized and below threshold, can only answer a volley with a volley", "one clock, outward from the hub: on, off, on, off"][ph]
        dr = chrome(im, "3 · The spread — from retrosplenial cortex outward", sub, "3 populations · sparse random synapses, 20–30 ms delays · Vesuna et al. 2020 for the RSP rhythm · Brad Caldwell", clock_units=(t + 0) * UPS, slow=SLOW3)
        # ring schematic on an RGBA layer
        lay = Image.new("RGBA", (W, H), (0, 0, 0, 0)); ld = ImageDraw.Draw(lay)
        for rr, col in ((228, TEAL), (392, VIOLET)): ld.ellipse([cx - rr, cy - rr, cx + rr, cy + rr], outline=hexrgb(col) + (60,), width=2)
        ld.ellipse([cx - 112, cy - 112, cx + 112, cy + 112], outline=hexrgb(GOLD) + (70,), width=2)
        coupled = t >= t_couple
        for (a_, b_) in spokes:
            g = min(1.0, glow[a_] / 1.5); al = int(28 + 150 * g)
            ld.line([tuple(pos[a_]), tuple(pos[b_])], fill=hexrgb(cols[a_]) + (al,), width=1 + (1 if g > 0.4 else 0))
        for i in range(N):
            g = min(1.0, glow[i] / 1.5); base = hexrgb(cols[i]); c = tuple(int(b + (255 - b) * 0.75 * g) for b in base)
            r = 5 + 5 * g; x, y = pos[i]
            if g > 0.25: ld.ellipse([x - r * 2.2, y - r * 2.2, x + r * 2.2, y + r * 2.2], fill=base + (int(90 * g),))
            ld.ellipse([x - r, y - r, x + r, y + r], fill=c + (255,))
        im = Image.alpha_composite(im.convert("RGBA"), lay).convert("RGB"); dr = ImageDraw.Draw(im)
        dr.text((cx, cy + 128), "retrosplenial L5 · HCN1 bursters", font=FONT("SemiBold", 17), fill=GOLD, anchor="ma")
        for j, lab in enumerate(RING_LABELS_DMN):
            th = j / len(RING_LABELS_DMN) * 2 * math.pi - math.pi / 2 + 0.52; x, y = cx + 262 * math.cos(th), cy + 262 * math.sin(th)
            dr.text((x, y), lab, font=FONT("Regular", 14), fill=TEAL, anchor="mm")
        dr.text((cx, cy - 300), "default mode network", font=FONT("SemiBold", 17), fill=TEAL, anchor="ma")
        for j, lab in enumerate(RING_LABELS_CTX):
            th = j / len(RING_LABELS_CTX) * 2 * math.pi - math.pi / 2 + 0.39; x, y = cx + 430 * math.cos(th), cy + 430 * math.sin(th)
            dr.text((x, y), lab, font=FONT("Regular", 14), fill=VIOLET, anchor="mm")
        dr.text((cx, cy - 470), "the rest of cortex", font=FONT("SemiBold", 17), fill=VIOLET, anchor="ma")
        if coupled and t < t_couple + 6: dr.text((cx, cy - 150), "coupling on", font=FONT("SemiBold", 22), fill=GOLD, anchor="ma")
        if t_sed <= t < t_sed + sed_dur: dr.text((cx, cy - 150), "sedation deepens…", font=FONT("SemiBold", 22), fill=CYAN, anchor="ma")
        im.save(f"{OUT}/s3/f_{f:05d}.png", compress_level=1)
        if f % 300 == 0: print("s3", f, nF, flush=True)
    json.dump(dict(frames=nF, dur=T, audio=[(k, float(s)) for k, s in zip(keys, starts)]), open(f"{OUT}/s3.json", "w"))

# ================================================================ intro / outro
def pulse_rings(im, t, cx, cy, col=GOLD, period=0.5, nring=5, rmax=520):
    lay = Image.new("RGBA", im.size, (0, 0, 0, 0)); ld = ImageDraw.Draw(lay)
    for k in range(nring):
        u = ((t / period) + k / nring) % 1.0; r = 40 + u * rmax; a = int(120 * (1 - u) ** 1.5)
        ld.ellipse([cx - r, cy - r * 0.62, cx + r, cy + r * 0.62], outline=hexrgb(col) + (a,), width=3)
    return Image.alpha_composite(im.convert("RGBA"), lay).convert("RGB")

def job_intro():
    La = wav_len("intro"); T = La + 1.5; nF = int(round(T * FPS)); os.makedirs(f"{OUT}/intro", exist_ok=True)
    steps = [("1", "the rate code", "awake: 10 → 60 Hz, graded"), ("2", "the burst code", "hyperpolarized + HCN1: 80 Hz or off"), ("3", "the spread", "one 2-hertz clock, hub outward")]
    t_three = 0.78 * La
    for f in range(nF):
        if f % SKIP: continue
        t = f / FPS; im = Image.new("RGB", (W, H), INK); im = pulse_rings(im, t, 1380, 560)
        d = ImageDraw.Draw(im)
        d.text((100, 250), "Rate code, burst code", font=FONT("SemiBold", 96), fill=PAPER)
        d.text((100, 380), "how one hub puts the whole brain on a two-hertz clock", font=FONT("Regular", 44), fill=GOLD)
        d.text((100, 460), "hyperpolarization, HCN1 (or T-channels), and the retrosplenial / default-mode spread", font=FONT("Regular", 26), fill=MUTED)
        for k, (num, name, sub) in enumerate(steps):
            if t > t_three + k * 1.1:
                y = 640 + k * 110; d.ellipse([100, y, 160, y + 60], outline=GOLD, width=3); d.text((130, y + 30), num, font=FONT("SemiBold", 30), fill=GOLD, anchor="mm")
                d.text((190, y - 2), name, font=FONT("SemiBold", 40), fill=PAPER); d.text((190, y + 50), sub, font=FONT("Regular", 22), fill=MUTED)
        d.text((W - 60, H - 26), "Brad Caldwell  ·  theoriesofconsciousness.com  ·  Journal of Ring Bank Phenomenology", font=FONT("Regular", 20), fill=MUTED, anchor="rd")
        im.save(f"{OUT}/intro/f_{f:05d}.png", compress_level=1)
    json.dump(dict(frames=nF, dur=T, audio=[("intro", 0.6)]), open(f"{OUT}/intro.json", "w"))

def job_outro():
    La = wav_len("outro"); T = La + 3.0; nF = int(round(T * FPS)); os.makedirs(f"{OUT}/outro", exist_ok=True)
    # a representative frame from each scene (they must already be rendered)
    pick = []
    for job, frac in (("s1", 0.33), ("s2", 0.42), ("s3", 0.9)):
        n = json.load(open(f"{OUT}/{job}.json"))["frames"]; pick.append(Image.open(f"{OUT}/{job}/f_{int(n * frac):05d}.png").convert("RGB").resize((560, 315), Image.LANCZOS))
    caps = [("awake · rate code", "an event = a change of rate: 10 → 60 Hz, graded, never silent", GOLD),
            ("hyperpolarized · burst code", "HCN1 (or T-channels): 80 Hz or off — a switch gated by a 2 Hz clock", CYAN),
            ("the spread", "the hub flips first and shouts loudest; its clock becomes the brain's", VIOLET)]
    for f in range(nF):
        if f % SKIP: continue
        t = f / FPS; im = Image.new("RGB", (W, H), INK); im = pulse_rings(im, t, 960, 1150, col=GOLD, rmax=900, nring=4)
        d = ImageDraw.Draw(im)
        d.text((60, 28), "Rate code → burst code → one clock for the whole brain", font=FONT("SemiBold", 34), fill=PAPER)
        d.text((60, 72), "Izhikevich I_Na,p + I_K neurons · SN fast subsystem + HCN1 gate · three sparsely coupled populations", font=FONT("Regular", 20), fill=MUTED)
        for k in range(3):
            if t > 0.22 * La * k + 0.5:
                x = 60 + k * 620; y = 180; im.paste(pick[k], (x, y)); d = ImageDraw.Draw(im)
                d.rectangle([x, y, x + 560, y + 315], outline=GREY, width=2)
                d.text((x, y + 335), caps[k][0], font=FONT("SemiBold", 28), fill=caps[k][2])
                # wrap caption
                words = caps[k][1].split(); lines = []; cur = ""
                for w_ in words:
                    if len(cur) + len(w_) > 46: lines.append(cur); cur = w_
                    else: cur = (cur + " " + w_).strip()
                lines.append(cur)
                for j, ln in enumerate(lines): d.text((x, y + 378 + j * 28), ln, font=FONT("Regular", 21), fill=PAPER)
        if t > La - 9:
            d.text((60, 800), "Vesuna, S. et al. Deep posteromedial cortical rhythm in dissociation. Nature 586, 87–94 (2020)", font=FONT("Regular", 20), fill=MUTED)
            d.text((60, 832), "Izhikevich, E. M. Dynamical Systems in Neuroscience (MIT Press, 2007) — the I_Na,p + I_K model and the four regimes", font=FONT("Regular", 20), fill=MUTED)
            d.text((60, 864), "Model and code in the description.", font=FONT("Regular", 20), fill=MUTED)
        if t > La - 6: d.text((W - 60, H - 26), "Brad Caldwell  ·  theoriesofconsciousness.com  ·  Journal of Ring Bank Phenomenology", font=FONT("Regular", 22), fill=MUTED, anchor="rd")
        im.save(f"{OUT}/outro/f_{f:05d}.png", compress_level=1)
    json.dump(dict(frames=nF, dur=T, audio=[("outro", 0.6)]), open(f"{OUT}/outro.json", "w"))

if __name__ == "__main__":
    job = sys.argv[1]
    {"intro": job_intro, "s1": job_s1, "s2": job_s2, "s3": job_s3, "outro": job_outro}[job]()
    print("done", job)
