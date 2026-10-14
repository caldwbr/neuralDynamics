"""Two YouTube thumbnails (1280x720, rendered at 2x) for the rate-code / burst-code video."""
import numpy as np, os, math
from PIL import Image, ImageDraw, ImageFont, ImageFilter
from rb_net import Net, burster_params, HCN, I_BURST, SCALE_MS
from model import REGIMES
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out_rb"; os.makedirs(OUT, exist_ok=True)
W, H = 2560, 1440
INK, PAPER, MUTED, GOLD, TEAL, VIOLET, RED, CYAN = "#0B0F19", "#F4F1EA", "#A9B1C1", "#E3B23C", "#35B8A6", "#9C86E6", "#E8453C", "#4CC9F0"
F = lambda w, s: ImageFont.truetype(f"/usr/share/fonts/opentype/inter/Inter-{w}.otf", s)
def hexrgb(h): return tuple(int(h[i:i + 2], 16) for i in (1, 3, 5))

def burst_trace(T=620.0):
    net = Net([dict(name="b", N=1, p=burster_params(), I=I_BURST, hcn=HCN, noise=0.0)], dt=0.01, seed=1); net.set_rest(-70.0); net.h[:] = 0.40
    rec = net.step(int(T / 0.01), record_every=5); t = np.array([r[0] for r in rec]); v = np.array([r[1][0] for r in rec])
    m = t > 300; return t[m] - 300, v[m]
def tonic_trace(T=250.0):
    net = Net([dict(name="c", N=1, p=REGIMES["SNIC"], I=4.63, noise=0.0)], dt=0.01, seed=1); net.set_rest(-60.0); net.V[:] = -58
    net.I[:] = 4.63; rec1 = net.step(int(T / 0.01), record_every=5)
    t = np.array([r[0] for r in rec1]); v = np.array([r[1][0] for r in rec1]); m = t > 40; return t[m] - 40, v[m]

def glow_line(im, pts, col, width, blur=18, alpha=180):
    lay = Image.new("RGBA", im.size, (0, 0, 0, 0)); d = ImageDraw.Draw(lay)
    d.line(pts, fill=hexrgb(col) + (alpha,), width=width * 3); lay = lay.filter(ImageFilter.GaussianBlur(blur))
    im.alpha_composite(lay); d = ImageDraw.Draw(im); d.line(pts, fill=hexrgb(col) + (255,), width=width); return im

def rings(im, cx, cy, n=6, r0=60, dr=95, col=GOLD):
    lay = Image.new("RGBA", im.size, (0, 0, 0, 0)); d = ImageDraw.Draw(lay)
    for k in range(n):
        r = r0 + k * dr; a = int(200 * (1 - k / n) ** 1.2) + 20
        d.ellipse([cx - r, cy - r * 0.62, cx + r, cy + r * 0.62], outline=hexrgb(col) + (a,), width=6 if k == 0 else 4)
    g = lay.filter(ImageFilter.GaussianBlur(14)); im.alpha_composite(g); im.alpha_composite(lay); return im

def sticker(d, cx, cy, w, h, l1, l2, rot=-7):
    w = max(w, int(F("Bold", 60 if l2 else 66).getlength(l1)) + 70, int(F("Bold", 36).getlength(l2)) + 70 if l2 else 0)
    tmp = Image.new("RGBA", (w + 80, h + 80), (0, 0, 0, 0)); td = ImageDraw.Draw(tmp)
    td.rounded_rectangle([40, 40, 40 + w, 40 + h], radius=26, fill=RED)
    td.text((40 + w / 2, 40 + (h / 2 - 12 if l2 else h / 2)), l1, font=F("Bold", 60 if l2 else 66), fill="#FFFFFF", anchor="mm")
    if l2: td.text((40 + w / 2, 40 + h / 2 + 44), l2, font=F("Bold", 36), fill="#FFE9E7", anchor="mm")
    tmp = tmp.rotate(rot, resample=Image.BICUBIC, expand=True); return tmp

def base():
    im = Image.new("RGBA", (W, H), INK)
    glow = Image.new("RGBA", (W, H), (0, 0, 0, 0)); gd = ImageDraw.Draw(glow)
    gd.ellipse([1500, 100, 2700, 1300], fill=hexrgb(TEAL) + (55,)); gd.ellipse([-300, 200, 900, 1300], fill=hexrgb(VIOLET) + (40,))
    im.alpha_composite(glow.filter(ImageFilter.GaussianBlur(220))); return im

def variantA():
    im = base()
    # right: the burst train as a glowing trace, with rings behind
    im = rings(im, 1960, 640)
    t, v = burst_trace(); x = 1480 + (t / t[-1]) * 1020; y = 980 - (v + 80) * 5.6
    im = glow_line(im, list(zip(x, y)), GOLD, 7)
    t2, v2 = tonic_trace(); x2 = 1480 + (t2 / t2[-1]) * 1020; y2 = 1300 - (v2 + 80) * 2.2
    im = glow_line(im, list(zip(x2, y2)), TEAL, 4, blur=10, alpha=120)
    d = ImageDraw.Draw(im)
    d.text((1480, 1010), "sedated: 80 Hz bursts, 2 a second", font=F("SemiBold", 34), fill=GOLD)
    d.text((1480, 1318), "awake: 10 Hz, graded", font=F("SemiBold", 30), fill=TEAL)
    d.text((110, 300), "80 Hz", font=F("Bold", 300), fill=GOLD)
    d.text((110, 620), "OR NOTHING", font=F("Bold", 200), fill=PAPER)
    d.text((116, 880), "how sedation puts your brain on a 2 Hz clock", font=F("Italic", 64), fill=TEAL)
    d.text((116, 1000), "hyperpolarization  ·  HCN1  ·  the retrosplenial hub", font=F("Regular", 44), fill=MUTED)
    st = sticker(d, 0, 0, 560, 150, "rate code → burst code", None, -6); im.alpha_composite(st, (W - st.size[0] - 40, 70))
    d = ImageDraw.Draw(im); d.text((116, 1310), "Brad Caldwell  ·  theoriesofconsciousness.com", font=F("Regular", 36), fill=MUTED)
    return im

def variantB():
    im = base()
    # right: the burst-loop 3-D panel from scene 2, if rendered; else rings
    panel = None
    try:
        import json
        n = json.load(open(f"{HERE}/rb_frames/s2.json"))["frames"]; fr = Image.open(f"{HERE}/rb_frames/s2/f_{int(n * 0.45):05d}.png").convert("RGB")
        panel = fr.crop((1080, 560, 1820, 1020)).resize((1110, 690), Image.LANCZOS)
    except Exception as e: print("no scene-2 frame yet:", e)
    if panel is not None:
        mask = Image.new("L", panel.size, 0); ImageDraw.Draw(mask).rounded_rectangle([0, 0, panel.size[0], panel.size[1]], radius=40, fill=255)
        glow = Image.new("RGBA", (W, H), (0, 0, 0, 0)); ImageDraw.Draw(glow).rounded_rectangle([1400, 300, 1400 + 1110, 300 + 690], radius=40, fill=hexrgb(GOLD) + (90,))
        im.alpha_composite(glow.filter(ImageFilter.GaussianBlur(50))); im.paste(panel, (1400, 300), mask)
    else: im = rings(im, 1960, 640)
    d = ImageDraw.Draw(im)
    d.text((1955, 1030), "one loop: fold → 80 Hz → homoclinic → silence", font=F("SemiBold", 34), fill=GOLD, anchor="ma")
    d.text((110, 250), "YOUR BRAIN", font=F("Bold", 210), fill=PAPER)
    d.text((110, 480), "HAS A SWITCH", font=F("Bold", 210), fill=GOLD)
    d.text((116, 760), "a sedative flips it, and the whole brain", font=F("Italic", 60), fill=TEAL)
    d.text((116, 840), "starts ticking at 2 hertz", font=F("Italic", 60), fill=TEAL)
    d.text((116, 980), "rate code  →  burst code  →  one clock, hub outward", font=F("Regular", 44), fill=MUTED)
    st = sticker(d, 0, 0, 430, 150, "HCN1", "retrosplenial cortex", -7); im.alpha_composite(st, (W - st.size[0] - 60, 1090))
    d = ImageDraw.Draw(im); d.text((116, 1310), "Brad Caldwell  ·  theoriesofconsciousness.com", font=F("Regular", 36), fill=MUTED)
    return im

if __name__ == "__main__":
    for name, fn in (("Rate_code_thumbnail_A", variantA), ("Rate_code_thumbnail_B", variantB)):
        fn().convert("RGB").resize((1280, 720), Image.LANCZOS).save(f"{OUT}/{name}.png"); print("wrote", name)
