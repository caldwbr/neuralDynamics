"""Segments -> mp4 with narration -> concat -> 1080p + 720p"""
import json, os, subprocess, numpy as np, soundfile as sf
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out"; os.makedirs(OUT, exist_ok=True)
FPS = 30; SR = 48000
ORDER = ["intro", "SNIC", "SN", "subH", "supH", "outro"]
GAIN = None
def load(k):
    a, sr = sf.read(f"{HERE}/tts/{k}.wav", dtype="float32")
    if a.ndim > 1: a = a.mean(1)
    from scipy.signal import resample_poly
    return resample_poly(a, SR // np.gcd(SR, sr), sr // np.gcd(SR, sr)).astype(np.float32)
peak = max(float(np.abs(load(k)).max()) for k in ["intro", "outro"] + [f"{r}_{c}" for r in ORDER[1:-1] for c in "ABC"])
GAIN = 10 ** (-1 / 20) / peak
def run(c): subprocess.run(c, check=True)
parts = []; chapters = []; t = 0.0
for seg in ORDER:
    meta = json.load(open(f"{HERE}/frames/{seg}.json")); dur = meta["frames"] / FPS
    mix = np.zeros(int(dur * SR) + 1, np.float32)
    for k, off in meta["audio"]:
        a = load(k) * GAIN; s0 = int(off * SR); e0 = min(len(mix), s0 + len(a)); mix[s0:e0] += a[: e0 - s0]
    sf.write(f"{OUT}/a_{seg}.wav", np.clip(mix[: int(dur * SR)], -1, 1), SR, subtype="PCM_16")
    vf = f"fade=t=in:st=0:d=0.5:color=0x0B0F19,fade=t=out:st={dur-0.5:.3f}:d=0.5:color=0x0B0F19,format=yuv420p"
    run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", f"{HERE}/frames/{seg}/f_%05d.png", "-i", f"{OUT}/a_{seg}.wav", "-vf", vf, "-c:v", "libx264", "-preset", "veryfast", "-crf", "20", "-pix_fmt", "yuv420p", "-r", str(FPS),
         "-profile:v", "high", "-level", "4.1", "-c:a", "aac", "-b:a", "128k", "-ar", str(SR), "-ac", "2", "-shortest", "-movflags", "+faststart", f"{OUT}/v_{seg}.mp4"])
    parts.append(f"v_{seg}.mp4"); chapters.append((t, seg)); t += dur; print("encoded", seg, round(dur, 1), flush=True)
open(f"{OUT}/list.txt", "w").write("".join(f"file '{p}'\n" for p in parts))
final = f"{OUT}/Four_ways_a_neuron_starts_to_fire.mp4"
run(["ffmpeg", "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", f"{OUT}/list.txt", "-c", "copy", "-movflags", "+faststart", final])
run(["ffmpeg", "-y", "-loglevel", "error", "-i", final, "-vf", "scale=1280:720:flags=lanczos", "-c:v", "libx264", "-preset", "veryfast", "-crf", "26", "-pix_fmt", "yuv420p", "-c:a", "aac", "-b:a", "96k", "-movflags", "+faststart", f"{OUT}/Four_ways_a_neuron_starts_to_fire_720p.mp4"])
def ts(s): m, s = divmod(int(s), 60); return f"{m}:{s:02d}"
names = {"intro": "The phase plane, and the current as a knob", "SNIC": "1 · Saddle-node on an invariant circle (class 1)", "SN": "2 · Saddle-node off the circle (bistable)", "subH": "3 · Subcritical Hopf (resonator, hysteresis)", "supH": "4 · Supercritical Hopf (oscillations that grow)", "outro": "Four doorways from rest to rhythm"}
open(f"{OUT}/chapters.txt", "w").write("\n".join(f"{ts(t0)} {names[s]}" for t0, s in chapters) + f"\n\ntotal {ts(t)}\n")
print(open(f"{OUT}/chapters.txt").read())
