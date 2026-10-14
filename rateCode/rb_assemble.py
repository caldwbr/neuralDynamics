"""Segments -> mp4 with narration -> concat -> 1080p + 720p (rate-code / burst-code video)"""
import json, os, subprocess, numpy as np, soundfile as sf, sys
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out_rb"; os.makedirs(OUT, exist_ok=True)
FPS = 30; SR = 48000
ORDER = ["intro", "s1", "s2", "s3", "outro"]
NAME = "Rate_code_burst_code"
only = sys.argv[1:]   # optional: only re-encode these segments (the concat is always redone)
def load(k):
    a, sr = sf.read(f"{HERE}/rb_tts/{k}.wav", dtype="float32")
    if a.ndim > 1: a = a.mean(1)
    from scipy.signal import resample_poly
    return resample_poly(a, SR // np.gcd(SR, sr), sr // np.gcd(SR, sr)).astype(np.float32)
keys = [k for seg in ORDER for k, _ in json.load(open(f"{HERE}/rb_frames/{seg}.json"))["audio"]]
peak = max(float(np.abs(load(k)).max()) for k in keys); GAIN = 10 ** (-1 / 20) / peak
def run(c): subprocess.run(c, check=True)
parts = []; chapters = []; t = 0.0
for seg in ORDER:
    meta = json.load(open(f"{HERE}/rb_frames/{seg}.json")); dur = meta["frames"] / FPS
    if not only or seg in only:
        mix = np.zeros(int(dur * SR) + 1, np.float32)
        for k, off in meta["audio"]:
            a = load(k) * GAIN; s0 = int(off * SR); e0 = min(len(mix), s0 + len(a)); mix[s0:e0] += a[: e0 - s0]
        sf.write(f"{OUT}/a_{seg}.wav", np.clip(mix[: int(dur * SR)], -1, 1), SR, subtype="PCM_16")
        vf = f"fade=t=in:st=0:d=0.5:color=0x0B0F19,fade=t=out:st={dur-0.5:.3f}:d=0.5:color=0x0B0F19,format=yuv420p"
        run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", f"{HERE}/rb_frames/{seg}/f_%05d.png", "-i", f"{OUT}/a_{seg}.wav", "-vf", vf, "-c:v", "libx264", "-preset", "veryfast", "-crf", "20", "-pix_fmt", "yuv420p", "-r", str(FPS),
             "-profile:v", "high", "-level", "4.1", "-c:a", "aac", "-b:a", "128k", "-ar", str(SR), "-ac", "2", "-shortest", "-movflags", "+faststart", f"{OUT}/v_{seg}.mp4"])
        print("encoded", seg, round(dur, 1), flush=True)
    parts.append(f"v_{seg}.mp4"); chapters.append((t, seg)); t += dur
open(f"{OUT}/list.txt", "w").write("".join(f"file '{p}'\n" for p in parts))
final = f"{OUT}/{NAME}.mp4"
run(["ffmpeg", "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", f"{OUT}/list.txt", "-c", "copy", "-movflags", "+faststart", final])
run(["ffmpeg", "-y", "-loglevel", "error", "-i", final, "-vf", "scale=1280:720:flags=lanczos", "-c:v", "libx264", "-preset", "veryfast", "-crf", "26", "-pix_fmt", "yuv420p", "-c:a", "aac", "-b:a", "96k", "-movflags", "+faststart", f"{OUT}/{NAME}_720p.mp4"])
def ts(s): m, s = divmod(int(s), 60); return f"{m}:{s:02d}"
names = {"intro": "Three steps: rate code, burst code, the spread", "s1": "1 · The rate code — an event is a change of firing rate (10 → 60 Hz)", "s2": "2 · The burst code — hyperpolarization + HCN1: 80 Hz or off, on a 2 Hz clock", "s3": "3 · The spread — retrosplenial cortex → default mode network → the rest of cortex", "outro": "Summary and references"}
open(f"{OUT}/chapters.txt", "w").write("\n".join(f"{ts(t0)} {names[s]}" for t0, s in chapters) + f"\n\ntotal {ts(t)}\n")
json.dump(dict(chapters=[(t0, s) for t0, s in chapters], total=t), open(f"{OUT}/timing.json", "w"))
print(open(f"{OUT}/chapters.txt").read())
