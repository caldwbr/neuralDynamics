import sys, time, soundfile as sf, os
sys.path.insert(0, "/tmp/claude-0/-home-claude/4a4affb1-6bc6-5c1d-8059-94aa177e5a0d/scratchpad/tts")
os.chdir("/tmp/claude-0/-home-claude/4a4affb1-6bc6-5c1d-8059-94aa177e5a0d/scratchpad/tts")
from tts import make
sys.path.insert(0, "/home/claude/neuron")
from rb_narration import N
tts = make(16, threads=1); os.makedirs("/home/claude/neuron/rb_tts", exist_ok=True)
for k, text in N.items():
    t0 = time.time(); a = tts.generate(text, sid=16, speed=1.1)
    sf.write(f"/home/claude/neuron/rb_tts/{k}.wav", a.samples, a.sample_rate)
    print(f"{k}: {len(a.samples)/a.sample_rate:.1f}s ({time.time()-t0:.0f}s)", flush=True)
print("done")
