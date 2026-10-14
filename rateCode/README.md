# Rate code, burst code — model and animation code

Python 3 · numpy · scipy · matplotlib · Pillow · soundfile (+ ffmpeg for assembly, sherpa-onnx/Kokoro for the narration)

## The model
`model.py` — Izhikevich's I_Na,p + I_K neuron (V, n): dV/dt = (I − g_L(V−E_L) − g_Na m∞(V)(V−E_Na) − g_K n (V−E_K))/C, dn/dt = (n∞(V) − n)/τ.
Parameter sets `REGIMES`: SNIC (class 1, used for the rate-coding cell and the downstream populations) and SN (bistable; the burster's fast subsystem).

`rb_net.py` — a small network class: populations of these neurons, an optional HCN gate h per cell
(I_h = g_h · h · (E_h − V), h∞ = 1/(1+exp((V−V_h)/k)), dh/dt = (h∞ − h)/τ_h), synaptic gating s (jumps by 1 at each spike, decays with τ_s),
arbitrary weight matrix W with per-population-pair delays, voltage noise. Constants at the top: SCALE_MS (3.6 ms per model unit),
STRETCH (the burster's C and τ_n are both ×2.1 — a pure time rescaling of the fast subsystem), HCN parameters, I_BURST.

`rb_tune.py`, `rb_search2.py` — the parameter searches (SNIC period vs current; burster grid search).
`rb_analyze.py` — the burster's fast-subsystem geometry parametrized by h (equilibria, limit cycles, separatrices) → data/burst_h.pkl;
`analyze.py`/`render3d.py` do the same job for the four regimes (from the companion video) and draw the 3-D stacks.
`rb_pop.py` — scene 3: 30 RSP bursters (weak all-to-all excitation, switched on at a cue) → 40 DMN cells → 60 cortical cells,
sparse random synapses (8 and 10 inputs per cell), delays 6 and 8 units (22 and 29 ms); the downstream bias ramps down as "sedation deepens".

## The video
`rb_narration.py` → `rb_gen_tts.py` (Kokoro am_michael) → `rb_frames.py <job>` (intro, s1, s2, s3, outro; 1920×1080 @ 30 fps) → `rb_assemble.py` → `rb_describe.py`, `rb_thumb.py`.
Scenes 1–2 run 10× slower than real time, scene 3 4× slower; the on-screen clock is real time.
