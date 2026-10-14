# Four ways a neuron starts to fire — animation code (Python port of the MATLAB four-regime script)

model.py      the I_Na,p + I_K model, equilibria + classification, vectorized RK4, limit-cycle / unstable-cycle / separatrix finders
analyze.py    precomputes, per regime and per current level, everything the frames need  -> data/<regime>.pkl
render3d.py   the 3D manifold panel (equilibrium branches, separatrix surface, limit-cycle manifolds, bifurcation plane, current plane) at 61 current levels
frames.py     the animation: 2D phase plane + voltage trace + 3D panel, per regime (intro / SNIC / SN / subH / supH / outro)
narration.py  the narration text (voiced with Kokoro via sherpa-onnx; any TTS works — frames.py only reads the wav lengths)
assemble.py   frames + narration -> per-segment mp4 -> concat -> 1080p and 720p

Run order:  python3 analyze.py ; for r in SNIC SN subH supH: python3 render3d.py $r 61 ; (voice narration.py to tts/<key>.wav) ; python3 frames.py <job> for each job ; python3 assemble.py
Needs numpy, scipy, matplotlib, pillow, soundfile, ffmpeg.
