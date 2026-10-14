"""YouTube title options + description for the rate-code / burst-code video, from out_rb/timing.json"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = f"{HERE}/out_rb"
T = json.load(open(f"{OUT}/timing.json")); start = {s: t for t, s in T["chapters"]}
def ts(s): m, s = divmod(int(s), 60); return f"{m}:{s:02d}"
titles = """TITLE OPTIONS (pick one; pair it with the thumbnail that says something different)

1. Rate Code, Burst Code: How One Hub Puts the Whole Brain on a 2 Hz Clock              [pair with thumbnail A: "80 Hz OR NOTHING"]
2. Why a Sedated Brain Ticks at 2 Hz — HCN1, Bursts, and the Retrosplenial Hub           [pair with thumbnail B: "YOUR BRAIN HAS A SWITCH"]
3. From 10 Hz to 80-Hz-or-Nothing: The Burst Code Behind Sedation and Dissociation
4. A Neuron's Dial Becomes a Switch: Hyperpolarization, HCN1 and the Two-Hertz Brain
5. How Sedation Rewrites the Neural Code (and Why the Whole Brain Follows One Hub)

"""
desc = f"""Awake, a neuron reports an event with its rate: ten spikes a second at rest, fifty or sixty for something strong, thirty for something weaker — a graded dial. Hyperpolarize the same kind of cell (THC through presynaptic CB1 receptors, ketamine through NMDA receptors) and two things happen: left alone it goes silent, and its HCN1 channels start to open. The slow h-current then drives the cell around a loop — fold, fast spiking at 70–80 Hz, homoclinic, silence — about twice a second. Events no longer change a rate; they throw a switch: the whole 80 Hz burst, or nothing, and only when the clock allows. Thalamic relay cells do the same thing with T-type calcium channels.

The retrosplenial cortex is dense with HCN1 and, in Vesuna et al.'s mice, is where ketamine produced exactly this 1–3 Hz rhythm. The last part of the video shows how thirty such bursters lock into one volley, and how populations downstream — hyperpolarized and below threshold — can only answer a volley with a volley, stage by stage, until the sum is a two-hertz wave. In their mice the rhythm stayed coupled to thalamus and decoupled from most of the brain; my claim, from my own 16-channel EEG under heavier sedation, is that it spreads to the whole cap — and that the two-hertz strobe is what the frames of experience are printed with.

CHAPTERS
{ts(start['intro'])} Three steps: rate code, burst code, the spread
{ts(start['s1'])} 1 · The rate code — an event is a change of firing rate (10 → 60 Hz)
{ts(start['s2'])} 2 · The burst code — hyperpolarization + HCN1: 80 Hz or off, on a 2 Hz clock
{ts(start['s3'])} 3 · The spread — retrosplenial cortex → default mode network → the rest of cortex
{ts(start['outro'])} Summary and references

THE MODEL
• Neurons: Izhikevich's I_Na,p + I_K model (persistent sodium + potassium; two variables, V and n).
• Rate-coding cell: the class-1 (SNIC) parameter set; bias 4.63 → ~10 Hz, 6.7 → 30 Hz, 22 → 56 Hz. Time scale 3.6 ms per model unit.
• Burster: the bistable (saddle-node off the circle) parameter set, slowed 2.1×, plus an HCN gate h with
  I_h = g_h · h · (E_h − V), h∞(V) = 1 / (1 + exp((V − V_h)/k)), V_h = −60, k = 7, g_h = 1.2, E_h = −30, τ_h = 100, bias −12.
  Result: 11 spikes per burst at ~75–80 Hz, bursts at 2.05 Hz, silent phase ~350 ms. Burst onset = fold of the fast subsystem (h ≈ 0.45); burst end = homoclinic (h ≈ 0.36). Events are instantaneous voltage shoves.
• The spread: 30 retrosplenial bursters (weak all-to-all excitation), 40 "default-mode" and 60 "cortex" class-1 cells with sparse random synapses (8 and 10 inputs each), 20–30 ms delays; the downstream bias ramps down 1.2 units as "sedation deepens".
• Scenes 1–2 are slowed 10×, scene 3 is slowed 4×. Code (Python/NumPy/Matplotlib) is in the zip linked below.

RELATED
Four ways a neuron starts to fire — the four excitability regimes this model is built on: [link]
Phenomenology of Sedation (TSC 2026 talk, 22 min): [link]

REFERENCES
Vesuna, S. et al. Deep posteromedial cortical rhythm in dissociation. Nature 586, 87–94 (2020). https://doi.org/10.1038/s41586-020-2731-9
Izhikevich, E. M. Dynamical Systems in Neuroscience: The Geometry of Excitability and Bursting (MIT Press, 2007).

THE BOOK
Journal of Ring Bank Phenomenology — hardcover and EPUB (IngramSpark) · https://theoriesofconsciousness.com
Notes and case files: https://bradcaldwell.substack.com · Twitter/X @Caldwbr · caldwbr@gmail.com

Narration is a synthetic voice reading my script. Runtime {ts(T['total'])}.

#neuroscience #HCN1 #retrosplenial #defaultmodenetwork #bursting #phasespace #dynamicalsystems #sedation #dissociation #ketamine #THC #consciousness #EEG
"""
open(f"{OUT}/Rate_code_burst_code_description.md", "w").write(titles + "DESCRIPTION\n\n" + desc)
print(titles); print(desc)
