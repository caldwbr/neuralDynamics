import pickle, os, numpy as np
from render3d import render, surfaces, branches
d = pickle.load(open("data/burst_h.pkl", "rb")); surf = surfaces(d); br = branches(d)
os.makedirs("p3d/burst_h", exist_ok=True)
for k, h in enumerate(d["Ilist"]):
    render(d, h, f"p3d/burst_h/L{k:03d}.png", size=(740, 460), surf=surf, br=br, zlabel="h")
    if k % 10 == 0: print(k, flush=True)
print("done")
