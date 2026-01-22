import meep as mp
import numpy as np

# =========================
# Units: um
# =========================
resolution = 80  # px/um (increase later: 100-150)
dpml = 1.0

# Background + beam
n_bg = 1.0                 # air background
n_beam = 2.0               # e.g., SiN-ish effective index (2D)
beam_mat = mp.Medium(index=n_beam)

# Nanobeam / PhC parameters (2D)
a = 0.25                   # lattice period (um) ~ 250 nm
r = 0.075                  # hole radius (um) ~ 75 nm
w = 0.45                   # beam width (um) ~ 450 nm
Nholes_each_side = 12      # total holes = 2*N + (defect region)

# Cell
sx = 2*dpml + (2*Nholes_each_side + 6)*a
sy = 2*dpml + 4.0
cell = mp.Vector3(sx, sy, 0)

# Build geometry: beam block + periodic holes with central defect
geometry = []
geometry.append(
    mp.Block(material=beam_mat,
             center=mp.Vector3(0, 0),
             size=mp.Vector3(mp.inf, w, mp.inf))
)

# Holes along x. Defect = skip hole at x=0
for m in range(-Nholes_each_side, Nholes_each_side + 1):
    x = m * a
    if m == 0:
        continue  # defect (missing hole)
    geometry.append(
        mp.Cylinder(radius=r,
                    height=mp.inf,
                    center=mp.Vector3(x, 0),
                    material=mp.air)
    )

# =========================
# Source (broadband pulse)
# =========================
# Center frequency guess (1/um). For visible ~ 550–650 nm:
# lam = 0.55 um => f ~ 1.818 1/um ; lam = 0.65 => f ~ 1.538
f0_guess = 1.75
df = 0.6  # wide to find modes

src = mp.Source(
    src=mp.GaussianSource(frequency=f0_guess, fwidth=df),
    center=mp.Vector3(0, 0),
    component=mp.Ez,     # TE-like in 2D
    amplitude=1.0
)

sim = mp.Simulation(
    cell_size=cell,
    boundary_layers=[mp.PML(dpml)],
    geometry=geometry,
    sources=[src],
    default_material=mp.Medium(index=n_bg),
    resolution=resolution
)

# Harminv monitor at cavity center
har = mp.Harminv(mp.Ez, mp.Vector3(0, 0), f0_guess, df)

# Run: pulse + ringdown
sim.run(mp.after_sources(har), until_after_sources=400)

# Print + save the best mode
modes = sorted(har.modes, key=lambda m: m.Q, reverse=True)
print("\n=== Harminv Modes (sorted by Q) ===")
for m in modes[:8]:
    lam_nm = (1.0 / m.freq) * 1000.0
    print(f"freq={m.freq:.6f}  Q={m.Q:.1f}  lambda={lam_nm:.1f} nm  decay={m.decay:.3e}")

if len(modes) == 0:
    raise RuntimeError("No cavity modes found. Adjust f0_guess/df or geometry.")

best = modes[0]
lam_nm = (1.0 / best.freq) * 1000.0

np.savetxt(
    "cavity_mode_best.txt",
    np.array([[best.freq, best.Q, lam_nm]]),
    header="freq(1/um)  Q  lambda(nm)"
)

print("\nSaved best mode to cavity_mode_best.txt")