import meep as mp
import numpy as np

# =========================
# Load best cavity mode
# =========================
mode = np.loadtxt("cavity_mode_best.txt")
f_mode = float(mode[0])  # 1/um
Q_mode = float(mode[1])
lam_nm = float(mode[2])

print(f"Using cavity mode: f={f_mode:.6f}  Q={Q_mode:.1f}  lambda={lam_nm:.1f} nm")

# =========================
# Common parameters
# =========================
resolution = 80
dpml = 1.0

n_bg = 1.0
n_beam = 2.0
beam_mat = mp.Medium(index=n_beam)

a = 0.25
r = 0.075
w = 0.45
Nholes_each_side = 12

sx = 2*dpml + (2*Nholes_each_side + 6)*a
sy = 2*dpml + 4.0
cell = mp.Vector3(sx, sy, 0)

# Dipole (emitter) placement: near cavity center
dip_pos = mp.Vector3(0.0, 0.0)   # move later to test position sensitivity
comp = mp.Ez

# Frequency window around the cavity mode for Purcell spectrum
df = f_mode / 8.0
nfreq = 250

# Flux box around the dipole/cavity region
box_half = 1.2  # um (increase if needed)

def build_flux_box(sim):
    regions = [
        mp.FluxRegion(center=mp.Vector3( box_half, 0), size=mp.Vector3(0, 2*box_half)),
        mp.FluxRegion(center=mp.Vector3(-box_half, 0), size=mp.Vector3(0, 2*box_half)),
        mp.FluxRegion(center=mp.Vector3(0,  box_half), size=mp.Vector3(2*box_half, 0)),
        mp.FluxRegion(center=mp.Vector3(0, -box_half), size=mp.Vector3(2*box_half, 0)),
    ]
    return sim.add_flux(f_mode, df, nfreq, *regions)

def run_power_spectrum(geometry, out_prefix):
    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=[mp.PML(dpml)],
        geometry=geometry,
        sources=[mp.Source(
            src=mp.GaussianSource(frequency=f_mode, fwidth=df),
            center=dip_pos,
            component=comp,
            amplitude=1.0
        )],
        default_material=mp.Medium(index=n_bg),
        resolution=resolution
    )

    flux = build_flux_box(sim)

    # Important: cavity ringdown can be long if Q is high.
    # This stop condition is safer than a fixed time.
    sim.run(until_after_sources=mp.stop_when_fields_decayed(
        50, comp, dip_pos, 1e-8
    ))

    freqs = np.array(mp.get_flux_freqs(flux))
    P = np.array(mp.get_fluxes(flux))

    np.savetxt(
        f"{out_prefix}_power.csv",
        np.column_stack([freqs, P]),
        delimiter=",",
        header="freq(1/um),P(a.u.)",
        comments=""
    )
    return freqs, P

# =========================
# Reference geometry (no cavity): homogeneous background
# =========================
freq_ref, P_ref = run_power_spectrum(geometry=[], out_prefix="ref")

# =========================
# Device geometry: nanobeam PhC cavity
# =========================
geometry = []
geometry.append(
    mp.Block(material=beam_mat,
             center=mp.Vector3(0, 0),
             size=mp.Vector3(mp.inf, w, mp.inf))
)
for m in range(-Nholes_each_side, Nholes_each_side + 1):
    x = m * a
    if m == 0:
        continue
    geometry.append(
        mp.Cylinder(radius=r, height=mp.inf, center=mp.Vector3(x, 0), material=mp.air)
    )

freq_dev, P_dev = run_power_spectrum(geometry=geometry, out_prefix="device")

# Same grid check
assert np.allclose(freq_ref, freq_dev)

Fp = np.where(P_ref > 0, P_dev / P_ref, 0.0)

# Save Purcell spectrum
np.savetxt(
    "purcell_spectrum.csv",
    np.column_stack([freq_ref, Fp]),
    delimiter=",",
    header="freq(1/um),Fp",
    comments=""
)

# Report peak near the cavity mode
idx_peak = np.argmax(Fp)
f_peak = freq_ref[idx_peak]
lam_peak_nm = (1.0 / f_peak) * 1000.0
Fp_peak = Fp[idx_peak]

print(f"\nPurcell peak: Fp={Fp_peak:.3f} at f={f_peak:.6f} (lambda={lam_peak_nm:.1f} nm)")
print("Wrote: ref_power.csv, device_power.csv, purcell_spectrum.csv")