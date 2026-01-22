import os
import numpy as np
import pandas as pd
from gpaw import GPAW

BASE = "."
rows = []

for d in sorted(os.listdir(BASE)):
    path = os.path.join(BASE, d)
    if not os.path.isdir(path):
        continue

    gpw = None
    for fn in os.listdir(path):
        if fn.endswith("_fd.gpw"):
            gpw = os.path.join(path, fn)
            break
    if gpw is None:
        continue

    calc = GPAW(gpw, txt=None)

    # =========================
    # Eigenvalues / KS gap
    # =========================
    eigs = calc.get_eigenvalues(spin=0)
    fermi = calc.get_fermi_level()

    vbm = np.max(eigs[eigs <= fermi])
    cbm = np.min(eigs[eigs > fermi])
    ks_gap = cbm - vbm

    # =========================
    # Magnetic moment
    # =========================
    try:
        magmom = calc.get_magnetic_moment()
    except:
        magmom = 0.0

    spin_active = abs(magmom) > 0.1

    # =========================
    # Defect level positions
    # =========================
    defect_levels = eigs[(eigs > vbm - 1.0) & (eigs < cbm + 1.0)]
    defect_levels_rel = defect_levels - vbm

    # =========================
    # Charge density localization proxy
    # =========================
    rho = calc.get_all_electron_density(gridrefinement=2)
    rho_flat = rho.ravel()

    rho_mean = np.mean(rho_flat)
    rho2_mean = np.mean(rho_flat**2)

    density_ipr = rho2_mean / (rho_mean**2 + 1e-12)
    rho_max = np.max(rho_flat)

    rows.append({
        "Defect": d,
        "KS_gap (eV)": ks_gap,
        "VBM (eV)": vbm,
        "CBM (eV)": cbm,
        "Defect_levels_rel_VBM (eV)": np.round(defect_levels_rel, 3).tolist(),
        "Magnetic_moment (μB)": magmom,
        "Spin_active": spin_active,
        "Density_IPR_proxy": density_ipr,
        "Max_charge_density": rho_max,
    })

df = pd.DataFrame(rows)
df.to_csv("electronic_defect_metrics.csv", index=False)
print(df)