import os
from gpaw import GPAW
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

DEFECT = "hBN_5x5_C-VN"
gpw_path = os.path.join(DEFECT, f"{DEFECT}_fd.gpw")

calc = GPAW(gpw_path, txt=None)

energies, dos = calc.get_dos(npts=2000, width=0.1)
fermi = calc.get_fermi_level()
energies -= fermi

plt.figure(figsize=(6,4))
plt.plot(energies, dos, label="Total DOS")
plt.axvline(0, color="k", linestyle="--", label="VBM")
plt.xlim(-3, 5)
plt.xlabel("Energy relative to VBM (eV)")
plt.ylabel("DOS (a.u.)")
plt.legend()
plt.tight_layout()
plt.savefig("dos_hBN_5x5_C-VN.png", dpi=300)
plt.show()

# Spin-resolved DOS (if spin-polarized)
energies, dos_up = calc.get_dos(spin=0, npts=2000, width=0.1)
energies, dos_dn = calc.get_dos(spin=1, npts=2000, width=0.1)
energies -= fermi

plt.figure(figsize=(6,4))
plt.plot(energies, dos_up, label="Spin ↑")
plt.plot(energies, -dos_dn, label="Spin ↓")
plt.axvline(0, color="k", linestyle="--")
plt.xlabel("Energy relative to VBM (eV)")
plt.ylabel("Spin DOS (a.u.)")
plt.legend()
plt.tight_layout()
plt.savefig("spin_dos_hBN_5x5_C-VN.png", dpi=300)
plt.show()

# =========================
# Defect level positions
# =========================

df = pd.read_csv("electronic_defect_metrics.csv")

plt.figure(figsize=(7,4))
plt.scatter(
    df["Defect"],
    df["Defect_levels_rel_VBM (eV)"],
    c=df["Spin_active"].astype(int),
    cmap="coolwarm",
    s=80
)

plt.axhline(0, color="k", linestyle="--", label="VBM")
plt.ylabel("Defect level energy (eV)")
plt.title("Defect Levels Relative to VBM in h-BN")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("defect_levels_rel_vbm.png", dpi=300)
plt.show()