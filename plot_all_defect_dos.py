import os
import numpy as np
import matplotlib.pyplot as plt
from gpaw import GPAW

BASE = "."
NPTS = 2000
WIDTH = 0.1

# Containers
dos_data = {}
spin_dos_data = {}

for d in sorted(os.listdir(BASE)):
    path = os.path.join(BASE, d)
    if not os.path.isdir(path):
        continue

    gpw = os.path.join(path, f"{d}_fd.gpw")
    if not os.path.exists(gpw):
        continue

    print(f"Reading DOS for {d}")
    calc = GPAW(gpw, txt=None)

    # --- Total DOS ---
    energies, dos = calc.get_dos(npts=NPTS, width=WIDTH)
    fermi = calc.get_fermi_level()
    energies -= fermi

    dos_data[d] = (energies, dos)

    # --- Spin DOS if spin-polarized ---
    if calc.get_spin_polarized():
        e_up, dos_up = calc.get_dos(spin=0, npts=NPTS, width=WIDTH)
        e_dn, dos_dn = calc.get_dos(spin=1, npts=NPTS, width=WIDTH)
        e_up -= fermi
        e_dn -= fermi
        spin_dos_data[d] = (e_up, dos_up, dos_dn)

# =========================
# Plot: Total DOS (ALL)
# =========================
plt.figure(figsize=(7,5))

for d, (E, D) in dos_data.items():
    plt.plot(E, D, label=d.replace("_", " "))

plt.axvline(0, color="k", linestyle="--", label="VBM")
plt.xlim(-3, 5)
plt.xlabel("Energy relative to VBM (eV)")
plt.ylabel("DOS (a.u.)")
plt.title("Total DOS of Defect-Engineered h-BN")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig("all_defects_total_dos.png", dpi=300)
plt.show()

# =========================
# Plot: Spin DOS (ALL)
# =========================
plt.figure(figsize=(7,5))

for d, (E, up, dn) in spin_dos_data.items():
    plt.plot(E, up, label=f"{d} ↑")
    plt.plot(E, -dn, linestyle="--", label=f"{d} ↓")

plt.axvline(0, color="k", linestyle="--")
plt.xlim(-3, 5)
plt.xlabel("Energy relative to VBM (eV)")
plt.ylabel("Spin DOS (a.u.)")
plt.title("Spin-Resolved DOS of Spin-Active h-BN Defects")
plt.legend(fontsize=7, ncol=2)
plt.tight_layout()
plt.savefig("all_defects_spin_dos.png", dpi=300)
plt.show()