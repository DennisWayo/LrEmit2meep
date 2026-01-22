from gpaw import GPAW
import numpy as np
import matplotlib.pyplot as plt
import os
import gc
import warnings

warnings.filterwarnings("ignore")

# ==========================================================
# ZPL energies (eV relative to VBM ≈ EF)
# ==========================================================
ZPL_EV = {
    "hBN_5x5_CB":   2.17,
    "hBN_5x5_CN":   2.51,
    "hBN_5x5_VB":   4.13,
    "hBN_5x5_C-VN": 4.44,
    "hBN_5x5_VN":   4.50,
}

# ==========================================================
# Systems
# ==========================================================
DEFECTS = [
    "hBN_5x5_C-VN",
    "hBN_5x5_CB",
    "hBN_5x5_CN",
    "hBN_5x5_VB",
    "hBN_5x5_VN",
    "hBN_5x5_pristine",
]

NPTS  = 2000
WIDTH = 0.10

# ==========================================================
# Main loop
# ==========================================================
for DEFECT in DEFECTS:
    GPW = os.path.join(DEFECT, f"{DEFECT}_fd.gpw")
    if not os.path.exists(GPW):
        print(f"Skipping {DEFECT} (missing GPW)")
        continue

    print(f"Processing {DEFECT}")
    calc = GPAW(GPW, txt=None)
    atoms = calc.get_atoms()
    symbols = atoms.get_chemical_symbols()

    fermi = calc.get_fermi_level()

    # ------------------------------------------------------
    # Element-resolved PDOS via orbital LDOS
    # ------------------------------------------------------
    energy_ref = None
    pdos = {"B": None, "C": None, "N": None}

    for i, sym in enumerate(symbols):
        if sym not in pdos:
            continue

        e, ldos = calc.get_orbital_ldos(
            a=i,
            npts=NPTS,
            width=WIDTH
        )

        if energy_ref is None:
            energy_ref = e.copy()
            for k in pdos:
                pdos[k] = np.zeros_like(ldos)

        pdos[sym] += ldos

    if energy_ref is None:
        print(f"  No B/C/N atoms found — skipping.")
        calc.close()
        calc = None
        gc.collect()
        continue

    energies = energy_ref - fermi

    # ------------------------------------------------------
    # Normalize PDOS safely
    # ------------------------------------------------------
    for k in pdos:
        m = np.max(pdos[k])
        if m > 0:
            pdos[k] /= m

    # ------------------------------------------------------
    # Plot
    # ------------------------------------------------------
    plt.figure(figsize=(7, 4))

    plt.plot(energies, pdos["B"], lw=2, label="B PDOS")
    plt.plot(energies, pdos["C"], lw=2, label="C PDOS")
    plt.plot(energies, pdos["N"], lw=2, label="N PDOS")

    plt.axvline(0, color="k", ls="--", lw=1.5, label="EF (≈VBM)")

    if DEFECT in ZPL_EV:
        plt.axvline(
            ZPL_EV[DEFECT],
            color="crimson",
            ls=":",
            lw=2.5,
            label="ZPL"
        )

    plt.xlim(-3, 5)
    plt.xlabel("Energy relative to EF (eV)")
    plt.ylabel("Normalized PDOS (a.u.)")
    plt.title(f"Element-projected DOS: {DEFECT}")

    # De-duplicate legend entries
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(dict(zip(labels, handles)).values(),
               dict(zip(labels, handles)).keys(),
               frameon=True)

    plt.tight_layout()
    plt.savefig(f"pdos_BCN_{DEFECT}.png", dpi=300)
    plt.show()

    # ------------------------------------------------------
    # Cleanup (prevents GPAW destructor noise)
    # ------------------------------------------------------
    try:
        calc.close()
    except Exception:
        pass

    calc = None
    gc.collect()

print("All PDOS + ZPL plots generated successfully.")