from gpaw import GPAW
import numpy as np
import matplotlib.pyplot as plt
import os
import gc
import warnings
warnings.filterwarnings("ignore")



# =====================
# ZPL energies (eV, relative to VBM)
# =====================
ZPL_EV = {
    "hBN_5x5_CB":  2.17,
    "hBN_5x5_CN":  2.51,
    "hBN_5x5_VB":  4.13,
    "hBN_5x5_C-VN": 4.44,
    "hBN_5x5_VN":  4.50,
}

DEFECTS = [
    "hBN_5x5_C-VN",
    "hBN_5x5_CB",
    "hBN_5x5_CN",
    "hBN_5x5_VB",
    "hBN_5x5_VN",
    "hBN_5x5_pristine",
]

NPTS = 2000
WIDTH = 0.1

for DEFECT in DEFECTS:
    GPW = os.path.join(DEFECT, f"{DEFECT}_fd.gpw")
    if not os.path.exists(GPW):
        print(f"Skipping {DEFECT} (no GPW)")
        continue

    print(f"Processing {DEFECT}")
    calc = GPAW(GPW, txt=None)
    atoms = calc.get_atoms()
    symbols = atoms.get_chemical_symbols()

    fermi = calc.get_fermi_level()

    # Build element-resolved PDOS arrays
    e_ref = None
    pdos = {"B": None, "C": None, "N": None}

    for i, sym in enumerate(symbols):
        if sym not in pdos:
            continue

        e, ldos = calc.get_orbital_ldos(a=i, npts=NPTS, width=WIDTH)

        if e_ref is None:
            e_ref = e.copy()
            for k in pdos:
                pdos[k] = np.zeros_like(ldos)

        pdos[sym] += ldos

    if e_ref is None:
        print(f"  (No B/C/N atoms found?) Skipping plot for {DEFECT}")
        # cleanup
        calc = None
        gc.collect()
        continue

    energies = e_ref - fermi  # energy relative to VBM-ish reference (fermi)

    # Normalize safely
    for k in ["B", "C", "N"]:
        m = float(np.max(pdos[k]))
        if m > 0:
            pdos[k] /= m

    # Plot
    plt.figure(figsize=(7, 4))
    plt.plot(energies, pdos["B"], label="B PDOS")
    plt.plot(energies, pdos["C"], label="C PDOS")
    plt.plot(energies, pdos["N"], label="N PDOS")
    plt.axvline(0, color="k", ls="--", label="VBM/EF ref")

    plt.xlim(-3, 5)
    plt.xlabel("Energy relative to reference (eV)")
    plt.ylabel("Normalized PDOS (a.u.)")
    plt.title(f"Element-projected DOS: {DEFECT}")
    plt.legend()

    # =====================
    # Overlay ZPL line (if available)
    # =====================
    if DEFECT in ZPL_EV:
        plt.axvline(
            ZPL_EV[DEFECT],
            color="crimson",
            linestyle=":",
            linewidth=2,
            label="ZPL"
        )
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())

    plt.tight_layout()
    plt.savefig(f"pdos_BCN_{DEFECT}.png", dpi=300)
    plt.show()

    # Hard cleanup to reduce GPAW destructor noise
    try:
        calc.close()
    except Exception:
        pass
    calc = None
    gc.collect()