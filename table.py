import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

# ==========================
# Constants
# ==========================
HC = 1239.84193  # eV·nm
SIGMA_EV = 0.10
EMIN, EMAX = 1.0, 4.5

def eV_to_nm(E):
    return HC / E

# ==========================
# Load TDDFT ZPL data
# ==========================
df = pd.read_csv("all_spectra_merged.csv")

Egrid = np.linspace(0.5, 6.0, 3000)
dE = Egrid[1] - Egrid[0]
sigma_pts = SIGMA_EV / dE

zpl_rows = []

for defect, group in df.groupby("Molecule"):
    y = np.zeros_like(Egrid)

    for E, f in zip(group["Energy(eV)"], group["Osc"]):
        idx = np.argmin(np.abs(Egrid - E))
        y[idx] += f

    y = gaussian_filter1d(y, sigma_pts)

    if y.max() == 0:
        continue

    y /= y.max()

    mask = (Egrid >= EMIN) & (Egrid <= EMAX)
    idx = np.argmax(y[mask])

    Ezpl = Egrid[mask][idx]
    lam_zpl = eV_to_nm(Ezpl)

    zpl_rows.append({
        "Defect": defect,
        "ZPL_eV": Ezpl,
        "ZPL_nm": lam_zpl,
        "In_560_590_nm": 560 <= lam_zpl <= 590
    })

zpl_df = pd.DataFrame(zpl_rows)

# ==========================
# Load Purcell spectrum
# ==========================
purcell = pd.read_csv("purcell_spectrum.csv")

# ---- Detect frequency column robustly
freq_col = None
for c in purcell.columns:
    if "freq" in c.lower():
        freq_col = c
        break

if freq_col is None:
    raise RuntimeError(f"Cannot identify frequency column in purcell_spectrum.csv. Columns = {purcell.columns.tolist()}")

# Convert frequency (1/µm) → wavelength (nm)
freq_um = purcell[freq_col].values
lambda_nm = 1000.0 / freq_um

Fp = purcell["Fp"].values

# Peak Purcell mode
idx_max = np.argmax(Fp)
Fp_max = Fp[idx_max]
lambda_p = lambda_nm[idx_max]

# ==========================
# Merge ZPL–Purcell table
# ==========================
zpl_df["lambda_p_nm"] = lambda_p
zpl_df["Fp_max"] = Fp_max
zpl_df["Delta_lambda_nm"] = np.abs(zpl_df["ZPL_nm"] - lambda_p)

zpl_df = zpl_df.sort_values("Delta_lambda_nm")

# ==========================
# Save outputs
# ==========================
zpl_df.to_csv("ZPL_Purcell_matching.csv", index=False)

latex = zpl_df.copy()
latex["ZPL_eV"] = latex["ZPL_eV"].map(lambda x: f"{x:.2f}")
latex["ZPL_nm"] = latex["ZPL_nm"].map(lambda x: f"{x:.1f}")
latex["lambda_p_nm"] = latex["lambda_p_nm"].map(lambda x: f"{x:.1f}")
latex["Fp_max"] = latex["Fp_max"].map(lambda x: f"{x:.2f}")
latex["Delta_lambda_nm"] = latex["Delta_lambda_nm"].map(lambda x: f"{x:.1f}")

latex = latex.rename(columns={
    "Defect": "Defect",
    "ZPL_eV": "ZPL (eV)",
    "ZPL_nm": "ZPL (nm)",
    "In_560_590_nm": "560--590 nm?",
    "lambda_p_nm": r"$\lambda_p$ (nm)",
    "Fp_max": r"$F_p^{\max}$",
    "Delta_lambda_nm": r"$|\Delta\lambda|$ (nm)"
})

with open("ZPL_Purcell_matching.tex", "w") as f:
    f.write(latex.to_latex(
        index=False,
        escape=False,
        caption="ZPL--Purcell matching between defect-induced emitters in h-BN and nanobeam cavity modes.",
        label="tab:zpl_purcell",
        column_format="lcccccc"
    ))

print("✔ ZPL–Purcell matching table generated successfully")
print(f"✔ Cavity peak: λp = {lambda_p:.1f} nm, Fp = {Fp_max:.2f}")