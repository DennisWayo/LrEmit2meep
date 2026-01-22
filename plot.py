import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HC = 1239.84193  # eV·nm

# Load Purcell spectrum
df = pd.read_csv("purcell_spectrum.csv")

# Expected columns: f, Fp (or similar)
# Robustly detect:
cols = [c.lower() for c in df.columns]
fcol = df.columns[cols.index("f")] if "f" in cols else df.columns[0]
fpcol = df.columns[cols.index("fp")] if "fp" in cols else df.columns[1]

f = df[fcol].to_numpy()          # 1/um
Fp = df[fpcol].to_numpy()

lam_nm = 1000.0 / f              # because f is in 1/um

# Peak
imax = np.argmax(Fp)
lam_pk = lam_nm[imax]
Fp_pk = Fp[imax]
f_pk = f[imax]

fig, ax = plt.subplots(figsize=(9,5))
ax.plot(lam_nm, Fp, lw=2)

# Experimental window: 560–590 nm
ax.axvspan(560, 590, alpha=0.25, label="Experimental SPE window (560–590 nm)")

# Annotate peak
ax.plot(lam_pk, Fp_pk, "o")
ax.annotate(f"Fp={Fp_pk:.2f}\nλ={lam_pk:.1f} nm",
            xy=(lam_pk, Fp_pk),
            xytext=(lam_pk+10, Fp_pk*0.8),
            arrowprops=dict(arrowstyle="->", lw=0.8),
            fontsize=10)

ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Purcell Factor $F_p$")
ax.set_title("Purcell enhancement spectrum (2D nanobeam cavity)")
ax.grid(alpha=0.3)
ax.legend(frameon=False)

plt.tight_layout()
plt.savefig("purcell_spectrum.png", dpi=300)
plt.show()

print(f"Peak: Fp={Fp_pk:.3f} at f={f_pk:.6f}  (λ={lam_pk:.2f} nm)")


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load tables
tab = pd.read_csv("zpl_purcell_matching_with_lifetime.csv")
pur = pd.read_csv("purcell_spectrum.csv")

# Identify freq and Fp columns robustly
freq_col = [c for c in pur.columns if "freq" in c.lower()][0]
fp_col = [c for c in pur.columns if "fp" in c.lower()][0]

pur["wavelength_nm"] = 1000.0 / pur[freq_col]
pur["Fp"] = pur[fp_col]
pur = pur[(pur["Fp"] > 0) & (pur["wavelength_nm"] > 0)]

# Tolerance window (nm)
delta = 10.0   # ±10 nm tolerance

Fp_eff = []
for lam_zpl in tab["ZPL_nm"]:
    win = pur[
        (pur["wavelength_nm"] >= lam_zpl - delta) &
        (pur["wavelength_nm"] <= lam_zpl + delta)
    ]
    Fp_eff.append(win["Fp"].max() if len(win) else 0.0)

tab["Fp_eff_pm10nm"] = Fp_eff

# Plot
plt.figure(figsize=(7,5))
plt.scatter(tab["delta_lambda_to_peak_nm"], tab["Fp_eff_pm10nm"], s=80)
plt.xlabel(r"$|\lambda_{\mathrm{ZPL}}-\lambda_p|$ (nm)")
plt.ylabel(r"$\max_{\pm 10\,\mathrm{nm}}\,F_p$")
plt.title("Design tolerance: effective Purcell enhancement")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("design_tolerance_windowed.png", dpi=300)
plt.show()