import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("electronic_defect_metrics.csv")

# ---------- Band gap ----------
plt.figure(figsize=(7,4))
plt.bar(df["Defect"], df["KS_gap (eV)"], color="tab:blue")
plt.xticks(rotation=45, ha="right")
plt.ylabel("KS Band Gap (eV)")
plt.title("DFT Band Gap Across h-BN Defects")
plt.tight_layout()
plt.savefig("bandgap_comparison.png", dpi=300)
plt.close()

# ---------- Magnetic moment ----------
plt.figure(figsize=(7,4))
plt.bar(df["Defect"], df["Magnetic_moment (μB)"], color="tab:red")
plt.xticks(rotation=45, ha="right")
plt.ylabel("Magnetic Moment (μB)")
plt.title("Spin Polarization of h-BN Defects")
plt.tight_layout()
plt.savefig("magnetic_moment.png", dpi=300)
plt.close()

# ---------- Charge localization (IPR proxy) ----------
plt.figure(figsize=(7,4))
plt.bar(df["Defect"], df["Density_IPR_proxy"], color="tab:purple")
plt.xticks(rotation=45, ha="right")
plt.ylabel("Charge Localization Proxy (IPR-like)")
plt.title("Defect Charge Localization in h-BN")
plt.tight_layout()
plt.savefig("ipr_localization.png", dpi=300)
plt.close()

# ---------- Spin activity summary ----------
spin_df = df[["Defect", "Spin_active"]]
spin_df.to_csv("spin_activity_summary.csv", index=False)