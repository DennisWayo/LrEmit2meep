#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HC = 1239.84193  # eV*nm  (lambda[nm] = HC/E[eV])

# =========================
# USER SETTINGS
# =========================
ZPL_WINDOW_EV = (1.0, 4.5)     # where you searched for "lowest-energy bright"
BRIGHT_THR = 1e-6             # oscillator strength threshold for "bright"
TAU0_NS = 1.0                 # assumed free-space lifetime in ns (change if you want)
USE_GLOBAL_PEAK = True        # if False, compute per-defect Fp at ZPL only
CLIP_NEGATIVE_FP = True       # Purcell should be >=0 physically; negatives come from numerical noise
# =========================

def ev_to_nm(E):
    return HC / E

def load_purcell_csv(path="purcell_spectrum.csv"):
    pur = pd.read_csv(path)

    # find freq column
    freq_col = None
    for c in pur.columns:
        if "freq" in c.lower():
            freq_col = c
            break
    if freq_col is None:
        raise RuntimeError(f"Cannot find freq column in purcell_spectrum.csv columns={pur.columns.tolist()}")

    # find Fp column
    fp_col = None
    for c in pur.columns:
        cl = c.strip().lower()
        if cl == "fp" or "purcell" in cl:
            fp_col = c
            break
    if fp_col is None:
        # fallback: any col containing 'fp'
        for c in pur.columns:
            if "fp" in c.lower():
                fp_col = c
                break
    if fp_col is None:
        raise RuntimeError(f"Cannot find Purcell column in purcell_spectrum.csv columns={pur.columns.tolist()}")

    pur = pur.rename(columns={freq_col: "freq_1_per_um", fp_col: "Fp"}).copy()
    pur["freq_1_per_um"] = pd.to_numeric(pur["freq_1_per_um"], errors="coerce")
    pur["Fp"] = pd.to_numeric(pur["Fp"], errors="coerce")

    # wavelength in nm: freq is 1/um => lambda(um)=1/freq => lambda(nm)=1000/freq
    pur["wavelength_nm"] = 1000.0 / pur["freq_1_per_um"]

    pur = pur.replace([np.inf, -np.inf], np.nan).dropna(subset=["wavelength_nm", "Fp"])
    pur = pur[pur["wavelength_nm"] > 0].copy()

    if CLIP_NEGATIVE_FP:
        pur.loc[pur["Fp"] < 0, "Fp"] = 0.0

    pur = pur.sort_values("wavelength_nm").reset_index(drop=True)
    return pur

def zpl_from_transitions(df, mol):
    g = df[df["Molecule"] == mol].copy()
    g = g.replace([np.inf, -np.inf], np.nan).dropna(subset=["Energy(eV)", "Osc"])

    # window
    gw = g[(g["Energy(eV)"] >= ZPL_WINDOW_EV[0]) & (g["Energy(eV)"] <= ZPL_WINDOW_EV[1])].copy()
    if len(gw) == 0:
        gw = g.copy()

    bright = gw[gw["Osc"] > BRIGHT_THR].sort_values("Energy(eV)")
    if len(bright) > 0:
        Ezpl = float(bright.iloc[0]["Energy(eV)"])
        Fosc = float(bright.iloc[0]["Osc"])
    else:
        # fallback: lowest energy even if dim
        gw2 = gw.sort_values("Energy(eV)")
        Ezpl = float(gw2.iloc[0]["Energy(eV)"])
        Fosc = float(gw2.iloc[0]["Osc"])

    lam = float(ev_to_nm(Ezpl))
    return Ezpl, lam, Fosc

def fp_at_lambda(pur, lam_nm):
    # nearest-neighbor lookup
    idx = int(np.argmin(np.abs(pur["wavelength_nm"].to_numpy() - lam_nm)))
    return float(pur.iloc[idx]["Fp"]), float(pur.iloc[idx]["wavelength_nm"])

def main():
    # --- Load data ---
    df = pd.read_csv("all_spectra_merged.csv")
    pur = load_purcell_csv("purcell_spectrum.csv")

    # global Purcell peak (for reporting + table column)
    idx_peak = int(np.argmax(pur["Fp"].to_numpy()))
    lam_p = float(pur.iloc[idx_peak]["wavelength_nm"])
    Fp_max = float(pur.iloc[idx_peak]["Fp"])

    # --- Build ZPL–Purcell matching table ---
    rows = []
    for mol in sorted(df["Molecule"].unique()):
        Ezpl, lam_zpl, osc = zpl_from_transitions(df, mol)

        # Purcell at ZPL wavelength (more meaningful than “global peak for all defects”)
        Fp_zpl, lam_used = fp_at_lambda(pur, lam_zpl)

        # “Matching” to the peak (detuning to cavity best mode)
        dlam_to_peak = abs(lam_zpl - lam_p)

        in_window = (560.0 <= lam_zpl <= 590.0)

        # Lifetime shortening:
        # tau_cav = tau0 / Fp  (report both the factor and ns if tau0 given)
        if Fp_zpl > 0:
            tau_factor = 1.0 / Fp_zpl         # tau_cav / tau0
            tau_ns = TAU0_NS / Fp_zpl
        else:
            tau_factor = np.inf
            tau_ns = np.inf

        rows.append({
            "Defect": mol,
            "ZPL_eV": Ezpl,
            "ZPL_nm": lam_zpl,
            "in_560_590": in_window,
            "Fp_at_ZPL": Fp_zpl,
            "tau_cav_over_tau0": tau_factor,
            f"tau_cav_ns_assuming_tau0_{TAU0_NS}ns": tau_ns,
            "lambda_p_nm": lam_p,
            "Fp_max": Fp_max,
            "delta_lambda_to_peak_nm": dlam_to_peak,
        })

    tab = pd.DataFrame(rows).sort_values("delta_lambda_to_peak_nm").reset_index(drop=True)

    tab.to_csv("zpl_purcell_matching_with_lifetime.csv", index=False)

    # --- Save LaTeX table (clean) ---
    def yesno(x): return "Yes" if bool(x) else "No"
    tex = []
    tex.append(r"\begin{table}[t]")
    tex.append(r"\centering")
    tex.append(r"\caption{ZPL--Purcell matching between defect-induced emitters in h-BN and a 2D nanobeam cavity. "
               r"$F_p(\lambda_{\mathrm{ZPL}})$ is read from the simulated Purcell spectrum at each ZPL wavelength. "
               r"Radiative lifetime shortening is estimated as $\tau_{\mathrm{cav}}=\tau_0/F_p$.}")
    tex.append(r"\label{tab:zpl_purcell}")
    tex.append(r"\begin{tabular}{lcccccc}")
    tex.append(r"\toprule")
    tex.append(r"Defect & ZPL (eV) & ZPL (nm) & 560--590? & $F_p(\lambda_{\mathrm{ZPL}})$ & $\tau_{\mathrm{cav}}/\tau_0$ & $|\Delta\lambda|$ (nm) \\")
    tex.append(r"\midrule")
    for _, r in tab.iterrows():
        tex.append(
            f"{r['Defect'].replace('_', r'\\_')} & "
            f"{r['ZPL_eV']:.2f} & {r['ZPL_nm']:.1f} & {yesno(r['in_560_590'])} & "
            f"{r['Fp_at_ZPL']:.2f} & "
            f"{(r['tau_cav_over_tau0'] if np.isfinite(r['tau_cav_over_tau0']) else 0):.3g} & "
            f"{r['delta_lambda_to_peak_nm']:.1f} \\\\"
        )
    tex.append(r"\bottomrule")
    tex.append(r"\end{tabular}")
    tex.append(r"\end{table}")

    with open("zpl_purcell_matching_with_lifetime.tex", "w") as f:
        f.write("\n".join(tex))

    print("\nSaved:")
    print(" - zpl_purcell_matching_with_lifetime.csv")
    print(" - zpl_purcell_matching_with_lifetime.tex")

    # --- Design tolerance plot: |Δλ| vs Fp(λ_ZPL) ---
    plt.figure(figsize=(7.2, 4.8))
    plt.scatter(tab["delta_lambda_to_peak_nm"], tab["Fp_at_ZPL"])
    plt.xlabel(r"$|\Delta\lambda| = |\lambda_{\mathrm{ZPL}}-\lambda_p|$ (nm)")
    plt.ylabel(r"$F_p(\lambda_{\mathrm{ZPL}})$")
    plt.title("Design tolerance: detuning vs Purcell enhancement at ZPL")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("design_tolerance_delta_lambda_vs_Fp.png", dpi=300)
    plt.show()

    print(" - design_tolerance_delta_lambda_vs_Fp.png")

    # quick summary line
    print(f"\nCavity peak: lambda_p = {lam_p:.1f} nm, Fp_max = {Fp_max:.2f}")

if __name__ == "__main__":
    main()