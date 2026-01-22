import os
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

    calc = None
    try:
        # IMPORTANT: avoid txt=None in legacy GPAW
        calc = GPAW(gpw, txt=os.devnull)

        params = calc.parameters

        occ = params.get("occupations", {})
        if isinstance(occ, dict):
            smearing = occ.get("width", None)
            occ_name = occ.get("name", "unknown")
        else:
            smearing = getattr(occ, "width", None)
            occ_name = occ.__class__.__name__

        rows.append({
            "Defect": d,
            "XC": params.get("xc"),
            "Mode": params.get("mode"),
            "Grid_spacing_h (Å)": params.get("h"),
            "kpts": params.get("kpts"),
            "Spin_polarized": params.get("spinpol"),
            "Occupation": occ_name,
            "Smearing (eV)": smearing,
        })

    finally:
        # Explicit close prevents __del__ warnings
        if calc is not None:
            try:
                calc.close()
            except Exception:
                pass
            del calc

df = pd.DataFrame(rows)
df.to_csv("electronic_structure_summary.csv", index=False)
print(df)