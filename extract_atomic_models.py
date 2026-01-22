import os
import pandas as pd
import numpy as np
from ase.io import read
from ase.neighborlist import NeighborList

BASE = "."
rows = []

# Covalent radii (Å) – used to build neighbor list
COV_RADII = {
    "B": 0.84,
    "N": 0.71,
    "C": 0.76,
}

def find_defect_atoms(symbols, pristine_counts):
    """Return indices of atoms differing from pristine composition."""
    idx = []
    for i, s in enumerate(symbols):
        if s not in pristine_counts or pristine_counts[s] < symbols.count(s):
            idx.append(i)
    return idx

for d in sorted(os.listdir(BASE)):
    path = os.path.join(BASE, d)
    if not os.path.isdir(path):
        continue

    # Find relaxed structure
    struct_file = None
    for fn in os.listdir(path):
        if "relaxed" in fn and fn.endswith((".xyz", ".cif", ".traj")):
            struct_file = os.path.join(path, fn)
            break

    if struct_file is None:
        continue

    atoms = read(struct_file)
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    cell = atoms.cell

    # Composition
    unique, counts = np.unique(symbols, return_counts=True)
    composition = dict(zip(unique, counts))

    # Define pristine reference
    pristine_counts = {"B": 25, "N": 25}

    defect_indices = find_defect_atoms(symbols, pristine_counts)
    if not defect_indices:
        # pristine: use center atom
        defect_indices = [len(atoms) // 2]

    # Build neighbor list
    cutoffs = [1.2 * COV_RADII.get(s, 0.8) for s in symbols]
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    bond_lengths = []
    bond_angles = []

    for i in defect_indices:
        indices, offsets = nl.get_neighbors(i)
        ri = positions[i]

        for j in indices:
            rj = positions[j]
            bond_lengths.append(np.linalg.norm(ri - rj))

        # Bond angles: j–i–k
        for j in indices:
            for k in indices:
                if j >= k:
                    continue
                v1 = positions[j] - ri
                v2 = positions[k] - ri
                cosang = np.dot(v1, v2) / (
                    np.linalg.norm(v1) * np.linalg.norm(v2)
                )
                angle = np.degrees(np.arccos(np.clip(cosang, -1, 1)))
                bond_angles.append(angle)

    rows.append({
        "Defect": d,
        "Total_atoms": len(atoms),
        "Cell_x (Å)": round(cell.lengths()[0], 2),
        "Cell_y (Å)": round(cell.lengths()[1], 2),
        "Vacuum_z (Å)": round(cell.lengths()[2], 1),
        "Composition": ", ".join([f"{k}{v}" for k, v in composition.items()]),
        "⟨Bond length⟩ (Å)": f"{np.mean(bond_lengths):.2f} ± {np.std(bond_lengths):.2f}",
        "⟨Bond angle⟩ (deg)": f"{np.mean(bond_angles):.1f} ± {np.std(bond_angles):.1f}",
    })

df = pd.DataFrame(rows)
df.to_csv("atomic_models_summary_with_geometry.csv", index=False)
print(df)