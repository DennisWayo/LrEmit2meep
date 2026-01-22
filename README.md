![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Language: Python](https://img.shields.io/badge/Language-Python%203.10+-brightgreen)
![DFT: GPAW](https://img.shields.io/badge/DFT-GPAW-orange)
![TDDFT: LR-TDDFT](https://img.shields.io/badge/TDDFT-Linear--Response-blueviolet)
![Photonics: MEEP](https://img.shields.io/badge/Photonics-MEEP-red)
![Platform: Quantum Photonics](https://img.shields.io/badge/Platform-Quantum%20Photonics-black)
![Status: Research Code](https://img.shields.io/badge/Status-Research%20Code-yellow)
![Reproducibility](https://img.shields.io/badge/Reproducibility-Documented-success)


## LrEmit2meep

LrEmit2meep is a first-principles computational framework that links atomic-scale defect engineering in monolayer hexagonal boron nitride (h-BN) to cavity-enhanced quantum emission in integrated photonic structures.

The repository accompanies the manuscript:

First-Principles Design of Quantum Emitters in PLD-Grown h-BN Monolayers for Scalable Photonic Integration

and provides all scripts and workflows used to:
- model defect-induced electronic states using spin-polarized DFT,
- identify optically active transitions using linear-response TDDFT (LR-TDDFT),
- correlate PDOS features with zero-phonon line (ZPL) energies, and
- evaluate Purcell enhancement via MEEP-based FDTD simulations.

### Motivation 

Solid-state single-photon emitters (SPEs) in h-BN have demonstrated room-temperature operation, chemical robustness, and compatibility with photonic platforms. Recent experiments on PLD-grown carbon-doped h-BN thin films report record brightness, narrow linewidths, and near-ideal photon purity. However, the microscopic electronic origin of these emitters and their compatibility with photonic cavities remain incompletely understood.

LrEmit2meep addresses this gap by providing a unified, first-principles workflow that connects:
- defect chemistry →
- localized electronic states →
- optical transitions (ZPLs) →
- cavity-enhanced emission dynamics.

### Defects Studied

The following defect configurations in monolayer h-BN are investigated:
- Pristine h-BN (reference)
- Boron vacancy (V_\mathrm{B})
- Nitrogen vacancy (V_\mathrm{N})
- Carbon substitution on boron site (C_\mathrm{B})
- Carbon substitution on nitrogen site (C_\mathrm{N})
- Carbon–nitrogen-vacancy complex (C–V_\mathrm{N})

All models are based on a 5x5 supercell with sufficient vacuum separation to eliminate spurious interactions.

### Computational Workflow 

The end-to-end workflow implemented in this repository is:
	1.	Atomic structure construction
Pristine and defected h-BN supercells
	2.	Spin-polarized DFT relaxation (GPAW)
Geometry optimization and ground-state electronic structure
	3.	Element-projected DOS (PDOS)
Identification of localized defect-derived states
	4.	Linear-response TDDFT (LR-TDDFT)
Calculation of neutral excitation energies and ZPL candidates
	5.	PDOS–ZPL correlation
Direct overlay of ZPL energies onto PDOS spectra
	6.	Photonic-crystal cavity modeling (MEEP)
Nanobeam cavity design, mode extraction, and Purcell enhancement
	7.	Emitter–cavity spectral matching
Quantitative assessment of radiative lifetime reduction

<img width="6228" height="2548" alt="workflow" src="https://github.com/user-attachments/assets/5d453100-ec1a-4d3f-93cc-faad73fae5d6" />


### Key Results

| Category | Quantity | Value | Notes |
|--------|--------|-------|-------|
| **Electronic Structure** | Pristine h-BN band gap (PBE) | ~5.1 eV | Wide-gap insulating host |
|  | Defect-induced gap (C\_B) | < 3.0 eV | Mid-gap defect levels |
|  | Spin polarization | Yes (vacancy-related) | Enables spin–photon interfaces |
| **Optical Properties (LR-TDDFT)** | ZPL (C\_B defect) | 2.17 eV (571.6 nm) | Visible SPE window |
|  | ZPL linewidth (exp.) | ~3 nm | PLD-grown h-BN benchmark |
| **Photonic Coupling (FDTD)** | Cavity resonance | 587.7 nm | Nanobeam cavity |
|  | Quality factor (Q) | ~176 | Moderate-Q regime |
|  | Peak Purcell factor | ~26 | Order-of-magnitude lifetime reduction |
| **Experimental Benchmark** | Brightness | 466 kcps @ 216 µW | PLD h-BN (Science Adv. 2025) |
|  | Single-photon purity | g²(0) = 0.015 (pulsed) | Near-ideal SPE |

<img width="2700" height="1500" alt="purcell_spectrum" src="https://github.com/user-attachments/assets/17c84e36-53f2-4a67-ab69-fe8d1e0ef9cd" />


### Software 
	•	GPAW – real-space DFT and LR-TDDFT
	•	ASE – atomic structure handling
	•	MEEP – FDTD photonic simulations
	•	NumPy / Matplotlib – data analysis and visualization


### Citation





