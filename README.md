# LrEmit2meep

This repository preserves an exploratory calculation workflow for point defects in monolayer hexagonal boron nitride (h-BN) and a separate two-dimensional nanobeam model. It accompanies the manuscript *Carbon Defects in Monolayer Hexagonal Boron Nitride for Quantum Light Sources: Optical Screening and Cavity Limits*.

## What the calculations show

The archived neutral-defect screen used a 5 × 5 h-BN cell, PBE ground states, and GPAW linear-response TDDFT. Its reported optical features are **vertical excitation proxies**, not relaxed-emission zero-phonon lines (ZPLs). In particular, the archived 2.167 eV C_B feature cannot identify the emitter in a pulsed-laser-deposited film merely because its wavelength falls in a measured emission window. The archived C_B and C_N cells were also spin paired despite odd electron counts, and the archived spectrum used a nonstandard norm of GPAW's oscillator-strength components.

A separate matched, spin-polarized comparison in the manuscript gives dominant broadened vertical features of 2.483 eV for neutral C_B and 3.808 eV for neutral C_B–C_N–C_B (CBC) under one PBE/LR-TDDFT protocol. These values are screening descriptors. Cell size, response settings, structural alternatives, charge state, and relaxed emission remain to be converged or calculated. The comparison does not overturn the published CBC assignment based on broader evidence for the experimental film.

The archived photonic calculation reports a 611.2 nm Harminv mode with Q ≈ 176 in a two-dimensional effective-index nanobeam. Its sampled spectrum is a device-to-**air** local-density-of-states (LDOS) ratio, not an isolated cavity-only Purcell enhancement for an h-BN emitter. The largest sampled ratio, 3.57, occurs at the 550 nm boundary and is not established as a resonance. The ratio near the archived C_B proxy wavelength of 572.1 nm is about 0.94. These outputs do not establish cavity-enhanced emission or a reduced emitter lifetime.

## Archive and reproduction

The repository retains the original calculation scripts, input structures, and generated artifacts where present. Some legacy file and column names contain `zpl` or `purcell`; those names describe the historical workflow and **do not validate** ZPL energies or a cavity-only Purcell factor. The manuscript and its supplementary comparison supply the interpretation and the separate matched-comparison inputs, logs, spectra, and scripts.

GPAW and ASE were used for the electronic calculations; Meep was used for the two-dimensional nanobeam model. Reproducing a historical output is distinct from validating its physical interpretation. A measured emitter assignment would require converged defect and excited-state calculations, and a cavity claim would require a matched-material reference, source polarization and position, and a mode-resolved three-dimensional study.

## Citation

Please cite the accompanying manuscript when using the analysis, and identify any archived result explicitly as a vertical or photonic screening quantity.
