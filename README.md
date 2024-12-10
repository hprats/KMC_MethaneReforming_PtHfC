# KMC_MethaneReforming_PtHfC

This repository contains a collection of Python scripts to reproduce the Kinetic Monte Carlo (KMC) simulations performed in the research paper titled:

**"First-principles Kinetic Monte Carlo simulations for single-cluster catalysis: Study of CO₂ and CH₄ conversion on Pt/HfC"**

## 🛠️ Requirements

- **ZacrosTools** ([GitHub - hprats/ZacrosTools](https://github.com/hprats/ZacrosTools)) version 2.1 or higher.

## 📂 Scripts

- **`write_input_files.py`**  
  Generates all Zacros input files for a (Px, Py) scan of KMC simulations for the following reactions on HfC or Pt/HfC:
  - Dry Reforming of Methane (DRM)
  - Steam Reforming of Methane (SRM)
  - Partial Oxidation of Methane (POM)
  - Water-Gas Shift (WGS)
  - Reverse Water-Gas Shift (RWGS)

- **`plot_heatmaps.py`**  
  Creates figures with heatmaps for Turnover Frequency (TOF), selectivity, coverage, and kinetic phase diagrams based on the results of a (Px, Py) scan of KMC simulations for the above reactions on HfC or Pt/HfC.
