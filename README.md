## LigKit: Ligand Preparation Pipeline for Molecular Docking

🧪 A fully parallelised, 11-stage Colab notebook (`LigKit.ipynb`) that takes a raw multi-ligand SDF (or ZIP of SDFs) through deduplication, standardization, tautomer and protonation-state enumeration, stereoisomer enumeration, 3D embedding and minimization, and export to both a custom RDKit-native SYBYL MOL2 writer and Meeko-based PDBQT - producing a docking-ready ligand library for AutoDock-GPU. Supports multi-library mode (queue and prepare successive SDFs in one Colab session) and includes a built-in self-test suite that validates atom typing and charge conservation before any library is committed to a production run.

## Pipeline Overview

| # | Stage | Tool | Purpose |
|---|-------|------|---------|
| 1 | Split / Extract | OpenBabel | Split multi-ligand SDF or extract ZIP of SDFs |
| 2 | Pre-Enum Dedup | RDKit | Remove duplicates before expensive enumeration |
| 3 | Standardize / Neutralize | RDKit | Converts the molecule into neutral form |
| 4 | Tautomer Enumeration | RDKit | Top-N tautomers on neutral form (no protonation locks) |
| 5 | Protonation States | Dimorphite-DL | pH-dependent ionisation states per tautomer |
| 6 | Post-Enum Dedup | RDKit | Collapse identical tautomer × protonation cross-products |
| 7 | Stereoisomer Enumeration | RDKit | Enumerate undefined stereocenters before 3D |
| 8 | 3D Generation | RDKit | Embed 3D coordinates (tiered ETKDGv3 → ETKDGv2 → ETDG fallback) |
| 9 | Energy Minimization | RDKit | Geometry relaxation (MMFF94s) |
| 10 | SDF → MOL2 | RDKit (universal SYBYL writer) | CGenFF-compatible Tripos MOL2 |
| 11 | SDF → PDBQT | Meeko | Vina-native torsion trees (OpenBabel fallback) |

## How to Run

1. **Step 1** — Install dependencies (once per Colab session)
2. **Step 2** — Upload your `.sdf` or `.zip` of SDFs
3. **Step 3** — Adjust pipeline settings, or leave defaults
4. **Step 4** — Select which stages to run
5. **Infrastructure cell** — Run once to load helper functions
6. *(Optional)* **Self-Test cell** — Validates core chemistry functions (SYBYL typing, charge conservation, MOL2 round-trip) against a fixed aspirin reference before you commit a real library
7. **Run Pipeline cell** — Executes all selected stages
8. *(Optional)* **Run Another File** — Queue the next library without re-running Steps 1–4

**Output:** `<input_name>_prepared/` containing every intermediate stage folder, a `pipeline_report.txt` provenance log, and a final combined ZIP.

## Running `ligkit.py` Locally on Linux

The same 11-stage pipeline ships as a standalone CLI script (`ligkit.py`) for anyone who wants to run it off-Colab — on a workstation, an HPC node, or inside a WSL/Docker environment. It carries no notebook-specific code; the stage logic, tool assignment, and output layout are identical to `LigKit.ipynb`. The steps below set it up cleanly on Ubuntu/Debian-based systems.

### 1. Install Open Babel (system package)

`ligkit.py` shells out to the `obabel` binary for Stage 1 (SDF splitting) and as the Stage 11 PDBQT fallback, so it must be resolvable on `$PATH` — a `pip install` alone will not provide it:

```bash
sudo apt update
sudo apt install -y openbabel
obabel -V   # confirms the binary is on PATH
```

### 2. Create and activate an isolated virtual environment

Keeping the pipeline's dependencies out of the system interpreter avoids version collisions with other cheminformatics projects on the same machine:

```bash
python3 -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip
```

### 3. Install the Python dependencies

RDKit (Stages 2–4 and 7–10) and Dimorphite-DL (Stage 5) cover the core chemistry. Meeko (Stage 11's Vina-native PDBQT export) additionally needs **Gemmi** and **SciPy** as pip-installable packages in their own right — install all five explicitly rather than relying on Meeko to resolve them on its own:

```bash
pip install rdkit dimorphite-dl meeko gemmi scipy
```

<details>
<summary>Expected environment footprint</summary>

Once installed, `venv/lib/python3.13/site-packages/` should contain, at minimum: `rdkit`, `dimorphite_dl`, `meeko`, `gemmi`, `scipy`, plus `numpy`, `loguru`, and `pillow` (`PIL`) pulled in alongside them — each paired with its own `*.dist-info` metadata folder. This is the same dependency set recorded by `pipeline_report.txt` under "Toolchain versions" at the end of every run.

</details>

### 4. Configure and run

1. Open `ligkit.py` in a text editor and edit the **USER SETTINGS** block near the top:
   - `INPUT_PATH` — full path to your input `.sdf` or `.zip` of SDFs (the only mandatory field).
   - `OUTPUT_DIR` — full path to where the `<input_name>_prepared/` folder should be written; leave `""` to place it next to the input file.
   - All remaining parameters (`MAX_WORKERS`, `DEDUP_LEVEL`, tautomer/protonation/stereoisomer limits, gen3D and MMFF minimization settings, etc.) are optional — change any of them to suit your library and hardware, or leave them at their defaults.
2. Make the script executable (one-time step):
   ```bash
   chmod +x ligkit.py
   ```
3. Run it from inside the activated virtual environment:
   ```bash
   ./ligkit.py
   # or, without chmod:
   python3 ligkit.py
   ```
4. Answer the `yes`/`no` prompts to choose which of the 11 stages to execute for this run.

Output lands in `<input_name>_prepared/` next to the input file (or under `OUTPUT_DIR` if set), with the same stage-folder layout, `pipeline_report.txt`, and combined ZIP described above — deactivate the environment with `deactivate` when finished.

## Citations

Tools invoked by this pipeline should be cited alongside LigKit itself when reporting results derived from it.

- **RDKit** (standardization, tautomer/stereoisomer enumeration, 3D embedding, MMFF94s minimization, SYBYL MOL2 writer): RDKit: Open-source cheminformatics. https://www.rdkit.org
- **Open Babel** (SDF splitting, PDBQT fallback): O'Boyle, N. M., Banck, M., James, C. A., Morley, C., Vandermeersch, T., & Hutchison, G. R. (2011). Open Babel: An open chemical toolbox. *Journal of Cheminformatics*, 3, 33. https://doi.org/10.1186/1758-2946-3-33
- **Dimorphite-DL** (pH-dependent protonation state enumeration): Ropp, P. J., Kaminsky, J. C., Yablonski, S., & Durrant, J. D. (2019). Dimorphite-DL: an open-source program for enumerating the ionization states of drug-like small molecules. *Journal of Cheminformatics*, 11, 14. https://doi.org/10.1186/s13321-019-0336-9
- **Meeko** (Vina-native PDBQT generation and torsion tree assignment, Forli Lab, Scripps Research): https://github.com/forlilab/Meeko
- **Gemmi** (structure I/O, used internally by Meeko): Wojdyr, M. (2022). GEMMI: A library for structural biology. *Journal of Open Source Software*, 7(73), 4200. https://doi.org/10.21105/joss.04200
- **SciPy** (numerical routines used internally by Meeko): Virtanen, P. et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272. https://doi.org/10.1038/s41592-019-0686-2

## Citing This Pipeline

If LigKit is used to prepare ligands reported in a publication, thesis, or preprint, cite the repository directly, in addition to the underlying tool citations above:

**Repository citation:**
LigKit: An 11-Stage Ligand Preparation Pipeline for Molecular Docking. GitHub repository: https://github.com/MANAV2502/LigKit

**BibTeX:**
```bibtex
@software{ligkit2026,
  title  = {LigKit: An 11-Stage Ligand Preparation Pipeline for Molecular Docking},
  author = {{Manav Patel}},
  year   = {2026},
  url    = {https://github.com/MANAV2502/LigKit}
}
```

Include the specific commit hash of the repository used at the time of ligand preparation (e.g., via `git rev-parse HEAD` or the commit visible in the GitHub UI), since the notebook may be revised over time and downstream reproducibility depends on pinning the exact version consulted, alongside the toolchain versions recorded in `pipeline_report.txt`.

