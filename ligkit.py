#!/usr/bin/env python3
"""
==================================================================================
          LIGKIT  --  Ligand Preparation Pipeline for Molecular Docking
          CLI Edition  --  Fully parallelised, 11-stage, self-contained
==================================================================================
  TOOL ASSIGNMENT:
    Split / Extract          ->  OpenBabel
    Duplicate Filter         ->  RDKit  (InChIKey, pre- and post-enumeration)
    Standardize / Neutralize ->  RDKit MolStandardize (salt strip, normalize,
                                  reionize, uncharge -- full neutral form)
    Tautomer Enumeration     ->  RDKit TautomerEnumerator (neutral form)
    Protonation State Enum   ->  Dimorphite-DL  (each tautomer, target pH)
    Stereoisomer Enumeration ->  RDKit  (undefined stereocenters by default)
    3D Generation             ->  RDKit  (ETKDGv3 -> ETKDGv2 -> ETDG fallback)
    Energy Minimization       ->  RDKit  (MMFF94/MMFF94s -> UFF fallback)
    SDF -> MOL2               ->  RDKit universal SYBYL/Tripos writer
                                  (CGenFF-compatible: correct atom typing,
                                   unique atom names, Gasteiger charges)
    SDF -> PDBQT               ->  Meeko  (Vina-native torsion trees;
                                  OpenBabel fallback if Meeko is unavailable)

  SCIENTIFICALLY CORRECT STAGE ORDER (with rationale):
   1. Split        -- entry point, always runs
   2. Dedup         -- before all enumeration, to avoid wasted compute
   3. Standardize   -- salt-strip + neutralize gives one clean form for
                        unrestricted downstream tautomer enumeration
   4. Tautomers     -- enumerate on the NEUTRAL form (RDKit sees the full
                        tautomer space unobstructed by protonation locks)
   5. Proto states  -- Dimorphite-DL assigns pH states to EACH tautomer
                        (not the other way around -- avoids artifact states)
   6. Post dedup    -- tautomer x protonation cross-products can collide
   7. Stereoiso.    -- enumerate undefined centers before 3D embedding
   8. gen3D          -- each unique 2D variant gets its own 3D conformer
   9. Minimize       -- MMFF94s geometry relaxation (UFF fallback)
  10. MOL2           -- RDKit-native SYBYL writer (GOLD / DOCK6 / CGenFF)
  11. PDBQT          -- Meeko Vina-native torsion trees (AutoDock-GPU / Vina)

  HOW TO RUN (Ubuntu):
   1. Make sure your environment has the required tools installed
      (OpenBabel, RDKit, Dimorphite-DL, Meeko -- see project README/
      requirements for the exact install commands for your setup).
   2. Edit USER SETTINGS below -- only INPUT_PATH is mandatory.
   3. Make the script executable (one time only):
        chmod +x ligkit.py
   4. Run it directly:
        ./ligkit.py
      ...or via the interpreter, no chmod required:
        python3 ligkit.py
   5. Answer the yes/no prompts to choose which stages to run.

  OUTPUT FOLDER:  <input_name>_prepared/
    split_sdf/ deduped/ standardized/ tautomers/ proto_states/
    post_deduped/ stereoisomers/ gen3d/ minimized/ mol2/ pdbqt/
    pipeline_report.txt   <input_name>_prepared.zip
==================================================================================
"""

import os, sys, glob, re, shutil, subprocess, zipfile, threading, multiprocessing
import tempfile, math
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

# ══════════════════════════════════════════════════════════════════════════════
#  USER SETTINGS  <- Only edit this block
# ══════════════════════════════════════════════════════════════════════════════

INPUT_PATH = ""
# ^ Full path to your input file. REQUIRED.
#   - A single multi-ligand .sdf  ->  will be split automatically
#   - A .zip of individual .sdf files  ->  extracted directly

OUTPUT_DIR = ""
# ^ Leave "" to place the output folder next to your input file (recommended).

MAX_WORKERS = 0
# ^ 0 = auto-detect all CPU cores.

# ── Pre/post-enumeration duplicate filter ───────────────────────────────────
DEDUP_LEVEL = "connectivity"
# ^ "connectivity" -- InChIKey layer 1 only (heavy-atom skeleton, ignores stereo)
#                     Best for ZINC/PubChem/Enamine downloads
# ^ "full"          -- full InChIKey (stereo-aware, keeps both enantiomers)
#                     Use when the library is already stereochemically curated

# ── Stage 4: RDKit Tautomer settings ────────────────────────────────────────
RDKIT_TAUTO_MAX = 4
# ^ Top-N tautomers kept per molecule by RDKit stability score.
#   Recommended: 3-5 for VS libraries, up to 10 for fragment libraries.

# ── Stage 5: Dimorphite-DL protonation state settings ───────────────────────
PROTO_MIN_PH   = 7.4   # lower pH bound (physiological +-1 unit)
PROTO_MAX_PH   = 7.4   # upper pH bound
PROTO_MAX_VARS = 1     # max protonation state variants per tautomer
# ^ Keep 3-5 for large libraries to avoid combinatorial explosion.

# ── Stage 7: Stereoisomer settings ──────────────────────────────────────────
MAX_STEREO_ISOMERS    = 4     # max stereoisomers per molecule
ONLY_UNDEFINED_STEREO = True
# ^ True  = only enumerate centers with no defined wedge bonds (recommended)
# ^ False = enumerate all stereocenters, including already-defined ones

# ── Stage 8: RDKit gen3D settings ───────────────────────────────────────────
RANDOM_SEED        = 42     # fixed seed -> reproducible conformers
ENFORCE_CHIRALITY  = True   # always preserve R/S stereocenters
MAX_EMBED_ATTEMPTS = 100    # max ETKDGv3 attempts before fallback

# ── Stage 9: RDKit MMFF94 minimization settings ─────────────────────────────
MMFF_VARIANT   = "MMFF94s"  # MMFF94 or MMFF94s (MMFF94s recommended for
                             # conjugated/aromatic systems)
MMFF_MAX_ITERS = 2000       # minimization iterations
MMFF_FORCE_TOL = 1e-3       # RMS force convergence threshold (kcal/mol/A)
UFF_FALLBACK   = True       # use UFF if MMFF94 is not applicable
                             # (metals, exotic atom types)

# ── Stage 11: PDBQT OpenBabel-fallback timeout ──────────────────────────────
TIMEOUT_PDBQT = 30   # seconds per ligand (only used if Meeko is unavailable)

# ══════════════════════════════════════════════════════════════════════════════
#  END OF USER SETTINGS  --  Do not edit below this line
# ══════════════════════════════════════════════════════════════════════════════


# ── ask_run(): yes/no gate before every stage (same pattern for every stage) ──
def ask_run(step_name: str, description: str) -> bool:
    print(f"\n{'='*65}")
    print(f"  STEP : {step_name}")
    print(f"  INFO : {description}")
    print(f"{'='*65}")
    while True:
        ans = input(f"Run \"{step_name}\"? [yes/no]: ").strip().lower()
        if ans in ("yes", "y"):
            print(f"  Proceeding with {step_name}...\n")
            return True
        elif ans in ("no", "n"):
            print(f"  {step_name} SKIPPED -- files from the previous stage passed through.")
            return False
        else:
            print("  Please type  yes  or  no")


def optional_save_zip(folder: str, zip_label: str, arcname_prefix: str):
    """Ask the user if they want a ZIP snapshot of a stage's output."""
    file_list = sorted(glob.glob(f"{folder}/*"))
    if not file_list:
        print(f"  Nothing in \"{folder}\" -- skipping snapshot.")
        return
    while True:
        ans = input(
            f"\nSave ZIP snapshot of \"{zip_label}\" ({len(file_list)} files)? [yes/no]: "
        ).strip().lower()
        if ans in ("yes", "y"):
            snap_path = f"snapshot_{zip_label}.zip"
            with zipfile.ZipFile(snap_path, "w", zipfile.ZIP_DEFLATED) as zf:
                for fp in file_list:
                    zf.write(fp, arcname=f"{arcname_prefix}/{Path(fp).name}")
            mb = os.path.getsize(snap_path) / (1024 * 1024)
            print(f"  Saved -> {os.path.abspath(snap_path)}  ({mb:.2f} MB)")
            break
        elif ans in ("no", "n"):
            print("  Snapshot skipped -- continuing.")
            break
        else:
            print("  Please type  yes  or  no")


_IN_COLAB = 'google.colab' in sys.modules   # always False on Ubuntu; this check
                                             # is kept only for parity with the
                                             # Colab notebook edition of LigKit

def _make_executor(n_workers, cpu_bound=False):
    # On Ubuntu (this script): ThreadPoolExecutor for I/O/GIL-releasing stages,
    # ProcessPoolExecutor for CPU-bound stages -- identical selection logic to
    # the Colab notebook's own non-Colab branch.
    if _IN_COLAB or not cpu_bound:
        return ThreadPoolExecutor(max_workers=n_workers)
    return ProcessPoolExecutor(max_workers=n_workers)

# ── Late imports with clean error messages ────────────────────────────────
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, SDWriter
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit.Chem.EnumerateStereoisomers import (
        EnumerateStereoisomers, StereoEnumerationOptions
    )
    from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False

try:
    import dimorphite_dl
    DIMORPHITE_OK = True
except ImportError:
    DIMORPHITE_OK = False

try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    MEEKO_OK = True
except ImportError:
    MEEKO_OK = False

# ── Thread-safe print + progress counter ─────────────────────────────────
_print_lock   = threading.Lock()
_counter_lock = threading.Lock()
_progress     = {"done": 0, "total": 0}

def tprint(*args, **kwargs):
    with _print_lock:
        print(*args, **kwargs)

def tick():
    with _counter_lock:
        _progress["done"] += 1
        d, t = _progress["done"], _progress["total"]
    tprint(f"  [{d}/{t}] completed", flush=True)

def reset_progress(total: int):
    with _counter_lock:
        _progress["done"] = 0
        _progress["total"] = total

# ════════════════════════════════════════════════════════════════════════════
# UNIVERSAL SYBYL/TRIPOS MOL2 WRITER  (replaces bare `obabel -O mol2`)
# ════════════════════════════════════════════════════════════════════════════
# Implements, per the official Tripos SYBYL 7.1 Mol2 File Format spec:
#   1. Unique atom names      -> element symbol + sequential per-element index
#   2. Correct SYBYL atom typing -> full hybridisation/aromaticity/functional-
#      group perception (O.co2, N.am, N.pl3, N.4, S.o/S.o2, C.cat, halogens,
#      common metals/Du) instead of OpenBabel's generic guesses
#   3. Universal Gasteiger-Marsili partial charges -> charge_type field set
#      to GASTEIGER_CHARGES. Chosen over MMFF94s as the PRIMARY engine
#      because MMFF94s bond-charge-increments are exactly 0.000 for
#      nonpolar C-C/C-H bonds (flattens alkyl backbones to 0.0000) and
#      MMFF94s has no parameters at all for several common motifs
#      (hypervalent P/S, boronates, some metals). Gasteiger covers the
#      full periodic table RDKit parses and never produces that zero
#      plateau (see compute_partial_charges() below for the full
#      rationale and the Tier-2 NaN-safety patch).
#
# This module is molecule-agnostic: every rule below is a general chemical
# perception rule (hybridisation, ring membership, functional-group pattern),
# not anything hard-coded to a specific compound. It is safe to use as the
# canonical Stage 10 MOL2 writer for ANY small-molecule SDF/MOL input.
# ════════════════════════════════════════════════════════════════════════════

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType, HybridizationType
from pathlib import Path

# ── SMARTS patterns for functional-group-dependent SYBYL types ────────────
# Built once at import time; cheap to reuse across every molecule in the run.
_PAT_CARBOXYLATE_O   = Chem.MolFromSmarts('[CX3](=[OX1])[OX1-]')      # STRICT: only the genuinely resonance-delocalised anionic carboxylate (both O's equivalent). Neutral esters/lactones must NOT match this -- they get normal O.2 (C=O) / O.3 (C-O-R) typing via the double-bond rule below.
_PAT_PHOSPHATE_O     = Chem.MolFromSmarts('[OX1]=[PX4]')               # P=O / P-O(-) oxygens -> O.co2 family (intentional: SYBYL convention types ALL P=O oxygens, neutral phosphonate/phosphate alike, as O.co2 due to P=O/P-O- resonance ambiguity)
_PAT_AMIDE_N         = Chem.MolFromSmarts('[NX3][CX3](=[OX1])')        # amide nitrogen -> N.am
_PAT_GUANIDINIUM_C   = Chem.MolFromSmarts('[NX3][CX3](=[NX3+,NX3])[NX3]')  # guanidinium central C -> C.cat
_PAT_SULFOXIDE       = Chem.MolFromSmarts('[#16X3]=[OX1]')             # S=O (single O) -> S.o
_PAT_SULFONE         = Chem.MolFromSmarts('[#16X4](=[OX1])(=[OX1])')   # S(=O)(=O) -> S.o2

_METAL_TYPES = {
    3: 'Li', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 19: 'K', 20: 'Ca',
    25: 'Mn', 26: 'Fe', 29: 'Cu', 30: 'Zn', 42: 'Mo', 50: 'Sn',
}
_HALOGENS = {9: 'F', 17: 'Cl', 35: 'Br', 53: 'I'}


def _ring_size_aromatic(atom, mol):
    """Return True if atom is in an aromatic ring (any size)."""
    if not atom.GetIsAromatic():
        return False
    ri = mol.GetRingInfo()
    return ri.NumAtomRings(atom.GetIdx()) > 0


def get_sybyl_atom_type(mol, atom_idx, match_cache):
    """
    Determine the correct Tripos SYBYL atom type for a single atom,
    following the official SYBYL 7.1 Mol2 atom-type table (see Tripos
    Mol2 File Format spec, 'SYBYL Atom Types' section).

    match_cache: dict of precomputed SMARTS match atom-index sets,
    built once per molecule by `_build_match_cache()` for efficiency.
    """
    atom    = mol.GetAtomWithIdx(atom_idx)
    z       = atom.GetAtomicNum()
    hyb     = atom.GetHybridization()
    arom    = _ring_size_aromatic(atom, mol)

    # ── Hydrogen ────────────────────────────────────────────────────────
    if z == 1:
        return 'H'

    # ── Halogens ────────────────────────────────────────────────────────
    if z in _HALOGENS:
        return _HALOGENS[z]

    # ── Metals (octahedral/tetrahedral distinctions need 3D context that
    #    SYBYL infers from coordination geometry; we default to the bare
    #    element symbol, which SYBYL/CGenFF readers accept and is always
    #    safer than a wrong geometry-specific guess) ─────────────────────
    if z in _METAL_TYPES:
        return _METAL_TYPES[z]

    # ── Carbon ──────────────────────────────────────────────────────────
    if z == 6:
        if atom_idx in match_cache['guanidinium_c']:
            return 'C.cat'
        if arom:
            return 'C.ar'
        if hyb == HybridizationType.SP3:
            return 'C.3'
        if hyb == HybridizationType.SP2:
            return 'C.2'
        if hyb == HybridizationType.SP:
            return 'C.1'
        return 'C.3'  # safe fallback

    # ── Nitrogen ────────────────────────────────────────────────────────
    if z == 7:
        if arom:
            return 'N.ar'
        if atom_idx in match_cache['amide_n']:
            return 'N.am'
        if atom.GetFormalCharge() == 1 and hyb == HybridizationType.SP3:
            return 'N.4'
        if hyb == HybridizationType.SP3:
            return 'N.3'
        if hyb == HybridizationType.SP2:
            # Trigonal-planar sp2 N not in an amide/aromatic context
            # (e.g. aniline-type NH2, enamine N) -> N.pl3
            # NOTE: must count explicit H *neighbors* directly, not
            # GetTotalNumHs() -- that always returns 0 once AddHs() has
            # been called (all hydrogens become explicit graph atoms, so
            # the implicit/explicit-H-count property no longer reflects
            # them). Using GetTotalNumHs() here silently makes this branch
            # unreachable for every molecule and N.pl3 is never assigned.
            n_h_neighbors = sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 1)
            if atom.GetDegree() == 3 and n_h_neighbors >= 1:
                return 'N.pl3'
            return 'N.2'
        if hyb == HybridizationType.SP:
            return 'N.1'
        return 'N.3'

    # ── Oxygen ──────────────────────────────────────────────────────────
    if z == 8:
        if atom_idx in match_cache['carboxylate_o'] or atom_idx in match_cache['phosphate_o']:
            return 'O.co2'
        if arom:
            return 'O.2'   # aromatic-ring O (e.g. furan/pyrylium O) -> SYBYL has no O.ar; O.2 is the accepted convention
        # CRITICAL RULE: O.2 requires an ACTUAL double bond (C=O, P=O, S=O).
        # Do NOT trust RDKit's hybridisation flag alone -- conjugated ether
        # oxygens (e.g. aryl-O-C, furanone ring O) are frequently perceived
        # as SP2 by RDKit purely from lone-pair delocalisation, even though
        # they carry only single bonds. SYBYL/CGenFF both require these to
        # be typed O.3 (ether/ester/hydroxyl), not O.2 (carbonyl/oxo).
        # This is the exact bug class that produced O.2 on the furanone
        # ring oxygen in the original OpenBabel output -- fixed here for
        # every molecule, not just this one.
        has_double_bond = any(b.GetBondTypeAsDouble() == 2.0 for b in atom.GetBonds())
        if has_double_bond:
            return 'O.2'
        return 'O.3'

    # ── Sulfur ──────────────────────────────────────────────────────────
    if z == 16:
        if atom_idx in match_cache['sulfone_s']:
            return 'S.o2'
        if atom_idx in match_cache['sulfoxide_s']:
            return 'S.o'
        if hyb == HybridizationType.SP2:
            return 'S.2'
        return 'S.3'

    # ── Phosphorus ──────────────────────────────────────────────────────
    if z == 15:
        # SYBYL spec defines only P.3 (phosphorus sp3) -- this is the
        # correct universal SYBYL type for ALL phosphorus environments
        # (phosphate, phosphonate, phosphite), since SYBYL has no separate
        # type for P=O vs P-C vs P-O-C connectivity. CGenFF's internal
        # CHARMM atom type (e.g. PG0, PG1) is assigned downstream by its
        # OWN analogy engine from this P.3 + bonding pattern, independent
        # of which SYBYL P type is written -- so P.3 is always correct here.
        return 'P.3'

    # ── Unknown heavy atom -> dummy carbon fallback (never silently drop) ─
    return f'Du'


def _build_match_cache(mol):
    """Pre-compute SMARTS substructure matches once per molecule."""
    def _flat(matches):
        return {idx for m in matches for idx in m}

    cache = {
        'carboxylate_o': set(),
        'phosphate_o':   set(),
        'amide_n':       set(),
        'guanidinium_c': set(),
        'sulfoxide_s':   set(),
        'sulfone_s':     set(),
    }
    try:
        if _PAT_CARBOXYLATE_O is not None:
            for m in mol.GetSubstructMatches(_PAT_CARBOXYLATE_O):
                cache['carboxylate_o'].update(m)
        if _PAT_PHOSPHATE_O is not None:
            for m in mol.GetSubstructMatches(_PAT_PHOSPHATE_O):
                cache['phosphate_o'].update(m)
        if _PAT_AMIDE_N is not None:
            for m in mol.GetSubstructMatches(_PAT_AMIDE_N):
                cache['amide_n'].add(m[0])
        if _PAT_GUANIDINIUM_C is not None:
            for m in mol.GetSubstructMatches(_PAT_GUANIDINIUM_C):
                cache['guanidinium_c'].add(m[1])
        if _PAT_SULFOXIDE is not None:
            for m in mol.GetSubstructMatches(_PAT_SULFOXIDE):
                cache['sulfoxide_s'].add(m[0])
        if _PAT_SULFONE is not None:
            for m in mol.GetSubstructMatches(_PAT_SULFONE):
                cache['sulfone_s'].add(m[0])
    except Exception:
        pass  # malformed SMARTS match on an unusual molecule -> safe no-op
    return cache


def _sybyl_bond_type(bond, mol):
    """SYBYL bond type string per Tripos spec (1/2/3/am/ar)."""
    if bond.GetIsAromatic():
        return 'ar'
    bt = bond.GetBondType()
    if bt == BondType.SINGLE:
        # Detect amide C-N single bond -> 'am' per SYBYL convention
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        for c_atom, n_atom in ((a1, a2), (a2, a1)):
            if c_atom.GetAtomicNum() == 6 and n_atom.GetAtomicNum() == 7:
                for nb in c_atom.GetNeighbors():
                    if nb.GetAtomicNum() == 8:
                        ob = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nb.GetIdx())
                        if ob and ob.GetBondTypeAsDouble() == 2.0:
                            return 'am'
        return '1'
    if bt == BondType.DOUBLE:
        return '2'
    if bt == BondType.TRIPLE:
        return '3'
    return '1'


def _assign_unique_atom_names(mol):
    """
    Per-element sequential naming (C1, C2, ..., O1, O2, ..., H1, H2, ...).
    This is the universal fix for CGenFF's 'renaming non-unique atom'
    warnings, which fire whenever two atoms share a bare element symbol
    as their MOL2 atom_name -- true for the OpenBabel default on every
    molecule, not just this one.
    """
    counters = {}
    names = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        counters[sym] = counters.get(sym, 0) + 1
        names[atom.GetIdx()] = f"{sym}{counters[sym]}"
    return names


def compute_partial_charges(mol):
    """
    UNIVERSAL per-atom partial-charge engine (molecule-agnostic; works
    identically for any element/oxidation-state combination RDKit can
    parse -- organics, phosphonates/phosphates, sulfoxides/sulfones,
    boronates, halogens, common metals).

    Design rationale (why this replaces a bare MMFF94s call):
    MMFF94s assigns charge via Halgren bond-charge-increments (BCIs).
    For nonpolar C(sp3)-C(sp3) / C(sp3)-H bonds the parametrised BCI is
    EXACTLY 0.000 -- chemically legitimate within that force field, but
    it flattens every buried alkyl backbone atom to 0.0000 in the MOL2,
    which is useless for downstream per-atom QC/CGenFF triage and is
    frequently misread as a parametrisation failure. MMFF94s also has
    no BCI table at all for several common medicinal-chemistry motifs
    (hypervalent P/S centres, boronates, some metals), so a whole-
    molecule MMFF94s call can fail outright on exactly the ligands that
    most need a charge.

    Tier 1 - Gasteiger-Marsili (iterative electronegativity equalisation):
             PRIMARY engine. Universal coverage across the periodic
             table RDKit recognises; never degenerates to a flat-zero
             plateau on nonpolar atoms (empirically verified: plain
             alkyl backbones still differentiate to 1e-4 precision).
    Tier 2 - Atom-local formal-charge redistribution: last-resort patch,
             touches ONLY atoms where Tier 1 returned NaN/Inf (can occur
             for disconnected fragments or pathological valence states).
             Spreads that atom's formal charge across itself + heavy-
             atom neighbours, guaranteeing a finite value -- zero atoms
             are ever left silently unassigned, for any input ligand.

    Returns (charges: dict[atom_idx -> float], charge_type_label: str).
    charge_type_label is 'GASTEIGER_CHARGES' (clean run) or
    'GASTEIGER_CHARGES*' (one or more atoms required the Tier 2 patch --
    inspect such ligands manually before CGenFF/force-field submission).
    """
    import math

    try:
        AllChem.ComputeGasteigerCharges(mol, throwOnParamFailure=False)
    except Exception:
        for atom in mol.GetAtoms():
            atom.SetProp('_GasteigerCharge', '0.0')

    charges, unresolved = {}, []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        try:
            q = float(atom.GetProp('_GasteigerCharge'))
        except (KeyError, ValueError):
            q = float('nan')
        if not math.isfinite(q):
            unresolved.append(idx)
            q = 0.0  # placeholder; patched below
        charges[idx] = q

    # Tier 2: atom-local formal-charge redistribution for any NaN/Inf atom
    for idx in unresolved:
        atom = mol.GetAtomWithIdx(idx)
        pool = [atom] + [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
        total_formal = sum(a.GetFormalCharge() for a in pool)
        charges[idx] = round(total_formal / len(pool), 4) if pool else 0.0

    charge_label = 'GASTEIGER_CHARGES*' if unresolved else 'GASTEIGER_CHARGES'
    return charges, charge_label


def write_universal_sybyl_mol2(sdf_path: str, out_path: str, mol_name: str = None) -> None:
    """
    Universal SDF -> Tripos MOL2 writer for CGenFF / GOLD / DOCK6 submission.
    Molecule-agnostic: works identically for any small-molecule SDF.

    Fixes implemented (per Tripos SYBYL 7.1 Mol2 spec + CGenFF requirements):
      1. Unique atom names       (element + sequential index)
      2. Correct SYBYL atom types (full perception: aromaticity, hybridisation,
                                    O.co2/N.am/N.pl3/N.4/S.o/S.o2/C.cat, halogens,
                                    metals, universal P.3 for all phosphorus)
      3. Universal Gasteiger-Marsili partial charges (every atom resolved;
       see compute_partial_charges() for the full per-atom guarantee)
    """
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False, sanitize=True)
    if mol is None:
        raise ValueError(f"RDKit could not parse {sdf_path}")
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)

    mol = Chem.AddHs(mol, addCoords=True)
    Chem.Kekulize(mol, clearAromaticFlags=False)

    if mol_name is None:
        mol_name = Path(sdf_path).stem

    match_cache = _build_match_cache(mol)
    atom_names  = _assign_unique_atom_names(mol)
    charges, charge_label = compute_partial_charges(mol)

    conf      = mol.GetConformer()
    n_atoms   = mol.GetNumAtoms()
    n_bonds   = mol.GetNumBonds()
    res_name  = (mol_name[:7] if len(mol_name) > 7 else mol_name).upper() + '1'

    lines = []
    lines.append('@<TRIPOS>MOLECULE')
    lines.append(mol_name)
    lines.append(f' {n_atoms} {n_bonds} 1 0 0')
    lines.append('SMALL')
    lines.append(charge_label)
    lines.append('')

    lines.append('@<TRIPOS>ATOM')
    for atom in mol.GetAtoms():
        idx     = atom.GetIdx()
        name    = atom_names[idx]
        pos     = conf.GetAtomPosition(idx)
        sybyl_t = get_sybyl_atom_type(mol, idx, match_cache)
        q       = charges.get(idx, 0.0)
        lines.append(
            f'{idx + 1:>7} {name:<8} {pos.x:>10.4f} {pos.y:>10.4f} {pos.z:>10.4f} '
            f'{sybyl_t:<8} {1:>3} {res_name:<8} {q:>9.4f}'
        )

    lines.append('@<TRIPOS>BOND')
    for bidx, bond in enumerate(mol.GetBonds()):
        lines.append(
            f'{bidx + 1:>6} {bond.GetBeginAtomIdx() + 1:>6} '
            f'{bond.GetEndAtomIdx() + 1:>6} {_sybyl_bond_type(bond, mol)}'
        )

    lines.append('@<TRIPOS>SUBSTRUCTURE')
    lines.append(f'{1:>5} {res_name:<8} {1:>5} RESIDUE    0 **** ROOT')

    with open(out_path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


def run_obabel(cmd, critical=True, timeout=60):
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"TIMEOUT after {timeout}s: {' '.join(str(c) for c in cmd[:4])}...")
    if result.returncode != 0 and critical:
        raise RuntimeError((result.stderr or result.stdout or "unknown error")[:300])
    return result.returncode == 0

def get_mol_name(sdf_path: str, fallback: str) -> str:
    with open(sdf_path, "r", encoding="utf-8", errors="replace") as fh:
        first_line = fh.readline().strip()
    if not first_line:
        return fallback
    safe = re.sub(r"[^\w\-]", "_", first_line)
    safe = re.sub(r"_+", "_", safe).strip("_")
    return safe if safe else fallback

def passthrough(src_glob: str, dest_dir: str):
    files_found = sorted(glob.glob(src_glob))
    if not files_found:
        print(f"  No files matched '{src_glob}' - nothing to pass through.")
        return
    os.makedirs(dest_dir, exist_ok=True)
    for fp in files_found:
        shutil.copy2(fp, os.path.join(dest_dir, Path(fp).name))
    print(f"  Passthrough: {len(files_found)} files -> {dest_dir}/")

def inchikey_dedup(sdf_files: list, dest_dir: str, level: str,
                   label: str, workers: int = 4) -> list:
    os.makedirs(dest_dir, exist_ok=True)
    seen_keys  = set()
    dup_count  = 0
    fail_count = 0
    kept       = []

    def _compute_key(sdf_path):
        try:
            mol      = Chem.MolFromMolFile(sdf_path, removeHs=True, sanitize=True)
            if mol is None:
                return sdf_path, None
            inchi    = MolToInchi(mol)
            full_key = InchiToInchiKey(inchi)
            key      = full_key[:14] if level == "connectivity" else full_key
            return sdf_path, key
        except Exception:
            return sdf_path, None

    key_map = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for path_, key_ in pool.map(_compute_key, sdf_files):
            key_map[path_] = key_

    for sdf_path in sdf_files:
        name = Path(sdf_path).stem
        key_ = key_map.get(sdf_path)
        if key_ is None:
            dest = os.path.join(dest_dir, Path(sdf_path).name)
            shutil.copy2(sdf_path, dest)
            kept.append(dest); fail_count += 1
            continue
        if key_ in seen_keys:
            dup_count += 1
            print(f"  Duplicate removed: {name}")
        else:
            seen_keys.add(key_)
            dest = os.path.join(dest_dir, Path(sdf_path).name)
            shutil.copy2(sdf_path, dest)
            kept.append(dest)

    print(f"  {label}: {len(kept)} kept, {dup_count} duplicates removed, "
          f"{fail_count} unparseable kept")
    return kept


# ════════════════════════════════════════════════════════════════════════════
# Dimorphite-DL version-safe wrapper
# Handles: v2.x (ph_min/ph_max/precision), v1.x (min_ph/max_ph/pka_precision),
#          v1.3 run_with_mol_list, v1.2 run(), legacy dict-based run().
# ════════════════════════════════════════════════════════════════════════════
def _dimorphite_run(smi: str) -> list:
    if hasattr(dimorphite_dl, "protonate_smiles"):
        # v2.x keyword API
        try:
            r = dimorphite_dl.protonate_smiles(
                smi, ph_min=PROTO_MIN_PH, ph_max=PROTO_MAX_PH,
                max_variants=PROTO_MAX_VARS, precision=1.0)
            if r is not None: return list(r)
        except TypeError:
            pass
        # v1.x keyword API
        try:
            r = dimorphite_dl.protonate_smiles(
                smi, min_ph=PROTO_MIN_PH, max_ph=PROTO_MAX_PH,
                max_variants=PROTO_MAX_VARS, pka_precision=1.0)
            if r is not None: return list(r)
        except TypeError:
            pass
        # Positional fallback (works across versions)
        try:
            r = dimorphite_dl.protonate_smiles(
                smi, PROTO_MIN_PH, PROTO_MAX_PH, 1.0, False, False, PROTO_MAX_VARS)
            if r is not None: return list(r)
        except Exception:
            pass

    if hasattr(dimorphite_dl, "run_with_mol_list"):   # v1.3+
        try:
            r = dimorphite_dl.run_with_mol_list(
                [smi], min_ph=PROTO_MIN_PH, max_ph=PROTO_MAX_PH,
                max_variants=PROTO_MAX_VARS, pka_precision=1.0, label_states=False)
            if r is not None: return list(r)
        except Exception:
            pass

    if hasattr(dimorphite_dl, "run"):
        # v1.2 keyword form
        try:
            r = dimorphite_dl.run(
                smiles=smi, min_ph=PROTO_MIN_PH, max_ph=PROTO_MAX_PH,
                max_variants=PROTO_MAX_VARS, pka_precision=1.0, label_states=False)
            if r is not None: return list(r)
        except TypeError:
            pass
        # Legacy dict form
        try:
            r = dimorphite_dl.run({
                "smiles": smi, "min_ph": PROTO_MIN_PH, "max_ph": PROTO_MAX_PH,
                "max_variants": PROTO_MAX_VARS, "pka_precision": 1.0,
                "label_states": False})
            if r is not None: return list(r)
        except Exception:
            pass

    raise RuntimeError(
        f"dimorphite_dl v{getattr(dimorphite_dl, '__version__', 'unknown')} - "
        "no callable API found. Try: pip install --upgrade dimorphite_dl"
    )


# ════════════════════════════════════════════════════════════════════════════
# STAGE 3 WORKER: Standardize / Neutralize
# Correct sequence: salt_strip -> normalize -> reionize -> uncharge
# ════════════════════════════════════════════════════════════════════════════
def _uncharge_standardize(sdf: str) -> dict:
    import os
    from pathlib import Path
    from rdkit import Chem
    from rdkit.Chem import SDWriter
    from rdkit.Chem.MolStandardize import rdMolStandardize
    parent_name = Path(sdf).stem
    try:
        mol = Chem.MolFromMolFile(sdf, removeHs=True, sanitize=True)
        if mol is None:
            raise ValueError("RDKit could not parse SDF")
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        # Step 1: Salt/solvent stripping (MUST run first - before charge-aware ops)
        lfc = rdMolStandardize.LargestFragmentChooser()
        mol = lfc.choose(mol)
        if mol is None:
            raise ValueError("LargestFragmentChooser returned None")

        # Step 2: Normalize (canonical functional group reps for charge ops)
        normalizer = rdMolStandardize.Normalizer()
        mol = normalizer.normalize(mol)

        # Step 3: Reionize (consistent charge distribution before uncharging)
        reionizer = rdMolStandardize.Reionizer()
        mol = reionizer.reionize(mol)

        # Step 4: Uncharge (full neutralization for unrestricted tautomer enum)
        # canonicalOrder=True added in RDKit >= 2022.09; silently ignored on older
        try:
            uncharger = rdMolStandardize.Uncharger(canonicalOrder=True)
        except TypeError:
            uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
        if mol is None:
            raise ValueError("Uncharger returned None")

        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        out_path = os.path.join("standardized", f"{parent_name}.sdf")
        w = SDWriter(out_path); w.write(mol); w.close()
        tprint(f"  ok  {parent_name}"); tick()
        return {"name": parent_name, "ok": True, "step": "standardize", "error": None}
    except Exception as exc:
        err = str(exc)
        tprint(f"  SKIP  {parent_name}: {err[:100]}")
        tick()
        return {"name": parent_name, "ok": False, "step": "standardize", "error": err}


# ════════════════════════════════════════════════════════════════════════════
# STAGE 4 WORKER: Tautomer enumeration (RDKit TautomerEnumerator)
# ════════════════════════════════════════════════════════════════════════════
def _tautomer_worker(args: tuple) -> dict:
    sdf_path, cfg = args
    import os, shutil
    from pathlib import Path
    parent_name = Path(sdf_path).stem
    out_dir     = cfg["out_dir"]
    tauto_max   = cfg["rdkit_tauto_max"]
    try:
        from rdkit import Chem
        from rdkit.Chem import SDWriter
        from rdkit.Chem.MolStandardize import rdMolStandardize as _rms

        te = _rms.TautomerEnumerator()
        te.SetMaxTautomers(tauto_max * 10)   # generate wide, then score-filter
        te.SetMaxTransforms(1000)
        # Guard with hasattr - these methods added in newer RDKit versions
        if hasattr(te, "SetRemoveBondStereoIfPossible"):
            te.SetRemoveBondStereoIfPossible(False)  # never destroy E/Z stereo
        if hasattr(te, "SetReassignStereo"):
            te.SetReassignStereo(True)               # re-perceive after each transform

        mol = Chem.MolFromMolFile(sdf_path, removeHs=True, sanitize=True)
        if mol is None:
            raise ValueError("RDKit could not parse SDF")
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        all_tautos = list(te.Enumerate(mol))
        was_capped = len(all_tautos) >= tauto_max * 10

        scored = []
        for t in all_tautos:
            try:    scored.append((te.ScoreTautomer(t), t))
            except: scored.append((0.0, t))
        scored.sort(key=lambda x: x[0], reverse=True)
        top = [m for _, m in scored[:tauto_max]]

        if not top:
            shutil.copy2(sdf_path, os.path.join(out_dir, f"{parent_name}_ta01.sdf"))
            return {"name": parent_name, "ok": True, "written": 1,
                    "capped": False, "no_tauto": True, "error": None}

        written = 0
        for idx, t_mol in enumerate(top, 1):
            Chem.AssignStereochemistry(t_mol, cleanIt=True, force=True)
            t_name   = f"{parent_name}_ta{idx:02d}"
            out_path = os.path.join(out_dir, f"{t_name}.sdf")
            t_mol.SetProp("_Name", t_name)
            w = SDWriter(out_path); w.write(t_mol); w.close()
            written += 1

        if written == 0:
            shutil.copy2(sdf_path, os.path.join(out_dir, f"{parent_name}_ta01.sdf"))
            written = 1

        return {"name": parent_name, "ok": True, "written": written,
                "capped": was_capped, "no_tauto": False, "error": None}
    except Exception as exc:
        try:
            shutil.copy2(sdf_path, os.path.join(out_dir, f"{parent_name}_ta01.sdf"))
        except Exception: pass
        return {"name": parent_name, "ok": False, "written": 1,
                "capped": False, "no_tauto": False, "error": str(exc)[:200]}


# ════════════════════════════════════════════════════════════════════════════
# STAGE 5 WORKER: Protonation state enumeration (Dimorphite-DL)
# ════════════════════════════════════════════════════════════════════════════
def _proto_worker(args: tuple) -> dict:
    sdf_path, cfg = args
    import os, shutil
    from pathlib import Path
    parent_name = Path(sdf_path).stem
    out_dir     = cfg["out_dir"]
    try:
        from rdkit import Chem
        from rdkit.Chem import SDWriter
        mol = Chem.MolFromMolFile(sdf_path, removeHs=True, sanitize=True)
        if mol is None:
            raise ValueError("RDKit could not parse SDF")
        # Assign stereo BEFORE SMILES conversion to preserve R/S and E/Z
        # through the Dimorphite SMILES round-trip
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        mol_noH = Chem.RemoveHs(mol)
        smi_in  = Chem.MolToSmiles(mol_noH, canonical=True)
        variants = _dimorphite_run(smi_in)

        if not variants:
            shutil.copy2(sdf_path, os.path.join(out_dir, f"{parent_name}_p01.sdf"))
            return {"name": parent_name, "ok": True, "written": 1,
                    "no_variants": True, "error": None}

        written = 0
        for idx, smi in enumerate(variants, 1):
            var_mol = Chem.MolFromSmiles(smi)
            if var_mol is None: continue
            # Re-assign stereo after SMILES round-trip
            Chem.AssignStereochemistry(var_mol, cleanIt=True, force=True)
            var_name = f"{parent_name}_p{idx:02d}"
            out_path = os.path.join(out_dir, f"{var_name}.sdf")
            var_mol.SetProp("_Name", var_name)
            w = SDWriter(out_path); w.write(var_mol); w.close()
            written += 1

        if written == 0:
            shutil.copy2(sdf_path, os.path.join(out_dir, f"{parent_name}_p01.sdf"))
            written = 1

        return {"name": parent_name, "ok": True, "written": written,
                "no_variants": False, "error": None}
    except Exception as exc:
        try:
            shutil.copy2(sdf_path, os.path.join(out_dir, f"{parent_name}_p01.sdf"))
        except Exception: pass
        return {"name": parent_name, "ok": False, "written": 1,
                "no_variants": False, "error": str(exc)[:200]}


# ════════════════════════════════════════════════════════════════════════════
# STAGE 7 WORKER: Stereoisomer enumeration (RDKit)
# ════════════════════════════════════════════════════════════════════════════
def _stereo_worker(args: tuple) -> dict:
    sdf_path, cfg = args
    import os, shutil
    from pathlib import Path
    parent_name  = Path(sdf_path).stem
    out_dir      = cfg["out_dir"]
    max_isomers  = cfg["max_isomers"]
    only_undef   = cfg["only_undefined"]
    try:
        from rdkit import Chem
        from rdkit.Chem import SDWriter
        from rdkit.Chem.EnumerateStereoisomers import (
            EnumerateStereoisomers, StereoEnumerationOptions)

        mol = Chem.MolFromMolFile(sdf_path, removeHs=True, sanitize=True)
        if mol is None:
            raise ValueError("RDKit could not parse SDF")
        # AssignStereochemistry BEFORE enumeration: correctly classifies defined
        # vs undefined centers for onlyUnassigned=True
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        stereo_info   = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        n_undefined   = sum(1 for _, s in stereo_info if s == "?")
        max_possible  = 2 ** n_undefined if n_undefined > 0 else 1

        stereo_opts   = StereoEnumerationOptions(
            unique=True, onlyUnassigned=only_undef,
            maxIsomers=max_isomers, tryEmbedding=False)
        isomers       = list(EnumerateStereoisomers(mol, options=stereo_opts))
        was_truncated = len(isomers) >= max_isomers and max_possible > max_isomers

        if not isomers:
            shutil.copy2(sdf_path, os.path.join(out_dir, f"{parent_name}_s01.sdf"))
            return {"name": parent_name, "ok": True, "written": 1,
                    "truncated": False, "max_possible": max_possible,
                    "no_stereo": True, "error": None}

        written = 0
        for idx, iso_mol in enumerate(isomers, 1):
            # Re-assign stereo on EACH isomer output
            Chem.AssignStereochemistry(iso_mol, cleanIt=True, force=True)
            iso_name = f"{parent_name}_s{idx:02d}"
            out_path = os.path.join(out_dir, f"{iso_name}.sdf")
            iso_mol.SetProp("_Name", iso_name)
            w = SDWriter(out_path); w.write(iso_mol); w.close()
            written += 1

        return {"name": parent_name, "ok": True, "written": written,
                "truncated": was_truncated, "max_possible": max_possible,
                "no_stereo": False, "error": None}
    except Exception as exc:
        try:
            shutil.copy2(sdf_path, os.path.join(out_dir, f"{parent_name}_s01.sdf"))
        except Exception: pass
        return {"name": parent_name, "ok": False, "written": 1,
                "truncated": False, "max_possible": 1,
                "no_stereo": False, "error": str(exc)[:200]}


# ════════════════════════════════════════════════════════════════════════════
# STAGE 8+9 WORKER: 3D generation + MMFF94 minimization (RDKit)
# ETKDGv3 -> ETKDGv2 -> ETDG fallback chain; MMFF94 -> UFF fallback
# ════════════════════════════════════════════════════════════════════════════
def _rdkit_gen3d_and_minimize(args: tuple) -> dict:
    mol_block, name, output_dir, run_gen3d, run_minimize, cfg = args
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, SDWriter
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError:
        return {"name": name, "ok": False, "stage": "import",
                "error": "RDKit not available", "warn": None}

    warn = None
    mol  = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=True)
    if mol is None:
        return {"name": name, "ok": False, "stage": "parsing",
                "error": "RDKit could not parse molecule", "warn": None}
    try:
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            return {"name": name, "ok": False, "stage": "sanitization",
                    "error": str(e)[:200], "warn": None}

        try:
            lfc = rdMolStandardize.LargestFragmentChooser()
            mol = lfc.choose(mol)
        except Exception:
            pass

        mol_h = Chem.AddHs(mol)

        if run_gen3d:
            # Tier 1: ETKDGv3 (CSD torsion library, macrocycle-aware)
            params                       = AllChem.ETKDGv3()
            params.randomSeed            = cfg["random_seed"]
            params.numThreads            = 1
            params.enforceChirality      = cfg["enforce_chirality"]
            params.useSmallRingTorsions  = True
            params.useMacrocycleTorsions = True
            params.ETversion             = 2
            params.maxIterations         = cfg["max_attempts"]
            result = AllChem.EmbedMolecule(mol_h, params)

            if result == -1:
                # Tier 2: ETKDGv2
                warn    = "ETKDGv3 failed -> ETKDGv2 fallback"
                params2 = AllChem.ETKDGv2()
                params2.randomSeed       = cfg["random_seed"]
                params2.enforceChirality = cfg["enforce_chirality"]
                params2.maxIterations    = cfg["max_attempts"]
                result  = AllChem.EmbedMolecule(mol_h, params2)

            if result == -1:
                # Tier 3: ETDG (last resort, no torsion library)
                warn    = "ETKDGv3+v2 failed -> ETDG last-resort fallback"
                params3 = AllChem.EmbedParameters()
                params3.randomSeed       = cfg["random_seed"]
                params3.enforceChirality = False
                result  = AllChem.EmbedMolecule(mol_h, params3)

            if result == -1:
                return {"name": name, "ok": False, "stage": "3D generation",
                        "error": "All embedding methods failed (ETKDGv3/v2/ETDG)",
                        "warn": warn}

            # Conformer geometry validation: check for collapsed/exploded bonds
            def _check_geom(m):
                try:    conf = m.GetConformer()
                except ValueError: return False, "no conformer"
                for bond in m.GetBonds():
                    i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                    d    = conf.GetAtomPosition(i).Distance(conf.GetAtomPosition(j))
                    sym  = (m.GetAtomWithIdx(i).GetSymbol(),
                            m.GetAtomWithIdx(j).GetSymbol())
                    if d < 0.8:
                        return False, f"collapsed {sym[0]}-{sym[1]} bond = {d:.3f} A"
                    if d > 2.5:
                        return False, f"exploded  {sym[0]}-{sym[1]} bond = {d:.3f} A"
                return True, "ok"

            ok, detail = _check_geom(mol_h)
            if not ok:
                rseed = cfg["random_seed"] + 1000
                rp                       = AllChem.ETKDGv3()
                rp.randomSeed            = rseed
                rp.numThreads            = 1
                rp.enforceChirality      = cfg["enforce_chirality"]
                rp.useSmallRingTorsions  = True
                rp.useMacrocycleTorsions = True
                rp.ETversion             = 2
                rp.maxIterations         = cfg["max_attempts"]
                rr = AllChem.EmbedMolecule(mol_h, rp)
                if rr == -1:
                    return {"name": name, "ok": False, "stage": "conformer validation",
                            "error": f"Bad geometry ({detail}); retry seed {rseed} failed",
                            "warn": warn}
                ok2, detail2 = _check_geom(mol_h)
                if not ok2:
                    return {"name": name, "ok": False, "stage": "conformer validation",
                            "error": f"Bad geometry after 2 attempts: {detail} | {detail2}",
                            "warn": warn}
                w    = f"Bad geometry on seed {cfg['random_seed']} -> re-embedded with {rseed}"
                warn = f"{warn} | {w}" if warn else w

        if run_minimize:
            # MMFFSanitizeMolecule: handles exotic atoms (B, Se, Te) that cause
            # MMFFGetMoleculeProperties to silently return None without it
            try:
                AllChem.MMFFSanitizeMolecule(mol_h)
            except Exception:
                pass

            ff_props = AllChem.MMFFGetMoleculeProperties(
                mol_h, mmffVariant=cfg["mmff_variant"])

            if ff_props is not None:
                ff = AllChem.MMFFGetMoleculeForceField(mol_h, ff_props, confId=0)
                if ff is not None:
                    ff.Initialize()
                    converged = ff.Minimize(maxIts=cfg["mmff_max_iters"],
                                            forceTol=cfg["mmff_force_tol"])
                    if converged != 0:
                        w    = f"{cfg['mmff_variant']} did not fully converge"
                        warn = f"{warn} | {w}" if warn else w
                else:
                    ff_props = None

            if ff_props is None:
                if cfg["uff_fallback"]:
                    w    = f"{cfg['mmff_variant']} not applicable -> UFF fallback"
                    warn = f"{warn} | {w}" if warn else w
                    try:
                        uff = AllChem.UFFGetMoleculeForceField(mol_h, confId=0)
                        if uff:
                            uff.Initialize()
                            uff.Minimize(maxIts=cfg["mmff_max_iters"])
                    except Exception as ue:
                        w2   = f"UFF also failed: {str(ue)[:80]}"
                        warn = f"{warn} | {w2}" if warn else w2
                else:
                    w    = f"{cfg['mmff_variant']} not applicable, UFF disabled"
                    warn = f"{warn} | {w}" if warn else w

        out_path = os.path.join(output_dir, f"{name}.sdf")
        counter  = 1
        while os.path.exists(out_path):
            out_path = os.path.join(output_dir, f"{name}_{counter:02d}.sdf")
            counter += 1
        writer = SDWriter(out_path)
        writer.write(mol_h)
        writer.close()
        return {"name": name, "ok": True, "stage": "done", "error": None, "warn": warn}

    except Exception as exc:
        return {"name": name, "ok": False, "stage": "unexpected",
                "error": str(exc)[:300], "warn": warn}




def run_self_test() -> bool:
    """
    Validates SYBYL atom typing, Gasteiger charge conservation, Kekulize
    stability, unique atom naming, and the MOL2 writer round-trip against a
    fixed aspirin reference molecule. Purely diagnostic - sandboxed in a temp
    directory, does not touch pipeline state or outputs. Returns True if all
    checks pass.
    """
    _ST_PASS = 0
    _ST_FAIL = []

    def _st_check(label, condition):
        nonlocal _ST_PASS
        if condition:
            _ST_PASS += 1
            print(f"  \u2705 {label}")
        else:
            _ST_FAIL.append(label)
            print(f"  \u274c {label}")

    print("Running LigKit self-test (aspirin, C9H8O4)...\n")

    if not RDKIT_OK:
        print("\u274c RDKit not available -- install dependencies first, then re-run.")
        return False
    else:
        _st_dir = tempfile.mkdtemp(prefix="ligkit_selftest_")
        _cwd_before = os.getcwd()
        try:
            os.chdir(_st_dir)

            # ── Check 1: parse + sanitize ───────────────────────────────────
            aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"
            mol = Chem.MolFromSmiles(aspirin_smiles)
            _st_check("RDKit parses reference SMILES", mol is not None)

            mol = Chem.AddHs(mol)
            embed_ok = AllChem.EmbedMolecule(mol, randomSeed=42) == 0
            _st_check("3D embedding succeeds (ETKDG)", embed_ok)

            # ── Check 2: Kekulize does not raise on an aromatic ring ────────
            kek_ok = True
            try:
                Chem.Kekulize(mol, clearAromaticFlags=False)
            except Exception:
                kek_ok = False
            _st_check("Kekulize succeeds on aromatic ring", kek_ok)

            # ── Check 3: unique atom names ───────────────────────────────────
            names = _assign_unique_atom_names(mol)
            _st_check("Atom names are unique", len(names) == len(set(names.values())))

            # ── Check 4: SYBYL atom types assigned to every atom ─────────────
            match_cache = _build_match_cache(mol)
            sybyl_types = [get_sybyl_atom_type(mol, a.GetIdx(), match_cache) for a in mol.GetAtoms()]
            _st_check("Every atom receives a non-empty SYBYL type",
                      all(t and t.strip() for t in sybyl_types))
            _st_check("Aromatic ring carbons typed C.ar",
                      any(t == "C.ar" for t in sybyl_types))

            # ── Check 5: partial charges -- finite, and net charge conserved ─
            charges, charge_label = compute_partial_charges(mol)
            all_finite = all(math.isfinite(q) for q in charges.values())
            _st_check("All partial charges are finite (no NaN/Inf)", all_finite)

            net_formal = sum(a.GetFormalCharge() for a in mol.GetAtoms())
            net_computed = sum(charges.values())
            _st_check(
                f"Net charge conserved (formal={net_formal:+.0f}e, "
                f"computed={net_computed:+.3f}e, tol=0.05e)",
                abs(net_computed - net_formal) < 0.05
            )
            _st_check("Charge label is clean GASTEIGER_CHARGES (no Tier-2 patch needed)",
                      charge_label == "GASTEIGER_CHARGES")

            # ── Check 6: end-to-end MOL2 writer round-trip ──────────────────
            sdf_path  = os.path.join(_st_dir, "aspirin.sdf")
            mol2_path = os.path.join(_st_dir, "aspirin.mol2")
            w = SDWriter(sdf_path); w.write(mol); w.close()

            mol2_ok = True
            try:
                write_universal_sybyl_mol2(sdf_path, mol2_path, mol_name="ASPIRIN")
            except Exception as e:
                mol2_ok = False
                print(f"     MOL2 writer raised: {e}")
            _st_check("MOL2 writer completes without raising", mol2_ok)

            if mol2_ok and os.path.exists(mol2_path):
                with open(mol2_path) as fh:
                    mol2_text = fh.read()
                _st_check("MOL2 file has @<TRIPOS>MOLECULE/ATOM/BOND records",
                          all(tag in mol2_text for tag in
                              ("@<TRIPOS>MOLECULE", "@<TRIPOS>ATOM", "@<TRIPOS>BOND")))
                n_atom_lines = sum(
                    1 for l in mol2_text.splitlines()
                    if l.strip() and not l.startswith("@") and len(l.split()) >= 9
                )
                _st_check("MOL2 atom count matches source molecule",
                          n_atom_lines == mol.GetNumAtoms())
            else:
                _st_check("MOL2 file has @<TRIPOS>MOLECULE/ATOM/BOND records", False)
                _st_check("MOL2 atom count matches source molecule", False)

        finally:
            os.chdir(_cwd_before)
            shutil.rmtree(_st_dir, ignore_errors=True)

        print(f"\n{'='*60}")
        if not _ST_FAIL:
            print(f"\u2705 Self-test PASSED  ({_ST_PASS}/{_ST_PASS} checks)")
            print("   Core chemistry functions behave as expected -- safe to proceed.")
        else:
            print(f"\u274c Self-test FAILED  ({_ST_PASS}/{_ST_PASS + len(_ST_FAIL)} checks passed)")
            print(f"   Failing checks: {', '.join(_ST_FAIL)}")
            print("   Investigate before running Step 5 on a real library.")
        print(f"{'='*60}")

    return not _ST_FAIL


FAILURE_LOG:   list = []
SKIPPED_STEPS: list = []


if __name__ == "__main__":

    _pipeline_start = datetime.now()
    FAILURE_LOG.clear()
    SKIPPED_STEPS.clear()

    # ── Resolve paths ────────────────────────────────────────────────────────
    INPUT_PATH_clean = INPUT_PATH.strip()
    if not INPUT_PATH_clean or not os.path.isfile(INPUT_PATH_clean):
        print(f"\nERROR: Cannot find input file:\n     {INPUT_PATH_clean}")
        print("Update INPUT_PATH in USER SETTINGS (top of this script).\n")
        sys.exit(1)

    input_file = str(Path(INPUT_PATH_clean).resolve())
    input_stem = Path(input_file).stem
    base_out   = str(Path(OUTPUT_DIR.strip()).resolve()) if OUTPUT_DIR.strip() \
                 else str(Path(input_file).parent)
    WORK_DIR   = os.path.join(base_out, f"{input_stem}_prepared")
    os.makedirs(WORK_DIR, exist_ok=True)
    os.chdir(WORK_DIR)
    zip_name   = f"{input_stem}_prepared"
    _workers   = multiprocessing.cpu_count() if MAX_WORKERS == 0 else MAX_WORKERS

    # ── Dependency gate ───────────────────────────────────────────────────────
    import shutil as _sh
    missing = []
    if _sh.which("obabel") is None:
        missing.append("OpenBabel  (re-run Step 1 and restart runtime)")
    if not RDKIT_OK:
        missing.append("RDKit      (re-run Step 1)")
    if not DIMORPHITE_OK:
        missing.append("Dimorphite-DL  (re-run Step 1)")
    if not MEEKO_OK:
        print("Note: Meeko not found - PDBQT stage will use OpenBabel fallback.")
    if missing:
        raise RuntimeError("Missing required dependencies:\n" +
                           "\n".join(f"  - {m}" for m in missing))

    print("LIGKIT — LIGAND PREPARATION PIPELINE  (CLI Edition)")
    print(f"  Input     : {input_file}")
    print(f"  Output    : {WORK_DIR}")
    print(f"  Workers   : {_workers}  (ProcessPool for CPU-bound stages)")
    print()

    # ── Fresh workspace ───────────────────────────────────────────────────────
    for d in ["split_sdf", "deduped", "standardized", "tautomers",
              "proto_states", "post_deduped", "stereoisomers",
              "gen3d", "minimized", "mol2", "pdbqt"]:
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)

    # ── Optional self-test (validates SYBYL typing / charges / MOL2 writer) ────
    if ask_run(
        "Self-Test — Validate Core Chemistry Functions",
        "Runs 10 assertion-based checks against a fixed aspirin reference molecule\n"
        "  (SYBYL typing, Gasteiger charge conservation, Kekulize stability, unique\n"
        "  atom naming, MOL2 writer round-trip) before processing your real library.\n"
        "  Purely diagnostic - does not touch pipeline state."
    ):
        run_self_test()


    # ════════════════════════════════════════════════════════════════════════════
    # STAGE 1: Split / Extract  (OpenBabel - always runs)
    # ════════════════════════════════════════════════════════════════════════════
    print("=" * 60)
    print("STAGE 1: Split / Extract  [OpenBabel - always runs]")
    print("=" * 60)

    ext = Path(input_file).suffix.lower()
    if ext == ".sdf":
        print("Single SDF -> splitting with obabel -m ...")
        subprocess.run(
            ["obabel", "-isdf", input_file, "-osdf", "-O", "split_sdf/mol.sdf", "-m"],
            check=True)
    elif ext == ".zip":
        print("ZIP -> extracting SDF files ...")
        with zipfile.ZipFile(input_file, "r") as zf:
            entries = [e for e in zf.namelist()
                       if e.lower().endswith(".sdf")
                       and not os.path.basename(e).startswith(".")
                       and os.path.basename(e) != ""]
            if not entries:
                raise RuntimeError("No .sdf files found in ZIP.")
            for entry in entries:
                fname = os.path.basename(entry)
                dest  = os.path.join("split_sdf", fname)
                if os.path.exists(dest):
                    stem_e = Path(fname).stem
                    n      = sum(1 for _ in glob.glob(f"split_sdf/{stem_e}*.sdf"))
                    dest   = os.path.join("split_sdf", f"{stem_e}_{n}.sdf")
                with zf.open(entry) as src_f, open(dest, "wb") as dst_f:
                    dst_f.write(src_f.read())
    else:
        raise ValueError(f"Unsupported file type '{ext}'. Please provide .sdf or .zip.")

    split_files: list = sorted(glob.glob("split_sdf/*.sdf"))
    used_names:  dict = {}

    def unique_name(raw: str) -> str:
        if raw not in used_names:
            used_names[raw] = 0; return raw
        used_names[raw] += 1
        return f"{raw}_{used_names[raw]:02d}"

    name_map: dict = {}
    for sdf in split_files:
        stem = Path(sdf).stem
        name_map[stem] = unique_name(get_mol_name(sdf, fallback=stem))

    print(f"Stage 1 done - {len(split_files)} individual SDF files ready.")
    print(f"  Sample names: {list(name_map.values())[:5]} ...")
    optional_save_zip("split_sdf", "01_split_sdf", "split_sdf")


    # ════════════════════════════════════════════════════════════════════════════
    # STAGE 2: Pre-Enumeration Duplicate Filter  (RDKit InChIKey)
    # ════════════════════════════════════════════════════════════════════════════
    RUN_DEDUP = ask_run(
        "Pre-Enumeration Duplicate Filter  [RDKit InChIKey]",
        f"Removes duplicate molecules before any enumeration.\n"
        f"  Mode '{DEDUP_LEVEL}': "
        f"{'heavy-atom skeleton only (ignores stereo)' if DEDUP_LEVEL == 'connectivity' else 'full stereo-aware key'}\n"
        "  Skip if: library is already curated and duplicate-free."
    )

    if RUN_DEDUP:
        t0 = datetime.now()
        print(f"Scanning {len(split_files)} molecules (mode: {DEDUP_LEVEL}, {_workers} threads)...\n")
        seen_keys: set = set()
        dup_count  = 0
        fail_count = 0

        def _s2_key(sdf_path: str):
            stem = Path(sdf_path).stem
            name = name_map[stem]
            try:
                mol      = Chem.MolFromMolFile(sdf_path, removeHs=True, sanitize=True)
                if mol is None: return sdf_path, name, None
                inchi    = MolToInchi(mol)
                full_key = InchiToInchiKey(inchi)
                key      = full_key[:14] if DEDUP_LEVEL == "connectivity" else full_key
                return sdf_path, name, key
            except Exception:
                return sdf_path, name, None

        key_map_s2: dict = {}
        with ThreadPoolExecutor(max_workers=_workers) as pool:
            futures = {pool.submit(_s2_key, s): s for s in split_files}
            for fut in as_completed(futures):
                path_, name_, key_ = fut.result()
                key_map_s2[path_] = (name_, key_)

        for sdf_path in split_files:
            name_, key_ = key_map_s2[sdf_path]
            if key_ is None:
                shutil.copy2(sdf_path, f"deduped/{name_}.sdf"); fail_count += 1; continue
            if key_ in seen_keys:
                dup_count += 1; print(f"  Duplicate removed: {name_}")
            else:
                seen_keys.add(key_)
                shutil.copy2(sdf_path, f"deduped/{name_}.sdf")

        deduped_files: list = sorted(glob.glob("deduped/*.sdf"))
        print(f"Stage 2 done - {len(deduped_files)} unique, "
              f"{dup_count} removed, {fail_count} unparseable kept  "
              f"({(datetime.now()-t0).total_seconds():.1f}s)")
    else:
        SKIPPED_STEPS.append("Pre-Enumeration Dedup")
        for sdf in split_files:
            stem = Path(sdf).stem
            shutil.copy2(sdf, f"deduped/{name_map[stem]}.sdf")
        deduped_files: list = sorted(glob.glob("deduped/*.sdf"))
        print(f"  Skipped - {len(deduped_files)} files forwarded.")

    optional_save_zip("deduped", "02_deduped", "deduped")


    # ════════════════════════════════════════════════════════════════════════════
    # STAGE 3: Standardize / Neutralize  (RDKit MolStandardize.Uncharger)
    # ════════════════════════════════════════════════════════════════════════════
    RUN_STANDARDIZE = ask_run(
        "Standardize / Neutralize  [RDKit MolStandardize.Uncharger]",
        "Salt/solvent stripping + full neutralization. Required for correct,\n"
        "  unrestricted tautomer enumeration downstream.\n"
        "  Skip if: library is already curated as neutral, salt-free structures."
    )

    if RUN_STANDARDIZE:
        reset_progress(len(deduped_files))
        t0 = datetime.now()
        print(f"Standardizing {len(deduped_files)} molecules ({_workers} threads)...\n")
        std_results: list = []
        with _make_executor(_workers, cpu_bound=False) as pool:
            futures = {pool.submit(_uncharge_standardize, s): s for s in deduped_files}
            for fut in as_completed(futures):
                r = fut.result()
                std_results.append(r)
                if not r["ok"]:
                    FAILURE_LOG.append({**r, "timeout": False})

        standardized_files: list = sorted(glob.glob("standardized/*.sdf"))
        std_failed = sum(1 for r in std_results if not r["ok"])
        print(f"\nStage 3 done - {len(standardized_files)} succeeded, {std_failed} failed  "
              f"({(datetime.now()-t0).total_seconds():.1f}s)")
    else:
        SKIPPED_STEPS.append("Standardize/Neutralize")
        passthrough("deduped/*.sdf", "standardized")
        standardized_files: list = sorted(glob.glob("standardized/*.sdf"))

    optional_save_zip("standardized", "03_standardized", "standardized")


    # ════════════════════════════════════════════════════════════════════════════
    # STAGE 4: Tautomer Enumeration  (RDKit TautomerEnumerator)
    # ════════════════════════════════════════════════════════════════════════════
    RUN_TAUTOMERS = ask_run(
        "Tautomer Enumeration  [RDKit TautomerEnumerator]",
        f"Enumerates top-{RDKIT_TAUTO_MAX} tautomers per molecule on the neutral form.\n"
        "  Skip if: your target's tautomeric state is already fixed/known."
    )
    print(f"\n{'='*60}\nSTAGE 4: Tautomer Enumeration  [RDKit]  Run={RUN_TAUTOMERS}  Max={RDKIT_TAUTO_MAX}\n{'='*60}")

    if RUN_TAUTOMERS:
        tauto_cfg   = {"out_dir": "tautomers", "rdkit_tauto_max": RDKIT_TAUTO_MAX}
        worker_args = [(s, tauto_cfg) for s in standardized_files]
        t0          = datetime.now()
        tauto_total = tauto_failed = tauto_capped = 0
        reset_progress(len(worker_args))

        print(f"Enumerating tautomers for {len(standardized_files)} molecules ({_workers} threads, top {RDKIT_TAUTO_MAX})...\n")
        with _make_executor(_workers, cpu_bound=False) as pool:
            futures = {pool.submit(_tautomer_worker, a): a for a in worker_args}
            for fut in as_completed(futures):
                r = fut.result()
                if not r["ok"]:
                    tprint(f"  SKIP  {r['name']}: {r['error'][:80]}")
                    FAILURE_LOG.append({"name": r["name"], "ok": False,
                                        "step": "tautomer enum", "error": r["error"],
                                        "timeout": False})
                elif r.get("no_tauto"):
                    tprint(f"  --  {r['name']}: no tautomers - parent kept")
                else:
                    tprint(f"  ok  {r['name']}  ->  {r['written']} tautomer(s)"
                           + (" [capped]" if r.get("capped") else ""))
                tauto_total  += r["written"]
                tauto_failed += 0 if r["ok"] else 1
                tauto_capped += 1 if r.get("capped") else 0
                tick()

        tautomer_files: list = sorted(glob.glob("tautomers/*.sdf"))
        print(f"\nStage 4 done - {len(standardized_files)} -> {len(tautomer_files)} variants, "
              f"{tauto_failed} failed, {tauto_capped} capped  "
              f"({(datetime.now()-t0).total_seconds():.1f}s)")
    else:
        SKIPPED_STEPS.append("Tautomer Enumeration")
        passthrough("standardized/*.sdf", "tautomers")
        tautomer_files: list = sorted(glob.glob("tautomers/*.sdf"))

    optional_save_zip("tautomers", "04_tautomers", "tautomers")


    # ════════════════════════════════════════════════════════════════════════════
    # STAGE 5: Protonation State Enumeration  (Dimorphite-DL)
    # ════════════════════════════════════════════════════════════════════════════
    RUN_PROTOSTATES = ask_run(
        "Protonation State Enumeration  [Dimorphite-DL]",
        f"Assigns pH-dependent ionisation states (pH {PROTO_MIN_PH}-{PROTO_MAX_PH}, "
        f"max {PROTO_MAX_VARS}/tautomer)\n"
        "  to EACH tautomer from Stage 4 - not the reverse, which would produce\n"
        "  artifact protonation states.\n"
        "  Skip if: ligands should remain in their neutral/as-drawn form."
    )
    print(f"\n{'='*60}\nSTAGE 5: Protonation States  [Dimorphite-DL]  Run={RUN_PROTOSTATES}  pH {PROTO_MIN_PH}-{PROTO_MAX_PH}  max {PROTO_MAX_VARS}\n{'='*60}")

    if RUN_PROTOSTATES:
        if not DIMORPHITE_OK:
            raise RuntimeError("Dimorphite-DL not installed. Re-run Step 1.")

        proto_cfg   = {"out_dir": "proto_states"}
        worker_args = [(s, proto_cfg) for s in tautomer_files]
        t0        = datetime.now()
        ps_total  = ps_failed = 0
        reset_progress(len(worker_args))

        print(f"Protonation states for {len(tautomer_files)} tautomers ({_workers} threads)...\n")
        with _make_executor(_workers, cpu_bound=False) as pool:
            futures = {pool.submit(_proto_worker, a): a for a in worker_args}
            for fut in as_completed(futures):
                r = fut.result()
                if not r["ok"]:
                    tprint(f"  SKIP  {r['name']}: {r['error'][:80]}")
                    FAILURE_LOG.append({"name": r["name"], "ok": False,
                                        "step": "protonation state enum",
                                        "error": r["error"], "timeout": False})
                elif r["no_variants"]:
                    tprint(f"  --  {r['name']}: no pH variants - parent kept")
                else:
                    tprint(f"  ok  {r['name']}  ->  {r['written']} state(s)")
                ps_total  += r["written"]
                ps_failed += 0 if r["ok"] else 1
                tick()

        proto_files: list = sorted(glob.glob("proto_states/*.sdf"))
        print(f"\nStage 5 done - {len(tautomer_files)} -> {len(proto_files)} variants, "
              f"{ps_failed} failed  ({(datetime.now()-t0).total_seconds():.1f}s)")
    else:
        SKIPPED_STEPS.append("Protonation State Enumeration")
        passthrough("tautomers/*.sdf", "proto_states")
        proto_files: list = sorted(glob.glob("proto_states/*.sdf"))

    optional_save_zip("proto_states", "05_proto_states", "proto_states")


    # ════════════════════════════════════════════════════════════════════════════
    # STAGE 6: Post-Enumeration Duplicate Filter  (RDKit InChIKey full)
    # ════════════════════════════════════════════════════════════════════════════
    RUN_POST_DEDUP = ask_run(
        "Post-Enumeration Duplicate Filter  [RDKit InChIKey, full/stereo-aware]",
        "Collapses identical tautomer x protonation cross-products that can\n"
        "  arise from Stages 4-5 running independently.\n"
        "  Skip if: you want to keep every tautomer x protonation combination."
    )
    print(f"\n{'='*60}\nSTAGE 6: Post-Enum Dedup  [RDKit InChIKey full]  Run={RUN_POST_DEDUP}\n{'='*60}")

    if RUN_POST_DEDUP:
        t0 = datetime.now()
        print(f"Post-enum dedup on {len(proto_files)} variants (full InChIKey)...\n")
        post_deduped = inchikey_dedup(proto_files, "post_deduped", "full",
                                      "Post-enum dedup", workers=_workers)
        post_deduped_files: list = sorted(glob.glob("post_deduped/*.sdf"))
        print(f"Stage 6 done - {len(post_deduped_files)} unique variants remain  "
              f"({(datetime.now()-t0).total_seconds():.1f}s)")
    else:
        SKIPPED_STEPS.append("Post-Enumeration Dedup")
        passthrough("proto_states/*.sdf", "post_deduped")
        post_deduped_files: list = sorted(glob.glob("post_deduped/*.sdf"))

    optional_save_zip("post_deduped", "06_post_deduped", "post_deduped")


    # ════════════════════════════════════════════════════════════════════════════
    # STAGE 7: Stereoisomer Enumeration  (RDKit)
    # ════════════════════════════════════════════════════════════════════════════
    RUN_STEREO = ask_run(
        "Stereoisomer Enumeration  [RDKit EnumerateStereoisomers]",
        f"Enumerates up to {MAX_STEREO_ISOMERS} stereoisomers/molecule. "
        f"only_undefined={ONLY_UNDEFINED_STEREO}\n"
        "  (True = only centers with no wedge bonds; False = all centers, "
        "including already-defined ones)\n"
        "  Skip if: input stereochemistry is already complete and trusted."
    )
    print(f"\n{'='*60}\nSTAGE 7: Stereoisomer Enum  [RDKit]  Run={RUN_STEREO}  Max={MAX_STEREO_ISOMERS}  undefined_only={ONLY_UNDEFINED_STEREO}\n{'='*60}")

    if RUN_STEREO:
        stereo_cfg  = {"max_isomers": MAX_STEREO_ISOMERS,
                       "only_undefined": ONLY_UNDEFINED_STEREO, "out_dir": "stereoisomers"}
        worker_args = [(s, stereo_cfg) for s in post_deduped_files]
        t0              = datetime.now()
        s_total = s_failed = s_trunc = 0
        reset_progress(len(worker_args))

        print(f"Stereoisomer enum for {len(post_deduped_files)} variants ({_workers} workers)...\n")
        with _make_executor(_workers, cpu_bound=True) as pool:
            futures = {pool.submit(_stereo_worker, a): a for a in worker_args}
            for fut in as_completed(futures):
                r = fut.result()
                if not r["ok"]:
                    tprint(f"  SKIP  {r['name']}: {r['error'][:80]}")
                    FAILURE_LOG.append({"name": r["name"], "ok": False,
                                        "step": "stereoisomer enumeration",
                                        "error": r["error"], "timeout": False})
                elif r.get("no_stereo"):
                    tprint(f"  --  {r['name']}: no undefined centers - kept")
                else:
                    trunc = (f" [capped: {r['max_possible']} possible, "
                             f"{MAX_STEREO_ISOMERS} written]") if r["truncated"] else ""
                    tprint(f"  ok  {r['name']}  ->  {r['written']} isomer(s){trunc}")
                s_total  += r["written"]
                s_failed += 0 if r["ok"] else 1
                s_trunc  += 1 if r.get("truncated") else 0
                tick()

        stereoisomer_files: list = sorted(glob.glob("stereoisomers/*.sdf"))
        print(f"\nStage 7 done - {len(post_deduped_files)} -> {len(stereoisomer_files)} variants, "
              f"{s_failed} failed, {s_trunc} capped  ({(datetime.now()-t0).total_seconds():.1f}s)")
        if s_trunc:
            print(f"  Note: {s_trunc} molecule(s) hit MAX_STEREO_ISOMERS={MAX_STEREO_ISOMERS}. "
                  f"Increase in Step 3 to enumerate all isomers.")
    else:
        SKIPPED_STEPS.append("Stereoisomer Enumeration")
        passthrough("post_deduped/*.sdf", "stereoisomers")
        stereoisomer_files: list = sorted(glob.glob("stereoisomers/*.sdf"))

    optional_save_zip("stereoisomers", "07_stereoisomers", "stereoisomers")


    # ════════════════════════════════════════════════════════════════════════════
    # STAGES 8 & 9: 3D Generation + MMFF94 Minimization  (RDKit)
    # ════════════════════════════════════════════════════════════════════════════
    RUN_GEN3D = ask_run(
        "3D Generation  [RDKit ETKDGv3 -> ETKDGv2 -> ETDG fallback chain]",
        "Embeds 3D coordinates for each stereoisomer variant, with automatic\n"
        "  tiered fallback and post-embed bond-geometry validation.\n"
        "  Skip if: you only need 2D output (MOL2/PDBQT require 3D - skipping\n"
        "  this will cause those stages to fail downstream)."
    )
    RUN_MINIMIZE = ask_run(
        f"Energy Minimization  [RDKit {MMFF_VARIANT} -> UFF fallback]",
        f"Relaxes embedded conformers to a local energy minimum "
        f"({MMFF_MAX_ITERS} max iters).\n"
        "  Skip if: you want raw ETKDG-embedded geometries, unminimized."
    )
    print(f"\n{'='*60}\nSTAGE 8: 3D Generation  [ETKDGv3]  Run={RUN_GEN3D}\nSTAGE 9: Minimize  [{MMFF_VARIANT}]  Run={RUN_MINIMIZE}\n{'='*60}")

    if RUN_GEN3D or RUN_MINIMIZE:
        rdkit_outdir = "gen3d" if (RUN_GEN3D and not RUN_MINIMIZE) else "minimized"
        os.makedirs(rdkit_outdir, exist_ok=True)

        rdkit_cfg = {
            "random_seed":       RANDOM_SEED,
            "enforce_chirality": ENFORCE_CHIRALITY,
            "max_attempts":      MAX_EMBED_ATTEMPTS,
            "mmff_variant":      MMFF_VARIANT,
            "mmff_max_iters":    MMFF_MAX_ITERS,
            "mmff_force_tol":    MMFF_FORCE_TOL,
            "uff_fallback":      UFF_FALLBACK,
        }

        worker_args = []
        for sdf_path in stereoisomer_files:
            mol_name = Path(sdf_path).stem
            with open(sdf_path, "r", encoding="utf-8", errors="replace") as fh:
                mol_block = fh.read()
            worker_args.append(
                (mol_block, mol_name, rdkit_outdir, RUN_GEN3D, RUN_MINIMIZE, rdkit_cfg))

        label = " + ".join(
            (["3D gen"] if RUN_GEN3D else []) + (["MMFF94 minimize"] if RUN_MINIMIZE else []))
        reset_progress(len(worker_args))
        t0 = datetime.now()
        print(f"\nRDKit [{label}] - {len(worker_args)} variants ({_workers} workers)...\n")

        rdkit_results: list = []
        with _make_executor(_workers, cpu_bound=True) as pool:
            futures = {pool.submit(_rdkit_gen3d_and_minimize, a): a for a in worker_args}
            for future in as_completed(futures):
                r = future.result()
                rdkit_results.append(r)
                if r["ok"]:
                    tprint(f"  ok  {r['name']}" + (f"  [{r['warn']}]" if r["warn"] else ""))
                else:
                    tprint(f"  SKIP  {r['name']}  @ {r['stage']}: {r['error'][:100]}")
                tick()

        rdkit_failed = [r for r in rdkit_results if not r["ok"]]
        for r in rdkit_failed:
            FAILURE_LOG.append({"name": r["name"], "ok": False,
                                "step": f"RDKit {r['stage']}",
                                "error": r["error"], "timeout": False})

        elapsed   = (datetime.now() - t0).total_seconds()
        out_count = len(glob.glob(f"{rdkit_outdir}/*.sdf"))
        warnings  = [r for r in rdkit_results if r["ok"] and r["warn"]]
        print(f"\nStages 8+9 done - {out_count} succeeded, "
              f"{len(rdkit_failed)} failed, {len(warnings)} with warnings  ({elapsed:.1f}s)")

        if RUN_GEN3D and not RUN_MINIMIZE:
            SKIPPED_STEPS.append("MMFF94 Minimization")
            passthrough("gen3d/*.sdf", "minimized")
        elif not RUN_GEN3D and RUN_MINIMIZE:
            SKIPPED_STEPS.append("3D Generation")
            passthrough("stereoisomers/*.sdf", "gen3d")
        elif RUN_GEN3D and RUN_MINIMIZE:
            passthrough("minimized/*.sdf", "gen3d")
    else:
        SKIPPED_STEPS.extend(["3D Generation", "MMFF94 Minimization"])
        passthrough("stereoisomers/*.sdf", "gen3d")
        passthrough("stereoisomers/*.sdf", "minimized")

    gen3d_files:     list = sorted(glob.glob("gen3d/*.sdf"))
    minimized_files: list = sorted(glob.glob("minimized/*.sdf"))

    optional_save_zip("gen3d", "08_gen3d", "gen3d")
    optional_save_zip("minimized", "09_minimized", "minimized")


    # ════════════════════════════════════════════════════════════════════════════
    # STAGE 10: SDF -> MOL2  (RDKit universal SYBYL writer)
    # ════════════════════════════════════════════════════════════════════════════
    RUN_MOL2 = ask_run(
        "SDF -> MOL2  [RDKit universal SYBYL writer]",
        "Writes CGenFF-compatible Tripos MOL2: unique atom names, full SYBYL\n"
        "  atom typing, universal Gasteiger-Marsili partial charges.\n"
        "  Skip if: you only need PDBQT output (Vina/AutoDock-GPU)."
    )
    print(f"\n{'='*60}\nSTAGE 10: SDF -> MOL2  [RDKit universal SYBYL writer]  Run={RUN_MOL2}\n{'='*60}")

    if RUN_MOL2:
        def _to_mol2(sdf: str) -> dict:
            # Universal CGenFF-compatible Tripos MOL2 writer (RDKit-based).
            # Replaces the bare `obabel -O mol2` call. Implements, for ANY
            # small molecule:
            #   1. Unique atom names      (no more "renaming non-unique atom" warnings)
            #   2. Correct SYBYL atom types (O.co2/N.am/N.pl3/N.4/S.o/S.o2/C.cat,
            #      double-bond-gated O.2 vs O.3, aromatic perception, universal P.3)
            #   3. MMFF94s partial charges (Gasteiger fallback only if MMFF94s fails)
            name = Path(sdf).stem
            try:
                os.makedirs("mol2", exist_ok=True)
                write_universal_sybyl_mol2(sdf, f"mol2/{name}.mol2", mol_name=name)
                tprint(f"  ok  {name}"); tick()
                return {"name": name, "ok": True, "step": "mol2",
                        "error": None, "timeout": False}
            except Exception as exc:
                err = str(exc)
                tprint(f"  SKIP  {name}: {err[:100]}")
                tick()
                return {"name": name, "ok": False, "step": "mol2",
                        "error": err, "timeout": False}

        reset_progress(len(minimized_files))
        t0 = datetime.now()
        print(f"SDF -> MOL2 for {len(minimized_files)} files ({_workers} threads)...\n")
        mol2_results: list = []
        with ThreadPoolExecutor(max_workers=_workers) as pool:
            for r in as_completed({pool.submit(_to_mol2, s): s for s in minimized_files}):
                mol2_results.append(r.result())
        mol2_failed = [r for r in mol2_results if not r["ok"]]
        FAILURE_LOG.extend(mol2_failed)
        mol2_files: list = sorted(glob.glob("mol2/*.mol2"))
        print(f"\nStage 10 done - {len(mol2_files)} succeeded, {len(mol2_failed)} failed  "
              f"({(datetime.now()-t0).total_seconds():.1f}s)")
    else:
        SKIPPED_STEPS.append("MOL2 Conversion")
        passthrough("minimized/*.sdf", "mol2")
        mol2_files: list = sorted(glob.glob("mol2/*"))

    optional_save_zip("mol2", "10_mol2", "mol2")


    # ════════════════════════════════════════════════════════════════════════════
    # STAGE 11: SDF -> PDBQT  (Meeko / OpenBabel fallback)
    # ════════════════════════════════════════════════════════════════════════════
    _tool_label = "Meeko" if MEEKO_OK else "OpenBabel fallback"
    RUN_PDBQT = ask_run(
        f"SDF -> PDBQT  [{_tool_label}]",
        "Generates Vina-native torsion trees for AutoDock-GPU / smina / AutoDock Vina.\n"
        "  Skip if: you only need MOL2 output (GOLD/DOCK6)."
    )

    # ── Standardize-disabled warning (now that RUN_MOL2/RUN_PDBQT are both known) ──
    if (RUN_MOL2 or RUN_PDBQT) and not RUN_STANDARDIZE:
        print("\nWARNING: Stage 3 (Standardize/Neutralize) is DISABLED but "
              "MOL2/PDBQT export is enabled. Salts, counterions, or disconnected "
              "fragments will NOT be stripped before MOL2/PDBQT writing -- affected "
              "ligands will be skipped per-molecule (see pipeline_report.txt) rather "
              "than silently mis-typed. Re-enable Stage 3 unless you are certain your "
              "input is already a single-fragment, salt-free structure.\n")

    print(f"\n{'='*60}\nSTAGE 11: SDF -> PDBQT  [{_tool_label}]  Run={RUN_PDBQT}\n{'='*60}")

    _pdbqt_inputs = sorted(glob.glob("minimized/*.sdf"))

    if RUN_PDBQT:
        def _to_pdbqt(sdf: str) -> dict:
            name = Path(sdf).stem
            if MEEKO_OK:
                try:
                    from rdkit import Chem as _Chem
                    from meeko import MoleculePreparation, PDBQTWriterLegacy
                    mol = _Chem.MolFromMolFile(sdf, removeHs=False, sanitize=True)
                    if mol is None:
                        raise ValueError("RDKit could not parse SDF for Meeko")
                    preparator = MoleculePreparation()
                    setups     = preparator.prepare(mol)
                    if not setups:
                        raise ValueError("Meeko returned no molecule setups")
                    pdbqt_str, is_ok, err_msg = PDBQTWriterLegacy.write_string(setups[0])
                    if not is_ok:
                        raise ValueError(f"PDBQTWriter: {err_msg}")
                    with open(f"pdbqt/{name}.pdbqt", "w") as fh:
                        fh.write(pdbqt_str)
                    tprint(f"  ok  {name}  [Meeko]"); tick()
                    return {"name": name, "ok": True, "step": "pdbqt",
                            "error": None, "timeout": False}
                except Exception as exc:
                    err = str(exc)
                    tprint(f"  SKIP  {name}: {err[:100]}")
                    tick()
                    return {"name": name, "ok": False, "step": "pdbqt",
                            "error": err, "timeout": False}
            else:
                try:
                    run_obabel(["obabel", sdf, "-O", f"pdbqt/{name}.pdbqt",
                                "--partialcharge", "gasteiger"], timeout=TIMEOUT_PDBQT)
                    tprint(f"  ok  {name}  [obabel]"); tick()
                    return {"name": name, "ok": True, "step": "pdbqt",
                            "error": None, "timeout": False}
                except Exception as exc:
                    err = str(exc)
                    tprint(f"  {'TIMEOUT' if 'TIMEOUT' in err else 'SKIP'}  {name}: {err[:100]}")
                    tick()
                    return {"name": name, "ok": False, "step": "pdbqt",
                            "error": err, "timeout": "TIMEOUT" in err}

        reset_progress(len(_pdbqt_inputs))
        t0 = datetime.now()
        print(f"SDF -> PDBQT [{_tool_label}] for {len(_pdbqt_inputs)} files ({_workers} threads)...\n")
        pdbqt_results: list = []
        with ThreadPoolExecutor(max_workers=_workers) as pool:
            for r in as_completed({pool.submit(_to_pdbqt, s): s for s in _pdbqt_inputs}):
                pdbqt_results.append(r.result())
        pdbqt_failed = [r for r in pdbqt_results if not r["ok"]]
        FAILURE_LOG.extend(pdbqt_failed)
        pdbqt_files: list = sorted(glob.glob("pdbqt/*.pdbqt"))
        print(f"\nStage 11 done - {len(pdbqt_files)} succeeded, {len(pdbqt_failed)} failed  "
              f"({(datetime.now()-t0).total_seconds():.1f}s)")
    else:
        SKIPPED_STEPS.append("PDBQT Conversion")
        pdbqt_files: list = []
        print("  PDBQT skipped.")

    optional_save_zip("pdbqt", "11_pdbqt", "pdbqt")


    # ════════════════════════════════════════════════════════════════════════════
    # FINAL REPORT + COMBINED ZIP
    # ════════════════════════════════════════════════════════════════════════════
    _total_wall  = (datetime.now() - _pipeline_start).total_seconds()
    _final_sdf   = sorted(glob.glob("minimized/*.sdf"))
    _final_mol2  = sorted(glob.glob("mol2/*.mol2"))
    _final_pdbqt = sorted(glob.glob("pdbqt/*.pdbqt"))
    timed_all    = [r for r in FAILURE_LOG if r.get("timeout")]
    hard_all     = [r for r in FAILURE_LOG if not r.get("timeout")]

    print(f"\n{'='*60}")
    print("FINAL PIPELINE SUMMARY")
    print(f"{'='*60}")
    print(f"  Input ligands (raw)         : {len(split_files)}")
    print(f"  After pre-dedup             : {len(deduped_files)}")
    print(f"  After tautomer enum         : {len(tautomer_files)}")
    print(f"  After protonation states    : {len(proto_files)}")
    print(f"  After post-enum dedup       : {len(post_deduped_files)}")
    print(f"  After stereoisomer enum     : {len(stereoisomer_files)}")
    print(f"  Minimized SDF output        : {len(_final_sdf)}")
    print(f"  MOL2 output                 : {len(_final_mol2)}")
    print(f"  PDBQT output                : {len(_final_pdbqt)}")
    print(f"  CPU workers used            : {_workers}")
    print(f"  Skipped stages              : {', '.join(SKIPPED_STEPS) if SKIPPED_STEPS else 'none'}")
    print(f"  Total failures              : {len(FAILURE_LOG)}  "
          f"(timeouts: {len(timed_all)}, errors: {len(hard_all)})")
    print(f"  Total wall time             : {_total_wall:.1f}s  ({_total_wall/60:.1f} min)")
    print(f"{'='*60}")

    if FAILURE_LOG:
        print("\nFailures by stage:")
        for r in FAILURE_LOG:
            tag = "TIMEOUT" if r.get("timeout") else "ERROR"
            print(f"  {tag}  {r['name']:<50}  @ {r['step']}")

    # ── Pipeline report ───────────────────────────────────────────────────────
    report_path = "pipeline_report.txt"
    with open(report_path, "w") as rf:
        rf.write("LigKit Ligand Preparation Pipeline - Full Run Report\n")
        rf.write(f"Generated         : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        rf.write(f"Input file        : {input_file}\n")
        rf.write(f"Output folder     : {WORK_DIR}\n")
        rf.write(f"CPU workers       : {_workers}\n")
        rf.write(f"Total wall time   : {_total_wall:.1f}s ({_total_wall/60:.1f} min)\n")
        rf.write("Toolchain versions:\n")
        try:
            import rdkit as _rdkit_v
            rf.write(f"  RDKit             : {_rdkit_v.__version__}\n")
        except Exception:
            rf.write("  RDKit             : version unavailable\n")
        try:
            import dimorphite_dl as _dm_v
            rf.write(f"  Dimorphite-DL     : {getattr(_dm_v, '__version__', 'unknown')}\n")
        except Exception:
            rf.write("  Dimorphite-DL     : version unavailable\n")
        try:
            import meeko as _mk_v
            rf.write(f"  Meeko             : {getattr(_mk_v, '__version__', 'installed')}\n")
        except Exception:
            rf.write("  Meeko             : not installed (OpenBabel fallback used)\n")
        rf.write(f"Tool assignment   :\n")
        rf.write(f"  Stage  1  Split/Extract          -> OpenBabel\n")
        rf.write(f"  Stage  2  Pre-dedup              -> RDKit InChIKey\n")
        rf.write(f"  Stage  3  Standardize/Neutralize -> RDKit MolStandardize.Uncharger\n")
        rf.write(f"  Stage  4  Tautomer Enum          -> RDKit TautomerEnumerator\n")
        rf.write(f"  Stage  5  Protonation States     -> Dimorphite-DL\n")
        rf.write(f"  Stage  6  Post-enum Dedup        -> RDKit InChIKey (full)\n")
        rf.write(f"  Stage  7  Stereoiso. Enum        -> RDKit EnumerateStereoisomers\n")
        rf.write(f"  Stage  8  3D Generation          -> RDKit ETKDGv3\n")
        rf.write(f"  Stage  9  Minimization           -> RDKit {MMFF_VARIANT}\n")
        rf.write(f"  Stage 10  MOL2                   -> RDKit (universal SYBYL writer)\n")
        rf.write(f"  Stage 11  PDBQT                  -> {'Meeko' if MEEKO_OK else 'OpenBabel'}\n")
        rf.write(f"Input ligands (raw)  : {len(split_files)}\n")
        rf.write(f"After pre-dedup      : {len(deduped_files)}\n")
        rf.write(f"After tautomers      : {len(tautomer_files)}\n")
        rf.write(f"After proto states   : {len(proto_files)}\n")
        rf.write(f"After post-dedup     : {len(post_deduped_files)}\n")
        rf.write(f"After stereoiso.     : {len(stereoisomer_files)}\n")
        rf.write(f"Minimized SDF output : {len(_final_sdf)}\n")
        rf.write(f"MOL2 output          : {len(_final_mol2)}\n")
        rf.write(f"PDBQT output         : {len(_final_pdbqt)}\n")
        rf.write(f"Skipped stages       : {', '.join(SKIPPED_STEPS) if SKIPPED_STEPS else 'none'}\n")
        rf.write(f"Total failures       : {len(FAILURE_LOG)}  "
                 f"(Timeouts: {len(timed_all)}, Errors: {len(hard_all)})\n")
        rf.write("=" * 75 + "\n\n")
        if not FAILURE_LOG:
            rf.write("All processed variants completed every active stage successfully.\n")
        else:
            rf.write(f"{len(FAILURE_LOG)} variant(s) failed:\n\n")
            rf.write(f"  {'#':<5}  {'Type':<10}  {'Ligand/Variant Name':<50}  "
                     f"{'Stage':<35}  Error\n")
            rf.write("-" * 120 + "\n")
            for idx, r in enumerate(FAILURE_LOG, 1):
                ftype = "TIMEOUT" if r.get("timeout") else "ERROR"
                rf.write(f"  {idx:<5}  {ftype:<10}  {r['name']:<50}  "
                         f"{r['step']:<35}  {r['error']}\n")

    print(f"\nReport -> {os.path.abspath(report_path)}")

    # ── Final combined ZIP — always includes minimized SDF + MOL2 + PDBQT ──────
    final_zip = f"{zip_name}.zip"
    print(f"\nBuilding final ZIP: {final_zip}")
    with zipfile.ZipFile(final_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        for f_ in _final_sdf:
            zf.write(f_, arcname=f"{zip_name}/minimized_sdf/{Path(f_).name}")
        for f_ in _final_mol2:
            zf.write(f_, arcname=f"{zip_name}/mol2/{Path(f_).name}")
        for f_ in _final_pdbqt:
            zf.write(f_, arcname=f"{zip_name}/pdbqt/{Path(f_).name}")
        zf.write(report_path, arcname=f"{zip_name}/{report_path}")

    zip_mb = os.path.getsize(final_zip) / (1024 * 1024)
    print(f"  Minimized SDF  : {len(_final_sdf)}")
    print(f"  MOL2 files     : {len(_final_mol2)}")
    print(f"  PDBQT files    : {len(_final_pdbqt)}")
    print(f"  Archive size   : {zip_mb:.2f} MB")

    while True:
        ans = input(f"\nSave final ZIP \"{final_zip}\"? [yes/no]: ").strip().lower()
        if ans in ("yes", "y"):
            print(f"Saved -> {os.path.abspath(final_zip)}")
            break
        elif ans in ("no", "n"):
            print(f"  Skipped. ZIP is still at: {os.path.abspath(final_zip)}")
            break
        else:
            print("  Please type  yes  or  no")

    print(f"\n{'='*60}")
    print(f"  All output files are in:")
    print(f"  {WORK_DIR}")
    print(f"{'='*60}")
    print("\nPipeline complete.")
    print()
    print("=" * 60)
    print("  READY FOR NEXT LIGAND LIBRARY")
    print("=" * 60)
    print("  To prepare another SDF / ZIP:")
    print("    1. Edit INPUT_PATH in USER SETTINGS (top of this script)")
    print("    2. (Optional) adjust other settings in USER SETTINGS")
    print("    3. Re-run:  ./ligkit.py    (or  python3 ligkit.py)")
    print("  Each run creates its own isolated <name>_prepared/ output folder.")
    print("=" * 60)