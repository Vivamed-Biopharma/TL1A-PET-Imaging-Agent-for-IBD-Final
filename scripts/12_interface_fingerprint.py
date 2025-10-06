#!/usr/bin/env python3
"""
Experiment 12: Interface Fingerprinting (3D)

Performs a simple 3D interface analysis on Fab–TL1A complexes using ProDy.
It expects a complex PDB with chains A/B (Fab) and C (antigen), and quantifies
interface residues within 5 Å along with a hydrophobic ratio.
"""

import logging
import os
from pathlib import Path
from typing import Optional, Dict, Any

import pandas as pd
import prody as pr

# Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def analyze_interface_fingerprint(fab_name: str, complex_pdb_path: str) -> Optional[Dict[str, Any]]:
    if not os.path.exists(complex_pdb_path):
        logger.warning(f"Complex PDB not found for {fab_name} at {complex_pdb_path}. Skipping.")
        return None

    complex_structure = pr.parsePDB(complex_pdb_path)

    # Assumptions: chain A/B are Fab, chain C is antigen
    fab_part = complex_structure.select('chain A B')
    antigen_part = complex_structure.select('chain C')

    if fab_part is None or antigen_part is None:
        logger.error("Could not select Fab and Antigen chains (A/B vs C).")
        return None

    # Neighbor-based interface within 5 Å
    # ProDy's findNeighbors signature is findNeighbors(atoms, radius, sel2=None)
    neighbors = pr.findNeighbors(fab_part, 5.0, antigen_part)
    if not neighbors:
        logger.warning(f"No interface detected for {fab_name}")
        return {"Fab_Name": fab_name, "Total_Interface_Residues": 0, "Hydrophobic_Ratio": 0.0}

    interface_residues = set()
    # neighbors is a ProDy Neighbors instance; iterate over pairs via getPairs()
    try:
        pairs = neighbors.getPairs()
    except Exception:
        pairs = []
    for atom_a, atom_b in pairs:
        res = atom_a.getResidue()
        if res is not None:
            interface_residues.add(res)

    hydrophobic = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'}
    total = len(interface_residues)
    num_hydrophobic = sum(1 for r in interface_residues if r.getResname() in hydrophobic)

    return {
        "Fab_Name": fab_name,
        "Total_Interface_Residues": total,
        "Hydrophobic_Ratio": (num_hydrophobic / total) if total else 0.0,
    }


def main():
    # Discover complexes from experiment 11 output if available; otherwise look for a default path
    results_dir = Path("results/biobetter")
    complexes_dir = results_dir / "11_complex_models"

    analyses = []

    # If we have models from experiment 11, iterate; otherwise, attempt a fallback example
    candidate_files = []
    if complexes_dir.exists():
        candidate_files.extend(sorted(str(p) for p in complexes_dir.glob("*.pdb")))

    # Fallback example file (user-provided or generated)
    fallback = os.environ.get("TL1A_COMPLEX_PDB", "data/fab_model.pdb")
    if os.path.exists(fallback):
        candidate_files.append(fallback)

    # Deduplicate while preserving order
    seen = set()
    unique_candidates = []
    for p in candidate_files:
        if p not in seen:
            unique_candidates.append(p)
            seen.add(p)

    for pdb_path in unique_candidates:
        fab_name = Path(pdb_path).stem.capitalize()
        logger.info(f"Analyzing interface fingerprint for {fab_name}")
        result = analyze_interface_fingerprint(fab_name, pdb_path)
        if result:
            analyses.append(result)

    if not analyses:
        logger.warning("No complexes analyzed. Ensure PDBs are present from Experiment 11.")
        return

    df = pd.DataFrame(analyses)
    out_dir = results_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "12_interface_fingerprint.csv"
    df.to_csv(out_path, index=False)
    logger.info(f"Results saved to {out_path}")

    print("Interface Fingerprint Analysis:")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
