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
    # If chain C doesn't exist, analyze all chains
    antigen_part = complex_structure.select('chain C')

    if antigen_part is None:
        # No separate antigen chain - analyze intra-Fab interfaces (A-B)
        logger.info(f"No chain C found. Analyzing A-B interface instead.")
        fab_part_a = complex_structure.select('chain A')
        fab_part_b = complex_structure.select('chain B')
        if fab_part_a is None or fab_part_b is None:
            logger.error("Could not select chains A and B")
            return None
        interface_analysis_1 = fab_part_a
        interface_analysis_2 = fab_part_b
    else:
        # Standard Fab-antigen interface
        fab_part = complex_structure.select('chain A B')
        if fab_part is None:
            logger.error("Could not select Fab chains (A/B)")
            return None
        interface_analysis_1 = fab_part
        interface_analysis_2 = antigen_part

    # Find interface residues using distance criterion
    interface_residues = set()
    interface_distance = 5.0  # Angstroms

    # Get all atoms from both parts
    for atom1 in interface_analysis_1.iterAtoms():
        coords1 = atom1.getCoords()
        for atom2 in interface_analysis_2.iterAtoms():
            coords2 = atom2.getCoords()
            distance = ((coords1[0]-coords2[0])**2 +
                       (coords1[1]-coords2[1])**2 +
                       (coords1[2]-coords2[2])**2)**0.5
            if distance <= interface_distance:
                # Use ProDy's attribute accessors
                chid = atom1.getChid()
                resnum = atom1.getResnum()
                resname = atom1.getResname()
                interface_residues.add((chid, resnum, resname))
                break  # Found one close atom, residue is at interface

    hydrophobic = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'}
    total = len(interface_residues)
    num_hydrophobic = sum(1 for r in interface_residues if r[2] in hydrophobic)  # r[2] is resname

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
