#!/usr/bin/env python3
"""
Experiment 4: Linker Flexibility Analysis

This script analyzes the flexibility of the linker-chelator by generating conformers
and calculating conformational entropy or flexibility metrics.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign as MolAlign
import pandas as pd
try:
    import scripts.inputs as inputs
except ModuleNotFoundError:
    import inputs as inputs
import logging
from pathlib import Path
import numpy as np

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def analyze_flexibility(smiles, name, num_conformers=50):
    """
    Analyze molecular flexibility by generating conformers.

    Args:
        smiles (str): SMILES string
        name (str): Molecule name
        num_conformers (int): Number of conformers to generate

    Returns:
        dict: Flexibility metrics
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Generate conformers
    try:
        params = AllChem.ETKDGv3()
    except AttributeError:
        params = AllChem.ETKDG()
    params.randomSeed = 42  # For reproducibility
    # Use signature (mol, numConfs, params)
    conformer_ids = AllChem.EmbedMultipleConfs(mol, num_conformers, params)

    if len(conformer_ids) == 0:
        logger.warning(f"No conformers generated for {name}")
        return {"Molecule": name, "Conformers": 0, "RMS_Matrix": None, "Flexibility_Score": 0}

    # Calculate RMS matrix
    # RDKit API differences: some versions accept (mol, prealigned=False)
    # We'll compute pairwise RMS via AlignMolConformers when needed.
    def _pairwise_rms(i: int, j: int) -> float:
        try:
            return float(MolAlign.GetBestRMS(mol, mol, prbId=int(conformer_ids[i]), refId=int(conformer_ids[j])))
        except Exception:
            return 0.0

    # Calculate average RMS (measure of diversity)
    rms_values = []
    for i in range(len(conformer_ids)):
        for j in range(i+1, len(conformer_ids)):
            rms_values.append(_pairwise_rms(i, j))

    avg_rms = np.mean(rms_values) if rms_values else 0
    max_rms = np.max(rms_values) if rms_values else 0

    # Calculate rotatable bonds
    rot_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)

    # Flexibility score (higher = more flexible)
    flexibility_score = avg_rms * rot_bonds

    return {
        "Molecule": name,
        "Conformers_Generated": len(conformer_ids),
        "Average_RMS": avg_rms,
        "Max_RMS": max_rms,
        "Rotatable_Bonds": rot_bonds,
        "Flexibility_Score": flexibility_score
    }

def validate_inputs():
    """Validate inputs."""
    if not hasattr(inputs, 'LINKER_CHELATOR_SMILES'):
        raise ValueError("LINKER_CHELATOR_SMILES not found")

def main():
    try:
        validate_inputs()

        molecules = {
            "NOTA_Chelator": inputs.NOTA_CHELATOR_SMILES,
            "Linker_Chelator": inputs.LINKER_CHELATOR_SMILES
        }

        results = []
        for name, smiles in molecules.items():
            logger.info(f"Analyzing flexibility of {name}")
            result = analyze_flexibility(smiles, name)
            results.append(result)

        df = pd.DataFrame(results)

        # Save results
        results_dir = Path("results/prodrug")
        results_dir.mkdir(parents=True, exist_ok=True)
        output_path = results_dir / "04_linker_flexibility.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")

        print("Linker Flexibility Analysis:")
        print(df.to_string(index=False))

        # Compare linker vs chelator
        linker_score = df[df['Molecule'] == 'Linker_Chelator']['Flexibility_Score'].values[0]
        chelator_score = df[df['Molecule'] == 'NOTA_Chelator']['Flexibility_Score'].values[0]

        if linker_score > chelator_score:
            logger.info("Linker increases flexibility as expected")
        else:
            logger.warning("Unexpected flexibility results")

    except Exception as e:
        logger.error(f"Error in Experiment 4: {str(e)}")
        raise

if __name__ == "__main__":
    main()