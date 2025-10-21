#!/usr/bin/env python3
"""
Experiment 1: Foundational Physicochemical Profiling

This script calculates basic drug-like properties for small molecules in the TL1A program.
Properties include molecular weight, logP, hydrogen bond donors/acceptors, rotatable bonds, and TPSA.
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
try:
    import scripts.inputs as inputs
except ModuleNotFoundError:
    import inputs as inputs
import logging
import os
from pathlib import Path
try:
    from scripts.error_handling import handle_errors
except ModuleNotFoundError:
    from error_handling import handle_errors
from typing import Dict

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@handle_errors
def calculate_properties(smiles_dict: Dict[str, str]) -> pd.DataFrame:
    """
    Calculate physicochemical properties for a dictionary of SMILES strings.

    Args:
        smiles_dict (dict): Dictionary with names as keys and SMILES as values.

    Returns:
        pd.DataFrame: DataFrame with calculated properties.
    """
    results = []
    for name, smiles in smiles_dict.items():
        logger.info(f"Processing molecule: {name}")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES for {name}: {smiles}")
            continue

        # Calculate properties
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)

        results.append({
            "Name": name,
            "SMILES": smiles,
            "MW": mw,
            "LogP": logp,
            "HBD": hbd,
            "HBA": hba,
            "Rotatable_Bonds": rot_bonds,
            "TPSA": tpsa
        })

    return pd.DataFrame(results)

def validate_inputs():
    """Validate that required inputs are present."""
    if not hasattr(inputs, 'NOTA_CHELATOR_SMILES'):
        raise ValueError("NOTA_CHELATOR_SMILES not found in inputs.py")
    if not hasattr(inputs, 'LINKER_CHELATOR_SMILES'):
        raise ValueError("LINKER_CHELATOR_SMILES not found in inputs.py")

@handle_errors
def main():
    try:
        validate_inputs()

        # Define molecules to analyze
        smiles_dict = {
            "NOTA_Chelator": inputs.NOTA_CHELATOR_SMILES,
            "Linker_Chelator": inputs.LINKER_CHELATOR_SMILES
        }

        # Calculate properties
        df = calculate_properties(smiles_dict)

        # Ensure results directory exists
        results_dir = Path("results/prodrug")
        results_dir.mkdir(parents=True, exist_ok=True)

        # Save results
        output_path = results_dir / "01_physchem_properties.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")

        # Print summary
        print("Physicochemical Properties Summary:")
        print(df.to_string(index=False))

        # Basic validation
        for _, row in df.iterrows():
            if row['MW'] > 1000:
                logger.warning(f"High MW detected for {row['Name']}: {row['MW']}")
            if row['LogP'] > 5 or row['LogP'] < -2:
                logger.warning(f"Extreme LogP for {row['Name']}: {row['LogP']}")

    except Exception as e:
        logger.error(f"Error in Experiment 1: {str(e)}")
        raise

if __name__ == "__main__":
    main()