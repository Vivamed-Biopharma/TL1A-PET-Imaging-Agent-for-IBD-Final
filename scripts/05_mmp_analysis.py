#!/usr/bin/env python3
"""
Experiment 5: MMP Analysis (Matched Molecular Pair)

This script compares the NOTA chelator and linker-chelator as a matched pair
to understand the impact of linker addition on properties.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
import pandas as pd
import scripts.inputs as inputs
import logging
from pathlib import Path

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def calculate_mmp_properties(smiles_dict):
    """
    Calculate properties for MMP analysis.

    Args:
        smiles_dict (dict): Name: SMILES pairs

    Returns:
        pd.DataFrame: Properties for each molecule
    """
    results = []
    for name, smiles in smiles_dict.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES: {smiles}")
            continue

        props = {
            "Molecule": name,
            "SMILES": smiles,
            "MW": Descriptors.MolWt(mol),
            "LogP": Crippen.MolLogP(mol),
            "TPSA": Descriptors.TPSA(mol),
            "HBD": Descriptors.NumHDonors(mol),
            "HBA": Descriptors.NumHAcceptors(mol),
            "Rotatable_Bonds": Descriptors.NumRotatableBonds(mol),
            "Heavy_Atoms": Descriptors.HeavyAtomCount(mol),
            "Ring_Count": rdMolDescriptors.CalcNumRings(mol),
            "Aromatic_Rings": rdMolDescriptors.CalcNumAromaticRings(mol)
        }
        results.append(props)

    return pd.DataFrame(results)

def analyze_mmp_changes(df):
    """
    Analyze changes between chelator and linker-chelator.

    Args:
        df (pd.DataFrame): Properties DataFrame

    Returns:
        pd.DataFrame: Changes analysis
    """
    chelator = df[df['Molecule'] == 'NOTA_Chelator'].iloc[0]
    linker = df[df['Molecule'] == 'Linker_Chelator'].iloc[0]

    changes = {}
    for col in df.columns:
        if col not in ['Molecule', 'SMILES'] and pd.api.types.is_numeric_dtype(df[col]):
            changes[f"Delta_{col}"] = linker[col] - chelator[col]
            changes[f"Percent_Change_{col}"] = (changes[f"Delta_{col}"] / chelator[col] * 100) if chelator[col] != 0 else 0

    return pd.DataFrame([changes])

def validate_inputs():
    """Validate inputs."""
    if not hasattr(inputs, 'NOTA_CHELATOR_SMILES') or not hasattr(inputs, 'LINKER_CHELATOR_SMILES'):
        raise ValueError("Required SMILES not found in inputs.py")

def main():
    try:
        validate_inputs()

        molecules = {
            "NOTA_Chelator": inputs.NOTA_CHELATOR_SMILES,
            "Linker_Chelator": inputs.LINKER_CHELATOR_SMILES
        }

        # Calculate properties
        props_df = calculate_mmp_properties(molecules)

        # Analyze changes
        changes_df = analyze_mmp_changes(props_df)

        # Save results
        results_dir = Path("results/prodrug")
        results_dir.mkdir(parents=True, exist_ok=True)

        props_path = results_dir / "05_mmp_properties.csv"
        props_df.to_csv(props_path, index=False)

        changes_path = results_dir / "05_mmp_changes.csv"
        changes_df.to_csv(changes_path, index=False)

        logger.info(f"Properties saved to {props_path}")
        logger.info(f"Changes saved to {changes_path}")

        print("MMP Properties:")
        print(props_df.to_string(index=False))
        print("\nMMP Changes (Linker - Chelator):")
        print(changes_df.to_string(index=False))

        # Key insights
        mw_change = changes_df['Delta_MW'].values[0]
        logp_change = changes_df['Delta_LogP'].values[0]

        logger.info(f"MW increase: {mw_change:.1f} Da")
        logger.info(f"LogP change: {logp_change:.2f}")

        if mw_change > 300:
            logger.warning("Significant MW increase - may affect PK")
        if logp_change < -1:
            logger.info("Linker decreases lipophilicity - good for solubility")

    except Exception as e:
        logger.error(f"Error in Experiment 5: {str(e)}")
        raise

if __name__ == "__main__":
    main()