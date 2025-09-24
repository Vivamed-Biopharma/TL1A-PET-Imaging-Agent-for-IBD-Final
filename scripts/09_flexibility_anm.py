#!/usr/bin/env python3
"""
Experiment 9: Local Flexibility Analysis with ANM

This script uses Anisotropic Network Model (ANM) to predict local flexibility
in Fab structures. Requires a PDB file.
"""

import prody as pr
import pandas as pd
import matplotlib.pyplot as plt
import scripts.inputs as inputs
import logging
import os
from pathlib import Path

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def analyze_flexibility_anm(pdb_path, name):
    """
    Analyze flexibility using ANM.

    Args:
        pdb_path (str): Path to PDB file
        name (str): Structure name

    Returns:
        dict: Flexibility analysis results
    """
    try:
        # Load structure
        protein = pr.parsePDB(pdb_path)
        calphas = protein.select('calpha')

        if calphas is None or len(calphas) == 0:
            raise ValueError("No CA atoms found in PDB")

        logger.info(f"Loaded {len(calphas)} CA atoms from {pdb_path}")

        # Build ANM
        anm = pr.ANM(name)
        anm.buildHessian(calphas)
        anm.calcModes(n_modes=10)

        # Calculate mean square fluctuations
        sqflucts = pr.calcSqFlucts(anm)

        # Calculate per-residue flexibility metrics
        residues = calphas.getResnums()
        flexibility_data = []

        for i, (resnum, fluct) in enumerate(zip(residues, sqflucts)):
            flexibility_data.append({
                "Residue_Number": resnum,
                "MSF": fluct,
                "Normalized_Flexibility": fluct / max(sqflucts) if max(sqflucts) > 0 else 0
            })

        # Summary statistics
        avg_flex = sum(sqflucts) / len(sqflucts)
        max_flex = max(sqflucts)
        flex_range = max_flex - min(sqflucts)

        result = {
            "Name": name,
            "Total_Residues": len(calphas),
            "Average_MSF": avg_flex,
            "Max_MSF": max_flex,
            "MSF_Range": flex_range,
            "Flexibility_Data": flexibility_data
        }

        return result

    except Exception as e:
        logger.error(f"Error in ANM analysis for {name}: {str(e)}")
        return {"Name": name, "Error": str(e)}

def plot_flexibility(flexibility_data, name, output_path):
    """
    Plot flexibility profile.

    Args:
        flexibility_data (list): Flexibility data
        name (str): Name for plot
        output_path (str): Output plot path
    """
    if not flexibility_data:
        return

    df = pd.DataFrame(flexibility_data)
    plt.figure(figsize=(12, 6))
    plt.plot(df['Residue_Number'], df['MSF'], 'b-', alpha=0.7)
    plt.xlabel('Residue Number')
    plt.ylabel('Mean Square Fluctuation')
    plt.title(f'Local Flexibility Profile: {name}')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', format='png', transparent=False)
    plt.close()
    logger.info(f"Flexibility plot saved to {output_path}")

def validate_inputs():
    """Validate that PDB path is available."""
    # Assume we have a placeholder or download a Fab PDB
    pass

def main():
    try:
        pdb_path = "data/fab_model.pdb"

        if not os.path.exists(pdb_path):
            logger.error(f"PDB file not found: {pdb_path}")
            logger.error("Please download a Fab structure from PDB (e.g., 1IGT) and save as data/fab_model.pdb")
            raise FileNotFoundError(f"Required PDB file missing: {pdb_path}")

        # Real ANM analysis using ProDy
        result = analyze_flexibility_anm(pdb_path, "Fab_Model")
        results = [result]

        # Process results
        for result in results:
            if "Error" in result:
                continue

            # Save detailed data
            results_dir = Path("results/formulation")
            results_dir.mkdir(parents=True, exist_ok=True)

            df_detailed = pd.DataFrame(result["Flexibility_Data"])
            detailed_path = results_dir / f"09_flexibility_{result['Name']}.csv"
            df_detailed.to_csv(detailed_path, index=False)

            # Plot
            plot_path = results_dir / f"09_flexibility_{result['Name']}.png"
            plot_flexibility(result["Flexibility_Data"], result["Name"], plot_path)

        # Summary
        summary_df = pd.DataFrame([{k: v for k, v in r.items() if k != "Flexibility_Data"} for r in results])
        summary_path = results_dir / "09_flexibility_summary.csv"
        summary_df.to_csv(summary_path, index=False)

        logger.info(f"Summary saved to {summary_path}")

        print("ANM Flexibility Analysis Summary:")
        print(summary_df.to_string(index=False))

    except Exception as e:
        logger.error(f"Error in Experiment 9: {str(e)}")
        raise

if __name__ == "__main__":
    main()