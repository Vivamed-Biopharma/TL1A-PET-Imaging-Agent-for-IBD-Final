#!/usr/bin/env python3
"""
Experiment 2: Activation & Metabolism Prediction

This script uses BioTransformer to predict metabolic transformations of the linker-chelator.
BioTransformer simulates human metabolism pathways.
"""

import subprocess
import pandas as pd
import scripts.inputs as inputs
import logging
import os
from pathlib import Path

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_biotransformer(smiles, output_csv, jar_path):
    """
    Run BioTransformer on a SMILES string.

    Args:
        smiles (str): Input SMILES
        output_csv (str): Output CSV path
        jar_path (str): Path to BioTransformer JAR

    Returns:
        bool: True if successful
    """
    cmd = [
        "java", "-jar", jar_path,
        "-k", "h_all",  # Human all-in-one metabolism
        "-ismi", smiles,
        "-ocsv", output_csv
    ]

    logger.info(f"Running BioTransformer: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            logger.error(f"BioTransformer failed: {result.stderr}")
            return False
        logger.info("BioTransformer completed successfully")
        return True
    except subprocess.TimeoutExpired:
        logger.error("BioTransformer timed out")
        return False
    except FileNotFoundError:
        logger.error("Java not found. Please install Java.")
        return False

def validate_inputs():
    """Validate inputs and dependencies."""
    if not hasattr(inputs, 'LINKER_CHELATOR_SMILES'):
        raise ValueError("LINKER_CHELATOR_SMILES not found in inputs.py")

def main():
    try:
        validate_inputs()

        prodrug_smiles = inputs.LINKER_CHELATOR_SMILES
        jar_path = Path("scripts/BioTransformer3.0.jar")
        results_dir = Path("results/prodrug")
        results_dir.mkdir(parents=True, exist_ok=True)
        output_csv = results_dir / "02_biotransformer_metabolites.csv"

        if not jar_path.exists():
            logger.error(f"BioTransformer JAR not found at {jar_path}")
            logger.error("Please download from https://github.com/BioTransformer/BioTransformer3.0-cli/releases")
            raise FileNotFoundError(f"BioTransformer JAR missing: {jar_path}")

        # Run BioTransformer
        success = run_biotransformer(prodrug_smiles, str(output_csv), str(jar_path))

        if success and output_csv.exists():
            # Load and display results
            metabolites_df = pd.read_csv(output_csv)
            logger.info(f"Metabolites predicted: {len(metabolites_df)}")

            print("Predicted Metabolites:")
            if len(metabolites_df) > 0:
                print(metabolites_df[['SMILES', 'Reaction', 'Enzymes']].head(10).to_string(index=False))
            else:
                print("No metabolites predicted - molecule may be stable")

            # Analyze results
            if 'Reaction' in metabolites_df.columns:
                hydrolysis_count = metabolites_df['Reaction'].str.contains('hydrolysis', case=False).sum()
                logger.info(f"Hydrolysis reactions predicted: {hydrolysis_count}")

        else:
            logger.error("BioTransformer did not produce output")

    except Exception as e:
        logger.error(f"Error in Experiment 2: {str(e)}")
        raise

if __name__ == "__main__":
    main()