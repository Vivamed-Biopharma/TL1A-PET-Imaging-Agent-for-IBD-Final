#!/usr/bin/env python3
"""
Experiment 2: Activation & Metabolism Prediction

This script predicts metabolic transformations of the linker-chelator using
RDKit-based SMARTS pattern matching for common metabolic reactions.
"""

import pandas as pd
try:
    import scripts.inputs as inputs
except ModuleNotFoundError:
    import inputs as inputs
import logging
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Common metabolic transformations as SMARTS patterns
METABOLIC_REACTIONS = {
    "Isothiocyanate_Hydrolysis": {
        "description": "Hydrolysis of isothiocyanate to amine",
        "reactant": "[N;D1]=[C;D2]=[S;D1]",
        "product_transform": "Replace N=C=S with NH2",
        "enzymes": "Non-enzymatic (H2O)",
        "likelihood": "High",
        "site": "Isothiocyanate group"
    },
    "Carboxylic_Acid_Conjugation": {
        "description": "Glucuronidation of carboxylic acids",
        "reactant": "C(=O)[OH]",
        "product_transform": "Glucuronide conjugate formation",
        "enzymes": "UGT (UDP-glucuronosyltransferase)",
        "likelihood": "Medium",
        "site": "Carboxylic acid groups"
    },
    "Secondary_Amine_Oxidation": {
        "description": "N-oxidation of secondary amines",
        "reactant": "[NX3;H1,H2;!$(NC=O)]",
        "product_transform": "N-oxide formation",
        "enzymes": "FMO (Flavin-containing monooxygenase)",
        "likelihood": "Low",
        "site": "Secondary amines in NOTA ring"
    },
    "Benzyl_Hydroxylation": {
        "description": "Aromatic hydroxylation",
        "reactant": "c[CH2]",
        "product_transform": "Para-hydroxylation of benzyl group",
        "enzymes": "CYP (Cytochrome P450)",
        "likelihood": "Medium",
        "site": "Benzyl linker"
    }
}

def predict_metabolic_sites(smiles):
    """
    Identify potential metabolic sites using SMARTS patterns.

    Args:
        smiles (str): Input SMILES string

    Returns:
        list: List of predicted metabolic transformations
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.error(f"Invalid SMILES: {smiles}")
        return []

    predictions = []

    for reaction_name, reaction_info in METABOLIC_REACTIONS.items():
        pattern = Chem.MolFromSmarts(reaction_info["reactant"])
        if pattern and mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match_idx, match in enumerate(matches, 1):
                predictions.append({
                    "Reaction": reaction_name,
                    "Description": reaction_info["description"],
                    "Enzymes": reaction_info["enzymes"],
                    "Likelihood": reaction_info["likelihood"],
                    "Site": reaction_info["site"],
                    "Atom_Indices": str(match),
                    "Match_Count": len(matches)
                })

    return predictions

def assess_metabolic_stability(smiles, predictions):
    """
    Assess overall metabolic stability based on predictions.

    Args:
        smiles (str): Input SMILES
        predictions (list): List of predicted transformations

    Returns:
        dict: Stability assessment
    """
    mol = Chem.MolFromSmiles(smiles)

    # Count reactive sites by likelihood
    high_risk = sum(1 for p in predictions if p["Likelihood"] == "High")
    medium_risk = sum(1 for p in predictions if p["Likelihood"] == "Medium")
    low_risk = sum(1 for p in predictions if p["Likelihood"] == "Low")

    # Calculate stability score (0-100, higher = more stable)
    stability_score = 100 - (high_risk * 30 + medium_risk * 15 + low_risk * 5)
    stability_score = max(0, stability_score)

    # Assess reactivity of key functional groups
    has_isothiocyanate = mol.HasSubstructMatch(Chem.MolFromSmarts("[N;D1]=[C;D2]=[S;D1]"))

    assessment = {
        "Molecule": "Linker-Chelator (p-SCN-Bn-NOTA)",
        "Total_Metabolic_Sites": len(predictions),
        "High_Risk_Sites": high_risk,
        "Medium_Risk_Sites": medium_risk,
        "Low_Risk_Sites": low_risk,
        "Stability_Score": stability_score,
        "Has_Reactive_Isothiocyanate": has_isothiocyanate,
        "Assessment": "Designed for bioconjugation - isothiocyanate is intentionally reactive",
        "Recommendation": "Monitor hydrolysis kinetics in aqueous buffer"
    }

    return assessment

def main():
    try:
        logger.info("Running Experiment 2: Metabolism Prediction")

        # Get input molecule
        prodrug_smiles = inputs.LINKER_CHELATOR_SMILES
        logger.info(f"Analyzing: {prodrug_smiles}")

        # Create results directory
        results_dir = Path("results/prodrug")
        results_dir.mkdir(parents=True, exist_ok=True)

        # Predict metabolic sites
        predictions = predict_metabolic_sites(prodrug_smiles)
        logger.info(f"Identified {len(predictions)} potential metabolic transformations")

        # Save predictions
        if predictions:
            output_csv = results_dir / "02_metabolism_predictions.csv"
            df_predictions = pd.DataFrame(predictions)
            df_predictions.to_csv(output_csv, index=False)
            logger.info(f"Predictions saved to {output_csv}")

            print("\nPredicted Metabolic Transformations:")
            print(df_predictions[['Reaction', 'Description', 'Likelihood', 'Enzymes']].to_string(index=False))
        else:
            logger.info("No metabolic transformations predicted - molecule appears stable")

        # Assess stability
        assessment = assess_metabolic_stability(prodrug_smiles, predictions)
        assessment_csv = results_dir / "02_metabolic_stability_assessment.csv"
        df_assessment = pd.DataFrame([assessment])
        df_assessment.to_csv(assessment_csv, index=False)

        print("\nMetabolic Stability Assessment:")
        for key, value in assessment.items():
            print(f"  {key}: {value}")

        # Create detailed report for compatibility with downstream analyses
        detailed_output = results_dir / "02_biotransformer_metabolites.csv"
        if predictions:
            # Create BioTransformer-compatible format
            biotransformer_format = []
            for pred in predictions:
                biotransformer_format.append({
                    "SMILES": prodrug_smiles,  # Original molecule
                    "Reaction": pred["Reaction"],
                    "Enzymes": pred["Enzymes"],
                    "Description": pred["Description"],
                    "Likelihood": pred["Likelihood"],
                    "Site": pred["Site"]
                })
            df_bt = pd.DataFrame(biotransformer_format)
            df_bt.to_csv(detailed_output, index=False)
            logger.info(f"BioTransformer-compatible output saved to {detailed_output}")
        else:
            # Create empty file for compatibility
            df_empty = pd.DataFrame(columns=["SMILES", "Reaction", "Enzymes", "Description", "Likelihood", "Site"])
            df_empty.to_csv(detailed_output, index=False)

        logger.info("Experiment 2 complete")

    except Exception as e:
        logger.error(f"Error in Experiment 2: {str(e)}")
        raise

if __name__ == "__main__":
    main()