#!/usr/bin/env python3
"""
Experiment 3: In-Silico Toxicity & Liability Flagging

This script scans molecules for known toxic or reactive chemical groups (PAINS, reactive groups).
"""

from rdkit import Chem
import pandas as pd
try:
    import scripts.inputs as inputs
except ModuleNotFoundError:
    import inputs as inputs
import logging
from pathlib import Path
try:
    from scripts.neurosnap_wrappers import predict_toxicity as ns_predict_toxicity
except ModuleNotFoundError:
    from neurosnap_wrappers import predict_toxicity as ns_predict_toxicity

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Define liability filters as SMARTS patterns
LIABILITY_FILTERS = {
    "Isothiocyanate": "[N;D1]=[C;D2]=[S;D1]",  # Reactive group in linker
    "Michael_Acceptor": "[C;D1]=[C;D1][C;D2]=[O;D1]",  # Potential for nucleophilic attack
    "Epoxide": "[O;D1]1[C;D2][C;D2]1",  # Reactive ring
    "Aromatic_Nitro": "[c;D2][N;D2](=[O;D1])=[O;D1]",  # Potential mutagen
    "Alkyl_Halide": "[C;D1][Cl,Br,I]",  # Reactive carbon-halogen
    "PAINS_A": "[a;D2][N;D2][C;D2](=[O;D1])[C;D2](=[O;D1])[N;D2][a;D2]",  # Example PAINS pattern
}

def scan_liabilities(smiles_dict):
    """
    Scan molecules for liabilities using both structural patterns and AI prediction.

    Args:
        smiles_dict (dict): Dictionary of name: SMILES

    Returns:
        pd.DataFrame: DataFrame with detected liabilities and predictions
    """
    results = []
    for name, smiles in smiles_dict.items():
        logger.info(f"Analyzing {name}")

        # Structural pattern matching (fast pre-screen)
        mol = Chem.MolFromSmiles(smiles)
        structural_hits = []
        if mol is not None:
            for liability_name, smarts in LIABILITY_FILTERS.items():
                patt = Chem.MolFromSmarts(smarts)
                if patt and mol.HasSubstructMatch(patt):
                    structural_hits.append(liability_name)
                    logger.warning(f"Structural liability '{liability_name}' found in {name}")

        # AI-based toxicity prediction
        try:
            ai_results = ns_predict_toxicity(smiles)
            # Normalize returned schema from NeuroSnap output
            prediction = ai_results.get("prediction")
            ai_toxic = ai_results.get("toxic")
            if ai_toxic is None and prediction is not None:
                ai_toxic = str(prediction).upper() in {"TOXIC", "YES", "POSITIVE"}
            ai_confidence = float(ai_results.get("probability", ai_results.get("confidence", 0.0)) or 0.0)
            ai_details = ai_results.get("details", prediction or "No details")
        except Exception as e:
            logger.warning(f"AI toxicity prediction failed for {name}: {str(e)}")
            ai_toxic = None
            ai_confidence = 0.0
            ai_details = f"Prediction failed: {str(e)}"

        # Combine results
        all_liabilities = structural_hits.copy()
        if ai_toxic:
            all_liabilities.append(f"AI_Predicted_Toxic (confidence: {ai_confidence:.2f})")

        results.append({
            "Molecule": name,
            "SMILES": smiles,
            "Structural_Liabilities": "; ".join(structural_hits) if structural_hits else "None detected",
            "AI_Toxic": ai_toxic,
            "AI_Confidence": ai_confidence,
            "AI_Details": ai_details,
            "All_Liabilities": "; ".join(all_liabilities) if all_liabilities else "None detected"
        })

    return pd.DataFrame(results)

def validate_inputs():
    """Validate required inputs."""
    if not hasattr(inputs, 'LINKER_CHELATOR_SMILES'):
        raise ValueError("LINKER_CHELATOR_SMILES not found in inputs.py")

def main():
    try:
        validate_inputs()

        smiles_dict = {
            "NOTA_Chelator": inputs.NOTA_CHELATOR_SMILES,
            "Linker_Chelator": inputs.LINKER_CHELATOR_SMILES
        }

        df = scan_liabilities(smiles_dict)

        # Save results
        results_dir = Path("results/prodrug")
        results_dir.mkdir(parents=True, exist_ok=True)
        output_path = results_dir / "03_liability_hits.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")

        print("Liability Scan Results:")
        print(df.to_string(index=False))

        # Check for expected isothiocyanate
        linker_row = df[df['Molecule'] == 'Linker_Chelator']
        if not linker_row.empty and 'Isothiocyanate' in (linker_row['Structural_Liabilities'].values[0] or ''):
            logger.info("Expected isothiocyanate group detected in linker")
        else:
            logger.warning("Isothiocyanate not detected - check SMILES")

    except Exception as e:
        logger.error(f"Error in Experiment 3: {str(e)}")
        raise

if __name__ == "__main__":
    main()