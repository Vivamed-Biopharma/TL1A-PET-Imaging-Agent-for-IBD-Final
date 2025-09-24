#!/usr/bin/env python3
"""
Experiment 14: Immunogenicity Prediction

This script predicts immunogenicity of Fab sequences.
In a real implementation, this would use DeepImmuno or similar tools.
"""

import pandas as pd
import scripts.inputs as inputs
import logging
from pathlib import Path
import random
from scripts.neurosnap_wrappers import predict_immunogenicity

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def predict_immunogenicity(fab_name, sequence):
    """
    Predict immunogenicity using both sequence analysis and AI prediction.

    Args:
        fab_name (str): Fab name
        sequence (str): Sequence

    Returns:
        dict: Immunogenicity analysis
    """
    # Sequence-based analysis (fast calculation)
    # Factors that increase immunogenicity:
    # - High hydrophobicity
    # - Certain sequence motifs
    # - Non-human sequences

    # Since these are human framework Fabs, immunogenicity should be low
    base_score = 0.15  # Low immunogenicity for human frameworks

    # Check for potential T-cell epitopes (simplified)
    t_cell_motifs = ['VW', 'FW', 'YW', 'FL', 'WL', 'IL']
    t_cell_count = sum(sequence.count(motif) for motif in t_cell_motifs)

    # B-cell epitopes (hydrophobic patches)
    import re
    hydrophobic_stretches = len(re.findall(r'[FILVWYM]{3,}', sequence))

    # Adjust score based on sequence features
    seq_immunogenicity_score = base_score + (t_cell_count * 0.05) + (hydrophobic_stretches * 0.02)

    # AI-based prediction using DeepImmuno
    try:
        ai_results = predict_immunogenicity(sequence)
        ai_score = ai_results.get("immunogenicity_score", 0)
        ai_confidence = ai_results.get("confidence", 0.0)
        ai_epitopes = ai_results.get("predicted_epitopes", [])
    except Exception as e:
        logger.warning(f"AI immunogenicity prediction failed for {fab_name}: {str(e)}")
        ai_score = 0
        ai_confidence = 0.0
        ai_epitopes = []

    # Combine results (weighted average)
    if ai_confidence > 0.5:
        combined_score = (seq_immunogenicity_score * 0.4) + (ai_score * 0.6)
    else:
        combined_score = seq_immunogenicity_score

    # Classify risk level
    if combined_score < 0.2:
        risk_level = "Low"
    elif combined_score < 0.4:
        risk_level = "Medium"
    else:
        risk_level = "High"

    result = {
        "Fab_Name": fab_name,
        "Seq_Immunogenicity_Score": seq_immunogenicity_score,
        "AI_Immunogenicity_Score": ai_score,
        "Combined_Immunogenicity_Score": combined_score,
        "Risk_Level": risk_level,
        "AI_Confidence": ai_confidence,
        "T_Cell_Motifs": t_cell_count,
        "Hydrophobic_Stretches": hydrophobic_stretches,
        "Predicted_Epitopes": len(ai_epitopes),
        "Sequence_Length": len(sequence)
    }

    return result

def validate_inputs():
    """Validate sequences."""
    if not hasattr(inputs, 'fab_sequences'):
        raise ValueError("fab_sequences not found")

def main():
    try:
        validate_inputs()

        results = []
        for name, seq in inputs.fab_sequences.items():
            logger.info(f"Predicting immunogenicity for {name}")
            result = predict_immunogenicity(name, seq)
            results.append(result)

        df = pd.DataFrame(results)

        # Save results
        results_dir = Path("results/biobetter")
        results_dir.mkdir(parents=True, exist_ok=True)
        output_path = results_dir / "14_immunogenicity.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")

        print("Immunogenicity Prediction:")
        print(df[['Fab_Name', 'Immunogenicity_Score', 'Risk_Level']].to_string(index=False))

        # Validation
        high_risk = df[df['Risk_Level'] == 'High']
        if len(high_risk) > 0:
            logger.warning(f"High immunogenicity risk detected for: {', '.join(high_risk['Fab_Name'].tolist())}")

    except Exception as e:
        logger.error(f"Error in Experiment 14: {str(e)}")
        raise

if __name__ == "__main__":
    main()