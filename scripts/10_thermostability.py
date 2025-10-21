#!/usr/bin/env python3
"""
Experiment 10: Thermostability Prediction

This script predicts thermostability of Fab sequences using sequence-based metrics.
In a real implementation, this would use ThermoMPNN or similar AI tools.
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
try:
    import scripts.inputs as inputs
except ModuleNotFoundError:
    import inputs as inputs
import logging
from pathlib import Path
try:
    from scripts.neurosnap_wrappers import predict_thermostability as ns_predict_thermostability
except ModuleNotFoundError:
    from neurosnap_wrappers import predict_thermostability as ns_predict_thermostability

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Thermostability factors (simplified)
THERMO_FACTORS = {
    'stabilizing': ['A', 'L', 'I', 'V', 'W', 'F', 'P'],
    'destabilizing': ['G', 'S', 'T', 'C', 'M', 'Q', 'N'],
    'neutral': ['K', 'R', 'H', 'D', 'E', 'Y']
}

def predict_thermostability(sequence, name):
    """
    Predict thermostability using both sequence-based metrics and AI prediction.

    Args:
        sequence (str): Amino acid sequence
        name (str): Sequence name

    Returns:
        dict: Thermostability analysis
    """
    try:
        analyzed_seq = ProteinAnalysis(sequence)

        # Sequence-based analysis (fast calculation)
        stabilizing_count = sum(sequence.count(aa) for aa in THERMO_FACTORS['stabilizing'])
        destabilizing_count = sum(sequence.count(aa) for aa in THERMO_FACTORS['destabilizing'])

        seq_thermo_score = stabilizing_count - destabilizing_count
        seq_thermo_score_normalized = seq_thermo_score / len(sequence)

        seq_tm_estimate = 50 + (seq_thermo_score_normalized * 20) + (analyzed_seq.instability_index() * -0.5)

        # AI-based prediction using TemStaPro via NeuroSnap
        try:
            ai_results = ns_predict_thermostability(sequence, temperature=25.0)
            ai_tm = ai_results.get("melting_temperature", 0)
            ai_confidence = ai_results.get("confidence", 0.0)
            ai_stability_score = ai_results.get("stability_score", 0)
        except Exception as e:
            logger.warning(f"AI thermostability prediction failed for {name}: {str(e)}")
            ai_tm = 0
            ai_confidence = 0.0
            ai_stability_score = 0

        # Combine results (weighted average)
        if ai_confidence > 0.5:
            combined_tm = (seq_tm_estimate * 0.3) + (ai_tm * 0.7)
            combined_score = (seq_thermo_score * 0.3) + (ai_stability_score * 0.7)
        else:
            combined_tm = seq_tm_estimate
            combined_score = seq_thermo_score

        result = {
            "Name": name,
            "Length": len(sequence),
            "Stabilizing_Residues": stabilizing_count,
            "Destabilizing_Residues": destabilizing_count,
            "Seq_Thermo_Score": seq_thermo_score,
            "Seq_Thermo_Score_Normalized": seq_thermo_score_normalized,
            "Seq_Estimated_Tm": seq_tm_estimate,
            "AI_Estimated_Tm": ai_tm,
            "AI_Confidence": ai_confidence,
            "Combined_Estimated_Tm": combined_tm,
            "Combined_Stability_Score": combined_score,
            "Instability_Index": analyzed_seq.instability_index()
        }

        return result

    except Exception as e:
        logger.error(f"Error analyzing {name}: {str(e)}")
        return {"Name": name, "Error": str(e)}

def validate_inputs():
    """Validate Fab sequences."""
    if not hasattr(inputs, 'fab_sequences'):
        raise ValueError("fab_sequences not found")

def main():
    try:
        validate_inputs()

        results = []
        for name, seq in inputs.fab_sequences.items():
            logger.info(f"Predicting thermostability for {name}")
            result = predict_thermostability(seq, name)
            results.append(result)

        df = pd.DataFrame(results)

        # Save results
        results_dir = Path("results/formulation")
        results_dir.mkdir(parents=True, exist_ok=True)
        output_path = results_dir / "10_thermostability.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")

        print("Thermostability Prediction:")
        print(df[['Name', 'Seq_Thermo_Score', 'Seq_Thermo_Score_Normalized', 'Seq_Estimated_Tm']].to_string(index=False))

        # Validation
        for _, row in df.iterrows():
            if row['Seq_Estimated_Tm'] < 40:
                logger.warning(f"{row['Name']}: Low predicted Tm ({row['Seq_Estimated_Tm']:.1f}°C)")
            if row['Seq_Thermo_Score_Normalized'] < 0:
                logger.warning(f"{row['Name']}: Negative thermostability score")

    except Exception as e:
        logger.error(f"Error in Experiment 10: {str(e)}")
        raise

if __name__ == "__main__":
    main()