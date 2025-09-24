#!/usr/bin/env python3
"""
Experiment 7: Aggregation Hotspot Prediction

This script predicts aggregation-prone regions in Fab sequences using pattern-based analysis.
In a real implementation, this would use Aggrescan3D or similar AI tools.
"""

import re
import pandas as pd
import scripts.inputs as inputs
import logging
from pathlib import Path
from scripts.neurosnap_wrappers import predict_aggregation

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Aggregation-prone patterns (simplified from literature)
AGGREGATION_PATTERNS = {
    "Amyloid_Like": r"V[VIL]V[VIL]V",  # VVV, VIV, etc.
    "Beta_Sheet": r"[VILM][VILM][VILM][VILM]",  # Four hydrophobic residues
    "Hydrophobic_Patch": r"[FILVWYM][FILVWYM][FILVWYM]",  # Three hydrophobic
    "Proline_Break": r"P[FILVWYM]{3,}",  # Proline followed by hydrophobic stretch
}

def scan_aggregation_hotspots(sequence, name):
    """
    Scan sequence for aggregation-prone regions using both pattern matching and AI prediction.

    Args:
        sequence (str): Amino acid sequence
        name (str): Sequence name

    Returns:
        dict: Analysis results
    """
    # Pattern-based analysis (fast pre-screen)
    pattern_hotspots = []
    pattern_score = 0

    for pattern_name, pattern in AGGREGATION_PATTERNS.items():
        matches = re.finditer(pattern, sequence)
        count = 0
        positions = []

        for match in matches:
            count += 1
            start = match.start()
            end = match.end()
            positions.append(f"{start+1}-{end}")  # 1-based indexing
            pattern_score += len(match.group())  # Score based on pattern length

        if count > 0:
            pattern_hotspots.append({
                "Pattern": pattern_name,
                "Count": count,
                "Positions": "; ".join(positions)
            })

    # AI-based prediction using Aggrescan3D
    try:
        ai_results = predict_aggregation(sequence)
        ai_score = ai_results.get("aggregation_score", 0)
        ai_hotspots = ai_results.get("hotspots", [])
        ai_confidence = ai_results.get("confidence", 0.0)
    except Exception as e:
        logger.warning(f"AI aggregation prediction failed for {name}: {str(e)}")
        ai_score = 0
        ai_hotspots = []
        ai_confidence = 0.0

    # Combine results
    combined_hotspots = pattern_hotspots + ai_hotspots
    total_score = pattern_score + ai_score

    result = {
        "Name": name,
        "Sequence_Length": len(sequence),
        "Pattern_Hotspots": sum(h['Count'] for h in pattern_hotspots),
        "AI_Hotspots": len(ai_hotspots),
        "Total_Hotspots": sum(h['Count'] for h in pattern_hotspots) + len(ai_hotspots),
        "Pattern_Score": pattern_score,
        "AI_Score": ai_score,
        "Aggregation_Score": total_score,
        "AI_Confidence": ai_confidence,
        "Hotspots_Details": combined_hotspots
    }

    return result

def validate_inputs():
    """Validate Fab sequences."""
    if not hasattr(inputs, 'fab_sequences'):
        raise ValueError("fab_sequences not found")

def main():
    try:
        validate_inputs()

        results = []
        for name, seq in inputs.fab_sequences.items():
            logger.info(f"Scanning aggregation hotspots in {name}")
            result = scan_aggregation_hotspots(seq, name)
            results.append(result)

        # Create summary DataFrame
        summary_data = []
        for r in results:
            summary_data.append({
                "Name": r["Name"],
                "Sequence_Length": r["Sequence_Length"],
                "Total_Hotspots": r["Total_Hotspots"],
                "Aggregation_Score": r["Aggregation_Score"]
            })

        summary_df = pd.DataFrame(summary_data)

        # Save summary
        results_dir = Path("results/formulation")
        results_dir.mkdir(parents=True, exist_ok=True)
        summary_path = results_dir / "07_agg_hotspots_summary.csv"
        summary_df.to_csv(summary_path, index=False)

        # Save detailed results
        detailed_data = []
        for r in results:
            for hotspot in r["Hotspots_Details"]:
                detailed_data.append({
                    "Name": r["Name"],
                    "Pattern": hotspot["Pattern"],
                    "Count": hotspot["Count"],
                    "Positions": hotspot["Positions"]
                })

        if detailed_data:
            detailed_df = pd.DataFrame(detailed_data)
            detailed_path = results_dir / "07_agg_hotspots_detailed.csv"
            detailed_df.to_csv(detailed_path, index=False)
            logger.info(f"Detailed results saved to {detailed_path}")

        logger.info(f"Summary saved to {summary_path}")

        print("Aggregation Hotspot Analysis Summary:")
        print(summary_df.to_string(index=False))

        # Validation
        for _, row in summary_df.iterrows():
            if row['Aggregation_Score'] > 20:
                logger.warning(f"{row['Name']}: High aggregation risk (score: {row['Aggregation_Score']})")

    except Exception as e:
        logger.error(f"Error in Experiment 7: {str(e)}")
        raise

if __name__ == "__main__":
    main()