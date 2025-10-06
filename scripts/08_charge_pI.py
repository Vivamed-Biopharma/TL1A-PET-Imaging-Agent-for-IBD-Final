#!/usr/bin/env python3
"""
Experiment 8: Charge Distribution & pI Analysis

This script analyzes charge distribution and isoelectric points for Fab sequences,
including per-residue charge analysis.
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
try:
    import scripts.inputs as inputs
except ModuleNotFoundError:
    import inputs as inputs
import logging
from pathlib import Path

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Amino acid charge properties
AA_CHARGES = {
    'D': -1, 'E': -1, 'K': +1, 'R': +1, 'H': +1,  # Charged
    'C': -1, 'Y': -1,  # Potentially charged
    # Others: 0
}

def analyze_charge_distribution(sequence, name):
    """
    Analyze charge distribution in sequence.

    Args:
        sequence (str): Amino acid sequence
        name (str): Sequence name

    Returns:
        dict: Charge analysis
    """
    try:
        analyzed_seq = ProteinAnalysis(sequence)

        # Count charged residues
        charge_counts = {'Positive': 0, 'Negative': 0, 'Neutral': 0}
        charge_positions = {'Positive': [], 'Negative': []}

        for i, aa in enumerate(sequence):
            if aa in ['K', 'R', 'H']:
                charge_counts['Positive'] += 1
                charge_positions['Positive'].append(i+1)  # 1-based
            elif aa in ['D', 'E']:
                charge_counts['Negative'] += 1
                charge_positions['Negative'].append(i+1)

        charge_counts['Neutral'] = len(sequence) - sum(charge_counts.values())

        # Calculate charge-related properties
        net_charge = charge_counts['Positive'] - charge_counts['Negative']
        charge_density = net_charge / len(sequence)

        result = {
            "Name": name,
            "Length": len(sequence),
            "pI": analyzed_seq.isoelectric_point(),
            "Net_Charge": net_charge,
            "Charge_Density": charge_density,
            "Positive_Residues": charge_counts['Positive'],
            "Negative_Residues": charge_counts['Negative'],
            "Neutral_Residues": charge_counts['Neutral'],
            "Positive_Positions": "; ".join(map(str, charge_positions['Positive'])),
            "Negative_Positions": "; ".join(map(str, charge_positions['Negative']))
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
            logger.info(f"Analyzing charge distribution in {name}")
            result = analyze_charge_distribution(seq, name)
            results.append(result)

        df = pd.DataFrame(results)

        # Save results
        results_dir = Path("results/formulation")
        results_dir.mkdir(parents=True, exist_ok=True)
        output_path = results_dir / "08_charge_distribution.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")

        print("Charge Distribution Analysis:")
        print(df[['Name', 'pI', 'Net_Charge', 'Charge_Density', 'Positive_Residues', 'Negative_Residues']].to_string(index=False))

        # Validation checks
        for _, row in df.iterrows():
            if 6.0 <= row['pI'] <= 8.0:
                logger.warning(f"{row['Name']}: pI in aggregation-prone range ({row['pI']:.2f})")
            if abs(row['Net_Charge']) < 2:
                logger.info(f"{row['Name']}: Low net charge - may affect solubility")

    except Exception as e:
        logger.error(f"Error in Experiment 8: {str(e)}")
        raise

if __name__ == "__main__":
    main()