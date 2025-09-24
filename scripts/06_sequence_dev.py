#!/usr/bin/env python3
"""
Experiment 6: Sequence-Based Developability Scoring

This script calculates developability metrics for Fab sequences including pI, instability index,
GRAVY score, and molecular weight.
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import scripts.inputs as inputs
import logging
from pathlib import Path

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def analyze_sequence(name, sequence):
    """
    Analyze a protein sequence for developability metrics.

    Args:
        name (str): Sequence name
        sequence (str): Amino acid sequence

    Returns:
        dict: Analysis results
    """
    try:
        analyzed_seq = ProteinAnalysis(sequence)

        result = {
            "Name": name,
            "Length": len(sequence),
            "MW": analyzed_seq.molecular_weight(),
            "pI": analyzed_seq.isoelectric_point(),
            "Instability_Index": analyzed_seq.instability_index(),
            "GRAVY": analyzed_seq.gravy(),
            "Aromaticity": analyzed_seq.aromaticity(),
            "Secondary_Structure_Fraction": analyzed_seq.secondary_structure_fraction()
        }

        # Add amino acid composition
        aa_comp = analyzed_seq.get_amino_acids_percent()
        result.update({f"Percent_{aa}": percent for aa, percent in aa_comp.items()})

        return result

    except Exception as e:
        logger.error(f"Error analyzing {name}: {str(e)}")
        return {"Name": name, "Error": str(e)}

def validate_inputs():
    """Validate that Fab sequences are present."""
    if not hasattr(inputs, 'fab_sequences'):
        raise ValueError("fab_sequences not found in inputs.py")

def main():
    try:
        validate_inputs()

        results = []
        for name, seq in inputs.fab_sequences.items():
            logger.info(f"Analyzing {name}")
            result = analyze_sequence(name, seq)
            results.append(result)

        df = pd.DataFrame(results)

        # Save results
        results_dir = Path("results/formulation")
        results_dir.mkdir(parents=True, exist_ok=True)
        output_path = results_dir / "06_sequence_developability.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")

        print("Sequence Developability Analysis:")
        print(df[['Name', 'Length', 'MW', 'pI', 'Instability_Index', 'GRAVY']].to_string(index=False))

        # Validation checks
        for _, row in df.iterrows():
            if pd.isna(row.get('Instability_Index')):
                continue
            if row['Instability_Index'] > 40:
                logger.warning(f"{row['Name']}: High instability index ({row['Instability_Index']:.1f})")
            if row['GRAVY'] > 0:
                logger.warning(f"{row['Name']}: Hydrophobic GRAVY ({row['GRAVY']:.2f})")
            if 6.0 <= row['pI'] <= 8.0:
                logger.warning(f"{row['Name']}: pI in neutral range ({row['pI']:.1f}) - potential aggregation")

    except Exception as e:
        logger.error(f"Error in Experiment 6: {str(e)}")
        raise

if __name__ == "__main__":
    main()