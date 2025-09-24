#!/usr/bin/env python3
"""
Experiment 12: Interface Fingerprinting

This script analyzes the binding interface fingerprints for Fab-TL1A complexes.
"""

import pandas as pd
import scripts.inputs as inputs
import logging
from pathlib import Path
import random

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def analyze_interface_fingerprint(fab_name, sequence):
    """
    Analyze interface fingerprint.

    Args:
        fab_name (str): Fab name
        sequence (str): Sequence

    Returns:
        dict: Fingerprint analysis
    """
    # Simulate interface analysis
    # In reality, this would analyze the modeled complex

    # CDR regions (approximate positions for Fab)
    cdr_h1 = sequence[26:35] if len(sequence) > 35 else ""
    cdr_h2 = sequence[50:65] if len(sequence) > 65 else ""
    cdr_h3 = sequence[95:110] if len(sequence) > 110 else ""
    cdr_l1 = sequence[24:34] if len(sequence) > 34 else ""
    cdr_l2 = sequence[50:56] if len(sequence) > 56 else ""
    cdr_l3 = sequence[89:97] if len(sequence) > 97 else ""

    cdrs = {
        "CDR_H1": cdr_h1,
        "CDR_H2": cdr_h2,
        "CDR_H3": cdr_h3,
        "CDR_L1": cdr_l1,
        "CDR_L2": cdr_l2,
        "CDR_L3": cdr_l3
    }

    # Analyze CDR composition
    interface_aa = ['Y', 'W', 'R', 'D', 'E', 'S', 'T', 'N', 'Q', 'G', 'A']
    fingerprint = {}

    for cdr_name, cdr_seq in cdrs.items():
        composition = {aa: cdr_seq.count(aa) for aa in interface_aa}
        fingerprint[f"{cdr_name}_Composition"] = composition
        fingerprint[f"{cdr_name}_Length"] = len(cdr_seq)
        fingerprint[f"{cdr_name}_Hydrophobic"] = sum(cdr_seq.count(aa) for aa in ['Y', 'W', 'F', 'I', 'L', 'V', 'M'])

    # Overall interface score
    total_interface_residues = sum(len(cdr) for cdr in cdrs.values())
    hydrophobic_ratio = sum(fingerprint[k] for k in fingerprint.keys() if 'Hydrophobic' in k) / total_interface_residues if total_interface_residues > 0 else 0

    result = {
        "Fab_Name": fab_name,
        "Total_Interface_Length": total_interface_residues,
        "Hydrophobic_Ratio": hydrophobic_ratio,
        "CDR_Counts": {k: v for k, v in fingerprint.items() if 'Length' in k}
    }
    result.update(fingerprint)

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
            logger.info(f"Analyzing interface fingerprint for {name}")
            result = analyze_interface_fingerprint(name, seq)
            results.append(result)

        # Create summary DataFrame
        summary_data = []
        for r in results:
            summary_data.append({
                "Fab_Name": r["Fab_Name"],
                "Total_Interface_Length": r["Total_Interface_Length"],
                "Hydrophobic_Ratio": r["Hydrophobic_Ratio"],
                "CDR_H3_Length": r.get("CDR_H3_Length", 0)
            })

        summary_df = pd.DataFrame(summary_data)

        # Save summary
        results_dir = Path("results/biobetter")
        results_dir.mkdir(parents=True, exist_ok=True)
        summary_path = results_dir / "12_interface_fingerprint.csv"
        summary_df.to_csv(summary_path, index=False)
        logger.info(f"Results saved to {summary_path}")

        print("Interface Fingerprint Analysis:")
        print(summary_df.to_string(index=False))

        # Validation
        for _, row in summary_df.iterrows():
            if row['Hydrophobic_Ratio'] > 0.6:
                logger.info(f"{row['Fab_Name']}: Hydrophobic-rich interface - may affect solubility")

    except Exception as e:
        logger.error(f"Error in Experiment 12: {str(e)}")
        raise

if __name__ == "__main__":
    main()