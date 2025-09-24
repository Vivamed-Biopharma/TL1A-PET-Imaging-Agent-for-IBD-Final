#!/usr/bin/env python3
"""
Experiment 13: Alanine Scanning

This script performs in-silico alanine scanning on Fab sequences to identify
key binding residues. In a real implementation, this would use StaB-ddG.
"""

import pandas as pd
import scripts.inputs as inputs
import logging
from pathlib import Path
from scripts.neurosnap_wrappers import predict_stability_change

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def perform_alanine_scanning(fab_name, sequence):
    """
    Perform alanine scanning on CDR regions using StaB-ddG.

    Args:
        fab_name (str): Fab name
        sequence (str): Sequence

    Returns:
        list: Scanning results
    """
    # Define CDR positions (approximate for Fab)
    cdr_regions = {
        "CDR_H1": (26, 35),
        "CDR_H2": (50, 65),
        "CDR_H3": (95, 110),
        "CDR_L1": (24, 34),
        "CDR_L2": (50, 56),
        "CDR_L3": (89, 97)
    }

    results = []
    all_mutations = []

    for cdr_name, (start, end) in cdr_regions.items():
        if end > len(sequence):
            continue

        cdr_seq = sequence[start:end]
        for i, aa in enumerate(cdr_seq):
            if aa == 'A':  # Skip if already alanine
                continue

            position = start + i + 1  # 1-based
            mutation = f"{aa}{position}A"
            all_mutations.append(mutation)

    if all_mutations:
        try:
            # Use StaB-ddG for real predictions
            stab_results = predict_stability_change(sequence, all_mutations, max_wait_time=1800)

            # Map back to individual mutations
            mutation_ddg = stab_results.get("ddg_values", {})

            idx = 0
            for cdr_name, (start, end) in cdr_regions.items():
                if end > len(sequence):
                    continue

                cdr_seq = sequence[start:end]
                for i, aa in enumerate(cdr_seq):
                    if aa == 'A':
                        continue

                    position = start + i + 1
                    mutation = f"{aa}{position}A"
                    ddG = mutation_ddg.get(mutation, 0.0)  # Default to 0 if not found

                    results.append({
                        "Fab_Name": fab_name,
                        "CDR": cdr_name,
                        "Position": position,
                        "Wild_Type": aa,
                        "Mutant": "A",
                        "Predicted_ddG": ddG,
                        "Binding_Impact": "Disrupting" if ddG > 0.5 else "Neutral" if ddG > -0.5 else "Stabilizing"
                    })
                    idx += 1

        except Exception as e:
            logger.error(f"StaB-ddG prediction failed for {fab_name}: {str(e)}")
            raise RuntimeError(f"Alanine scanning failed for {fab_name}: {str(e)}")
    else:
        logger.info(f"No mutations to scan for {fab_name}")

    return results

def validate_inputs():
    """Validate sequences."""
    if not hasattr(inputs, 'fab_sequences'):
        raise ValueError("fab_sequences not found")

def main():
    try:
        validate_inputs()

        all_results = []
        for name, seq in inputs.fab_sequences.items():
            logger.info(f"Performing alanine scanning for {name}")
            results = perform_alanine_scanning(name, seq)
            all_results.extend(results)

        df = pd.DataFrame(all_results)

        # Save results
        results_dir = Path("results/biobetter")
        results_dir.mkdir(parents=True, exist_ok=True)
        output_path = results_dir / "13_alanine_scanning.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")

        # Summary by Fab
        summary = df.groupby('Fab_Name').agg({
            'Predicted_ddG': ['mean', 'min', 'max', 'count'],
            'Binding_Impact': lambda x: (x == 'Disrupting').sum()
        }).round(2)

        summary.columns = ['Mean_ddG', 'Min_ddG', 'Max_ddG', 'Total_Mutations', 'Disrupting_Mutations']
        summary_path = results_dir / "13_scanning_summary.csv"
        summary.to_csv(summary_path)
        logger.info(f"Summary saved to {summary_path}")

        print("Alanine Scanning Summary:")
        print(summary.to_string())

        # Identify key residues
        disrupting = df[df['Binding_Impact'] == 'Disrupting']
        if len(disrupting) > 0:
            print(f"\nTop disrupting mutations:")
            print(disrupting.nlargest(5, 'Predicted_ddG')[['Fab_Name', 'Position', 'Wild_Type', 'Predicted_ddG']].to_string(index=False))

    except Exception as e:
        logger.error(f"Error in Experiment 13: {str(e)}")
        raise

if __name__ == "__main__":
    main()