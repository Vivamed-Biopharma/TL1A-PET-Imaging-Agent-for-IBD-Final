#!/usr/bin/env python3
"""
Experiment 15: Integrated Decision Scorecard

This script creates a comprehensive scorecard integrating all experimental results
to guide decision-making for Fab selection and optimization.
"""

import pandas as pd
try:
    import scripts.inputs as inputs
except ModuleNotFoundError:
    import inputs as inputs
import logging
from pathlib import Path
import os

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_experiment_results():
    """
    Load results from all experiments.

    Returns:
        dict: Results by experiment
    """
    results_dir = Path("results")
    experiment_results = {}

    # Define expected result files
    expected_files = {
        "physchem": "prodrug/01_physchem_properties.csv",
        "metabolism": "prodrug/02_biotransformer_metabolites.csv",
        "toxicity": "prodrug/03_liability_hits.csv",
        "flexibility": "prodrug/04_linker_flexibility.csv",
        "mmp": "prodrug/05_mmp_changes.csv",
        "developability": "formulation/06_sequence_developability.csv",
        "aggregation": "formulation/07_agg_hotspots_summary.csv",
        "charge": "formulation/08_charge_distribution.csv",
        "anm_flex": "formulation/09_flexibility_summary.csv",
        "thermostability": "formulation/10_thermostability.csv",
        "complex": "biobetter/11_complex_modeling.csv",
        "interface": "biobetter/12_interface_fingerprint.csv",
        "scanning": "biobetter/13_scanning_summary.csv",
        "immunogenicity": "biobetter/14_immunogenicity.csv"
    }

    for exp_name, file_path in expected_files.items():
        full_path = results_dir / file_path
        if full_path.exists():
            try:
                df = pd.read_csv(full_path)
                experiment_results[exp_name] = df
                logger.info(f"Loaded {exp_name} results: {len(df)} rows")
            except Exception as e:
                logger.warning(f"Failed to load {exp_name}: {str(e)}")
        else:
            logger.warning(f"Missing results file: {full_path}")

    return experiment_results

def calculate_decision_scores(experiment_results):
    """
    Calculate integrated decision scores for each Fab.

    Args:
        experiment_results (dict): Results from all experiments

    Returns:
        pd.DataFrame: Decision scorecard
    """
    fab_names = list(inputs.fab_sequences.keys())
    scorecard_data = []

    for fab in fab_names:
        scores = {"Fab_Name": fab}

        # Physicochemical score (for linker, but we'll adapt)
        if "physchem" in experiment_results:
            physchem = experiment_results["physchem"]
            linker_row = physchem[physchem["Name"] == "Linker_Chelator"]
            if len(linker_row) > 0:
                scores["Linker_LogP"] = linker_row["LogP"].values[0]
                scores["Linker_TPSA"] = linker_row["TPSA"].values[0]

        # Developability score
        if "developability" in experiment_results:
            dev = experiment_results["developability"]
            fab_row = dev[dev["Name"] == fab]
            if len(fab_row) > 0:
                scores["Instability_Index"] = fab_row["Instability_Index"].values[0]
                scores["GRAVY"] = fab_row["GRAVY"].values[0]
                scores["pI"] = fab_row["pI"].values[0]

        # Aggregation score
        if "aggregation" in experiment_results:
            agg = experiment_results["aggregation"]
            fab_row = agg[agg["Name"] == fab]
            if len(fab_row) > 0:
                scores["Aggregation_Score"] = fab_row["Aggregation_Score"].values[0]

        # Charge score
        if "charge" in experiment_results:
            charge = experiment_results["charge"]
            fab_row = charge[charge["Name"] == fab]
            if len(fab_row) > 0:
                scores["Net_Charge"] = fab_row["Net_Charge"].values[0]

        # Thermostability
        if "thermostability" in experiment_results:
            thermo = experiment_results["thermostability"]
            fab_row = thermo[thermo["Name"] == fab]
            if len(fab_row) > 0:
                scores["Combined_Stability_Score"] = fab_row["Combined_Stability_Score"].values[0]
                scores["Combined_Estimated_Tm"] = fab_row["Combined_Estimated_Tm"].values[0]

        # Immunogenicity
        if "immunogenicity" in experiment_results:
            immuno = experiment_results["immunogenicity"]
            fab_row = immuno[immuno["Fab_Name"] == fab]
            if len(fab_row) > 0:
                scores["Combined_Immunogenicity_Score"] = fab_row["Combined_Immunogenicity_Score"].values[0]
                scores["Risk_Level"] = fab_row["Risk_Level"].values[0]

        # Calculate overall score (lower is better)
        overall_score = 0
        if "Instability_Index" in scores:
            overall_score += max(0, scores["Instability_Index"] - 30)  # Penalty for >30
        if "GRAVY" in scores:
            overall_score += max(0, scores["GRAVY"]) * 10  # Penalty for hydrophobicity
        if "Aggregation_Score" in scores:
            overall_score += scores["Aggregation_Score"] * 0.1
        if "Combined_Immunogenicity_Score" in scores:
            overall_score += scores["Combined_Immunogenicity_Score"] * 100
        if "Combined_Stability_Score" in scores:
            overall_score += abs(scores["Combined_Stability_Score"]) * 0.5  # Penalty for instability

        scores["Overall_Score"] = overall_score

        # Recommendation
        if overall_score < 5:
            scores["Recommendation"] = "Advance"
        elif overall_score < 15:
            scores["Recommendation"] = "Review"
        else:
            scores["Recommendation"] = "Optimize"

        scorecard_data.append(scores)

    return pd.DataFrame(scorecard_data)

def create_html_report(df):
    """
    Create a styled HTML report.

    Args:
        df (pd.DataFrame): Scorecard DataFrame
    """
    # Create styled HTML
    html_content = f"""
    <html>
    <head>
        <title>TL1A PET Imaging Agent - Decision Scorecard</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
            .advance {{ background-color: #d4edda; }}
            .review {{ background-color: #fff3cd; }}
            .optimize {{ background-color: #f8d7da; }}
        </style>
    </head>
    <body>
        <h1>TL1A PET Imaging Agent - Integrated Decision Scorecard</h1>
        <p>This scorecard integrates results from all 14 experiments to guide Fab selection.</p>
        {df.to_html(index=False, classes='table table-striped')}
    </body>
    </html>
    """

    output_path = Path("results/biobetter/15_decision_scorecard.html")
    with open(output_path, 'w') as f:
        f.write(html_content)
    logger.info(f"HTML report saved to {output_path}")

def main():
    try:
        logger.info("Generating integrated decision scorecard")

        # Load all results
        experiment_results = load_experiment_results()

        # Calculate scores
        scorecard_df = calculate_decision_scores(experiment_results)

        # Save CSV
        results_dir = Path("results/biobetter")
        results_dir.mkdir(parents=True, exist_ok=True)
        csv_path = results_dir / "15_decision_scorecard.csv"
        scorecard_df.to_csv(csv_path, index=False)
        logger.info(f"Scorecard saved to {csv_path}")

        # Create HTML report
        create_html_report(scorecard_df)

        print("Decision Scorecard:")
        print(scorecard_df.to_string(index=False))

        # Final recommendation
        best_fab = scorecard_df.loc[scorecard_df["Overall_Score"].idxmin()]
        logger.info(f"Recommended Fab: {best_fab['Fab_Name']} (Score: {best_fab['Overall_Score']:.2f})")

    except Exception as e:
        logger.error(f"Error in Experiment 15: {str(e)}")
        raise

if __name__ == "__main__":
    main()