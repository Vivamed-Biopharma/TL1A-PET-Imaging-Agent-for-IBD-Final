#!/usr/bin/env python3
"""
Run all TL1A PET Imaging Agent experiments in sequence.

This script executes all 15 experiments with proper error handling and progress tracking.
"""

import subprocess
import sys
import logging
from pathlib import Path
from typing import List, Tuple

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('experiment_run.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

EXPERIMENTS = [
    ("01_physchem", "Foundational Physicochemical Profiling"),
    ("02_metabolism", "Activation & Metabolism Prediction"),
    ("03_toxicity_flags", "In-Silico Toxicity & Liability Flagging"),
    ("04_linker_flex", "Linker Flexibility Analysis"),
    ("05_mmp_analysis", "MMP Analysis"),
    ("06_sequence_dev", "Sequence-Based Developability Scoring"),
    ("07_agg_hotspots", "Aggregation Hotspot Prediction"),
    ("08_charge_pI", "Charge Distribution & pI Analysis"),
    ("09_flexibility_anm", "Local Flexibility Analysis with ANM"),
    ("10_thermostability", "Thermostability Prediction"),
    ("11_complex_model", "Complex Modeling"),
    ("12_interface_fingerprint", "Interface Fingerprinting"),
    ("13_ala_scanning", "Alanine Scanning"),
    ("14_immunogenicity", "Immunogenicity Prediction"),
    ("15_decision_scorecard", "Integrated Decision Scorecard")
]

def run_experiment(exp_script: str, exp_name: str) -> Tuple[bool, str]:
    """
    Run a single experiment.

    Args:
        exp_script: Script filename (without .py)
        exp_name: Human-readable experiment name

    Returns:
        Tuple of (success, error_message)
    """
    script_path = f"scripts/{exp_script}.py"

    if not Path(script_path).exists():
        error_msg = f"Script not found: {script_path}"
        logger.error(error_msg)
        return False, error_msg

    logger.info(f"Starting experiment: {exp_name}")
    logger.info(f"Running script: {script_path}")

    try:
        result = subprocess.run(
            [sys.executable, script_path],
            capture_output=True,
            text=True,
            timeout=1800,  # 30 minutes timeout
            cwd="."
        )

        if result.returncode == 0:
            logger.info(f"✅ Experiment {exp_script} completed successfully")
            if result.stdout:
                logger.debug(f"Output: {result.stdout}")
            return True, ""
        else:
            error_msg = f"Script failed with return code {result.returncode}"
            if result.stderr:
                error_msg += f"\nStderr: {result.stderr}"
            logger.error(error_msg)
            return False, error_msg

    except subprocess.TimeoutExpired:
        error_msg = f"Experiment {exp_script} timed out after 30 minutes"
        logger.error(error_msg)
        return False, error_msg
    except Exception as e:
        error_msg = f"Unexpected error running {exp_script}: {str(e)}"
        logger.error(error_msg)
        return False, error_msg

def create_summary_report(results: List[Tuple[str, str, bool, str]]) -> None:
    """
    Create a summary report of all experiment runs.

    Args:
        results: List of (script, name, success, error) tuples
    """
    report_path = "experiment_summary_report.md"
    logger.info(f"Creating summary report: {report_path}")

    successful = sum(1 for _, _, success, _ in results if success)
    total = len(results)

    with open(report_path, 'w') as f:
        f.write("# TL1A PET Imaging Agent - Experiment Run Summary\n\n")
        f.write(f"**Total Experiments:** {total}\n")
        f.write(f"**Successful:** {successful}\n")
        f.write(f"**Failed:** {total - successful}\n\n")

        f.write("## Results\n\n")
        f.write("| Experiment | Status | Details |\n")
        f.write("|------------|--------|---------|\n")

        for script, name, success, error in results:
            status = "✅ Pass" if success else "❌ Fail"
            details = "Completed successfully" if success else error.replace('\n', '<br>')
            f.write(f"| {script} | {status} | {details} |\n")

        f.write("\n## Failed Experiments\n\n")
        failed = [(script, name, error) for script, name, success, error in results if not success]
        if failed:
            for script, name, error in failed:
                f.write(f"### {script}: {name}\n")
                f.write(f"**Error:** {error}\n\n")
        else:
            f.write("All experiments completed successfully! 🎉\n")

    logger.info(f"Summary report saved to {report_path}")

def main():
    """Run all experiments."""
    logger.info("Starting TL1A PET Imaging Agent experiment suite")
    logger.info(f"Total experiments to run: {len(EXPERIMENTS)}")

    results = []

    for exp_script, exp_name in EXPERIMENTS:
        success, error = run_experiment(exp_script, exp_name)
        results.append((exp_script, exp_name, success, error))

        # Log progress
        completed = len(results)
        successful = sum(1 for _, _, s, _ in results if s)
        logger.info(f"Progress: {completed}/{len(EXPERIMENTS)} experiments completed ({successful} successful)")

    # Create summary report
    create_summary_report(results)

    # Final summary
    successful = sum(1 for _, _, s, _ in results if s)
    if successful == len(EXPERIMENTS):
        logger.info("🎉 All experiments completed successfully!")
        return 0
    else:
        logger.error(f"❌ {len(EXPERIMENTS) - successful} experiment(s) failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())